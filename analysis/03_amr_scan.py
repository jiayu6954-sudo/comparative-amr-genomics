"""
Step 3 — Species-aware AMR gene detection + CARD sequence verification
───────────────────────────────────────────────────────────────────────
Tier 1: GFF3 regex (species-specific patterns)
Tier 2: phmmer against CARD v3.3 protein homolog models

Output: data/processed/{short}_amr_hits.tsv

Usage:
  python analysis/03_amr_scan.py [--species eco|aba|both]
"""
import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, DATA_VAL, LOGS, SPECIES_CONFIG

LOGS.mkdir(parents=True, exist_ok=True)
DATA_PROC.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'amr_scan.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('amr_scan')

# CARD database path (reuse from KPN project if available)
CARD_DB_CANDIDATES = [
    Path(__file__).resolve().parents[2] / 'amr_project' / 'tools' /
    'card' / 'protein_fasta_protein_homolog_model.fasta',
    Path(__file__).resolve().parents[1] / 'tools' /
    'card' / 'protein_fasta_protein_homolog_model.fasta',
]
CARD_DB = next((p for p in CARD_DB_CANDIDATES if p.exists()), None)


# ── GFF3 parsing ──────────────────────────────────────────────────────────────

def _parse_attr(attr_str: str) -> dict:
    attrs = {}
    for part in attr_str.strip().split(';'):
        if '=' in part:
            k, _, v = part.partition('=')
            attrs[k.strip()] = v.strip()
    return attrs


def _normalize_gene_name(raw: str) -> str:
    """Extract canonical allele name from raw GFF product/gene string."""
    raw = raw.strip()
    # bla-prefix priority: blaKPC-2, blaNDM-1 etc.
    m = re.search(r'bla([A-Z]{2,5}-\d+)', raw, re.IGNORECASE)
    if m:
        return m.group(1).upper()
    # NDM-\d, KPC-\d, OXA-\d, IMP-\d, VIM-\d
    m = re.search(r'\b(NDM|KPC|OXA|IMP|VIM)-(\d+)\b', raw, re.IGNORECASE)
    if m:
        return f'{m.group(1).upper()}-{m.group(2)}'
    return raw[:60]


def scan_gff_amr(gff_path: Path, patterns: list, gene_set: set) -> list:
    """Return list of hit dicts from a single GFF3 file."""
    compiled = [re.compile(p, re.IGNORECASE) for p in patterns]
    gene_lower = {g.lower() for g in gene_set}
    hits = []
    try:
        with gzip.open(gff_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue
                attrs = _parse_attr(parts[8])
                product = attrs.get('product', '')
                gene    = attrs.get('gene', '')
                note    = attrs.get('Note', '')
                search_text = f'{product} {gene} {note}'

                matched = any(pat.search(search_text) for pat in compiled)
                if not matched:
                    continue

                gene_name = _normalize_gene_name(product or gene)
                if gene_name.lower() not in gene_lower and \
                        not any(pat.search(gene_name) for pat in compiled):
                    # keep anyway — pattern matched, may be novel allele
                    pass

                hits.append({
                    'contig':     parts[0],
                    'start':      int(parts[3]),
                    'end':        int(parts[4]),
                    'strand':     parts[6],
                    'gene_name':  gene_name,
                    'product':    product,
                    'confidence': 'GFF_KEYWORD',
                    'card_status': 'PENDING',
                })
    except Exception as e:
        log.warning(f'GFF parse error {gff_path}: {e}')
    return hits


# ── CARD sequence verification ────────────────────────────────────────────────

def _extract_protein(fna_path: Path, contig: str, start: int, end: int,
                     strand: str) -> str | None:
    """Extract and translate CDS from genome FASTA."""
    try:
        seqs: dict[str, list] = {}
        current = None
        with gzip.open(fna_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                line = line.rstrip()
                if line.startswith('>'):
                    current = line[1:].split()[0]
                    seqs[current] = []
                elif current:
                    seqs[current].append(line)

        seq_str = ''.join(seqs.get(contig, []))
        if not seq_str:
            return None

        # 0-based slicing
        s, e = start - 1, end
        subseq = seq_str[s:e].upper()
        if strand == '-':
            subseq = subseq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]

        codon_table = {
            'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L',
            'CTA':'L','CTG':'L','ATT':'I','ATC':'I','ATA':'I','ATG':'M',
            'GTT':'V','GTC':'V','GTA':'V','GTG':'V','TCT':'S','TCC':'S',
            'TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
            'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A',
            'GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
            'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','AAT':'N','AAC':'N',
            'AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
            'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R',
            'CGA':'R','CGG':'R','AGT':'S','AGC':'S','AGA':'R','AGG':'R',
            'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
        }
        aa = []
        for i in range(0, len(subseq) - 2, 3):
            codon = subseq[i:i+3]
            aa.append(codon_table.get(codon, 'X'))
        protein = ''.join(aa).rstrip('*')
        return protein if len(protein) >= 50 else None
    except Exception:
        return None


def verify_with_card(hits: list, fna_path: Path, card_db: Path) -> list:
    """Run phmmer against CARD for each hit. Returns updated hits."""
    if card_db is None or not card_db.exists():
        log.warning('CARD database not found — skipping sequence verification')
        for h in hits:
            h['card_status'] = 'NO_CARD_DB'
        return hits

    try:
        import pyhmmer
    except ImportError:
        log.warning('pyhmmer not installed — skipping CARD verification')
        for h in hits:
            h['card_status'] = 'NO_PYHMMER'
        return hits

    # Load CARD sequences once
    card_seqs = {}
    try:
        with pyhmmer.easel.SequenceFile(str(card_db), digital=True) as sf:
            for seq in sf:
                name = seq.name.decode() if isinstance(seq.name, bytes) else seq.name
                card_seqs[name] = seq
    except Exception as e:
        log.warning(f'CARD DB load error: {e}')
        for h in hits:
            h['card_status'] = 'CARD_LOAD_ERROR'
        return hits

    alphabet = pyhmmer.easel.Alphabet.amino()

    for h in hits:
        protein = _extract_protein(fna_path, h['contig'],
                                   h['start'], h['end'], h['strand'])
        if not protein:
            h['card_status'] = 'EXTRACT_FAIL'
            continue

        try:
            query_seq = pyhmmer.easel.TextSequence(
                name=b'query', sequence=protein).digitize(alphabet)
            targets_list = list(card_seqs.values())
            all_hits = list(pyhmmer.hmmer.phmmer(
                query_seq, targets_list, E=1e-5))
            if all_hits and len(all_hits) > 0:
                top = all_hits[0]
                hit_name = top.name.decode() if isinstance(top.name, bytes) else str(top.name)
                gene_upper = h['gene_name'].upper().replace('-', '').replace('BLA', '')
                hit_upper  = hit_name.upper().replace('-', '').replace('BLA', '')
                if gene_upper[:3] in hit_upper or hit_upper[:3] in gene_upper:
                    h['card_status'] = 'CONFIRMED'
                else:
                    h['card_status'] = 'NAME_MISMATCH'
                h['card_hit'] = hit_name
            else:
                h['card_status'] = 'NO_HIT'
        except Exception as e:
            h['card_status'] = f'ERROR:{e}'

    return hits


# ── main ──────────────────────────────────────────────────────────────────────

def run_species(species_name: str, cfg: dict):
    short = cfg['short']
    manifest_path = DATA_VAL / f'{short}_manifest.tsv'
    if not manifest_path.exists():
        log.warning(f'No manifest for {species_name} — run 02_validate.py first')
        return

    manifest = pd.read_csv(manifest_path, sep='\t', dtype=str)
    log.info(f'{species_name}: scanning {len(manifest)} genomes for AMR genes…')

    all_hits = []
    patterns  = cfg['amr_gff_keywords']
    gene_set  = cfg['carbapenem_genes']

    for _, row in manifest.iterrows():
        acc = row['accession']
        gff = Path(row['gff_path'])
        fna = Path(row['fna_path'])
        if not gff.exists():
            continue

        hits = scan_gff_amr(gff, patterns, gene_set)
        if hits and fna.exists() and CARD_DB:
            hits = verify_with_card(hits, fna, CARD_DB)
        for h in hits:
            h['accession'] = acc
            h['species']   = species_name
        all_hits.extend(hits)

    df = pd.DataFrame(all_hits) if all_hits else pd.DataFrame()
    out = DATA_PROC / f'{short}_amr_hits.tsv'
    df.to_csv(out, sep='\t', index=False)

    n_confirmed = (df['card_status'] == 'CONFIRMED').sum() if not df.empty else 0
    log.info(f'{species_name}: {len(df)} hits, '
             f'{df["accession"].nunique() if not df.empty else 0} resistant genomes, '
             f'{n_confirmed} CARD-confirmed → {out}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', choices=['eco', 'aba', 'both'],
                        default='both')
    args = parser.parse_args()

    for name, cfg in SPECIES_CONFIG.items():
        if args.species == 'both' or cfg['short'] == args.species:
            run_species(name, cfg)


if __name__ == '__main__':
    main()
