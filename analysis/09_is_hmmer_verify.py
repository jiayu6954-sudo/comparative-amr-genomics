"""
Step 9 — IS element family verification using PFAM HMM (PF01527, IS6 family)
─────────────────────────────────────────────────────────────────────────────
Addresses annotation bias: IS26/IS6-family transposases (ECO) and ISAba1
(ABA, also IS6-family) may be annotated only as 'transposase' or 'insertion
sequence' without a specific family name, causing their counts to be missed.

Method:
  1. For each species (eco/aba), scan manifested GFF3 files for IS/transposase
     features, separating family-named vs IS_unknown entries
  2. Extract protein sequences for IS_unknown features (FNA + GFF3 coords)
  3. Run pyhmmer.hmmsearch — PF01527 (IS6/ISAba1 transposase) vs IS_unknown proteins
  4. Reclassify IS_unknown → IS6_PFAM where E ≤ 1e-5
  5. Recount focal IS per genome (name-based + PFAM-reclassified)
  6. Re-run Mann-Whitney U, Cliff's delta, AUC with corrected counts
  7. Report annotation bias direction per species

Note: ISAba1 belongs to the IS6 superfamily (PF01527), same as IS26/IS257.
      Both ECO (IS26) and ABA (ISAba1) focal IS elements use PF01527.

Outputs (per species):
  data/processed/{short}_is_hmmer_results.tsv    per-feature HMMER hits
  data/processed/{short}_is_burden_corrected.tsv  corrected focal IS counts
  data/processed/{short}_is_burden_corrected_stats.json  AUC etc.
  logs/is_hmmer_{short}.log

Usage:
  python analysis/09_is_hmmer_verify.py [--species eco|aba|both]
"""
import argparse
import gzip
import json
import logging
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests

try:
    import pyhmmer.easel
    import pyhmmer.hmmer
    import pyhmmer.plan7
except ImportError:
    print('ERROR: pyhmmer not installed. Run: pip install pyhmmer', file=sys.stderr)
    sys.exit(1)

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_PROC, DATA_VAL, LOGS, SPECIES_CONFIG,
    PFAM_IS6_ID, HMMER_EVALUE, HMMER_MIN_IDENTITY,
)

# ── constants ─────────────────────────────────────────────────────────────────
DB_DIR   = Path(__file__).resolve().parents[1] / 'data' / 'db' / 'pfam'
HMM_PATH = DB_DIR / f'{PFAM_IS6_ID}.hmm'
MIN_PROT_LEN = 50   # residues; skip fragments

_PFAM_URLS = [
    f'https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{PFAM_IS6_ID}/?annotation=hmm',
    f'https://www.ebi.ac.uk/interpro/api/entry/pfam/{PFAM_IS6_ID}/?annotation=hmm&download',
]

_ATTRS_RE = re.compile(r'(\w+)=([^;]+)')

_CODON = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}
_RC = str.maketrans('ACGTacgt', 'TGCAtgca')


# ── sequence utilities ────────────────────────────────────────────────────────

def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]


def translate(nuc: str) -> str:
    aa = []
    for i in range(0, len(nuc) - 2, 3):
        codon = nuc[i:i + 3].upper()
        aa.append(_CODON.get(codon, 'X'))
    return ''.join(aa).rstrip('*')


def load_contigs(fna_gz: Path, log) -> dict[str, str]:
    contigs: dict[str, str] = {}
    cid, parts = None, []
    try:
        with gzip.open(fna_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('>'):
                    if cid:
                        contigs[cid] = ''.join(parts)
                    cid = line[1:].split()[0]
                    parts = []
                elif cid:
                    parts.append(line)
        if cid:
            contigs[cid] = ''.join(parts)
    except Exception as e:
        log.error(f'FNA load error ({fna_gz.name}): {e}')
    return contigs


def extract_protein(contigs: dict, contig: str, start: int, stop: int,
                    strand: str) -> str | None:
    seq = contigs.get(contig)
    if not seq or stop > len(seq):
        return None
    nuc = seq[start:stop]
    if strand == '-':
        nuc = revcomp(nuc)
    if len(nuc) < MIN_PROT_LEN * 3:
        return None
    prot = translate(nuc)
    return prot if len(prot) >= MIN_PROT_LEN else None


def _domain_identity(dom) -> float:
    aln = dom.alignment
    total = matched = 0
    for h, t in zip(aln.hmm_sequence, aln.target_sequence):
        if h != '-' and t != '-':
            total += 1
            if h.upper() == t.upper():
                matched += 1
    return matched / total if total else 0.0


# ── GFF3 parsing ──────────────────────────────────────────────────────────────

def _parse_attrs(attr_str: str) -> dict:
    return dict(m.groups() for m in _ATTRS_RE.finditer(attr_str))


def _is_family_from_text(text: str, focal_pattern: re.Pattern) -> str:
    """Extract IS family name; prefer focal-pattern match, then generic."""
    m = focal_pattern.search(text)
    if m:
        return m.group(0)
    m = re.search(r'\b(IS[A-Za-z0-9]+)\b', text, re.IGNORECASE)
    if m:
        return m.group(1)
    m = re.search(r'\b(Tn\d+)\b', text, re.IGNORECASE)
    if m:
        return m.group(1)
    return 'IS_unknown'


def scan_is_features(gff_gz: Path, is_kw_patterns: list,
                     focal_pattern: re.Pattern, log) -> list[dict]:
    features = []
    try:
        with gzip.open(gff_gz, 'rt', encoding='utf-8', errors='replace') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 9:
                    continue
                feat_type = cols[2]
                if feat_type not in ('CDS', 'gene', 'mobile_element', 'repeat_region'):
                    continue
                attrs   = _parse_attrs(cols[8])
                product = attrs.get('product', '')
                gene    = attrs.get('gene', '')
                note    = attrs.get('Note', '')
                combined = f'{product} {gene} {note}'
                if not any(p.search(combined) for p in is_kw_patterns):
                    continue
                try:
                    start = int(cols[3]) - 1
                    stop  = int(cols[4])
                except ValueError:
                    continue
                features.append({
                    'contig':    cols[0],
                    'start':     start,
                    'stop':      stop,
                    'strand':    cols[6],
                    'is_family': _is_family_from_text(combined, focal_pattern),
                    'product':   product[:80],
                    'feat_type': feat_type,
                })
    except Exception as e:
        log.error(f'GFF scan error ({gff_gz.name}): {e}')
    return features


# ── PFAM HMM management ───────────────────────────────────────────────────────

def ensure_pfam_hmm(log) -> list:
    import gzip as _gzip
    if not HMM_PATH.exists():
        log.info(f'{HMM_PATH} not found — downloading from EBI …')
        downloaded = False
        for url in _PFAM_URLS:
            try:
                resp = requests.get(url, timeout=120,
                                    headers={'User-Agent': 'ComparativeAMR/1.0'})
                if not resp.ok:
                    log.warning(f'HTTP {resp.status_code} from {url}')
                    continue
                data = resp.content
                if data[:2] == b'\x1f\x8b':
                    data = _gzip.decompress(data)
                if b'HMMER3' not in data[:20]:
                    log.warning(f'{url}: not a HMMER3 HMM file')
                    continue
                HMM_PATH.write_bytes(data)
                log.info(f'Saved {HMM_PATH} ({len(data):,} bytes)')
                downloaded = True
                break
            except Exception as e:
                log.warning(f'Download failed ({url}): {e}')
        if not downloaded:
            log.error(
                f'Could not auto-download {PFAM_IS6_ID}.hmm.\n'
                f'Manual: python -c "import requests,gzip; '
                f'r=requests.get(\'{_PFAM_URLS[0]}\'); '
                f'open(\'{HMM_PATH}\',\'wb\').write(gzip.decompress(r.content))"'
            )
            sys.exit(1)
    try:
        with pyhmmer.plan7.HMMFile(str(HMM_PATH)) as hf:
            hmms = list(hf)
        log.info(f'Loaded {len(hmms)} HMM profile(s) from {HMM_PATH.name}')
        return hmms
    except Exception as e:
        log.error(f'HMM load error: {e}')
        sys.exit(1)


# ── statistics ────────────────────────────────────────────────────────────────

def mann_whitney(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    from scipy.stats import mannwhitneyu
    stat, p = mannwhitneyu(a, b, alternative='two-sided')
    return float(stat), float(p)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sign(a[:, None] - b[None, :]).mean())


def auc_via_mwu(scores: np.ndarray, labels: np.ndarray) -> float:
    from scipy.stats import mannwhitneyu
    pos = scores[labels == 1]
    neg = scores[labels == 0]
    if len(pos) == 0 or len(neg) == 0:
        return float('nan')
    u, _ = mannwhitneyu(pos, neg, alternative='greater')
    return float(u) / (len(pos) * len(neg))


# ── per-species processing ────────────────────────────────────────────────────

def run_species(short: str, cfg: dict, hmms: list, log) -> dict:
    manifest_path = DATA_VAL / f'{short}_manifest.tsv'
    amr_path      = DATA_PROC / f'{short}_amr_hits.tsv'

    if not manifest_path.exists():
        log.error(f'Manifest not found: {manifest_path}'); return {}
    if not amr_path.exists():
        log.error(f'AMR hits not found: {amr_path}'); return {}

    manifest = pd.read_csv(manifest_path, sep='\t', dtype=str)
    amr_df   = pd.read_csv(amr_path, sep='\t', dtype=str)

    # Build set of carbapenem-resistant accessions
    # AMR hits = those with any carbapenem gene detected; all entries = resistant
    resistant_acc = set(amr_df['accession'].unique())
    log.info(f'{short.upper()}: {len(manifest)} genomes | {len(resistant_acc)} resistant')

    # Compile IS keyword patterns from species config
    is_kw_patterns = [re.compile(p, re.IGNORECASE)
                      for p in cfg.get('is_gff_keywords', [])]
    focal_re = re.compile(cfg['focal_is_pattern'], re.IGNORECASE)
    # Pattern to detect focal IS family by name (for counting)
    focal_count_re = re.compile(cfg['focal_is_pattern'], re.IGNORECASE)

    alphabet = pyhmmer.easel.Alphabet.amino()
    all_features: list[dict] = []
    query_seqs:   list       = []
    seq_to_idx:   dict[str, int] = {}

    log.info(f'Scanning {len(manifest)} GFF3 files …')
    for i, row in manifest.iterrows():
        acc      = row['accession']
        loc      = Path(row.get('local_dir', row.get('gff_path', '')))
        is_res   = acc in resistant_acc

        # Resolve GFF/FNA paths — local_dir may be directory or direct file
        if loc.is_dir():
            gff_files = sorted(loc.glob('*_genomic.gff.gz'))
            fna_files = sorted(loc.glob('*_genomic.fna.gz'))
        else:
            # Fallback: look in data/raw/{short}/{acc}/
            d = Path(__file__).resolve().parents[1] / 'data' / 'raw' / short / acc
            gff_files = sorted(d.glob('*_genomic.gff.gz'))
            fna_files = sorted(d.glob('*_genomic.fna.gz'))

        if not gff_files:
            log.debug(f'{acc}: no GFF, skip')
            continue

        features = scan_is_features(gff_files[0], is_kw_patterns, focal_re, log)

        contigs: dict[str, str] = {}
        if fna_files and any(f['is_family'] == 'IS_unknown' for f in features):
            contigs = load_contigs(fna_files[0], log)

        for feat in features:
            feat['accession']        = acc
            feat['is_resistant']     = is_res
            feat['hmm_evalue']       = None
            feat['hmm_bitscore']     = None
            feat['hmm_identity']     = None
            feat['corrected_family'] = feat['is_family']
            feat['was_reclassified'] = False
            idx = len(all_features)
            all_features.append(feat)

            if feat['is_family'] != 'IS_unknown':
                continue
            if not contigs:
                continue

            prot = extract_protein(contigs, feat['contig'],
                                   feat['start'], feat['stop'], feat['strand'])
            if prot is None:
                continue

            seq_name = f'{acc}|{feat["contig"]}|{feat["start"]}|{feat["stop"]}'
            try:
                ts = pyhmmer.easel.TextSequence(name=seq_name.encode(), sequence=prot)
                query_seqs.append(ts.digitize(alphabet))
                seq_to_idx[seq_name] = idx
            except Exception as e:
                log.debug(f'{seq_name}: digitize error: {e}')

        if (i + 1) % 50 == 0:
            log.info(f'  {i + 1}/{len(manifest)} …')

    n_unknown = sum(1 for f in all_features if f['is_family'] == 'IS_unknown')
    log.info(f'{short.upper()}: IS features {len(all_features):,} | '
             f'IS_unknown {n_unknown:,} | queries {len(query_seqs):,}')

    # ── HMMER search ──────────────────────────────────────────────────────────
    n_reclassified = 0
    label = f'IS6_PFAM_{PFAM_IS6_ID}' if short == 'eco' else f'ISAba_PFAM_{PFAM_IS6_ID}'

    if query_seqs:
        log.info(f'Running hmmsearch ({PFAM_IS6_ID}) …')
        target_block = pyhmmer.easel.DigitalSequenceBlock(alphabet, query_seqs)
        for hits in pyhmmer.hmmer.hmmsearch(hmms, target_block, cpus=0):
            for hit in hits:
                seq_name = hit.name
                if isinstance(seq_name, bytes):
                    seq_name = seq_name.decode()
                idx = seq_to_idx.get(seq_name)
                if idx is None:
                    continue
                evalue   = float(hit.evalue)
                bitscore = float(hit.score)
                identity = 0.0
                for dom in hit.domains:
                    if dom.included:
                        d = _domain_identity(dom)
                        if d > identity:
                            identity = d
                all_features[idx]['hmm_evalue']   = evalue
                all_features[idx]['hmm_bitscore'] = bitscore
                all_features[idx]['hmm_identity'] = round(identity, 4)
                if evalue <= HMMER_EVALUE and identity >= HMMER_MIN_IDENTITY:
                    all_features[idx]['corrected_family'] = label
                    all_features[idx]['was_reclassified'] = True
                    n_reclassified += 1
        log.info(f'Reclassified IS_unknown → {label}: {n_reclassified:,}')
    else:
        log.warning(f'{short.upper()}: No IS_unknown proteins to search')

    # ── Save per-feature results ───────────────────────────────────────────────
    feat_df = pd.DataFrame(all_features)
    feat_df.to_csv(DATA_PROC / f'{short}_is_hmmer_results.tsv', sep='\t', index=False)

    # ── Compute corrected focal IS counts per genome ───────────────────────────
    burden_rows = []
    for _, man_row in manifest.iterrows():
        acc   = man_row['accession']
        is_res = acc in resistant_acc
        gf    = feat_df[feat_df['accession'] == acc]

        n_focal_orig = int(gf['is_family'].apply(
            lambda x: bool(focal_count_re.search(str(x)))).sum())
        n_new        = int(gf['was_reclassified'].sum())
        n_focal_corr = n_focal_orig + n_new

        burden_rows.append({
            'accession':          acc,
            'is_resistant':       is_res,
            'n_is_total':         len(gf),
            'n_focal_original':   n_focal_orig,
            'n_focal_pfam_new':   n_new,
            'n_focal_corrected':  n_focal_corr,
        })

    burden_df = pd.DataFrame(burden_rows)
    burden_df.to_csv(DATA_PROC / f'{short}_is_burden_corrected.tsv',
                     sep='\t', index=False)

    # ── Statistics ────────────────────────────────────────────────────────────
    res_df  = burden_df[burden_df['is_resistant'] == True]
    susc_df = burden_df[burden_df['is_resistant'] == False]
    y       = burden_df['is_resistant'].values.astype(int)

    def _stats(col: str) -> dict:
        a = res_df[col].values.astype(float)
        b = susc_df[col].values.astype(float)
        _, p = mann_whitney(a, b)
        d    = cliffs_delta(a, b)
        auc  = auc_via_mwu(burden_df[col].values.astype(float), y)
        return {
            'n_resistant':        len(a),
            'n_susceptible':      len(b),
            'resistant_median':   float(np.median(a)),
            'resistant_iqr':      [float(np.percentile(a, 25)),
                                   float(np.percentile(a, 75))],
            'susceptible_median': float(np.median(b)),
            'susceptible_iqr':    [float(np.percentile(b, 25)),
                                   float(np.percentile(b, 75))],
            'cliffs_delta':       round(d, 4),
            'mannwhitney_p':      float(p),
            'auc':                round(auc, 4),
        }

    orig = _stats('n_focal_original')
    corr = _stats('n_focal_corrected')

    n_reclass_res  = int(res_df['n_focal_pfam_new'].sum())
    n_reclass_susc = int(susc_df['n_focal_pfam_new'].sum())
    bias_confirmed = (len(susc_df) > 0 and n_reclass_susc > n_reclass_res * 2)

    auc_change   = round(corr['auc'] - orig['auc'], 4)
    delta_change = round(corr['cliffs_delta'] - orig['cliffs_delta'], 4)

    if bias_confirmed:
        verdict = (f'BIAS CONFIRMED — susceptible genomes gained '
                   f'{n_reclass_susc} vs resistant {n_reclass_res} '
                   f'reclassifications. AUC {orig["auc"]:.3f} → {corr["auc"]:.3f} '
                   f'({auc_change:+.4f}). Corrected AUC is publishable figure.')
    else:
        verdict = (f'No disproportionate annotation bias detected — '
                   f'reclassifications: resistant {n_reclass_res} / '
                   f'susceptible {n_reclass_susc}. '
                   f'AUC {orig["auc"]:.3f} → {corr["auc"]:.3f} ({auc_change:+.4f}).')

    stats_out = {
        'species':    short,
        'focal_is':   cfg['focal_is_family'],
        'pfam_id':    PFAM_IS6_ID,
        'original':   orig,
        'corrected':  {**corr,
                       'n_reclassified_total':       n_reclassified,
                       'n_reclassified_resistant':   n_reclass_res,
                       'n_reclassified_susceptible': n_reclass_susc},
        'bias_assessment': {
            'annotation_bias_confirmed':   bias_confirmed,
            'reclassification_toward':     'susceptible' if n_reclass_susc > n_reclass_res else 'resistant',
            'auc_change':                  auc_change,
            'cliffs_delta_change':         delta_change,
            'verdict':                     verdict,
        },
    }
    out_json = DATA_PROC / f'{short}_is_burden_corrected_stats.json'
    out_json.write_text(json.dumps(stats_out, indent=2, default=str))

    # ── Summary log ───────────────────────────────────────────────────────────
    log.info('═' * 65)
    log.info(f'{short.upper()} IS HMMER VERIFICATION COMPLETE')
    log.info(f'  Focal IS family   : {cfg["focal_is_family"]} (PF01527)')
    log.info(f'  IS features total : {len(all_features):,}')
    log.info(f'  IS_unknown queries: {len(query_seqs):,}')
    log.info(f'  Reclassified      : {n_reclassified:,}  '
             f'(res {n_reclass_res} / susc {n_reclass_susc})')
    log.info('─' * 65)
    log.info(f'  ORIGINAL  focal IS resistant median : {orig["resistant_median"]:.1f}')
    log.info(f'  ORIGINAL  focal IS susceptible median: {orig["susceptible_median"]:.1f}')
    log.info(f"  ORIGINAL  Cliff's delta : {orig['cliffs_delta']:.3f}")
    log.info(f'  ORIGINAL  AUC           : {orig["auc"]:.3f}')
    log.info('─' * 65)
    log.info(f'  CORRECTED focal IS resistant median : {corr["resistant_median"]:.1f}')
    log.info(f'  CORRECTED focal IS susceptible median: {corr["susceptible_median"]:.1f}')
    log.info(f"  CORRECTED Cliff's delta : {corr['cliffs_delta']:.3f}")
    log.info(f'  CORRECTED AUC           : {corr["auc"]:.3f}')
    log.info('─' * 65)
    log.info(f'  BIAS ASSESSMENT: {verdict}')
    log.info('═' * 65)

    return stats_out


# ── entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description='PFAM IS HMM verification for ECO and ABA IS elements')
    ap.add_argument('--species', choices=['eco', 'aba', 'both'], default='both')
    args = ap.parse_args()

    LOGS.mkdir(parents=True, exist_ok=True)
    DB_DIR.mkdir(parents=True, exist_ok=True)

    species_list = (['eco', 'aba'] if args.species == 'both'
                    else [args.species])
    short_to_name = {'eco': 'Escherichia coli', 'aba': 'Acinetobacter baumannii'}

    # Set up combined log
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s  %(levelname)-8s  %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        handlers=[
            logging.FileHandler(LOGS / 'is_hmmer.log', encoding='utf-8'),
            logging.StreamHandler(sys.stdout),
        ],
    )
    log = logging.getLogger('is_hmmer')

    hmms = ensure_pfam_hmm(log)

    all_results = {}
    for short in species_list:
        name = short_to_name[short]
        cfg  = SPECIES_CONFIG[name]
        log.info(f'\n{"=" * 65}')
        log.info(f'Processing {name} ({short})')
        log.info(f'{"=" * 65}')
        result = run_species(short, cfg, hmms, log)
        if result:
            all_results[short] = result

    # ── Write combined summary ────────────────────────────────────────────────
    if len(all_results) > 1:
        summary_path = DATA_PROC / 'is_hmmer_combined_summary.json'
        summary_path.write_text(json.dumps(all_results, indent=2, default=str))
        log.info(f'\nCombined summary → {summary_path}')

    log.info('\nAll done.')


if __name__ == '__main__':
    main()
