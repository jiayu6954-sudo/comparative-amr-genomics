"""
Step 4 — IS element flanking analysis (species-aware)
───────────────────────────────────────────────────────
For each AMR locus: find IS elements within ±10 kb.
Classifies: COMPOSITE_TRANSPOSON / SINGLE_IS_UP / SINGLE_IS_DOWN / NO_IS

Output: data/processed/{short}_is_context.tsv

Usage:
  python analysis/04_is_context.py [--species eco|aba|both]
"""
import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, DATA_VAL, IS_FLANK_WINDOW_BP, LOGS, SPECIES_CONFIG

LOGS.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'is_context.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('is_context')


# ── IS family classification ──────────────────────────────────────────────────

_IS_FAMILY_PATTERNS = [
    (re.compile(r'ISAba\d*|ISaba', re.I),  'ISAba'),
    (re.compile(r'IS26\b|IS257\b|IS1353\b|IS1006\b', re.I), 'IS6'),
    (re.compile(r'IS6\b', re.I),           'IS6'),
    (re.compile(r'Tn4401', re.I),          'Tn4401'),
    (re.compile(r'Tn3\b|Tn\d+', re.I),    'Tn3'),
    (re.compile(r'IS5\b|IS5\d', re.I),     'IS5'),
    (re.compile(r'IS1\b',  re.I),          'IS1'),
    (re.compile(r'IS3\b',  re.I),          'IS3'),
    (re.compile(r'IS91\b', re.I),          'IS91'),
    (re.compile(r'IS481\b',re.I),          'IS481'),
    (re.compile(r'IS110\b',re.I),          'IS110'),
    (re.compile(r'IS1182\b',re.I),         'IS1182'),
    (re.compile(r'IS\d+',  re.I),          'IS_other'),
]


def _classify_is_family(annotation: str) -> str:
    for pat, family in _IS_FAMILY_PATTERNS:
        if pat.search(annotation):
            return family
    return 'IS_unknown'


def _load_is_features(gff_path: Path, is_patterns: list) -> list:
    """Return list of IS feature dicts from GFF3."""
    compiled = [re.compile(p, re.IGNORECASE) for p in is_patterns]
    features = []
    try:
        with gzip.open(gff_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue
                attrs_str = parts[8]
                product = ''
                gene    = ''
                for seg in attrs_str.split(';'):
                    if seg.lower().startswith('product='):
                        product = seg.split('=', 1)[1]
                    elif seg.lower().startswith('gene='):
                        gene = seg.split('=', 1)[1]
                search_text = f'{product} {gene}'
                if not any(pat.search(search_text) for pat in compiled):
                    continue
                family = _classify_is_family(product or gene)
                features.append({
                    'contig':  parts[0],
                    'start':   int(parts[3]),
                    'end':     int(parts[4]),
                    'strand':  parts[6],
                    'product': product,
                    'family':  family,
                })
    except Exception as e:
        log.debug(f'IS parse error {gff_path}: {e}')
    return features


def _classify_flanking(up_count: int, down_count: int) -> str:
    if up_count > 0 and down_count > 0:
        return 'COMPOSITE_TRANSPOSON'
    elif up_count > 0:
        return 'SINGLE_IS_UPSTREAM'
    elif down_count > 0:
        return 'SINGLE_IS_DOWNSTREAM'
    return 'NO_IS'


def analyze_is_context(amr_hit: dict, is_features: list) -> list:
    """Generate IS–AMR context pairs for one AMR locus."""
    contig  = amr_hit['contig']
    a_start = int(amr_hit['start'])
    a_end   = int(amr_hit['end'])
    pairs   = []

    nearby = [f for f in is_features
              if f['contig'] == contig
              and f['end']   >= a_start - IS_FLANK_WINDOW_BP
              and f['start'] <= a_end   + IS_FLANK_WINDOW_BP]

    for is_f in nearby:
        # gap distance
        if is_f['end'] < a_start:
            gap = a_start - is_f['end']
            position = 'UPSTREAM'
        elif is_f['start'] > a_end:
            gap = is_f['start'] - a_end
            position = 'DOWNSTREAM'
        else:
            gap = 0
            position = 'OVERLAPPING'

        strand_concordant = (is_f['strand'] == amr_hit['strand'])

        pairs.append({
            'accession':     amr_hit['accession'],
            'species':       amr_hit['species'],
            'amr_gene':      amr_hit['gene_name'],
            'amr_contig':    contig,
            'amr_start':     a_start,
            'amr_end':       a_end,
            'is_family':     is_f['family'],
            'is_start':      is_f['start'],
            'is_end':        is_f['end'],
            'gap_bp':        gap,
            'position':      position,
            'same_strand':   strand_concordant,
        })

    return pairs


def run_species(species_name: str, cfg: dict):
    short = cfg['short']
    amr_path  = DATA_PROC / f'{short}_amr_hits.tsv'
    mani_path = DATA_VAL   / f'{short}_manifest.tsv'

    if not amr_path.exists():
        log.warning(f'No AMR hits for {species_name}')
        return
    if not mani_path.exists():
        log.warning(f'No manifest for {species_name}')
        return

    amr_df  = pd.read_csv(amr_path,  sep='\t', dtype=str)
    mani_df = pd.read_csv(mani_path, sep='\t', dtype=str)

    if amr_df.empty:
        log.info(f'{species_name}: no AMR hits — skipping IS context')
        return

    gff_map = dict(zip(mani_df['accession'], mani_df['gff_path']))
    is_patterns = cfg['is_gff_keywords']

    all_pairs = []
    resistant_accs = amr_df['accession'].unique()
    log.info(f'{species_name}: IS context for {len(resistant_accs)} resistant genomes…')

    for acc in resistant_accs:
        gff_path = Path(gff_map.get(acc, ''))
        if not gff_path.exists():
            continue
        is_features = _load_is_features(gff_path, is_patterns)
        hits = amr_df[amr_df['accession'] == acc].to_dict('records')
        for hit in hits:
            pairs = analyze_is_context(hit, is_features)
            all_pairs.extend(pairs)

    df = pd.DataFrame(all_pairs) if all_pairs else pd.DataFrame()

    # Locus-level flanking classification
    if not df.empty:
        locus_records = []
        for (acc, gene, contig, start, end), grp in df.groupby(
                ['accession', 'amr_gene', 'amr_contig', 'amr_start', 'amr_end']):
            up   = (grp['position'] == 'UPSTREAM').sum()
            down = (grp['position'] == 'DOWNSTREAM').sum()
            locus_records.append({
                'accession': acc, 'gene': gene, 'contig': contig,
                'start': start, 'end': end,
                'n_upstream_is': up, 'n_downstream_is': down,
                'flanking_class': _classify_flanking(up, down),
                'species': species_name,
            })
        locus_df = pd.DataFrame(locus_records)
        locus_df.to_csv(DATA_PROC / f'{short}_locus_classification.tsv',
                        sep='\t', index=False)

    out = DATA_PROC / f'{short}_is_context.tsv'
    df.to_csv(out, sep='\t', index=False)

    n_pairs = len(df)
    composite = 0
    if not df.empty and 'flanking_class' in locals().get('locus_df', pd.DataFrame()).columns:
        composite = (locus_df['flanking_class'] == 'COMPOSITE_TRANSPOSON').sum()
    log.info(f'{species_name}: {n_pairs} IS–AMR pairs → {out}')


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
