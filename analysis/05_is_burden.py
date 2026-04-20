"""
Step 5 — Focal IS burden across all genomes + AUC calculation
──────────────────────────────────────────────────────────────
Counts focal IS copies (IS6 for E. coli; ISAba for A. baumannii)
per genome, then calculates Mann-Whitney AUC for predicting resistance.

Also loads K. pneumoniae results from existing project for comparison.

Output: data/processed/cross_species_burden.tsv
        data/processed/cross_species_stats.json

Usage:
  python analysis/05_is_burden.py
"""
import gzip
import json
import logging
import math
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_PROC, DATA_VAL, KPN_AMR_HITS,
    KPN_IS_CONTEXT, KPN_MANIFEST, KPN_PROJECT, LOGS, SPECIES_CONFIG,
)

LOGS.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'is_burden.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('is_burden')


# ── statistics helpers ────────────────────────────────────────────────────────

def mann_whitney_auc(group1: list, group2: list) -> tuple[float, float]:
    """
    One-sided Mann-Whitney U: P(group1 > group2).
    Returns (AUC, p_value) using normal approximation.
    """
    n1, n2 = len(group1), len(group2)
    if n1 == 0 or n2 == 0:
        return 0.0, 1.0

    u = sum(1 for x in group1 for y in group2 if x > y) + \
        0.5 * sum(1 for x in group1 for y in group2 if x == y)
    auc = u / (n1 * n2)

    # Normal approximation
    mean_u = n1 * n2 / 2
    std_u  = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    if std_u == 0:
        return auc, 1.0
    z = (u - mean_u) / std_u

    # Two-sided p via standard normal CDF approximation
    def norm_cdf(z):
        t = 1 / (1 + 0.2316419 * abs(z))
        poly = t * (0.319381530 + t * (-0.356563782 + t * (
               1.781477937 + t * (-1.821255978 + t * 1.330274429))))
        p_one = 1 - (1 / math.sqrt(2 * math.pi)) * math.exp(-z**2 / 2) * poly
        return p_one if z >= 0 else 1 - p_one

    p = 2 * min(norm_cdf(z), 1 - norm_cdf(z))
    return auc, p


def cliffs_delta(group1: list, group2: list) -> float:
    """Cliff's delta: proportion of (g1>g2) minus (g1<g2) pairs."""
    n1, n2 = len(group1), len(group2)
    if n1 == 0 or n2 == 0:
        return 0.0
    gt = sum(1 for x in group1 for y in group2 if x > y)
    lt = sum(1 for x in group1 for y in group2 if x < y)
    return (gt - lt) / (n1 * n2)


def median(vals: list) -> float:
    if not vals:
        return 0.0
    s = sorted(vals)
    n = len(s)
    return (s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2)


# ── IS burden counting ────────────────────────────────────────────────────────

def count_focal_is(gff_path: Path, focal_pattern: str) -> int:
    """Count IS elements matching focal_pattern in a GFF3 file."""
    pat = re.compile(focal_pattern, re.IGNORECASE)
    count = 0
    try:
        with gzip.open(gff_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 9 or parts[2] != 'CDS':
                    continue
                attrs = parts[8]
                product = ''
                gene    = ''
                for seg in attrs.split(';'):
                    if seg.lower().startswith('product='):
                        product = seg.split('=', 1)[1]
                    elif seg.lower().startswith('gene='):
                        gene = seg.split('=', 1)[1]
                if pat.search(product) or pat.search(gene):
                    count += 1
    except Exception:
        pass
    return count


def compute_burden_for_species(species_name: str, cfg: dict) -> pd.DataFrame:
    short       = cfg['short']
    focal_pat   = cfg['focal_is_pattern']
    mani_path   = DATA_VAL / f'{short}_manifest.tsv'
    amr_path    = DATA_PROC / f'{short}_amr_hits.tsv'

    if not mani_path.exists():
        log.warning(f'No manifest for {species_name}')
        return pd.DataFrame()

    manifest = pd.read_csv(mani_path, sep='\t', dtype=str)
    resistant_accs = set()
    if amr_path.exists():
        amr_df = pd.read_csv(amr_path, sep='\t', dtype=str)
        if not amr_df.empty:
            resistant_accs = set(amr_df['accession'].unique())

    records = []
    for _, row in manifest.iterrows():
        acc      = row['accession']
        gff_path = Path(row['gff_path'])
        if not gff_path.exists():
            continue
        focal_count = count_focal_is(gff_path, focal_pat)
        records.append({
            'accession':    acc,
            'species':      species_name,
            'focal_is_family': cfg['focal_is_family'],
            'focal_is_count': focal_count,
            'resistant':    acc in resistant_accs,
        })

    log.info(f'{species_name}: {len(records)} genomes, '
             f'{len(resistant_accs)} resistant')
    return pd.DataFrame(records)


def compute_kpn_burden() -> pd.DataFrame:
    """Load K. pneumoniae burden from existing project."""
    # Use the processed is_burden_corrected_stats.json if available,
    # otherwise rebuild from manifests
    mani = KPN_MANIFEST
    amr  = KPN_AMR_HITS

    if not mani.exists():
        log.warning('KPN manifest not found')
        return pd.DataFrame()

    manifest = pd.read_csv(mani, sep='\t', dtype=str)
    resistant_accs = set()
    if amr.exists():
        amr_df = pd.read_csv(amr, sep='\t', dtype=str)
        if not amr_df.empty:
            resistant_accs = set(amr_df['accession'].unique())

    # IS6 pattern for K. pneumoniae (same as focal_is_pattern would be)
    focal_pat = r'IS26|IS257|IS1353|IS1006|ISKpn\d*|IS6\b'
    records = []
    for _, row in manifest.iterrows():
        acc = row['accession']
        # KPN manifest uses 'local_dir' column; find the GFF inside it
        local_dir = Path(row.get('local_dir', row.get('gff_path', '')))
        gff_path  = Path('')
        if local_dir.is_dir():
            candidates = sorted(local_dir.glob('*_genomic.gff.gz'))
            if candidates:
                gff_path = candidates[0]
        elif str(local_dir).endswith('.gz') and local_dir.exists():
            gff_path = local_dir
        if not gff_path.exists():
            continue
        focal_count = count_focal_is(gff_path, focal_pat)
        records.append({
            'accession':      acc,
            'species':        'Klebsiella pneumoniae',
            'focal_is_family':'IS6',
            'focal_is_count': focal_count,
            'resistant':      acc in resistant_accs,
        })

    log.info(f'Klebsiella pneumoniae: {len(records)} genomes loaded from existing project')
    return pd.DataFrame(records)


# ── statistics per species ────────────────────────────────────────────────────

def compute_stats(df: pd.DataFrame, species: str) -> dict:
    res = df[df['resistant'] == True]['focal_is_count'].tolist()
    sus = df[df['resistant'] == False]['focal_is_count'].tolist()
    if not res or not sus:
        return {'species': species, 'error': 'insufficient_data'}

    auc, pval = mann_whitney_auc(res, sus)
    delta      = cliffs_delta(res, sus)

    return {
        'species':          species,
        'n_total':          len(df),
        'n_resistant':      len(res),
        'n_susceptible':    len(sus),
        'resistance_rate':  round(len(res) / len(df) * 100, 1),
        'focal_is_family':  df['focal_is_family'].iloc[0],
        'median_resistant': median(res),
        'median_susceptible': median(sus),
        'auc':              round(auc, 4),
        'cliffs_delta':     round(delta, 4),
        'mwu_p':            f'{pval:.2e}',
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    all_dfs = []

    # K. pneumoniae from existing project
    kpn_df = compute_kpn_burden()
    if not kpn_df.empty:
        all_dfs.append(kpn_df)

    # E. coli and A. baumannii from new downloads
    for name, cfg in SPECIES_CONFIG.items():
        sp_df = compute_burden_for_species(name, cfg)
        if not sp_df.empty:
            all_dfs.append(sp_df)

    if not all_dfs:
        log.error('No burden data found for any species')
        return

    combined = pd.concat(all_dfs, ignore_index=True)
    combined.to_csv(DATA_PROC / 'cross_species_burden.tsv', sep='\t', index=False)

    # Per-species statistics
    stats_list = []
    for species in combined['species'].unique():
        sp_df = combined[combined['species'] == species]
        stats = compute_stats(sp_df, species)
        stats_list.append(stats)
        log.info(f'\n{species}:')
        for k, v in stats.items():
            if k != 'species':
                log.info(f'  {k}: {v}')

    with open(DATA_PROC / 'cross_species_stats.json', 'w') as fh:
        json.dump(stats_list, fh, indent=2)

    log.info(f'\nResults saved to {DATA_PROC}/cross_species_*.tsv/json')


if __name__ == '__main__':
    main()
