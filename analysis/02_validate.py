"""
Step 2 — Species-aware quality gate
────────────────────────────────────
For each downloaded genome: checks size, N50, CDS count, GFF integrity.
Writes data/validated/{species_short}_manifest.tsv

Usage:
  python analysis/02_validate.py [--species eco|aba|both]
"""
import argparse
import gzip
import logging
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_RAW, DATA_VAL, LOGS, SPECIES_CONFIG

LOGS.mkdir(parents=True, exist_ok=True)
DATA_VAL.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'validate.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('validate')


class ValidationError(Exception):
    pass


def _find_file(d: Path, suffix: str) -> Path | None:
    hits = list(d.glob(f'*{suffix}'))
    return hits[0] if hits else None


def _check_gz(path: Path, n: int = 65536):
    try:
        with gzip.open(path, 'rb') as fh:
            data = fh.read(n)
        if not data:
            raise ValidationError(f'Empty gz: {path}')
    except Exception as e:
        raise ValidationError(str(e))


def _parse_assembly_stats(path: Path) -> dict:
    stats = {}
    try:
        text = path.read_text(encoding='utf-8', errors='replace')
    except Exception:
        return stats
    patterns = {
        'total_length': r'total-length\s+all\s+(\d+)',
        'scaffold_n50': r'scaffold-N50\s+(\d+)',
        'contig_n50':   r'contig-N50\s+(\d+)',
        'cds_count':    r'Gene\s+Coding\s+all\s+(\d+)',
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, re.IGNORECASE)
        if m:
            stats[key] = int(m.group(1))
    return stats


def _count_cds_in_gff(gff_path: Path) -> int:
    count = 0
    try:
        with gzip.open(gff_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3 and parts[2] == 'CDS':
                    count += 1
    except Exception:
        pass
    return count


def validate_genome(accession: str, genome_dir: Path, cfg: dict) -> dict:
    result = {
        'accession': accession,
        'status': 'PASS',
        'gff_path': '',
        'fna_path': '',
        'genome_size': 0,
        'n50': 0,
        'cds_count': 0,
        'reason': '',
    }

    gff = _find_file(genome_dir, '_genomic.gff.gz')
    fna = _find_file(genome_dir, '_genomic.fna.gz')
    stats_file = _find_file(genome_dir, '_assembly_stats.txt')

    if not gff:
        result.update(status='FAIL', reason='missing_gff')
        return result
    if not fna:
        result.update(status='FAIL', reason='missing_fna')
        return result

    # GFF readable
    try:
        _check_gz(gff)
    except ValidationError as e:
        result.update(status='FAIL', reason=f'gff_unreadable:{e}')
        return result

    # FASTA readable
    try:
        _check_gz(fna)
    except ValidationError as e:
        result.update(status='FAIL', reason=f'fna_unreadable:{e}')
        return result

    # Assembly stats
    if stats_file:
        stats = _parse_assembly_stats(stats_file)
    else:
        stats = {}

    genome_size = stats.get('total_length', 0)
    n50 = stats.get('scaffold_n50') or stats.get('contig_n50', 0)
    cds_count = stats.get('cds_count', 0)

    # If stats file missing CDS, count from GFF
    if cds_count == 0:
        cds_count = _count_cds_in_gff(gff)

    result['gff_path'] = str(gff)
    result['fna_path'] = str(fna)
    result['genome_size'] = genome_size
    result['n50'] = n50
    result['cds_count'] = cds_count

    # Size check
    if genome_size > 0:
        if genome_size < cfg['genome_size_min_bp']:
            result.update(status='FAIL',
                          reason=f'size_too_small:{genome_size}')
            return result
        if genome_size > cfg['genome_size_max_bp']:
            result.update(status='FAIL',
                          reason=f'size_too_large:{genome_size}')
            return result

    # N50 check
    if n50 > 0 and n50 < cfg['min_n50']:
        result.update(status='FAIL', reason=f'n50_too_low:{n50}')
        return result

    # CDS check
    if cds_count > 0 and cds_count < cfg['min_cds']:
        result.update(status='FAIL', reason=f'cds_too_few:{cds_count}')
        return result

    return result


def run_species(species_name: str, cfg: dict):
    short = cfg['short']
    log.info(f'Validating {species_name} ({short})…')
    raw_dir = DATA_RAW / short
    if not raw_dir.exists():
        log.warning(f'{raw_dir} not found — run 01_download.py first')
        return

    status_tsv = raw_dir / 'download_status.tsv'
    if not status_tsv.exists():
        log.warning(f'No download_status.tsv in {raw_dir}')
        return

    dl_df = pd.read_csv(status_tsv, sep='\t', dtype=str)
    ok_rows = dl_df[dl_df['status'].isin(['OK', 'PARTIAL'])]
    log.info(f'  {len(ok_rows)} assemblies to validate')

    records = []
    for _, row in ok_rows.iterrows():
        acc = str(row['accession'])
        genome_dir = raw_dir / acc
        rec = validate_genome(acc, genome_dir, cfg)
        rec['species'] = species_name
        records.append(rec)

    if not records:
        log.warning(f'No records for {species_name}')
        return

    df = pd.DataFrame(records)
    passed = df[df['status'] == 'PASS']
    out = DATA_VAL / f'{short}_manifest.tsv'
    passed.to_csv(out, sep='\t', index=False)

    n_pass = len(passed)
    n_fail = len(df) - n_pass
    log.info(f'{species_name}: PASS={n_pass}  FAIL={n_fail} → {out}')

    # Failure summary
    if n_fail:
        fail_summary = df[df['status'] == 'FAIL']['reason'].value_counts()
        log.info(f'Failure reasons:\n{fail_summary.to_string()}')


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
