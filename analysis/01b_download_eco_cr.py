"""
Step 1b — Targeted download of carbapenem-resistant E. coli genomes
───────────────────────────────────────────────────────────────────
Supplements the general ECO cohort (01_download.py) with a targeted
search that enriches for NDM/KPC-carrying strains from China.

Strategy:
  - Multi-query Entrez search combining organism + resistance keywords
    (NDM, KPC, carbapenem, OXA-48) as [All Fields] to capture strains
    annotated in assembly name, strain name, or notes fields
  - Downloads into data/raw/eco_cr/ (separate from general eco/)
  - Deduplicates against existing eco/ downloads

Usage:
  python analysis/01b_download_eco_cr.py [--limit N] [--dry-run]
"""
import argparse
import logging
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

sys.path.insert(0, str(Path(__file__).parent))
from config import (
    DATA_RAW, DOWNLOAD_FILES, DOWNLOAD_TIMEOUT_S, LOGS, MAX_RETRIES,
    NCBI_API_KEY, NCBI_EMAIL, RATE_LIMIT_DELAY_S, RETRY_DELAY_S,
)

LOGS.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'download_eco_cr.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('eco_cr')

SESSION = requests.Session()
SESSION.headers['User-Agent'] = f'ComparativeAMR/1.0 ({NCBI_EMAIL})'
ENTREZ = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'

ECO_CR_QUERIES = [
    ('"Escherichia coli"[Organism] AND "China"[Country] AND "NDM"[All Fields] '
     'AND latest[filter] AND ("Complete Genome"[Assembly Level] OR '
     '"Chromosome"[Assembly Level] OR "Scaffold"[Assembly Level])'),
    ('"Escherichia coli"[Organism] AND "China"[Country] AND "KPC"[All Fields] '
     'AND latest[filter] AND ("Complete Genome"[Assembly Level] OR '
     '"Chromosome"[Assembly Level] OR "Scaffold"[Assembly Level])'),
    ('"Escherichia coli"[Organism] AND "China"[Country] AND "carbapenem"[All Fields] '
     'AND latest[filter] AND ("Complete Genome"[Assembly Level] OR '
     '"Chromosome"[Assembly Level] OR "Scaffold"[Assembly Level])'),
    ('"Escherichia coli"[Organism] AND "China"[Country] AND "OXA-48"[All Fields] '
     'AND latest[filter] AND ("Complete Genome"[Assembly Level] OR '
     '"Chromosome"[Assembly Level] OR "Scaffold"[Assembly Level])'),
]


def _get(url, params=None, stream=False):
    if NCBI_API_KEY:
        params = params or {}
        params['api_key'] = NCBI_API_KEY
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = SESSION.get(url, params=params, stream=stream,
                            timeout=DOWNLOAD_TIMEOUT_S)
            r.raise_for_status()
            return r
        except requests.HTTPError as e:
            if getattr(e.response, 'status_code', 0) == 404:
                raise
            if attempt == MAX_RETRIES:
                raise
            time.sleep(RETRY_DELAY_S * attempt)
        except requests.RequestException:
            if attempt == MAX_RETRIES:
                raise
            time.sleep(RETRY_DELAY_S * attempt)


def _search(query, retmax=500):
    r    = _get(f'{ENTREZ}/esearch.fcgi',
                {'db': 'assembly', 'term': query, 'retmax': min(retmax, 5000),
                 'retmode': 'xml', 'email': NCBI_EMAIL})
    root = ET.fromstring(r.content)
    count = int(root.findtext('Count') or '0')
    ids   = [el.text for el in root.findall('.//IdList/Id')]
    log.info(f'  {count} total hits → {len(ids)} UIDs retrieved')
    time.sleep(RATE_LIMIT_DELAY_S)
    return ids[:retmax]


def _summary(uids):
    records = []
    for i in range(0, len(uids), 200):
        batch  = uids[i:i + 200]
        r      = _get(f'{ENTREZ}/esummary.fcgi',
                      {'db': 'assembly', 'id': ','.join(batch),
                       'retmode': 'xml', 'email': NCBI_EMAIL})
        root   = ET.fromstring(r.content)
        time.sleep(RATE_LIMIT_DELAY_S)
        for doc in root.findall('.//DocumentSummary'):
            def v(tag):
                el = doc.find(tag)
                return el.text.strip() if el is not None and el.text else ''
            ftp = v('FtpPath_GenBank') or v('FtpPath_RefSeq')
            records.append({'assembly_accession': v('AssemblyAccession'),
                             'assembly_level':     v('AssemblyStatus'),
                             'ftp_path':           ftp})
    return records


def _dl_genome(acc, ftp_base, out_root, dry_run):
    if ftp_base.startswith('ftp://'):
        ftp_base = 'https://' + ftp_base[6:]
    ftp_base   = ftp_base.rstrip('/')
    fname_base = ftp_base.split('/')[-1]
    out_dir    = out_root / acc
    res        = {'accession': acc, 'status': 'OK', 'files': []}
    if dry_run:
        res['status'] = 'DRY_RUN'
        return res
    out_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    for suffix in DOWNLOAD_FILES:
        url  = f'{ftp_base}/{fname_base}{suffix}'
        dest = out_dir / f'{fname_base}{suffix}'
        if dest.exists() and dest.stat().st_size > 0:
            res['files'].append(dest.name)
            continue
        try:
            r = _get(url, stream=True)
            with open(dest, 'wb') as fh:
                for chunk in r.iter_content(1 << 20):
                    fh.write(chunk)
            res['files'].append(dest.name)
        except Exception as e:
            log.debug(f'Failed {url}: {e}')
            if dest.exists():
                dest.unlink()
            failed.append(suffix)
        time.sleep(RATE_LIMIT_DELAY_S)
    if failed:
        res['status'] = 'FAILED' if any('gff' in f or 'fna' in f
                                        for f in failed) else 'PARTIAL'
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--limit', type=int, default=300)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    out_root = DATA_RAW / 'eco_cr'
    out_root.mkdir(parents=True, exist_ok=True)

    # Gather UIDs from all targeted queries
    all_uids: set[str] = set()
    for i, q in enumerate(ECO_CR_QUERIES, 1):
        log.info(f'\nQuery {i}/{len(ECO_CR_QUERIES)}: {q[:70]}...')
        all_uids.update(_search(q, retmax=args.limit))
    log.info(f'\nUnique UIDs across all queries: {len(all_uids)}')

    # Deduplicate against existing downloads
    def _existing(subdir):
        p = DATA_RAW / subdir
        return {d.name for d in p.iterdir() if d.is_dir()} if p.exists() else set()

    already = _existing('eco') | _existing('eco_cr')
    records  = _summary(list(all_uids))
    log.info(f'Metadata retrieved: {len(records)}')

    ACCEPTED = {'Complete Genome', 'Chromosome', 'Scaffold'}
    records = [r for r in records
               if r['assembly_level'] in ACCEPTED
               and r['ftp_path'] and r['ftp_path'] != 'na'
               and r['assembly_accession'] not in already]
    log.info(f'New candidates (not yet downloaded): {len(records)}')

    if not records:
        log.info('Nothing new to download — all candidates already present.')
        return

    results = []
    for rec in tqdm(records, desc='eco_cr'):
        r = _dl_genome(rec['assembly_accession'], rec['ftp_path'],
                       out_root, args.dry_run)
        r['assembly_level'] = rec['assembly_level']
        results.append(r)
        time.sleep(RATE_LIMIT_DELAY_S)

    df = pd.DataFrame(results)
    df.to_csv(out_root / 'download_status_cr.tsv', sep='\t', index=False)
    ok = (df['status'] == 'OK').sum()
    log.info(f'Result: OK={ok}  PARTIAL={(df["status"]=="PARTIAL").sum()}'
             f'  FAILED={(df["status"]=="FAILED").sum()}')


if __name__ == '__main__':
    main()
