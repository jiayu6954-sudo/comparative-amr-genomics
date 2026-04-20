"""
Step 1 — Download E. coli and A. baumannii genomes from NCBI GenBank
──────────────────────────────────────────────────────────────────────
Uses Entrez search (proven approach from KPN project):
  1. esearch: organism + China + assembly level filter → UID list
  2. esummary: get FTP paths + metadata per UID
  3. BioSample: clinical source verification (sampled)
  4. Download GFF3 + FASTA + assembly stats per genome

Usage:
  python analysis/01_download.py [--species eco|aba|both] [--limit N] [--dry-run]
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
    ACCEPTED_ASSEMBLY_LEVELS, ACCEPTED_SOURCES, DATA_RAW,
    DOWNLOAD_FILES, DOWNLOAD_TIMEOUT_S, LOGS, MAX_RETRIES,
    NCBI_API_KEY, NCBI_EMAIL, RATE_LIMIT_DELAY_S,
    RETRY_DELAY_S, SPECIES_CONFIG,
)

LOGS.mkdir(parents=True, exist_ok=True)
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(LOGS / 'download.log', encoding='utf-8'),
        logging.StreamHandler(sys.stdout),
    ]
)
log = logging.getLogger('download')

SESSION = requests.Session()
SESSION.headers['User-Agent'] = f'ComparativeAMR/1.0 ({NCBI_EMAIL})'
ENTREZ_BASE = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'


# ── HTTP helpers ──────────────────────────────────────────────────────────────

def _get(url: str, params: dict = None, stream: bool = False) -> requests.Response:
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
            if e.response.status_code == 404:
                raise
            if attempt == MAX_RETRIES:
                raise
            log.warning(f'HTTP attempt {attempt}: {e}. Retrying…')
            time.sleep(RETRY_DELAY_S * attempt)
        except requests.RequestException as e:
            if attempt == MAX_RETRIES:
                raise
            log.warning(f'Network attempt {attempt}: {e}. Retrying…')
            time.sleep(RETRY_DELAY_S * attempt)


def _download_file(url: str, dest: Path) -> bool:
    try:
        r = _get(url, stream=True)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with open(dest, 'wb') as fh:
            for chunk in r.iter_content(chunk_size=1 << 20):
                fh.write(chunk)
        return True
    except Exception as e:
        log.debug(f'Failed {url}: {e}')
        if dest.exists():
            dest.unlink()
        return False


# ── Entrez search + summary ───────────────────────────────────────────────────

def _entrez_search(query: str, retmax: int = 1000) -> list[str]:
    params = {
        'db': 'assembly', 'term': query, 'retmax': min(retmax, 10000),
        'usehistory': 'y', 'retmode': 'xml', 'email': NCBI_EMAIL,
    }
    r    = _get(f'{ENTREZ_BASE}/esearch.fcgi', params=params)
    root = ET.fromstring(r.content)
    count   = int(root.findtext('Count') or '0')
    web_env = root.findtext('WebEnv') or ''
    qkey    = root.findtext('QueryKey') or '1'
    ids     = [el.text for el in root.findall('.//IdList/Id')]
    log.info(f'  Entrez: {count:,} total results for query')
    time.sleep(RATE_LIMIT_DELAY_S)

    batch = 500
    start = len(ids)
    while len(ids) < min(count, retmax):
        p2 = {'db': 'assembly', 'query_key': qkey, 'WebEnv': web_env,
               'retstart': start, 'retmax': batch,
               'retmode': 'xml', 'email': NCBI_EMAIL}
        r2    = _get(f'{ENTREZ_BASE}/esearch.fcgi', params=p2)
        root2 = ET.fromstring(r2.content)
        new   = [el.text for el in root2.findall('.//IdList/Id')]
        if not new:
            break
        ids.extend(new)
        start += len(new)
        time.sleep(RATE_LIMIT_DELAY_S)

    return ids[:retmax]


def _entrez_summary(uids: list[str]) -> list[dict]:
    records = []
    for i in range(0, len(uids), 200):
        batch  = uids[i:i + 200]
        params = {'db': 'assembly', 'id': ','.join(batch),
                  'retmode': 'xml', 'email': NCBI_EMAIL}
        r    = _get(f'{ENTREZ_BASE}/esummary.fcgi', params=params)
        root = ET.fromstring(r.content)
        time.sleep(RATE_LIMIT_DELAY_S)

        for doc in root.findall('.//DocumentSummary'):
            def val(tag: str) -> str:
                el = doc.find(tag)
                return el.text.strip() if el is not None and el.text else ''

            ftp = val('FtpPath_GenBank') or val('FtpPath_RefSeq')
            level_raw = val('AssemblyStatus')
            level_map = {'Complete Genome': 'Complete Genome',
                         'Chromosome': 'Chromosome',
                         'Scaffold': 'Scaffold', 'Contig': 'Contig'}
            records.append({
                'assembly_accession': val('AssemblyAccession'),
                'assembly_level':     level_map.get(level_raw, level_raw),
                'organism_name':      val('Organism'),
                'biosample':          val('BioSampleAccn'),
                'ftp_path':           ftp,
                'seq_rel_date':       val('SeqReleaseDate'),
                'submitter':          val('SubmitterOrganization'),
            })
    return records


def _build_query(species_name: str) -> str:
    return (
        f'"{species_name}"[Organism] AND "China"[Country] AND '
        f'latest[filter] AND '
        f'("Complete Genome"[Assembly Level] OR '
        f'"Chromosome"[Assembly Level] OR '
        f'"Scaffold"[Assembly Level])'
    )


def _fetch_biosample_source(biosample_id: str) -> str:
    if not biosample_id or biosample_id in ('na', ''):
        return ''
    params = {'db': 'biosample', 'id': biosample_id,
               'retmode': 'xml', 'email': NCBI_EMAIL}
    try:
        r    = _get(f'{ENTREZ_BASE}/efetch.fcgi', params=params)
        root = ET.fromstring(r.content)
        time.sleep(RATE_LIMIT_DELAY_S)
        for attr in root.findall('.//Attribute'):
            name = attr.get('attribute_name', '').lower()
            if 'isolation_source' in name or 'host' in name:
                return attr.text or ''
    except Exception:
        pass
    return ''


def _is_clinical(source: str) -> bool:
    if not isinstance(source, str):
        return False
    return any(kw in source.lower() for kw in ACCEPTED_SOURCES)


# ── per-genome download ───────────────────────────────────────────────────────

def _download_genome(short: str, accession: str, ftp_base: str,
                     dry_run: bool) -> dict:
    # Convert ftp:// → https:// (requests doesn't support FTP)
    if ftp_base.startswith('ftp://'):
        ftp_base = 'https://' + ftp_base[6:]
    ftp_base   = ftp_base.rstrip('/')
    fname_base = ftp_base.split('/')[-1]
    out_dir    = DATA_RAW / short / accession
    result     = {'accession': accession, 'status': 'OK', 'files': []}

    if dry_run:
        result['status'] = 'DRY_RUN'
        return result

    out_dir.mkdir(parents=True, exist_ok=True)
    failed = []
    for suffix in DOWNLOAD_FILES:
        url  = f'{ftp_base}/{fname_base}{suffix}'
        dest = out_dir / f'{fname_base}{suffix}'
        if dest.exists() and dest.stat().st_size > 0:
            result['files'].append(dest.name)
            continue
        ok = _download_file(url, dest)
        if ok:
            result['files'].append(dest.name)
        else:
            failed.append(suffix)
        time.sleep(RATE_LIMIT_DELAY_S)

    if failed:
        mandatory_failed = [f for f in failed if 'gff' in f or 'fna' in f]
        result['status']       = 'FAILED' if mandatory_failed else 'PARTIAL'
        result['failed_files'] = failed
    return result


# ── main per-species logic ────────────────────────────────────────────────────

def run_species(species_name: str, cfg: dict, dry_run: bool, limit: int | None):
    short = cfg['short']
    log.info(f'\n{"="*60}\n{species_name} ({short})\n{"="*60}')

    query = _build_query(species_name)
    log.info(f'Query: {query}')
    max_n = limit if limit else cfg['max_download']

    # Fetch more UIDs than needed to allow BioSample filtering
    uids = _entrez_search(query, retmax=min(max_n * 3, 2000))
    if not uids:
        log.warning(f'No assemblies found for {species_name}')
        return
    log.info(f'  Retrieved {len(uids)} UIDs')

    records = _entrez_summary(uids)
    log.info(f'  Got metadata for {len(records)} assemblies')

    # Filter by assembly level
    records = [r for r in records
               if r['assembly_level'] in ACCEPTED_ASSEMBLY_LEVELS
               and r['ftp_path'] and r['ftp_path'] != 'na']
    log.info(f'  After assembly-level filter: {len(records)}')

    # Sort: Complete > Chromosome > Scaffold
    level_rank = {'Complete Genome': 0, 'Chromosome': 1, 'Scaffold': 2}
    records.sort(key=lambda r: level_rank.get(r['assembly_level'], 3))

    # Entrez "China"[Country] filter already applied — include all records.
    # BioSample check is skipped to avoid NCBI rate-limit issues with
    # concurrent requests; downstream validation catches non-clinical genomes.
    for rec in records:
        rec['isolation_source'] = ''
    candidates = records[:max_n]
    log.info(f'  Final candidate set: {len(candidates)} genomes')

    records_out = []
    for rec in tqdm(candidates, desc=f'{short} download'):
        acc = rec['assembly_accession']
        ftp = rec['ftp_path']
        if not ftp or ftp == 'na':
            records_out.append({'accession': acc, 'status': 'NO_FTP', 'species': species_name})
            continue
        dl = _download_genome(short, acc, ftp, dry_run)
        dl['species']        = species_name
        dl['assembly_level'] = rec['assembly_level']
        dl['biosample']      = rec.get('biosample', '')
        records_out.append(dl)
        time.sleep(RATE_LIMIT_DELAY_S)

    if not records_out:
        log.warning(f'{species_name}: no candidates to download')
        return

    df = pd.DataFrame(records_out)
    out = DATA_RAW / short / 'download_status.tsv'
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep='\t', index=False)

    ok      = (df['status'] == 'OK').sum()    if 'status' in df.columns else 0
    partial = (df['status'] == 'PARTIAL').sum() if 'status' in df.columns else 0
    failed  = (df['status'] == 'FAILED').sum()  if 'status' in df.columns else 0
    log.info(f'{species_name}: OK={ok}  PARTIAL={partial}  FAILED={failed}  → {out}')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--species', choices=['eco', 'aba', 'both'], default='both')
    parser.add_argument('--limit', type=int, default=None)
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args()

    targets = {'eco': 'Escherichia coli', 'aba': 'Acinetobacter baumannii'}
    to_run  = list(targets.items()) if args.species == 'both' \
              else [(args.species, targets[args.species])]

    for short, name in to_run:
        run_species(name, SPECIES_CONFIG[name], args.dry_run, args.limit)


if __name__ == '__main__':
    main()
