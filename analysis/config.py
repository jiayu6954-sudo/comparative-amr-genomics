"""
Comparative AMR Project — Cross-species IS26/ISAba analysis
Species: Klebsiella pneumoniae · Escherichia coli · Acinetobacter baumannii
"""
from pathlib import Path

# ── project paths ─────────────────────────────────────────────────────────────
ROOT       = Path(__file__).resolve().parents[1]
DATA_RAW   = ROOT / 'data' / 'raw'
DATA_PROC  = ROOT / 'data' / 'processed'
DATA_VAL   = ROOT / 'data' / 'validated'
LOGS       = ROOT / 'logs'
FIGURES    = ROOT / 'figures'
REPORTS    = ROOT / 'reports'

# K. pneumoniae results from existing project (no re-download needed)
KPN_PROJECT = Path(__file__).resolve().parents[2] / 'amr_project'
KPN_MANIFEST    = KPN_PROJECT / 'data' / 'validated' / 'genome_manifest.tsv'
KPN_AMR_HITS    = KPN_PROJECT / 'data' / 'processed' / 'amr_hits.tsv'
KPN_IS_CONTEXT  = KPN_PROJECT / 'data' / 'processed' / 'is_context.tsv'
KPN_IS_BURDEN   = KPN_PROJECT / 'data' / 'processed' / 'is_burden_corrected_stats.json'

# ── NCBI credentials ──────────────────────────────────────────────────────────
NCBI_EMAIL   = 'jiayu6954@gmail.com'
NCBI_API_KEY = ''

# ── species-specific parameters ───────────────────────────────────────────────
SPECIES_CONFIG = {
    'Escherichia coli': {
        'short':              'eco',
        'genome_size_min_bp': 4_200_000,
        'genome_size_max_bp': 5_800_000,
        'min_cds':            3_000,
        'min_n50':            50_000,
        'max_download':       350,          # cap: China clinical complete/scaffold
        'assembly_summary_url':
            'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/'
            'Escherichia_coli/assembly_summary.txt',
        # IS family marker for this species (IS6 = IS26 clade)
        'focal_is_family':    'IS6',
        'focal_is_pattern':   r'IS26|IS257|IS1353|IS1006|IS6\b',
        # carbapenemases dominant in clinical E. coli China
        'carbapenem_genes': {
            'NDM-1','NDM-2','NDM-4','NDM-5','NDM-6','NDM-7','NDM-9',
            'NDM-12','NDM-13','NDM-14','NDM-15','NDM-16',
            'KPC-2','KPC-3',
            'OXA-48','OXA-181','OXA-232','OXA-162',
            'IMP-1','IMP-4','IMP-6','IMP-8','IMP-26',
            'VIM-1','VIM-2','VIM-4',
        },
        'amr_gff_keywords': [
            r'\bNDM-\d', r'\bKPC-\d|blaKPC',
            r'OXA-48', r'OXA-181', r'OXA-232',
            r'\bIMP-\d|blaIMP', r'\bVIM-\d|blaVIM',
            r'carbapenem.{0,30}(resistance|beta-lactamase)',
        ],
        'is_gff_keywords': [
            r'IS\d+', r'transposase', r'insertion sequence',
            r'Tn\d+', r'Tn4401',
        ],
    },

    'Acinetobacter baumannii': {
        'short':              'aba',
        'genome_size_min_bp': 3_700_000,
        'genome_size_max_bp': 4_500_000,
        'min_cds':            2_500,        # A. baumannii has fewer CDS
        'min_n50':            50_000,
        'max_download':       350,
        'assembly_summary_url':
            'https://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/'
            'Acinetobacter_baumannii/assembly_summary.txt',
        # ISAba1 is the canonical IS element flanking OXA carbapenemases
        'focal_is_family':    'ISAba',
        'focal_is_pattern':   r'ISAba\d*|IS\d+Aba|ISaba',
        'carbapenem_genes': {
            # OXA-type (dominant in A. baumannii)
            'OXA-23','OXA-24','OXA-25','OXA-26','OXA-40',
            'OXA-58','OXA-72','OXA-143','OXA-182',
            # NDM
            'NDM-1','NDM-2','NDM-9',
            # IMP
            'IMP-1','IMP-4','IMP-5','IMP-6',
            # VIM
            'VIM-1','VIM-2',
        },
        'amr_gff_keywords': [
            r'OXA-23', r'OXA-24', r'OXA-40', r'OXA-58', r'OXA-72',
            r'OXA-143', r'OXA-182', r'blaOXA',
            r'\bNDM-\d', r'\bIMP-\d|blaIMP', r'\bVIM-\d|blaVIM',
            r'carbapenem.{0,30}(resistance|beta-lactamase)',
        ],
        'is_gff_keywords': [
            r'IS\d+', r'ISAba\d*', r'transposase', r'insertion sequence',
            r'Tn\d+',
        ],
    },
}

# ── shared download settings ──────────────────────────────────────────────────
REQUIRED_ISOLATION_COUNTRY  = 'China'
ACCEPTED_ASSEMBLY_LEVELS    = {'Complete Genome', 'Chromosome', 'Scaffold'}
ACCEPTED_SOURCES = [
    'clinical', 'patient', 'blood', 'urine', 'sputum', 'wound',
    'hospital', 'human', 'infant', 'neonatal', 'respiratory',
    'cerebrospinal', 'abscess', 'bile', 'throat', 'burn',
]
DOWNLOAD_TIMEOUT_S  = 60
MAX_RETRIES         = 3
RETRY_DELAY_S       = 5
RATE_LIMIT_DELAY_S  = 0.40

DOWNLOAD_FILES = [
    '_genomic.gff.gz',
    '_genomic.fna.gz',
    '_assembly_stats.txt',
]

# IS flanking window
IS_FLANK_WINDOW_BP = 10_000

# HMMER / PFAM for IS6-family correction (same as KPN project)
PFAM_IS6_ID  = 'PF01527'   # IS6 family transposase
HMMER_EVALUE = 1e-5
HMMER_MIN_IDENTITY = 0.25
HMMER_MIN_COVERAGE = 0.80
