"""
Step 6 — Cross-species publication figures
────────────────────────────────────────────
Figure 1: Composite transposon rate by species + carbapenemase type (3-panel)
Figure 2: IS focal-family landscape donut charts (3-species comparison)
Figure 3: AUC comparison across species (ROC curves or bar)
Figure 4: Resistance prevalence + IS burden violin plots (3-species)
Figure 5: Temporal submission bias comparison (K. pn. vs E. coli)

Output: figures/fig{1-5}_{name}.pdf + .png

Usage:
  python analysis/06_cross_species_figures.py
"""
import json
import logging
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from config import DATA_PROC, FIGURES, LOGS, SPECIES_CONFIG

LOGS.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s  %(levelname)-8s  %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger('figures')

# ── matplotlib setup ──────────────────────────────────────────────────────────
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

SPECIES_COLORS = {
    'Klebsiella pneumoniae':  '#1b9e77',  # ColorBrewer Dark2 — teal
    'Escherichia coli':       '#d95f02',  # ColorBrewer Dark2 — orange
    'Acinetobacter baumannii':'#7570b3',  # ColorBrewer Dark2 — purple
}
SPECIES_SHORT = {
    'Klebsiella pneumoniae':  'K. pneumoniae',
    'Escherichia coli':       'E. coli',
    'Acinetobacter baumannii':'A. baumannii',
}

FIG_DPI    = 300
FIG_STYLE  = {'font.family': 'DejaVu Sans', 'font.size': 9,
               'axes.linewidth': 0.8, 'axes.spines.top': False,
               'axes.spines.right': False}


def _save(fig, name: str):
    for ext in ('pdf', 'png'):
        out = FIGURES / f'{name}.{ext}'
        fig.savefig(out, dpi=FIG_DPI, bbox_inches='tight')
    log.info(f'Saved {name}.pdf/.png')
    plt.close(fig)


# ── Figure 1: Composite transposon rates ─────────────────────────────────────

def fig1_composite_rates():
    species_order = ['Klebsiella pneumoniae', 'Escherichia coli',
                     'Acinetobacter baumannii']
    rates = {}
    for sp in species_order:
        cfg = next((c for n, c in SPECIES_CONFIG.items() if n == sp), None)
        if cfg is None:
            # K. pneumoniae: 47/50 loci = 94.0%, Clopper-Pearson 95%CI [83.5–98.8%]
            rates[sp] = {'locus': 94.0, 'ci_low': 83.5, 'ci_high': 98.8}
            continue
        short = cfg['short']
        locus_path = DATA_PROC / f'{short}_locus_classification.tsv'
        amr_path   = DATA_PROC / f'{short}_amr_hits.tsv'
        if locus_path.exists():
            ldf    = pd.read_csv(locus_path, sep='\t')
            n_comp = (ldf['flanking_class'] == 'COMPOSITE_TRANSPOSON').sum()
            # Use total AMR loci (rows in amr_hits) as denominator
            n_total = len(ldf)
            if amr_path.exists():
                amr_df = pd.read_csv(amr_path, sep='\t')
                if not amr_df.empty:
                    # ABA: exclude intrinsic OXA-51 from denominator
                    if short == 'aba' and 'gene_name' in amr_df.columns:
                        amr_df = amr_df[~amr_df['gene_name'].str.contains(
                            r'OXA-51', na=False)]
                    n_total = len(amr_df)
            rate = round(n_comp / n_total * 100, 1) if n_total else 0
            from math import sqrt
            # Wilson confidence interval
            z = 1.96; p = n_comp / n_total if n_total else 0; n = n_total
            ci_lo = max(0, round((p + z**2/(2*n) - z*sqrt(p*(1-p)/n + z**2/(4*n**2))) /
                                 (1 + z**2/n) * 100, 1)) if n else 0
            ci_hi = min(100, round((p + z**2/(2*n) + z*sqrt(p*(1-p)/n + z**2/(4*n**2))) /
                                   (1 + z**2/n) * 100, 1)) if n else 0
            rates[sp] = {'locus': rate, 'ci_low': ci_lo, 'ci_high': ci_hi}
        else:
            rates[sp] = {'locus': 0, 'ci_low': 0, 'ci_high': 0}

    plt.rcParams.update(FIG_STYLE)
    fig, ax = plt.subplots(figsize=(5.5, 3.5))

    x     = np.arange(len(species_order))
    vals  = [rates[sp]['locus'] for sp in species_order]
    lows  = [rates[sp]['locus'] - rates[sp]['ci_low'] for sp in species_order]
    highs = [rates[sp]['ci_high'] - rates[sp]['locus'] for sp in species_order]
    cols  = [SPECIES_COLORS[sp] for sp in species_order]

    bars = ax.bar(x, vals, color=cols, alpha=0.85, width=0.55,
                  edgecolor='white', linewidth=0.5)
    ax.errorbar(x, vals, yerr=[lows, highs], fmt='none',
                ecolor='#333333', capsize=4, linewidth=1.2)

    for bar, v in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, v + 1.5,
                f'{v:.1f}%', ha='center', va='bottom', fontsize=8, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([SPECIES_SHORT[sp] for sp in species_order],
                       fontsize=8.5, style='italic')
    ax.set_ylabel('Composite transposon rate (%, locus level)', fontsize=9)
    ax.set_ylim(0, 110)
    ax.axhline(90, color='gray', linestyle='--', linewidth=0.7, alpha=0.5)
    ax.set_title('IS-mediated composite transposon frequency diverges\n'
                 'across carbapenem-resistant Gram-negative pathogens',
                 fontsize=9, fontweight='bold', pad=8)
    _save(fig, 'fig1_composite_transposon_rates')


# ── Figure 2: IS focal burden violin plots ───────────────────────────────────

def fig2_burden_violins():
    burden_path = DATA_PROC / 'cross_species_burden.tsv'
    if not burden_path.exists():
        log.warning('cross_species_burden.tsv not found — skipping Figure 2')
        return

    df = pd.read_csv(burden_path, sep='\t')
    species_list = df['species'].unique().tolist()

    plt.rcParams.update(FIG_STYLE)
    n_sp = len(species_list)
    fig, axes = plt.subplots(1, n_sp, figsize=(3.5 * n_sp, 4), sharey=False)
    if n_sp == 1:
        axes = [axes]

    for ax, sp in zip(axes, species_list):
        res = df[(df['species'] == sp) & (df['resistant'] == True)]['focal_is_count'].tolist()
        sus = df[(df['species'] == sp) & (df['resistant'] == False)]['focal_is_count'].tolist()
        color = SPECIES_COLORS.get(sp, '#888888')
        short_name = SPECIES_SHORT.get(sp, sp)

        data   = [sus, res]
        labels = [f'Susceptible\n(n={len(sus)})', f'Resistant\n(n={len(res)})']
        colors_vp = ['#AAAAAA', color]

        parts = ax.violinplot(data, positions=[1, 2], showmedians=True,
                              showextrema=True)
        for pc, col in zip(parts['bodies'], colors_vp):
            pc.set_facecolor(col)
            pc.set_alpha(0.65)
        parts['cmedians'].set_color('#222222')
        parts['cmedians'].set_linewidth(1.5)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel(f'{df[df["species"]==sp]["focal_is_family"].iloc[0]} copies/genome',
                      fontsize=8)
        italic_name = short_name.replace(" ", "\\ ")
        ax.set_title(f'$\\it{{{italic_name}}}$',
                     fontsize=9, fontweight='bold')

        # AUC annotation
        stats_path = DATA_PROC / 'cross_species_stats.json'
        if stats_path.exists():
            with open(stats_path) as fh:
                stats = json.load(fh)
            sp_stat = next((s for s in stats if s['species'] == sp), None)
            if sp_stat:
                ax.text(0.97, 0.96,
                        f'AUC = {sp_stat["auc"]:.3f}\nδ = {sp_stat["cliffs_delta"]:.3f}',
                        transform=ax.transAxes, ha='right', va='top',
                        fontsize=7.5, bbox=dict(boxstyle='round,pad=0.3',
                        facecolor='white', edgecolor='#CCCCCC', alpha=0.9))

    fig.suptitle('Focal IS family genomic copy number predicts carbapenem resistance\n'
                 'across three Gram-negative pathogens',
                 fontsize=9.5, fontweight='bold', y=1.02)
    plt.tight_layout()
    _save(fig, 'fig2_is_burden_violins')


# ── Figure 3: AUC comparison bar chart ───────────────────────────────────────

def fig3_auc_comparison():
    stats_path = DATA_PROC / 'cross_species_stats.json'
    if not stats_path.exists():
        log.warning('cross_species_stats.json not found — skipping Figure 3')
        return

    with open(stats_path) as fh:
        stats = json.load(fh)

    plt.rcParams.update(FIG_STYLE)
    fig, ax = plt.subplots(figsize=(5, 3.5))

    species = [s['species'] for s in stats]
    aucs    = [s.get('auc', 0) for s in stats]
    cols    = [SPECIES_COLORS.get(sp, '#888888') for sp in species]
    labels  = [SPECIES_SHORT.get(sp, sp) for sp in species]

    x = np.arange(len(species))
    bars = ax.bar(x, aucs, color=cols, alpha=0.85, width=0.55,
                  edgecolor='white', linewidth=0.5)

    for bar, auc in zip(bars, aucs):
        ax.text(bar.get_x() + bar.get_width() / 2,
                auc + 0.008, f'{auc:.3f}',
                ha='center', va='bottom', fontsize=8.5, fontweight='bold')

    ax.axhline(0.5, color='gray', linestyle='--', linewidth=0.8, alpha=0.6,
               label='Random classifier (AUC=0.5)')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5, style='italic')
    ax.set_ylabel('AUC (focal IS copy number → carbapenem resistance)', fontsize=9)
    ax.set_ylim(0, 1.10)
    ax.legend(fontsize=7.5, framealpha=0.7)
    ax.set_title('IS genomic copy number as a cross-species\ncarbapenem resistance predictor',
                 fontsize=9, fontweight='bold', pad=8)
    _save(fig, 'fig3_auc_comparison')


# ── Figure 4: IS family composition donut ────────────────────────────────────

def fig4_is_family_landscape():
    plt.rcParams.update(FIG_STYLE)
    n_sp = len(SPECIES_CONFIG) + 1  # +1 for K. pn.
    fig, axes = plt.subplots(1, n_sp, figsize=(4 * n_sp, 4))
    if n_sp == 1:
        axes = [axes]

    species_list = [('Klebsiella pneumoniae', 'kpn'),
                    ('Escherichia coli', 'eco'),
                    ('Acinetobacter baumannii', 'aba')]

    for ax, (sp, short) in zip(axes, species_list):
        ctx_path = DATA_PROC / f'{short}_is_context.tsv'
        if sp == 'Klebsiella pneumoniae':
            ctx_path = DATA_PROC.parent.parent.parent / 'amr_project' / \
                       'data' / 'processed' / 'is_context.tsv'

        if not ctx_path.exists():
            ax.text(0.5, 0.5, 'Data pending', ha='center', va='center',
                    transform=ax.transAxes, fontsize=9, color='gray')
            ax.set_title(SPECIES_SHORT.get(sp, sp), fontsize=9, style='italic')
            continue

        df = pd.read_csv(ctx_path, sep='\t', dtype=str)
        if 'is_family' not in df.columns or df.empty:
            ax.text(0.5, 0.5, 'No IS data', ha='center', va='center',
                    transform=ax.transAxes)
            continue

        family_counts = df['is_family'].value_counts()
        top_n = 6
        top    = family_counts.head(top_n)
        other  = family_counts.iloc[top_n:].sum()
        if other > 0:
            top['Other'] = other

        colors = plt.cm.Set2(np.linspace(0, 1, len(top)))
        wedges, texts, autotexts = ax.pie(
            top.values, labels=top.index, colors=colors,
            autopct='%1.1f%%', pctdistance=0.8,
            wedgeprops={'width': 0.5, 'edgecolor': 'white', 'linewidth': 1},
            textprops={'fontsize': 7})
        for at in autotexts:
            at.set_fontsize(6.5)
        sp_italic = SPECIES_SHORT.get(sp, sp).replace(" ", "\\ ")
        ax.set_title(f'$\\it{{{sp_italic}}}$\n(n={len(df)} IS elements)',
                     fontsize=8.5, fontweight='bold')

    fig.suptitle('IS element family composition in carbapenem resistance gene contexts',
                 fontsize=10, fontweight='bold')
    plt.tight_layout()
    _save(fig, 'fig4_is_family_landscape')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    log.info('Generating cross-species figures…')
    fig1_composite_rates()
    fig2_burden_violins()
    fig3_auc_comparison()
    fig4_is_family_landscape()
    log.info('All figures complete.')


if __name__ == '__main__':
    main()
