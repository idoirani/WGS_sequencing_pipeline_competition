"""Figure: receptor loss, insertion-junction map, activator, tail-fibre recombinants.

Reproduces fig6_receptor_loss.png. Inputs: handoff/lamb_final.json,
handoff/junctions.json, handoff/recfreq2.json, handoff/tfa.json (stages 09-12).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from figstyle import apply_figure_style, panel_letter, set_frame
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# figure-style skill helpers (auto-injected)
def apply_figure_style(sizes=(8,7,6)):
    mpl.rcParams.update({
        'font.size': sizes[1],
        'axes.labelsize': sizes[1],
        'axes.titlesize': sizes[1],
        'xtick.labelsize': sizes[2],
        'ytick.labelsize': sizes[2],
        'legend.fontsize': sizes[2],
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'axes.linewidth': 0.6,
        'xtick.major.width': 0.6,
        'ytick.major.width': 0.6,
        'xtick.major.size': 2.5,
        'ytick.major.size': 2.5,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })

def set_frame(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

def panel_letter(ax, letter):
    ax.text(-0.13, 1.06, letter, transform=ax.transAxes,
            fontsize=9, fontweight='bold', va='top', ha='left')

apply_figure_style(sizes=(8,7,6))

LAM = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/febc290e-ec97-41f8-a55a-ff6e2d09ccfb/v293d3b83_lamb_final.json'))
JX  = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/bcf09069-2ea2-4432-bf6c-c4d192500ff4/v0e478774_junctions.json'))
MAL = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/76345acc-cca9-4396-b633-faea4e003692/v9d5c4581_recfreq2.json'))
TFA = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/27458e41-87bd-4647-8726-a92f670dca5a/v03b637ee_tfa.json'))
S   = C.SAMPLES
T37 = set(S[6:])                                   # H7..H11
COND = {'H1':'pure kan$^R$ lys$^+$','H2':'pure cm$^R$ lys$^+$','H3':'pure kan$^R$ lys$^-$',
        'H4':'50/50 +/−','H5':'20/80 +/−','H6':'double lysogen',
        'H7':'50/50 +/+','H8':'20/80 +/−','H9':'50/50 +/−','H10':'20/80 +/+','H11':'double lysogen'}

# IS-family split of the IS read pool, from per-junction read counts
C_IS1, C_IS2, C_IS5, C_NS, C_FS, C_WT = '#2b6cb0', '#7b3fa0', '#1b9e77', '#c0392b', '#d98c00', '#d9d9d9'
comp = {}
for s in S:
    v = LAM[s]; n = v['n']
    jr = {f:0 for f in ('IS1','IS2','IS5')}
    for j in JX:
        if j['sample']==s and j['gene']=='lamB' and j['family'] in jr: jr[j['family']] += j['reads']
    tj = sum(jr.values())
    isf = {f: (v['IS']*jr[f]/tj if tj else 0)/n for f in jr}
    if tj == 0 and v['IS']: isf['IS1'] = v['IS']/n
    comp[s] = dict(isf, NS=v['NS']/n, FS=v['FS']/n, WT=v['WT']/n, n=n,
                   f=v['f_dis'], lo=v['ci'][0], hi=v['ci'][1])
FSBG = sum(LAM[s]['FS'] for s in S[:6]) / sum(LAM[s]['n'] for s in S[:6])

fig = plt.figure(figsize=(7.2, 6.6))
gs  = fig.add_gridspec(2, 2, height_ratios=[1, 1.08], hspace=0.52, wspace=0.30,
                       left=0.085, right=0.985, top=0.90, bottom=0.085)
axA, axB = fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1])
axC, axD = fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1])
x = np.arange(11)

# ---- a: lamB allele composition -------------------------------------------
order = [('IS1',C_IS1,'IS$1$'),('IS2',C_IS2,'IS$2$'),('IS5',C_IS5,'IS$5$'),
         ('NS',C_NS,'nonsense'),('FS',C_FS,'frameshift'),('WT',C_WT,'intact')]
bot = np.zeros(11)
for k,c,lab in order:
    h = np.array([comp[s][k] for s in S])
    axA.bar(x, h, 0.74, bottom=bot, color=c, edgecolor='white', linewidth=0.4, label=lab)
    bot += h
axA.axhline(FSBG, color='0.35', lw=0.8, ls=(0,(4,2)), zorder=5)
axA.set_ylim(0, 1.48); axA.set_yticks([0,0.25,0.5,0.75,1.0])
axA.set_ylabel('fraction of $lamB$ alleles')
axA.set_title('Receptor is lost only at 37 °C', loc='left')
for i,(k,c,lab) in enumerate(order):
    col, row = i//3, i%3
    axA.add_patch(plt.Rectangle((0.15+col*5.0, 1.40-row*0.095), 0.42, 0.062,
                                color=c, clip_on=False))
    axA.text(0.72+col*5.0, 1.431-row*0.095, lab, va='center', fontsize=6)

# error-floor legend entry
axA.plot([0.15, 0.57], [1.146, 1.146], color='0.35', lw=0.9, ls=(0,(3,2)), clip_on=False)
axA.text(0.72, 1.146, 'sequencing error floor (4.9%)', va='center', fontsize=6, color='0.3')

# ---- b: junction map along lamB --------------------------------------------
GLO, GHI = 4324753, 4326093
ypos = {s:i for i,s in enumerate(S[6:][::-1])}
axB.add_patch(plt.Rectangle((GLO, -0.75), GHI-GLO, 0.30, color='0.86', zorder=1))
axB.text((GLO+GHI)/2, -0.60, '$lamB$  (446 codons, ←)', ha='center', va='center',
         fontsize=6, color='0.25', zorder=3)
for j in JX:
    if j['gene']!='lamB': continue
    c = {'IS1':C_IS1,'IS2':C_IS2,'IS5':C_IS5}[j['family']]
    axB.scatter(j['pos'], ypos[j['sample']], s=14+5.5*j['reads'], color=c,
                alpha=0.85, edgecolor='white', linewidth=0.5, zorder=4)
for p,lab in [(4325749,'frameshift\n(H11)'),(4325797,'nonsense\n(H7)')]:
    axB.axvline(p, color=C_FS if 'frame' in lab else C_NS, lw=0.9, ls=(0,(3,2)), zorder=2)
axB.set_ylim(-0.95, 5.45)
axB.text(4325749-18, 4.95, 'frameshift', ha='right', va='center', fontsize=6, color=C_FS)
axB.text(4325797+18, 4.95, 'nonsense',  ha='left',  va='center', fontsize=6, color=C_NS)
axB.set_yticks(range(5)); axB.set_yticklabels(S[6:][::-1])
axB.set_xlim(GLO-40, GHI+40)
axB.set_xticks([4324800, 4325300, 4325800, 4326093])
axB.set_xticklabels(['446', '280', '113', '1'])
axB.set_xlabel('codon position in $lamB$')
axB.set_title('15 independent insertion junctions', loc='left')
axB.text(GLO+60, 4.42, 'marker area ∝ reads', fontsize=6, color='0.35')

# ---- c: malT -----------------------------------------------------------------
fm  = np.array([MAL['%s|malT'%s]['f'] or 0 for s in S])
mlo = np.array([MAL['%s|malT'%s]['lo'] or 0 for s in S])
mhi = np.array([MAL['%s|malT'%s]['hi'] or 0 for s in S])
axC.bar(x, fm, 0.74, color=C_IS5, edgecolor='none')
axC.errorbar(x, fm, yerr=[np.clip(fm-mlo, 0, None), np.clip(mhi-fm, 0, None)],
             fmt='none', ecolor='0.25', elinewidth=0.8, capsize=1.6)
axC.scatter(x[fm==0], np.zeros((fm==0).sum()), s=5, color='0.45', zorder=4)
axC.set_ylim(0, 0.42); axC.set_yticks([0, 0.1, 0.2, 0.3])
axC.set_ylabel('fraction of $malT$ alleles')
axC.set_title('Activator disrupted at a single recurrent site', loc='left')
axC.text(-0.35, 0.385, 'IS$5$ at the identical base pair (330,623)\nin H7–H10: one ancestral allele, not four events',
         fontsize=6, color='0.3', linespacing=1.3, va='top')

# ---- d: tail-fibre recombinants ---------------------------------------------
ft  = np.array([TFA[s]['f'] or 0 for s in S])
tlo = np.array([TFA[s]['ci'][0] or 0 for s in S])
thi = np.array([TFA[s]['ci'][1] or 0 for s in S])
axD.bar(x, ft, 0.74, color='#b4484f', edgecolor='none')
axD.errorbar(x, ft, yerr=[np.clip(ft-tlo, 0, None), np.clip(thi-ft, 0, None)],
             fmt='none', ecolor='0.25', elinewidth=0.8, capsize=1.6)
axD.scatter(x[ft==0], np.zeros((ft==0).sum()), s=5, color='0.45', zorder=4)
axD.set_ylim(0, 0.275)
axD.set_ylabel('fraction of prophage $tfa$ alleles')
axD.set_title('Prophage tail fibre recombines with Rac', loc='left')
axD.text(6.45, 0.215, '16%', ha='left', fontsize=6.5, color='#b4484f')

for ax in (axA, axC, axD):
    ax.set_xticks(x); ax.set_xticklabels(S, fontsize=6)
    ax.axvline(5.5, color='0.55', lw=0.8, ls=(0,(2,2)))
    ax.text(2.5, -0.175, '30 °C', ha='center', transform=ax.get_xaxis_transform(),
            fontsize=6.5, color='0.25')
    ax.text(8.5, -0.175, '37 °C', ha='center', transform=ax.get_xaxis_transform(),
            fontsize=6.5, color='0.25')
for ax, L in ((axA,'a'), (axB,'b'), (axC,'c'), (axD,'d')):
    set_frame(ax); panel_letter(ax, L)

fig.savefig(os.path.join(C.ROOT, 'fig6_receptor_loss.png'), dpi=300, bbox_inches='tight')
# ---------------------------------------------------------------- corrections
# Applied after inspecting the first render. Each fixes a label that collided
# with data; none changes a plotted value.

# (a) the error-floor note sat across the line it describes -> legend entry
for t in list(axA.texts):
    if 'error' in t.get_text():
        t.remove()
axA.plot([0.15, 0.57], [1.146, 1.146], color='0.35', lw=0.9, ls=(0, (3, 2)), clip_on=False)
axA.text(0.72, 1.146, 'sequencing error floor (4.9%)', va='center', fontsize=6, color='0.3')

# (b) point-mutation labels overlapped the junction markers -> lifted clear
for t in list(axB.texts):
    if t.get_text() in ('frameshift', 'nonsense'):
        t.remove()
axB.set_ylim(-0.95, 5.45)
axB.text(4325749 - 18, 4.95, 'frameshift', ha='right', va='center', fontsize=6, color=C_FS)
axB.text(4325797 + 18, 4.95, 'nonsense', ha='left', va='center', fontsize=6, color=C_NS)

# (c) annotation crossed the divider and a confidence interval -> headroom
for t in list(axC.texts):
    if 'base pair' in t.get_text():
        t.remove()
axC.set_ylim(0, 0.42)
axC.set_yticks([0, 0.1, 0.2, 0.3])
axC.set_title('Activator disrupted at a single recurrent site', loc='left')
axC.text(-0.35, 0.385,
         'IS$5$ at the identical base pair (330,623)\nin H7-H10: one ancestral allele, not four events',
         fontsize=6, color='0.3', linespacing=1.3, va='top')

# (d) headline number sat on the divider -> shifted right
for t in list(axD.texts):
    if t.get_text() == '16%':
        t.remove()
axD.set_ylim(0, 0.275)
axD.text(6.45, 0.215, '16%', ha='left', fontsize=6.5, color='#b4484f')

for ax in (axA, axC, axD):
    ax.text(2.5, -0.175, '30 \u00b0C', ha='center', transform=ax.get_xaxis_transform(),
            fontsize=6.5, color='0.25')
    ax.text(8.5, -0.175, '37 \u00b0C', ha='center', transform=ax.get_xaxis_transform(),
            fontsize=6.5, color='0.25')

# the loop above re-added labels that already existed -> drop the duplicates
for ax in (axA, axC, axD):
    seen = set()
    for t in list(ax.texts):
        k = t.get_text()
        if k in ('30 \u00b0C', '37 \u00b0C'):
            if k in seen:
                t.remove()
            else:
                seen.add(k)

fig.savefig(os.path.join(C.ROOT, 'fig6_receptor_loss.png'), dpi=300, bbox_inches='tight')
