"""Figure: marker-allele composition, lysis outcome and prophage copy number.

Reproduces fig5_lysis_outcome.png. Inputs: handoff/final.json, handoff/locus.json,
handoff/boot.json (stages 04-06).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C
from figstyle import apply_figure_style, panel_letter, set_frame
import json
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


META_GREY = "#888888"


def apply_figure_style(*, frame="open", font=None, sizes=(8, 7, 6), grid=False):
    import matplotlib as mpl
    if frame not in ("open", "boxed", "none"):
        raise ValueError(f"frame must be 'open'|'boxed'|'none', got {frame!r}")
    try:
        import os, sys, glob, matplotlib.font_manager as fm
        fdir = os.path.join(os.environ.get("CONDA_PREFIX") or sys.prefix, "fonts")
        if os.path.isdir(fdir):
            known = {f.fname for f in fm.fontManager.ttflist}
            for f in glob.glob(os.path.join(fdir, "*.ttf")):
                if f not in known:
                    fm.fontManager.addfont(f)
    except Exception:
        pass
    base, secondary, tick = sizes
    boxed = (frame == "boxed")
    rc = {
        "font.family": "sans-serif",
        "font.size": base,
        "axes.labelsize": base,
        "axes.titlesize": base,
        "legend.fontsize": secondary,
        "xtick.labelsize": tick,
        "ytick.labelsize": tick,
        "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": boxed, "axes.spines.right": boxed,
        "axes.spines.left": frame != "none", "axes.spines.bottom": frame != "none",
        "axes.grid": bool(grid),
        "legend.frameon": False,
        "figure.dpi": 200,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelweight": "normal",
        "lines.linewidth": 1.2,
        "patch.linewidth": 0.6,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
    mpl.rcParams.update(rc)


def set_frame(ax, style="open"):
    show = {"open": (False, False, True, True),
            "boxed": (True, True, True, True),
            "none": (False, False, False, False)}[style]
    for side, vis in zip(("top", "right", "bottom", "left"), show):
        ax.spines[side].set_visible(vis)
        if vis:
            ax.spines[side].set_linewidth(0.6)
    ax.tick_params(direction="out", length=0 if style == "none" else 3, width=0.6)


def panel_letter(ax, letter, dx=-0.18, dy=1.02, case="lower", fontsize=None):
    import matplotlib.pyplot as plt
    if fontsize is None:
        fontsize = plt.rcParams.get("font.size", 8) + 1
    s = letter.lower() if case == "lower" else letter.upper()
    ax.text(dx, dy, s, transform=ax.transAxes,
            fontweight="bold", fontsize=fontsize, va="bottom", ha="left")


apply_figure_style(sizes=(8,7,6))

F  = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/506176e8-8886-4aa8-9465-fdcb600e0182/v53ffa67c_final.json'))
LOC= json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/eeb2728d-f5be-43ca-bd4a-dba76e84574a/v2f401975_locus.json'))
CN = json.load(open('/Users/ido/.claude-science/orgs/5fa4700b-f122-43a2-8399-05ad798f8dbe/artifacts/proj_56352a1d239a/715e108f-a9eb-430c-b382-144814974347/ve8733ba6_boot.json'))['cn']
S  = ['H%d'%i for i in range(1,12)]
T30= set(S[:6])

C_IA, C_CAT, C_IIA = '#2b6cb0', '#d98c00', '#8d4bb0'   # kanR-Ia / cmR-cat / kanR-IIa(lysis-)
C30, C37 = '#4a7fb5', '#c0392b'

fig = plt.figure(figsize=(7.2,5.6))
gs  = fig.add_gridspec(2,2, hspace=0.62, wspace=0.34,
                       left=0.085, right=0.975, top=0.90, bottom=0.095)
axA, axB, axC, axD = (fig.add_subplot(gs[0,0]), fig.add_subplot(gs[0,1]),
                      fig.add_subplot(gs[1,0]), fig.add_subplot(gs[1,1]))

# -- a: allele composition ------------------------------------------------
x = np.arange(11)
pIa  = np.array([F[s]['p_Ia']['f']  for s in S])
pCat = np.array([F[s]['p_cat']['f'] for s in S])
pIIa = np.array([F[s]['p_IIa']['f'] for s in S])
axA.bar(x, pIa,  .72, color=C_IA,  lw=.4, ec='w')
axA.bar(x, pCat, .72, bottom=pIa, color=C_CAT, lw=.4, ec='w')
axA.bar(x, pIIa, .72, bottom=pIa+pCat, color=C_IIA, lw=.4, ec='w')
axA.errorbar(x[pIIa>.01], (pIa+pCat+pIIa)[pIIa>.01],
             yerr=[ (pIIa-np.array([F[s]['p_IIa']['lo'] for s in S]))[pIIa>.01],
                    (np.array([F[s]['p_IIa']['hi'] for s in S])-pIIa)[pIIa>.01] ],
             fmt='none', ecolor='0.25', elinewidth=.8, capsize=1.6)
axA.axhline(1, color=META_GREY, lw=.6, ls=':', zorder=0)
axA.set_xticks(x); axA.set_xticklabels(S, rotation=90)
axA.set_ylim(0,1.14); axA.set_ylabel('fraction of prophages')
axA.set_title('Every mixture retains both competitors', loc='left')
for c,lab,yy in [(C_IA,"kan$^R$ $aph(3')$-$Ia$  lysis$^+$",.93),
                 (C_CAT,'cm$^R$ $catA1$  lysis$^+$',.845),
                 (C_IIA,"kan$^R$ $aph(3')$-$IIa$  lysis$^-$",.76)]:
    axA.add_patch(plt.Rectangle((6.45,yy),.42,.062,color=c,transform=axA.transData,clip_on=False))
    axA.text(7.0, yy+.031, lab, va='center', fontsize=6)

# -- b: the two lysis assays agree ---------------------------------------
pred = pIa + pCat
obs  = np.array([F[s]['lysis']['f'] for s in S])
elo  = obs - np.array([F[s]['lysis']['lo'] for s in S])
ehi  = np.array([F[s]['lysis']['hi'] for s in S]) - obs
axB.plot([0,1],[0,1], color=META_GREY, lw=.7, ls='--', zorder=0)
for i,s in enumerate(S):
    c = C30 if s in T30 else C37
    axB.errorbar(pred[i], obs[i], yerr=[[elo[i]],[ehi[i]]], fmt='o', ms=4.2,
                 color=c, mec='w', mew=.5, ecolor=c, elinewidth=.9, capsize=1.6)
for s,dx,dy,ha in [('H3',.03,-.015,'left'),('H9',.012,-.075,'left'),('H7',-.03,.02,'right')]:
    i=S.index(s)
    axB.annotate(s, (pred[i],obs[i]), (pred[i]+dx, obs[i]+dy), fontsize=6,
                 ha=ha, va='center', color=META_GREY,
                 arrowprops=dict(arrowstyle='-', lw=.5, color=META_GREY, shrinkA=0, shrinkB=2))
axB.set_xlabel('predicted from cassette identity'); axB.set_ylabel('measured at lysis junction')
axB.set_title('Lysis genotype follows the marker', loc='left')
axB.set_xlim(-.06,1.06); axB.set_ylim(-.06,1.06)
axB.text(.97,.06,'30 °C', color=C30, fontsize=6, ha='right'); axB.text(.97,.005,'37 °C', color=C37, fontsize=6, ha='right')

# -- c: who won -----------------------------------------------------------
MIX = [('H4',.50,30),('H5',.80,30),('H9',.50,37),('H8',.80,37)]
for s,st,t in MIX:
    c = C30 if t==30 else C37
    y1 = 1-F[s]['lysis']['f']
    lo = 1-F[s]['lysis']['hi']; hi = 1-F[s]['lysis']['lo']
    axC.plot([0,1],[st,y1], '-', color=c, lw=1.2, zorder=2)
    axC.errorbar([1],[y1], yerr=[[y1-lo],[hi-y1]], fmt='o', ms=4.6, color=c,
                 mec='w', mew=.5, ecolor=c, elinewidth=.9, capsize=2, zorder=3)
    axC.plot([0],[st],'o', ms=3.4, color=c, mec='w', mew=.5, zorder=3)
    axC.text(1.045, y1, '%s  %.2f'%(s,y1), fontsize=6, va='center', color=c)
axC.axhline(.5, color=META_GREY, lw=.6, ls=':', zorder=0)
axC.set_xticks([0,1]); axC.set_xticklabels(['inoculum','day 8'])
axC.set_xlim(-.12,1.38); axC.set_ylim(.33,.92)
axC.set_ylabel('fraction lysis-deficient')
axC.set_title('Lysis-deficient strain converges to parity', loc='left')

# -- d: prophage copy number ---------------------------------------------
cn  = np.array([CN[s]['cn'] for s in S]); sd = np.array([CN[s]['null_sd'] for s in S])
cols= [C30 if s in T30 else C37 for s in S]
axD.errorbar(x, cn, yerr=sd, fmt='o', ms=4, mec='w', mew=.5, elinewidth=.9, capsize=1.8,
             ecolor='0.35', color='none', zorder=3)
for i in range(11):
    axD.plot(x[i], cn[i], 'o', ms=4, color=cols[i], mec='w', mew=.5, zorder=4)
axD.axhline(1, color=META_GREY, lw=.6, ls=':', zorder=0)
axD.axvspan(5.5, 10.5, color=C37, alpha=.055, zorder=0)
axD.set_xticks(x); axD.set_xticklabels(S, rotation=90)
axD.set_ylabel('prophage copies per chromosome')
axD.set_title('Induction at 37 °C raises prophage copy number', loc='left')
axD.text(8.0, 3.15, '37 °C', color=C37, fontsize=6, ha='center')
axD.text(2.5, 3.15, '30 °C', color=C30, fontsize=6, ha='center')
axD.set_ylim(.45,3.45)

for a,l in [(axA,'a'),(axB,'b'),(axC,'c'),(axD,'d')]:
    set_frame(a); panel_letter(a, l)
fig.savefig(os.path.join(C.ROOT, 'fig5_lysis_outcome.png'), dpi=300)

for t in list(axA.texts): t.remove()
for p in list(axA.patches)[33:]: p.remove()
axA.set_ylim(0,1.40)
for c,lab,yy in [(C_IA, "$aph(3')$-$Ia$   (kan$^R$, lysis$^+$)", 1.29),
                 (C_CAT,"$catA1$   (cm$^R$, lysis$^+$)",        1.20),
                 (C_IIA,"$aph(3')$-$IIa$   (kan$^R$, lysis$^-$)",1.11)]:
    axA.add_patch(plt.Rectangle((-0.35,yy),.45,.055,color=c,clip_on=False))
    axA.text(0.30, yy+.028, lab, va='center', fontsize=6)

for t in list(axC.texts): t.remove()
for s,dy in [('H4',.014),('H5',0),('H9',-.016),('H8',.013)]:
    y1 = 1-F[s]['lysis']['f']
    axC.text(1.055, y1+dy, '%s  %.2f'%(s,y1), fontsize=6, va='center',
             color=C30 if s in T30 else C37)
fig.savefig(os.path.join(C.ROOT, 'fig5_lysis_outcome.png'), dpi=300)

for t in list(axA.texts): t.remove()
for p in list(axA.patches)[33:]: p.remove()
axA.set_ylim(0,1.40)
for c,lab,yy in [(C_IA, "$aph(3')$-$Ia$   (kan$^R$, lysis$^+$)", 1.29),
                 (C_CAT,"$catA1$   (cm$^R$, lysis$^+$)",        1.20),
                 (C_IIA,"$aph(3')$-$IIa$   (kan$^R$, lysis$^-$)",1.11)]:
    axA.add_patch(plt.Rectangle((-0.35,yy),.45,.055,color=c,clip_on=False))
    axA.text(0.30, yy+.028, lab, va='center', fontsize=6)

for t in list(axC.texts): t.remove()
for s,dy in [('H4',.014),('H5',0),('H9',-.016),('H8',.013)]:
    y1 = 1-F[s]['lysis']['f']
    axC.text(1.055, y1+dy, '%s  %.2f'%(s,y1), fontsize=6, va='center',
             color=C30 if s in T30 else C37)
fig.savefig(os.path.join(C.ROOT, 'fig5_lysis_outcome.png'), dpi=300)

panel_letter(axA,'a'); panel_letter(axC,'c')
fig.savefig(os.path.join(C.ROOT, 'fig5_lysis_outcome.png'), dpi=300)