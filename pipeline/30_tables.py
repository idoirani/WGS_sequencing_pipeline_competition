"""Assemble the two deliverable tables from the per-stage JSON."""
import os, sys, csv, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config as C

def L(n): return json.load(open(os.path.join(C.OUT, n + '.json')))
F, LOC, CN = L('final'), L('locus'), L('boot')['cn']
LAM, JX, MAL, TFA = L('lamb_final'), L('junctions'), L('recfreq2'), L('tfa')
OUT = C.ROOT

rows = []
for s in C.SAMPLES:
    v, m, t = LAM[s], MAL['%s|malT' % s], TFA[s]
    jx = [j for j in JX if j['sample'] == s and j['gene'] == 'lamB']
    n = v['n']
    rows.append([s, C.TEMP[s], n, round(v['f_dis'], 4), round(v['ci'][0], 4), round(v['ci'][1], 4),
                 round(v['IS'] / n, 4), round(v['NS'] / n, 4), round(v['FS'] / n, 4),
                 round(v['WT'] / n, 4), len(jx),
                 '+'.join(sorted({j['family'] for j in jx})) or '-',
                 m['n'], round(m['f'] or 0, 4), round(m['lo'] or 0, 4), round(m['hi'] or 0, 4),
                 t['n'], round(t['f'] or 0, 4), round(t['ci'][0] or 0, 4), round(t['ci'][1] or 0, 4)])
hdr = ['sample','temp_C','lamB_n_reads','lamB_f_disrupted','lamB_ci_lo','lamB_ci_hi',
       'lamB_f_IS','lamB_f_nonsense','lamB_f_frameshift','lamB_f_intact','lamB_n_junctions',
       'lamB_IS_families','malT_n_reads','malT_f_disrupted','malT_ci_lo','malT_ci_hi',
       'tfa_n_reads','tfa_f_recombinant','tfa_ci_lo','tfa_ci_hi']
with open(os.path.join(OUT, 'receptor_mutations.csv'), 'w', newline='') as fh:
    w = csv.writer(fh); w.writerow(hdr); w.writerows(rows)

jr = [[j['sample'], j['gene'], j['pos'],
       C.codon(j['gene'], j['pos']) if j['gene'] in C.GENES else '',
       j['reads'], j['len'], j['family']]
      for j in sorted(JX, key=lambda d: (d['gene'], d['sample'], d['pos']))]
with open(os.path.join(OUT, 'is_junctions.csv'), 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['sample','gene','position','codon','reads','insert_len_bp','IS_family'])
    w.writerows(jr)
print('receptor_mutations.csv: %d rows | is_junctions.csv: %d rows' % (len(rows), len(jr)))
