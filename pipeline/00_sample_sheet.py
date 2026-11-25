#!/usr/bin/env python
"""Build sample_sheet.csv from the lab's own plate record + the read archives.

The lab records conditions as a plate map: one sheet block per plate, rows =
serial-passage day, columns = condition.  A sequenced library is named
v<plate><row><col> (e.g. v1B3 = plate 1, row B, column 3), so the library name
alone determines replicate, day and condition.  Nothing here is hand-entered.

Edit the INPUTS block for a new experiment; everything else follows from it.
"""
import os, sys, csv, re, glob, zipfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import config_inputs as I          # see pipeline/config_inputs.py

ROW_DAY = {c: i + 1 for i, c in enumerate('ABCDEFGH')}


def plate_map(xlsx):
    """-> {plate_letter: {row_letter: {col_int: condition}}} from the plate blocks."""
    import openpyxl
    ws = openpyxl.load_workbook(xlsx, data_only=True).worksheets[0]
    grid = [list(r) for r in ws.iter_rows(values_only=True)]
    plates, cur = {}, None
    for r in grid:
        a = (str(r[0]) if r[0] is not None else '').strip()
        m = re.match(r'plate\s*(\d+)', a, re.I)
        if m:
            cur = 'v' + m.group(1)
            plates[cur] = {}
            continue
        if cur and a in ROW_DAY:
            cols = {i: str(v).strip() for i, v in enumerate(r[1:12], start=1)
                    if v is not None and str(v).strip()}
            if cols:
                plates[cur][a] = cols
    return plates


def group_key(cond):
    """Stable key for 'the same condition across days', from the lab's own text."""
    c = cond.lower()
    c = re.sub(r'ci\s*857', '', c)
    c = re.sub(r'\s+', ' ', c.replace('%', ' pct ')).strip()
    return re.sub(r'[^a-z0-9]+', '_', c).strip('_')


def classify(cond):
    """-> (temp_C, control_allele, lysis_plus_only) read off the condition text."""
    c = cond.lower()
    temp = 37 if '37' in c else 30
    ctrl = ''
    if '100' in c:
        if 'cmr' in c:
            ctrl = 'allele_cat'
        elif 'lysis-' in c or 'lysis -' in c:
            ctrl = 'allele_IIa'
        else:
            ctrl = 'allele_Ia'
    # a library contains the lysis-deficient strain only if the label says +/-
    # or says 100% lysis- ; everything else starts lysis-competent only
    has_minus = ('+/-' in c) or ('lysis-' in c) or ('lysis -' in c)
    return temp, ctrl, ('' if has_minus else '1')


def main():
    plates = plate_map(I.CONDITIONS_XLSX)
    members = {}
    if I.READ_ZIP:
        with zipfile.ZipFile(I.READ_ZIP) as z:
            for n in z.namelist():
                m = re.search(r'_(v\d+[A-H]\d+)\.fastq\.gz$', os.path.basename(n))
                if m:
                    members[m.group(1)] = (n, z.getinfo(n).file_size)

    rows = []
    for lib, (member, size) in sorted(members.items()):
        plate, rowl, col = re.match(r'(v\d+)([A-H])(\d+)', lib).groups()
        col = int(col)
        cond = plates.get(plate, {}).get(rowl, {}).get(col)
        if not cond:
            print('  ! no condition for %s' % lib); continue
        temp, ctrl, lpo = classify(cond)
        rows.append({'sample': lib, 'plate': plate, 'rep': plate[1:],
                     'day': ROW_DAY[rowl], 'well_row': rowl, 'column': col,
                     'condition': cond, 'group': group_key(cond), 'temp_C': temp,
                     'control_allele': ctrl, 'lysis_plus_only': lpo,
                     'fastq': 'zip:%s!%s' % (I.READ_ZIP, member),
                     'bytes': size,
                     'use': '0' if size < I.MIN_FASTQ_BYTES else '1',
                     'note': 'library too small (%d B)' % size if size < I.MIN_FASTQ_BYTES else ''})

    # libraries sequenced in an earlier run, given as <plate><row> + a glob
    for plate, rowl, pattern in I.EXTRA_RUNS:
        for p in sorted(glob.glob(pattern)):
            m = re.search(I.EXTRA_COL_RE, os.path.basename(p))
            if not m:
                continue
            col = int(m.group(1))
            cond = plates.get(plate, {}).get(rowl, {}).get(col)
            if not cond:
                print('  ! no condition for %s%s%d' % (plate, rowl, col)); continue
            temp, ctrl, lpo = classify(cond)
            rows.append({'sample': '%s%s%d' % (plate, rowl, col), 'plate': plate,
                         'rep': plate[1:], 'day': ROW_DAY[rowl], 'well_row': rowl,
                         'column': col, 'condition': cond, 'group': group_key(cond),
                         'temp_C': temp, 'control_allele': ctrl, 'lysis_plus_only': lpo,
                         'fastq': p, 'bytes': os.path.getsize(p), 'use': '1',
                         'note': 'original barcode %s' % re.sub(r'\.fastq\.gz$', '',
                                                                os.path.basename(p))})

    # flag columns whose condition differs between days -- the lab's record is
    # authoritative, but a difference must be visible, not silently absorbed
    bycol = {}
    for r in rows:
        bycol.setdefault(r['column'], set()).add(r['condition'])
    for r in rows:
        if len(bycol[r['column']]) > 1:
            r['note'] = (r['note'] + '; ' if r['note'] else '') + \
                        'column %d condition differs between days (record as written)' % r['column']

    rows.sort(key=lambda d: (d['rep'], d['day'], d['column']))
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       'sample_sheet.csv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    print('wrote %s: %d libraries (%d excluded)'
          % (out, len(rows), sum(1 for r in rows if r['use'] == '0')))
    days = sorted({(r['rep'], r['day']) for r in rows})
    print('rep/day coverage:', ', '.join('rep%s d%d' % d for d in days))
    print('controls:', ', '.join(sorted(r['sample'] for r in rows if r['control_allele'])))


if __name__ == '__main__':
    main()
