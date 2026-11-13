#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import sys

# parse shared strings
with open('PCR_Cornus_SDH/data/housekeeping/external/sharedStrings.xml', 'r', encoding='utf-8') as f:
    sst = ET.fromstring(f.read())
ns = {'ns': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
shared = []
for si in sst.findall('ns:si', ns):
    t = si.find('ns:t', ns)
    if t is not None:
        shared.append(t.text)
    else:
        # handle rich text? skip
        shared.append('')
print(f"Loaded {len(shared)} shared strings")

# parse sheet7
with open('PCR_Cornus_SDH/data/housekeeping/external/sheet7.xml', 'r', encoding='utf-8') as f:
    sheet = ET.fromstring(f.read())

ns_sheet = {'ns': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main',
            'x': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}

rows = {}
for row in sheet.findall('ns:sheetData/ns:row', ns_sheet):
    row_num = row.get('r')
    for cell in row.findall('ns:c', ns_sheet):
        col = cell.get('r')
        val = cell.find('ns:v', ns_sheet)
        if val is not None:
            v = val.text
            # check if cell type is 's' (shared string)
            t = cell.get('t')
            if t == 's':
                idx = int(v)
                if idx < len(shared):
                    v = shared[idx]
                else:
                    v = f'[shared:{idx}]'
            rows.setdefault(row_num, {})[col] = v

# print rows as table
max_col = 0
for r in rows:
    for c in rows[r]:
        col_idx = ord(c[0]) - ord('A') + 1
        if col_idx > max_col:
            max_col = col_idx

for r in sorted(rows.keys(), key=lambda x: int(x)):
    line = []
    for ci in range(max_col):
        col = chr(ord('A') + ci)
        key = col + r
        val = rows.get(r, {}).get(key, '')
        line.append(str(val))
    print('\t'.join(line))
