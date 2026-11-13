#!/usr/bin/env python3
import pandas as pd
import sys

filename = 'ijms-4064141-supplementary_rev2.xlsx'
try:
    xl = pd.ExcelFile(filename)
except Exception as e:
    print(f"Error opening Excel file: {e}")
    sys.exit(1)

print("Sheet names:", xl.sheet_names)
for sheet in xl.sheet_names:
    df = pd.read_excel(xl, sheet_name=sheet, nrows=5)
    print(f"\n--- Sheet: {sheet} ---")
    print(df.head())
    print(f"Shape: {df.shape}")