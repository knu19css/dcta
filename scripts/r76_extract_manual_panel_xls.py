"""r76_extract_manual_panel_xls.py - 최초 수기 마커 패널 XLS (final_unique_gene.xls) 의 형식 추출만 한다 (R76 입력 변환).

두 sheet (final_unique_gene = 양성 +1, Negative markers = 음성 -1) 의 비어 있지 않은 모든 칸을 원문 그대로 꺼낸다:
sheet, Excel 행 · 열 · 열 문자, 열 머리글 원문, 기호 원문, 앞뒤 공백만 뺀 기호, 부호, xlrd 칸 형식.
이름 교정 · 중복 제거 · 마커 선택 변경은 하지 않는다 (R76 이 따로 기록하며 한다). 결과 폴더가 비어 있지 않으면 쓰지 않는다.

실행 (W 에서, 기존 환경의 xlrd 만 사용):
  /Users/edemhw/miniforge3/envs/bonemarrow_scrna/bin/python -I scripts/r76_extract_manual_panel_xls.py \
      --xls /Users/edemhw/Honglab/CS/04_Claude/14_dcta/data/priors/final_unique_gene.xls \
      --out-dir results/neutrophil_atlas_rna_2026-10-19/gse253355_manual_panel_input_01
"""
import argparse
import csv
import datetime
import hashlib
import json
import os
import sys

import xlrd

SHEETS = [("final_unique_gene", 1), ("Negative markers", -1)]


def col_letter(j):
    s, j = "", j + 1
    while j:
        j, r = divmod(j - 1, 26)
        s = chr(65 + r) + s
    return s


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xls", required=True)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    if os.path.isdir(a.out_dir) and os.listdir(a.out_dir):
        sys.exit("결과 폴더가 비어 있지 않다(덮어쓰지 않는다): " + a.out_dir)
    os.makedirs(a.out_dir, exist_ok=True)
    book = xlrd.open_workbook(a.xls)
    rows, sheets = [], []
    for name, sign in SHEETS:
        if name not in book.sheet_names():
            sys.exit("sheet 가 없다: " + name)
        sh = book.sheet_by_name(name)
        header = [sh.cell_value(0, j) for j in range(sh.ncols)]
        n_cells = n_blank_text = n_nontext = 0
        for j in range(sh.ncols):
            for i in range(1, sh.nrows):
                c = sh.cell(i, j)
                if c.ctype == xlrd.XL_CELL_EMPTY:
                    continue
                if c.ctype == xlrd.XL_CELL_TEXT and c.value.strip() == "":
                    n_blank_text += 1
                    continue
                if c.ctype != xlrd.XL_CELL_TEXT:
                    n_nontext += 1
                raw = c.value if c.ctype == xlrd.XL_CELL_TEXT else repr(c.value)
                n_cells += 1
                rows.append({"sheet": name, "sign": sign, "excel_row": i + 1, "excel_col": j + 1, "col_letter": col_letter(j),
                             "header_raw": header[j], "symbol_raw": raw, "symbol": raw.strip(), "cell_type": c.ctype})
        sheets.append({"sheet": name, "sign": sign, "nrows": sh.nrows, "ncols": sh.ncols, "header": header,
                       "n_nonempty_cells": n_cells, "n_blank_text_cells": n_blank_text, "n_nontext_cells": n_nontext})
    out_csv = os.path.join(a.out_dir, "manual_panel_extracted.csv")
    fields = ["sheet", "sign", "excel_row", "excel_col", "col_letter", "header_raw", "symbol_raw", "symbol", "cell_type"]
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, quoting=csv.QUOTE_ALL)
        w.writeheader()
        w.writerows(rows)
    prov = {"source_xls": os.path.abspath(a.xls), "source_sha256": sha256(a.xls), "source_bytes": os.path.getsize(a.xls),
            "python": sys.executable, "python_version": sys.version.split()[0], "xlrd_version": xlrd.__version__,
            "command": [sys.executable, "-I"] + sys.argv, "written_at": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
            "sheets": sheets, "n_rows_extracted": len(rows), "outputs": {"manual_panel_extracted.csv": sha256(out_csv)},
            "note": "format extraction only; no name correction, deduplication or marker selection change"}
    with open(os.path.join(a.out_dir, "manual_panel_extract_provenance.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, ensure_ascii=False, indent=2)
    print("extracted", len(rows), "cells:", ", ".join("%s %d" % (s["sheet"], s["n_nonempty_cells"]) for s in sheets))


if __name__ == "__main__":
    main()
