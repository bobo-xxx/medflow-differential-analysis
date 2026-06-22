#!/usr/bin/env python3
import sys
import pandas as pd
from pathlib import Path

in_map = sys.argv[1]
incorrect_proportion = sys.argv[2].upper()
out_file = sys.argv[3]

group = pd.read_csv(in_map)
counts = group.iloc[:, 1].value_counts()
c_num, t_num = counts.iloc[0], counts.iloc[1]
ratio = max(c_num, t_num) / min(c_num, t_num)

if ratio > 10 and incorrect_proportion != "TRUE":
    g1, g2 = counts.index[0], counts.index[1]
    raise ValueError(
        f"\n样本比例异常，{g1}（{c_num} 样本）、{g2}（{t_num} 样本），"
        f"样本比例 {ratio:.1f}:1 超过 10:1，不建议进行差异分析。\n"
        "若要进行差异分析请将 Incorrect proportion 改为 TRUE。"
    )
Path(out_file).touch()