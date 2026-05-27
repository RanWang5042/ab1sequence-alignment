# AB1 测序文件比对分析

## 功能概述

批量读取 Sanger 测序 AB1 文件，**先将同一克隆的正向(F)和反向(R)读序合并为共识序列，再与参考序列进行比对**，最终输出详细的 Excel 报告。

适用于克隆鉴定的全流程分析。

## 核心工作流

```
AB1文件对 (F + R)
    ↓
1. 质量过滤 (Phred ≥ 15)
    ↓
2. 引物序列 trimming (pScbb23-F / pScbb46-F)
    ↓
3. R读序反向互补 (5'→3')
    ↓
4. F + R_RC 拼接 → 共识序列
    ↓
5. 与参考序列局部比对 (smith-waterman style)
    ↓
Excel 报告
```

## 使用方法

### 标准分析命令

```bash
python3 analyze_ab1.py \
  --ab1-dir ./ab1_files \
  --ref-xlsx ./reference.xlsx \
  --ref-sheet "Part" \
  --primer-f "GTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAATAAATGGTCTTC" \
  --output report.xlsx
```

### 参数说明

| 参数 | 说明 | 示例 |
|------|------|------|
| `--ab1-dir` | AB1 文件所在目录 | `./ab1_files` |
| `--ref-xlsx` | 参考序列 Excel 文件 | `./parts.xlsx` |
| `--ref-sheet` | Excel 工作表名称 | `Part` |
| `--ref-col-name` | 参考序列名称列(默认A) | `Part Name` |
| `--ref-col-seq` | 参考序列列(默认D) | `Sequence` |
| `--primer-f` | 正向引物序列(可选) | `GTTACTGCT...` |
| `--min-identity` | 最小 identity 阈值(默认 0.7) | `0.85` |
| `--output` | 输出报告路径 | `report.xlsx` |

### 从 Excel 读取参考序列

Excel 文件应包含：
- **Part Name 列**：克隆名称（如 PRO16、TER84）
- **Sequence 列**：对应的 DNA 序列

脚本会自动解析 T→TER、P→PRO 的命名约定。

## 输出报告格式

| 列名 | 说明 |
|------|------|
| 克隆ID | 从文件名提取的克隆编号 |
| 目标Part | Excel 中标注的目标序列 |
| 组装方式 | F_ONLY / CONCAT / OVERLAP |
| F清洗长度 | Forward 读序清洗引物后长度 |
| R_RC清洗长度 | Reverse 读序反向互补清洗后长度 |
| 共识序列长度 | F + R_RC 拼接后总长度 |
| Identity(%) | 最佳匹配 identity |
| 覆盖度(%) | 匹配区域占参考序列的比例 |
| 最佳匹配 | 比对到的 Part 名称 |
| 状态 | ✅ 正确 / ⚠️ 低 / ❌ 不匹配 |

## 判定标准

| Identity | 覆盖度 | 状态 |
|----------|--------|------|
| ≥ 95% | ≥ 70% | ✅ 正确 |
| 80-94% | ≥ 50% | ⚠️ 部分匹配 |
| < 80% | — | ❌ 不匹配 |

## 依赖

- Python 3.7+
- Biopython (`pip install biopython`)
- openpyxl (`pip install openpyxl`)

## 参考资料

- AB1 文件格式详见 [references/ab1_format.md](references/ab1_format.md)
- 引物序列详见 [references/primer_sequences.md](references/primer_sequences.md)
