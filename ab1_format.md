# AB1 文件格式参考

## 文件结构

AB1 (Applied Biosystem Chromatogram) 是 Sanger 测序仪输出的二进制格式，包含：

1. **序列数据** (nucleotide sequence)
2. **质量分数** (Phred quality scores)
3. **色谱数据** (trace data: A/T/G/C 荧光强度)
4. **测序信息** (引物、样品名、运行参数等)

## Biopython 读取示例

```python
from Bio import SeqIO

# 读取 AB1 文件
record = SeqIO.read("sample.ab1", "abi")

# 提取序列
sequence = str(record.seq)
print(f"序列长度: {len(sequence)} bp")
print(f"序列前50bp: {sequence[:50]}")

# 提取质量分数 (Phred scale)
qual_scores = record.letter_annotations['phred_quality']
print(f"质量分数数量: {len(qual_scores)}")
print(f"平均质量: {sum(qual_scores)/len(qual_scores):.1f}")
print(f"最低质量: {min(qual_scores)}")
print(f"最高质量: {max(qual_scores)}")

# 提取样本名
print(f"样本名: {record.name}")
print(f"描述: {record.description}")
```

## Phred 质量分数

Phred quality score Q 定义为：
```
Q = -10 * log10(P)
```
其中 P 是碱基调用错误的概率。

| Q 值 | 错误率 | 准确率 |
|------|--------|--------|
| 10 | 10% | 90% |
| 20 | 1% | 99% |
| 30 | 0.1% | 99.9% |
| 40 | 0.01% | 99.99% |

**推荐阈值**：
- 高质量序列：Q ≥ 20
- 清洗低端：Q ≥ 15

## 质量过滤

```python
def trim_by_quality(sequence, qual_scores, min_q=15):
    """从两端去除低质量碱基"""
    seq = list(sequence)
    qual = list(qual_scores)
    
    # 去除5'端低质量
    while qual and qual[0] < min_q:
        qual.pop(0); seq.pop(0)
    
    # 去除3'端低质量
    while qual and qual[-1] < min_q:
        qual.pop(); seq.pop()
    
    return ''.join(seq), qual
```

## 引物序列识别

测序文件中的引物序列通常位于序列开头（5'端），需要 trimming：

```
[引物序列][目标插入片段][poly-A尾巴]
```

常见引物：
- M13 Forward: `TGTAAAACGACGGCCAGT`
- M13 Reverse: `CAGGAAACAGCTATGACC`
- pUC/M13F: `CCCAGTCACGACGTTGTAAAACG`

 trimming 方法：
1. 精确匹配（完全匹配引物开头）
2. 模糊匹配（允许少量错配，通常 ≤ 2 个）
