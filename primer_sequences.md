# 常用测序引物序列参考

## pScbb23 载体引物

### pScbb23-F (正向测序引物)
```
GTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAATAAATGGTCTTC
```
- 长度: 61 bp
- 用途: PRO 克隆正向测序 (pScbb23 载体)
- 位置: 插入片段 5' 端上游

### pScbb23-R (反向测序引物)
```
位置待确认
```
- 用途: PRO 克隆反向测序
- 注意: 测序结果需要**反向互补**后才能与目标序列比对

## pScbb46 载体引物

### pScbb46-F (正向测序引物)
```
GTTACTGCTGCTGGTATTACCCATGGTATGGATGAATTGTACAAATAATAAATGGTCTTC
```
- 长度: 61 bp
- 用途: TER 克隆正向测序 (pScbb46 载体)
- 序列与 pScbb23-F **完全相同**

### pScbb46-R (反向测序引物)
```
位置待确认
```
- 用途: TER 克隆反向测序
- 注意: 测序结果需要**反向互补**后才能与目标序列比对

## 通用测序引物

### M13 Forward (-47)
```
CGCCAGGGTTTTCCCAGTCACGAC
```
- 长度: 23 bp
- 用途: 通用载体引物

### M13 Reverse (-48)
```
CTCACATTGATTGTTTGAGAGGG
```
- 长度: 23 bp
- 用途: 通用载体引物

### T7 Promoter
```
TAATACGACTCACTATAGGG
```
- 长度: 20 bp
- 用途: 表达载体测序

## 引物处理注意事项

### 1. 序列方向
```
原始序列 (5'→3'): ATGCGATCGATCG
                      ↓ SeqIO.read()
序列方向 (5'→3'): ATGCGATCGATCG  ← AB1 文件中序列的方向

反向互补 (3'→5' → 5'→3'): TAGCTAGCATGCA
```

### 2. R 引物处理流程
```python
# Step 1: 读取 R 向 AB1 文件
r_seq = str(record.seq)  # 方向: 3'→5'

# Step 2: 序列反向互补
r_rc = reverse_complement(r_seq)  # 方向: 5'→3'

# Step 3: Trim 引物
r_rc_clean = trim_primer(r_rc, primer_seq)
```

### 3. F + R 合并策略
```python
# 方案1: 简单拼接 (推荐用于短插入片段)
consensus = f_clean + r_rc_clean

# 方案2: Overlap 组装 (适用于长插入片段)
# 在重叠区域取较高质量的碱基
consensus = overlap_assemble(f_clean, r_rc_clean, min_overlap=20)
```

## 常见问题

| 问题 | 原因 | 解决方案 |
|------|------|----------|
| R 序列比对失败 | 未做反向互补 | 确认使用 `reverse_complement()` |
| 序列开头无引物 | 引物太短或测序失败 | 检查 quality score |
| F 和 R 拼接后重复 | 两端引物结合太近 | 使用 overlap 组装或只保留 F |
| 目标序列在载体附近 | 引物位置不理想 | 考虑更换测序引物 |
