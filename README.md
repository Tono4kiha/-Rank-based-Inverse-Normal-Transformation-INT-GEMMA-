# pheno_int.R

基于排秩的逆正态变换（Rank-based Inverse Normal Transformation, INT），将表型数据转换为 GEMMA 可直接使用的表型文件。

采用 **Blom transformation**：

```
z = qnorm((r - 3/8) / (N - 2×3/8 + 1))
```

## 功能

- 自动判断输入文件是否有表头
- 支持指定要转换的表型列
- 自动识别缺失值（`-9`, `NA`, `NaN`, `.` 等）
- 表型质量检查（样本数、缺失数、唯一值数）并给出警告
- 每个表型输出一个 `pheno_<trait>_GEMMA.txt`（`FID IID Pheno` 三列，无表头）
- 可选生成 INT 前后分布检查图

## 依赖

- R（无第三方包，仅用 base R）

## 用法

```bash
# 转换全部表型
Rscript pheno_int.R pheno.txt

# 只转换指定表型
Rscript pheno_int.R pheno.txt trait1 trait2

# 转换全部表型并生成检查图
Rscript pheno_int.R pheno.txt yes

# 指定表型 + 检查图
Rscript pheno_int.R pheno.txt trait1 trait2 yes
```

## 参数

| 参数 | 说明 |
|------|------|
| `<文件>` | 输入表型文件，第一列为 ID，其余为表型（空格或制表符分隔） |
| `trait1 trait2 ...` | 可选，指定要转换的表型列名（有表头时按列名匹配；无表头时按顺序命名，数量须与表型列数一致） |
| `yes` / `no` | 可选，是否输出检查图 `INT_check.jpg`（默认 `no`） |

## 输出

- `pheno_<trait>_GEMMA.txt`：GEMMA 表型文件（无表头，`NA` 表示缺失）
- `INT_check.jpg`：INT 后直方图 + QQ 图（`plot_flag = yes` 时生成）

## 示例

输入 `pheno.txt`：

```
id  height  weight
1   170     65
2   165     NA
3   180     80
```

运行：

```bash
Rscript pheno_int.R pheno.txt height weight
```

输出 `pheno_height_GEMMA.txt`、`pheno_weight_GEMMA.txt`。

Note:
```
该脚本只是使用了已有的成熟算法，在deepseek和chatgpt共同帮助下编写的快捷工具。INT转换只能在一定程度上减少假阳性，但是不能扭转本身表型记录非随机抽样所造成的偏差
```
