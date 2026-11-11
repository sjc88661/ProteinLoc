# PTS1 SVM 复现权重

本目录提供按 [oKoch/pts1-prediction-tool](https://github.com/oKoch/pts1-prediction-tool/tree/0c1041905d12b96b72063aec5dccfd6073f31435) 固定代码和数据重新拟合、导出的 SVM 参数，采用 BSD-3-Clause。上游构造预测器时使用随附 FASTA 训练 SVC，本项目使用已导出的 NPZ 进行 NumPy 推理。

| 文件 | 用途 |
| --- | --- |
| `svm_weights.npz` | 22,543 字节的复现参数 |
| `predict_pts1_numpy.py` | NumPy 预测器与 CLI |
| `model_meta.json` | 原导出环境、训练配方与元数据 |
| `LICENSE` | 原始 BSD-3-Clause 声明 |

权重 SHA256 为 `e699cf76224d06adb7c259cc9640c011e2315b25b82efa005d0248d7fa252327`。

## 核验结果

独立环境使用 scikit-learn 1.6.1，对固定源码重新训练。有效样本为 505 正样本和 11,269 负样本，C14 one-hot 为 308 维，得到 1,308 个支持向量，与 NPZ 重建的向量一致。系数最大差约 `1.92e-14`。

对 1,014 个原始训练记录进行数值对照，预测无差异，最大决策值误差约 `1.25e-13`。这些记录来自训练文件，数值一致性用于验证导出推理。

内部五折复核得到 Sensitivity `0.833663`、Specificity `0.998491`、Precision `0.961187`。将相同 C14 分组后的结果为 `0.841584`、`0.998491`、`0.961538`。独立泛化能力还需控制同源序列、训练筛选偏差和物种分布影响后评价。

适配器按上游类别编码处理特殊残基，未知类别形成全零块，保留 U / X。CSV 导出 17 位有效数字。应用分析真实 C 端 14 aa，保留完整输入和分析坐标，短序列与未建模的 C 端残基逐条报告。

## 使用

```bash
bash scripts/python.sh tools/pts1/predict_pts1_numpy.py \
  --sequence MMMMKLSKMLLLSLSKLSKLSKLSKL
bash scripts/python.sh tools/pts1/predict_pts1_numpy.py \
  --fasta proteins.fasta --out predictions.csv
```

输出包含候选分类与原始 `decision_function`。它是有符号决策值，正值表示该 SVM 的 PTS1 候选。科研引用应注明复现来源和内部验证含义。

本模型复现 oKoch 1.0.2 的固定分类器，ProteinLoc 的 PTS1 / PTS2 基序另行报告。
