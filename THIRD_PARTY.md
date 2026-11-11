# 第三方来源与引用

ProteinLoc 自有代码采用 MIT。第三方代码、数据、发行包和模型各自保留原授权，安装时需遵循官方条款。

| 项目 | 用途与来源 | 授权说明 |
| --- | --- | --- |
| [metaLoc](https://github.com/scottc-bio/metaLoc) | 核心工作流，commit `c10d8035603134a0c9916dcbcebdcf51d0559e95` | MIT |
| [PredGPI](https://github.com/BolognaBiocomp/predgpi) | GPI 预测，commit `35ac13788771aaf06b7174852476baea9321ae36` | GPL-3.0-or-later，保留上游 LICENSE |
| [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/) | 本地 fast 发行包 | 按官方发行与使用条款获取 |
| [Phobius](https://phobius.sbc.su.se/download) | 1.01 本地学术发行包 | 按作者学术分发条款使用 |
| [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) | 真核定位模型 | 按官方发行与使用条款获取 |
| [DeepLocPro](https://github.com/Jaimomar99/deeplocpro) | 原核定位，commit `32c0b37c7a31a0686ce022a697f48998fc2997b1` | CC BY-NC-SA 4.0 |
| [TargetP 2.0](https://github.com/JJAlmagro/TargetP-2.0) | 作者的 Linux 发行包 | 按上游发行包条款使用 |
| [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM/) | 本地 1.0 Academic 发行包 | 按随包学术授权使用 |
| [oKoch PTS1 SVM](https://github.com/oKoch/pts1-prediction-tool) | 复现 NPZ 与 NumPy 适配器 | 保留 [BSD-3-Clause](tools/pts1/LICENSE)，引用时说明复现来源 |
| [ESM](https://github.com/facebookresearch/esm) | 基础蛋白语言模型 | 按原代码、参数及模型文档条款获取 |
| [UniProt](https://www.uniprot.org/) | 酵母参考示例与注释 | [CC BY 4.0](https://www.uniprot.org/help/license) |

论文或报告应引用实际运行的工具及对应论文。PredGPI 对应 Pierleoni et al., BMC Bioinformatics 2008，其余论文入口见各工具官方页面。实际版本、配置、输入变更和分数写入每个任务结果。

ER / PTS 候选规则的参考来源随结果保存，包括 [ER 回收信号](https://pmc.ncbi.nlm.nih.gov/articles/PMC2120172/)、[ER 双赖氨酸信号](https://pubmed.ncbi.nlm.nih.gov/23481256/)、[PTS1 基序](https://prosite.expasy.org/PDOC00299) 和 [PTS2 规则](https://pubmed.ncbi.nlm.nih.gov/8125978/)。
