# 模型安装说明

先运行 `bash scripts/setup.sh` 准备 GUI、固定的 metaLoc 与 PredGPI 源码。模型环境位于 Linux 用户目录，项目使用 `.runtime/` 中的链接访问。安装需要网络与足够磁盘空间，ESM1b 文件约 7.8 GB，ESM-2 文件约 2.6 GB。

## 获取模型包

按工具官方条款获取以下文件，放入 `tools/assets/`，目录可自行创建。

| 文件 | 来源 |
| --- | --- |
| `signalp-6.0i.fast.tar.gz` | [SignalP 6.0](https://services.healthtech.dtu.dk/services/SignalP-6.0/)，fast 发行包 |
| `deeploc-2.1.All.tar.gz` | [DeepLoc 2.1](https://services.healthtech.dtu.dk/services/DeepLoc-2.1/) |
| `DeepTMHMM-Academic-License-v1.0.zip` | [DeepTMHMM](https://dtu.biolib.com/DeepTMHMM/)，本地学术发行包 |
| `targetp-2.0.Linux.tar.gz` | [TargetP 作者仓库](https://github.com/JJAlmagro/TargetP-2.0) |
| `phobius101_linux.tgz` | [Phobius 官方学术下载](https://phobius.sbc.su.se/download) |

TargetP、Phobius 和相关公开源码可通过现有脚本获取。下载前确认各工具的使用条款。

```bash
bash scripts/python.sh scripts/fetch_github.py
bash scripts/python.sh scripts/setup_support.py
```

第二条命令准备发行包中的 TargetP、Phobius 程序和运行链接。下载脚本核对固定包的 SHA256，记录在本地 `validation/` 中。包名、大小与校验值见 [upstream.lock.json](upstream.lock.json)，其他版本需要调整配置并验证适配。

## SignalP、DeepLoc 与 DeepTMHMM

准备前述三个授权发行包后依次运行。

```bash
bash scripts/python.sh scripts/install_core.py
bash scripts/python.sh scripts/install_deeptmhmm.py
```

核心环境使用 Python 3.10、PyTorch 1.13.1 / CUDA 11.7、NumPy 1.26.4、fair-esm 2.0.0 和 transformers 4.30.2。约束位于 `tools/core-environment.yml`、`tools/core-pip-requirements.txt` 和 `tools/core-constraints.txt`。DeepTMHMM 使用独立的旧版 NumPy、Biopython、fair-esm 环境。

DeepLoc Fast 另需官方 ESM1b 参数，可在首次推理时由其依赖获取，也可预先缓存。

```bash
mkdir -p .runtime/model-cache/torch/hub/checkpoints
curl --fail --location --retry 3 \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt \
  --output .runtime/model-cache/torch/hub/checkpoints/esm1b_t33_650M_UR50S.pt
```

GPU 核心环境面向 NVIDIA CUDA；DeepLocPro 与 TargetP 显式使用 CPU。其他设备可在高级配置中使用已安装并核验的本地工具命令。

## DeepLocPro

固定源码 commit `32c0b37c7a31a0686ce022a697f48998fc2997b1`，使用原始分类器及 ESM-2 参数。将以下文件下载至 `tools/assets/`，再运行安装器。

```bash
mkdir -p tools/assets
curl --fail --location --retry 3 \
  https://codeload.github.com/Jaimomar99/deeplocpro/zip/32c0b37c7a31a0686ce022a697f48998fc2997b1 \
  --output tools/assets/deeplocpro-32c0b37.zip
curl --fail --location --retry 3 \
  https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt \
  --output tools/assets/esm2_t33_650M_UR50D.pt
curl --fail --location --retry 3 \
  https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt \
  --output tools/assets/esm2_t33_650M_UR50D-contact-regression.pt
bash scripts/python.sh scripts/install_deeplocpro.py
```

安装器检查 SHA256。任务按类群选择 `negative`、`positive`、`archaea` 或 `any`，保存调用参数。当前范围为 10–6,000 aa，已有一条 1,077 aa 全长输入的本机验证。[官方服务说明](https://services.healthtech.dtu.dk/services/DeepLocPro-1.0/)

## PTS1 复现模型

`tools/pts1/` 包含复现的 `svm_weights.npz`、NumPy 预测器、导出元数据和 BSD-3-Clause 声明，日常推理使用 GUI 环境中的 NumPy。

```bash
bash scripts/python.sh tools/pts1/predict_pts1_numpy.py \
  --fasta examples/signals.fasta --out predictions.csv
```

模型按 oKoch 固定源码和数据重建，分析 C 端最后 14 aa。来源、数值核验和内部交叉验证见 [PTS1 说明](tools/pts1/README.md)。

## 检查与首次运行

```bash
bash scripts/python.sh pipeline.py doctor
bash run.sh
```

环境检查核对文件与命令。可使用 `examples/signals.fasta` 检查少量输入，真实预测状态以各模块的完整性检查、输出和日志为准。
