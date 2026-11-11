# 验证说明

截至 2026-11-11，本机环境为 Ubuntu 24.04 / WSL2 / WSLg，GUI 使用 Python 3.12.13、Tk 8.6、NumPy 1.26.4 和 openpyxl 3.1.5。模型使用独立运行环境。

| 验证项 | 观察结果 |
| --- | --- |
| 发布包回归 | 全新源码快照执行 73 项，71 项通过，2 项需要完整模型包而跳过，约 91.8 秒 |
| 真核联合流程 | 3 条公开酵母蛋白，九项各 3/3 完成 |
| 原核联合流程 | 50 条合成工程输入，四个适用工具各 50/50 完成 |
| PTS1 数值一致性 | 1,014 个记录预测无差异，最大决策值误差约 `1.25e-13` |
| Nextflow 续跑缓存 | 真实字节统计流程依次观察到 COMPLETED、CACHED、修改输入后 COMPLETED |
| 同名轨迹续跑 | 从 `_metaloc` 入口验证缓存命中与上一轮轨迹归档 |
| Phobius 短序列 | 官方程序对 5 aa 与 50 aa 各一条输入返回完整输出 |

发布包验证复用已安装的 GUI 基础环境，由安装脚本重新获取固定版本的 metaLoc 与 PredGPI 源码，并运行 PTS1 CLI。Phobius、TargetP 的真实发行包用例在此次源码快照中按依赖条件跳过。

回归覆盖 FASTA、重复 ID、输入边界、结果对齐、拓扑坐标、缺失输出、模型路由、缓存指纹、取消、孤儿进程回收、跨目录队列、历史任务、筛选、排序和导出。取消后的模型输出恢复使用受控工程样例验证。

真核示例为 UniProt `P00330`、`P16474` 和 `P19882`。原核工程输入为固定随机种子生成的 64–65 aa 序列，用于检查批量执行与 ID 对齐。模型预测精度需要独立标注数据与训练同源性控制后评价。

PTS1 还完成原作者四个单元测试和内部五折复核，数据和模型解释见 [PTS1 说明](../tools/pts1/README.md)。

```bash
bash scripts/python.sh -m unittest discover -s tests -v
```

基础回归复用 `tests/fixture_tool.py` 的受控输出，真实 Phobius、TargetP 和 PredGPI 用例在依赖可用时运行。完整模型流程可使用 `scripts/validate_real.py` 与 `scripts/validate_prok.py` 小规模验收。
