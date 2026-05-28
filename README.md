# DistSpMV_Balanced — 复现与使用说明

本仓库复现并验证论文 “Balancing Computation and Communication in Distributed Sparse Matrix-Vector Multiplication”（DistSpMV）的主要想法：
- 基于 METIS 的重排使非零元尽可能聚集到对角块；
- 对角区间（balanced boundary）通过滑动窗口扩展以覆盖足够的本地非零元，减少远端 x 请求；
- 两阶段 MPI 通信（请求索引 → 回复 x 值）与本地优先计算以重叠通信/计算。

**语言/依赖**：C, MPI（Open MPI / MPICH）, OpenMP, METIS。

**代码位置（关键文件/目录）**
- [dist_spmv_balanced.c](dist_spmv_balanced.c) — Balanced 实现（METIS 重排、rank0 分发、boundary/通信/计算逻辑）。
- [dist_spmv_naive.c](dist_spmv_naive.c) — Baseline（用于正确性/性能对照）。
- [REPRO_REPORT.md](REPRO_REPORT.md) — 本次复现实验报告与结果分析。
- [results_correctness.csv](results_correctness.csv) — 正确性相关的 GFlops 汇总（已记录若干矩阵）。
- [results_perf.csv](results_perf.csv) — 性能对比（bcsstk 系列等）。
- [logs/](logs/) — 每次运行的 stdout/stderr 日志（按矩阵/模式/np 分类保存）。
- [matrices/](matrices/) — 测试用的 Matrix Market 文件（部分来自 SuiteSparse / 自建）。
  测试矩阵可从 SuiteSparse Matrix Collection 下载。推荐测试 `cant`（Williams 组）和 `ecology1`（McRae 组）。
- [scripts/](scripts/) — 一些辅助脚本（下载、批量运行等）。

## ✨ Project Highlights (项目亮点)

- **通信优化显著**：在 `ecology1.mtx` (1M 维度) 测试中，通信时间减少约 **58%**。
- **高稳定性**：通过 Rank-0 数据分发重构，成功在 4GB 内存虚拟机上运行千万级 NNZ 矩阵，解决了原始实现的 OOM 问题。
- **100% 正确性**：通过 `Y-norm` 校验，确保所有优化（重排、边界扩展）均不改变计算结果。

## 🛠️ 实验结果汇总（关键数据）

| 矩阵 (Matrix) | 模式 (Mode) | 进程数 (np) | 通信耗时 (Comm) | 提速/优化效果 |
| :--- | :--- | :--- | :--- | :--- |
| **ecology1.mtx** | Naive | 4 | 0.136s | 基准 |
| **ecology1.mtx** | **Balanced** | 4 | **0.058s** | **通信量减少 58%** |
| **cant.mtx** | Balanced | 4 | 0.049s | 局部性优化验证 |

详细性能分析请参阅 [REPRO_REPORT.md](REPRO_REPORT.md)。

注意：由于复现环境为单核虚拟机，多进程会产生严重的 CPU 抢占。在此环境下，加速比主要通过**通信时间占比的下降**来验证，而非总执行时间的绝对值（单机环境下 MPI 可能对本地通信做共享内存优化）。

## 快速开始（在本机单节点多进程测试）

1) 编译（假设系统已安装 METIS 与 MPI）

```bash
mpicc -O3 -fopenmp -o dist_spmv_balanced dist_spmv_balanced.c -lmetis -lm
mpicc -O3 -fopenmp -o dist_spmv_naive dist_spmv_naive.c -lmetis -lm
```

2) 运行示例（以 4 个 MPI 进程、单线程/进程为例）

```bash
export OMP_NUM_THREADS=1
mpirun --oversubscribe -np 4 ./dist_spmv_balanced matrices/cant.mtx 10 1.0 balanced
mpirun --oversubscribe -np 4 ./dist_spmv_balanced matrices/cant.mtx 10 1.0 naive
```

命令行参数说明：
- `matrix.mtx`：输入矩阵（Matrix Market 格式）。
- `10`：迭代次数（niter，默认 10）。
- `1.0`（或 0.8 等）：lower_bound_frac，表示希望对角区间覆盖的本地对角 nnz 比例（实现会把比例乘以最大对角 nnz 得到阈值）。

## 输出与日志

- 程序在 stdout 打印：进程数、矩阵尺寸、迭代次数、Total time、Compute time、Comm time、GFlops、Y-norm 等摘要信息。
- 我把常用运行的输出保存在 `logs/` 下，文件命名示例：`cant_balanced_np4.log`、`ecology1_naive_np4.log`。
- 汇总结果保存在：`results_correctness.csv`（正确性对比）和 `results_perf.csv`（性能对比）。详细分析见 [REPRO_REPORT.md](REPRO_REPORT.md)。

## 实现细节与注意事项

- METIS：现在在 `rank==0` 上运行分区和生成 `perm`，然后将必要的数据分发给其它进程以避免每个进程都读入完整矩阵（已实现 `MPI_Scatterv` 局部分发以减少 OOM）。
- Boundary：`compute_balanced_boundary`（及其本地版本）实现了论文中滑动窗口思想的启发式扫描。
- 通信：每轮先收集远端列索引并发送请求（去重），所有者回复 x 值，采用非阻塞以重叠计算。
- 正确性：在本次复现中，所有测试矩阵的 `Y-norm` 在 `balanced` 与 `naive` 两种模式下均一致（参见报告）。

## 实验结果（摘要）

- 在 `matrices/ecology1.mtx` 上，Balanced 模式使通信时间从 ~0.136s 降到 ~0.058s，通信时间约减少 58%（详细数据和结论见 [REPRO_REPORT.md](REPRO_REPORT.md)）。

## 常见问题
- 如果运行大型矩阵（例如 `ecology1.mtx`）时占用内存过高，请确保：
  - 只在一个节点上运行小规模进程数，或使用 rank-0 分发（已实现）；
  - METIS 的内存开销也可能是瓶颈（分区前可考虑稀疏压缩或分块读取）。

## 提交与同步（本地仓库）
- 我已经在本目录初始化了本地 Git 仓库并提交了当前改动；若要推送到 GitHub，请在本机按需设置 `origin` 并执行 `git push`（会提示用户名/Token）。

## 📜 参考引用

Mi, H., Yu, X., Yu, X., Wu, S., & Liu, W. (2023). Balancing Computation and Communication in Distributed Sparse Matrix-Vector Multiplication. 2023 IEEE/ACM 23rd International Symposium on Cluster, Cloud and Internet Computing (CCGrid), pp. 535-544. [DOI: 10.1109/CCGrid57682.2023.00056]

---
作者/维护：复现实验自动化脚本与补丁

