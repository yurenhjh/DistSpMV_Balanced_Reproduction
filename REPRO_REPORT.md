# 论文复现报告 — "Balancing Computation and Communication in Distributed SpMV"（CCGrid 2023）

## 1. 论文技术分析

### 1.1 问题背景
- SpMV 在分布式系统上的性能瓶颈主要来自：跨进程通信开销与计算负载不均衡。
- 传统图划分（如 METIS）可改善通信局部性，但不能直接改变每行的非零元分布，因而对行级负载不均的改善有限。

### 1.2 两个核心优化策略
- 策略一（对角区间扩展，Algorithm 1）：基于每个对角块的非零元统计，按列单调扩展对角区间，直至覆盖足够多的本地非零元，从而减少对远端 x 值的请求。
- 策略二（远端行重划 + 计算/通信重叠）：在扩展对角区间后，对远端行做重划以均衡各进程的计算量；通信采用“请求索引 → 回复 x 值”的两阶段流程，并尽量重叠本地计算与远端值传输。

### 1.3 四步算法流程
1. 使用 METIS（rank 0）生成初始重排（perm）；
2. 对角区间扩展 + 基于扩展结果的行重划（Algorithm 1 与后续重划）；
3. 两阶段 MPI 通信（请求索引 → 回复 x 值），并在等待期间优先计算本地对角区块；
4. 汇总结果并做正确性校验（Y-norm）。

### 1.4 与现有方法的核心区别
- 论文不止依赖图划分，还主动对行的非零元进行局部重组/调度以减少远端访问量，同时兼顾线程级负载均衡。

## 2. 实现说明

### 2.1 代码与论文的对应关系
- Algorithm 1（boundary scan）  : `compute_balanced_boundary` / `compute_balanced_boundary_local` （`dist_spmv_balanced.c`）
- Algorithm 2/3（两阶段通信）    : Phase A/B/C/D/E 主循环实现（`dist_spmv_balanced.c`）
- Algorithm 4（线程级负载均衡） : `compute_thread_boundaries`（`dist_spmv_balanced.c`）
- 基线实现                        : `dist_spmv_naive.c`（Allgatherv）
- METIS+基线                       : `dist_spmv_metis_naive.c`

### 2.2 已知差异
- 为简化实现，当前仓库采用所有进程读取完整矩阵并在 rank 0 使用 `MPI_Scatterv` 将行分发给各进程的混合策略（见代码注释）；真实集群中应由 rank 0 精细分发以节省内存。
- 测试向量 `x` 在本次复现中固定为全 1（与论文测试设置略有差别）。
- 受限于本地虚拟机资源（2 逻辑 CPU，4GB 内存），实验规模远小于论文（论文使用高核数/高速网络）。

## 3. 修复记录
（参见代码注释；此处列出关键修复）
- 修复 `lower_bound` 的语义：命令行参数为比例系数 `lower_bound_frac`，阈值由 `max_diag_nnz * frac` 得到。
- 避免未初始化访问：`diag_nnz` 在所有进程上 `calloc` 并由 rank 0 填充后广播。
- 释放全矩阵以节省内存：rank 0 在分发完成后释放全矩阵 `A`。
- 初始化扩展区间的 `x_buf`：在迭代前做一次 `MPI_Allgatherv` 保证扩展列的 x 值正确。
- 预分配并复用 `send_lists` 减少每轮 malloc/realloc 开销。
- 合并 Phase C 与 Phase E 为单遍计算，避免对 nnz 的二次遍历与分支开销。

## 4. 实验设计与脚本

- 已实现脚本：`run_experiments.sh`（编译、运行、解析日志并输出 CSV）和 `plot_results.py`（生成 `figures/` 中的 PNG）。
- 主要 CSV：`results_main.csv`（主对比）、`results_ablation.csv`（阈值消融）、`results_balance.csv`（负载均衡度）。
- 日志目录：`logs/`（每次运行的 stdout/stderr 保存）

## 5. 实验结果占位（TBD）
### 5.1 正确性验证
矩阵 | np | Naive Y-norm | METIS Y-norm | Balanced Y-norm | 是否一致
--- | --- | --- | --- | --- | ---
cant.mtx | 2 | 4.170292046439e+04 | 4.170292046439e+04 | 4.170292046439e+04 | 是
ecology1.mtx | 2 | 1.309715002014e-13 | 1.309715002014e-13 | 1.309715002014e-13 | 是
bcsstk30.mtx | 2 | 1.317123426259e+04 | 1.317123426259e+04 | 1.317123426259e+04 | 是

### 5.2 主实验性能对比（TBD，插入图）
- 插图：`fig1_gflops_comparison.png`, `fig2_comm_time.png`, `fig3_speedup.png`

### 5.3 阈值消融（TBD，插图）
- 插图：`fig4_ablation.png`

### 5.4 负载均衡度（TBD，插图）
- 插图：`fig5_load_balance.png`, `fig6_remote_nnz.png`

## 6. 与论文结果对比分析（占位）
- 论文在高并发/高速网络下得到较高加速比，本复现受限环境下重点验证通信量下降趋势与负载均衡效果。

## 7. 使用说明（快速）
1. 编译并运行全部实验：

```bash
chmod +x run_experiments.sh
./run_experiments.sh
```

2. 生成图表：

```bash
python3 plot_results.py
```

3. 结果文件：`results_main.csv`, `results_ablation.csv`, `results_balance.csv`。

## 8. 结论与后续工作
- 本次复现实现了论文的核心思想并补充了正确性与负载统计输出，能够在本地环境中重现通信量下降的趋势。
- 后续可在多节点集群（真实网络）上运行更大规模实验，并将 rank-0 分发进一步替换为按需分发以提升可扩展性。

---

*注：本报告内容为当前代码状态的说明与占位，数值表格将在用户要求重新运行全部实验后填充。*

**自动分析（脚本生成）**

- Ecology1 通信缩减比（CommTime_Metis / CommTime_Balanced，np=2）: 1.259573

- Cant 性能下降说明（np=2）: 在低核心数（np=2）下，复杂负载均衡逻辑带来的软开销超过了微量通信缩减的收益。
