# DistSpMV_Balanced — 复现与使用说明（更新版）

本仓库复现并实现论文 “Balancing Computation and Communication in Distributed Sparse Matrix-Vector Multiplication” (CCGrid 2023) 的核心思想。

主要目标
- 实现论文提出的对角区间扩展与远端重划策略以减少通信开销；
- 提供基线（naive）、METIS+naive 与 Balanced 三种实现，便于比较；
- 提供自动化实验脚本与绘图脚本，便于复现实验与生成图表。

依赖
- C 编译器：gcc (>= 9.4.0 已测试)
- MPI：OpenMPI 或 MPICH
- OpenMP
- METIS（用于 METIS+naive 与 Balanced 的预处理）
- Python3 与 matplotlib（用于绘图）

快速安装（Ubuntu 示例）
```bash
sudo apt update
sudo apt install -y build-essential libopenmpi-dev openmpi-bin libmetis-dev python3-pip
pip3 install matplotlib
```

文件说明
- `dist_spmv_naive.c`        — Baseline Allgatherv 实现（已输出 Y-norm 与 Per-process nnz）
- `dist_spmv_metis_naive.c`  — METIS 重排 + Baseline（已输出 Y-norm 与 Per-process nnz）
- `dist_spmv_balanced.c`    — Balanced 实现（预处理索引交换 + 单次值交换 + 单遍计算），输出 local/remote nnz 与 Y-norm
- `run_experiments.sh`      — 一键编译、运行实验并将结果写入 CSV，日志写入 `logs/`
- `plot_results.py`         — 读取 CSV 并生成 `figures/` 中的 PNG 图表
- `REPRO_REPORT.md`         — 复现报告骨架（含实验占位）
- `matrices/`               — 测试矩阵（Matrix Market 格式）

运行实验
1. 确保已安装依赖并位于仓库根目录。
2. 赋予脚本执行权限并运行：

```bash
chmod +x run_experiments.sh
./run_experiments.sh
```

脚本说明（简要）
- `run_experiments.sh` 会：
  - 编译三份可执行文件（`dist_spmv_naive`、`dist_spmv_metis_naive`、`dist_spmv_balanced`）
  - 将 `OMP_NUM_THREADS=4` 并在 np=1,2 下对 `cant.mtx`、`ecology1.mtx`、`bcsstk30.mtx` 运行主实验
  - 对 `cant.mtx` 在多个 `lower_bound_frac` 下运行消融实验（np=2）
  - 运行负载均衡度测量（np=2）并计算 LIR
  - 所有程序 stdout/stderr 会写到 `logs/`（文件名按矩阵_算法_np.log）
  - 结果保存为：`results_main.csv`、`results_ablation.csv`、`results_balance.csv`

3. 运行后生成图表：

```bash
python3 plot_results.py
# 生成的图表位于 figures/ 目录
```

输出说明（CSV 字段）
- `results_main.csv`：Matrix,Algorithm,np,GFlops,TotalTime,ComputeTime,CommTime,Ynorm
- `results_ablation.csv`：frac,np,GFlops,TotalTime,ComputeTime,CommTime,LocalNNZ_p0,LocalNNZ_p1,RemoteNNZ_p0,RemoteNNZ_p1
- `results_balance.csv`：Matrix,Algorithm,np,NNZ_p0,NNZ_p1,LIR

后续建议
- 在真实集群（多节点、真实网络）上重跑实验以验证大规模的可伸缩性。
- 将 rank-0 的分发替换为更细粒度的分发以节省内存并提高可扩展性。

问题与贡献
- 如需对实验集、参数或画图方式做调整，我可以帮助修改脚本或添加新的指标。

---
作者/维护：复现实验自动化脚本与补丁
