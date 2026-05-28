# DistSpMV_Balanced — 论文复现说明

项目目标
- 复现论文 “Balancing Computation and Communication in Distributed Sparse Matrix-Vector Multiplication” 的核心思路 DistSpMV_Balanced。
- 使用纯 C，MPI（进程间通信）和 OpenMP（线程并行），并调用 METIS 做行划分与重排。

主要实现说明
- 输入格式：Matrix Market (.mtx)，支持 coordinate 格式（general / symmetric）。
- 存储格式：CSR（Compressed Sparse Row）。
- 划分与重排：在 rank 0 上调用 METIS（METIS_PartGraphKway）对行进行划分，随后把行和列重排，使非零元尽可能聚集到对角块附近。
- Balanced Boundary：使用基于覆盖本地非零元比例的滑动窗口（sliding window）启发式确定每个进程的对角区间（left_bound / right_bound），以尽可能把大部分计算放在本地，减少远程 x 元素请求。
- 通信（Algorithm 2 & 3）：每轮先统计本进程需要的远端列索引，向对应所有者发送请求（MPI_Isend / MPI_Irecv），所有者回复所需的 x 值。采用非阻塞发送/接收以便与计算重叠。
- 计算（Algorithm 4）：节点内使用 OpenMP 对行进行并行计算。首先计算仅依赖本地 x 值的部分，通信完成后再完成剩余部分。

性能计量
- 输出总运行时间、计算时间（compute）和通信时间（communication）。
- 计算 GFlops：采用 2*nnz（乘加）作为每次 SpMV 的浮点运算数，除以最大的计算时间得到全局 GFlops。

编译与运行
1. 依赖
   - MPI（Open MPI / MPICH）
   - OpenMP（编译器支持 `-fopenmp`）
   - METIS（开发头文件与库，通常通过 `-lmetis` 链接）

2. 编译示例
```
mpicc -O3 -fopenmp -I/usr/include -L/usr/lib -lmetis -o dist_spmv_balanced dist_spmv_balanced.c
```

3. 运行示例
```
export OMP_NUM_THREADS=4
mpirun -np 4 ./dist_spmv_balanced matrix.mtx 10 0.8
```
参数说明：
- `matrix.mtx`：输入矩阵（Matrix Market 格式）
- `10`：迭代次数（可选，默认 10）
- `0.8`：lower_bound（可选），表示对角区间应覆盖本地非零的比例阈值（在实现中用于启发式选择对角区间）。

实现与工程说明
- 本实现为可复现的研究原型，注重清晰与可读性，某些部分（如列所有者投票统计）为简化实现，可能不是最优的工程实现。
- METIS 调用与重排仅在 rank 0 上执行，然后把重排后的 CSR 广播到所有进程。真实大规模场景应当考虑并行划分/分布式输入。
- 通信采用两阶段：先交换请求列表（索引），再交换具体的 x 值；请求与回复使用非阻塞通信以便与本地计算重叠。

调优建议
- 根据本机内核数设置 `OMP_NUM_THREADS`。
- 对于大矩阵，注意进程数与 METIS 分区数一致。METIS 内存开销会随矩阵尺寸上升。
- 若要进一步优化：减少请求列表大小（去重与压缩）、用哈希/位图快速判重、采用 MPI_Alltoallv 或持久化通信通道以减少开销。

下一步
- 我可以帮助：
  - 调试/运行示例矩阵并收集性能数据；
  - 优化通信（合并消息、压缩索引）；
  - 增加单元测试或生成小型示例矩阵用于验证正确性。

作者：自动生成的复现实验原型（供教学与科研参考）。
