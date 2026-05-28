/*
 * 编译：
 *   mpicc -O3 -fopenmp -o dist_spmv_balanced dist_spmv_balanced.c -lm
 * 运行：
 *   export OMP_NUM_THREADS=4
 *   mpirun -np 4 ./dist_spmv_balanced matrix.mtx 10 0.8
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <metis.h>

/* ─────────────────────────────────────────────
 * CSR 数据结构
 * ───────────────────────────────────────────── */
typedef struct {
    int    nrows, ncols;
    long long nnz;
    int   *rowptr;   /* 长度 nrows+1 */
    int   *colidx;   /* 长度 nnz     */
    double *vals;    /* 长度 nnz     */
} CSR;

static void csr_free(CSR *A) {
    free(A->rowptr); free(A->colidx); free(A->vals);
    A->rowptr = NULL; A->colidx = NULL; A->vals = NULL;
}

/* 局部 CSR 版本的 compute_balanced_boundary：对每个进程接收到的本地行段运行 scanline 扩展边界 */
static void compute_balanced_boundary_local(
        const int *lrowptr, const int *lcolidx,
        int global_rstart, int local_nrows, int ncols,
        int cstart, int cend,
        long long lower_bound,
        int *leftBound, int *rightBound,
        int *rowptr_start_out, int *rowptr_end_out)
{
    if (local_nrows <= 0) {
        *leftBound = cstart;
        *rightBound = cend - 1;
        return;
    }
    int *rptr_start = (int*)malloc((size_t)local_nrows * sizeof(int));
    int *rptr_end   = (int*)malloc((size_t)local_nrows * sizeof(int));
    long long diag_nnz = 0;

    for (int lr = 0; lr < local_nrows; ++lr) {
        rptr_start[lr] = -1;
        rptr_end[lr] = cstart - 1;
        for (int k = lrowptr[lr]; k < lrowptr[lr+1]; ++k) {
            int c = lcolidx[k];
            if (c >= cstart && c < cend) {
                if (rptr_start[lr] == -1) rptr_start[lr] = c;
                if (c > rptr_end[lr]) rptr_end[lr] = c;
                diag_nnz++;
            }
        }
        if (rptr_start[lr] == -1) rptr_start[lr] = cstart;
    }

    int *ptr = (int*)malloc((size_t)local_nrows * sizeof(int));
    for (int lr = 0; lr < local_nrows; ++lr) {
        int k = lrowptr[lr];
        while (k < lrowptr[lr+1] && lcolidx[k] < cend) k++;
        ptr[lr] = k;
    }

    long long cur_nnz = diag_nnz;
    int cur_right = cend - 1;

    while (cur_nnz < lower_bound && cur_right + 1 < ncols) {
        int next_col = ncols;
        for (int lr = 0; lr < local_nrows; ++lr) {
            while (ptr[lr] < lrowptr[lr+1] && lcolidx[ptr[lr]] <= cur_right) ptr[lr]++;
            if (ptr[lr] < lrowptr[lr+1] && lcolidx[ptr[lr]] < next_col)
                next_col = lcolidx[ptr[lr]];
        }
        if (next_col == ncols) break;
        cur_right = next_col;
        for (int lr = 0; lr < local_nrows; ++lr) {
            if (ptr[lr] < lrowptr[lr+1] && lcolidx[ptr[lr]] == next_col) {
                if (next_col > rptr_end[lr]) rptr_end[lr] = next_col;
                cur_nnz++;
            }
        }
    }

    if (rowptr_start_out && rowptr_end_out) {
        for (int lr = 0; lr < local_nrows; ++lr) {
            int gr = global_rstart + lr;
            rowptr_start_out[gr] = rptr_start[lr];
            rowptr_end_out[gr]   = rptr_end[lr];
        }
    }
    *leftBound = cstart;
    *rightBound = cur_right;

    free(ptr);
    free(rptr_start); free(rptr_end);
}

/* ─────────────────────────────────────────────
 * Matrix Market 读取（coordinate real/integer/pattern, general/symmetric）
 * ───────────────────────────────────────────── */
static int int_cmp(const void *a, const void *b) {
    return (*(const int*)a > *(const int*)b) - (*(const int*)a < *(const int*)b);
}

/* 按 (row, col) 排序辅助 */
typedef struct { int row, col; double val; } Triplet;
static int triplet_cmp(const void *a, const void *b) {
    const Triplet *ta = (const Triplet*)a, *tb = (const Triplet*)b;
    if (ta->row != tb->row) return ta->row - tb->row;
    return ta->col - tb->col;
}

/*
 * 返回 0 成功；负数失败。
 * 支持 symmetric：自动展开下三角 → 同时写上三角（跳过对角线重复）。
 */
int read_matrix_market(const char *path, CSR *A) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n无法打开文件：%s\n", path, path); return -1; }

    /* ── 解析 banner ── */
    char line[256];
    int is_sym = 0, is_pattern = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
        if (strstr(line, "symmetric") || strstr(line, "Symmetric")) is_sym = 1;
        if (strstr(line, "pattern")   || strstr(line, "Pattern"))   is_pattern = 1;
    }

    /* ── 读尺寸行（已被上面循环读入 line） ── */
    int nrows, ncols; long long nnz_file;
    if (sscanf(line, "%d %d %lld", &nrows, &ncols, &nnz_file) != 3) {
        fclose(f); return -2;
    }

    /* 预分配（symmetric 最多 2× 非零元）*/
    long long cap = is_sym ? nnz_file * 2 : nnz_file;
    Triplet *trips = (Triplet*)malloc((size_t)cap * sizeof(Triplet));
    if (!trips) { fclose(f); return -3; }

    long long cnt = 0;
    for (long long k = 0; k < nnz_file; ++k) {
        int i, j; double v = 1.0;
        if (is_pattern) {
            if (fscanf(f, "%d %d",      &i, &j)    != 2) break;
        } else {
            if (fscanf(f, "%d %d %lf", &i, &j, &v) != 3) break;
        }
        i--; j--;   /* 转为 0-based */
        if (i < 0 || i >= nrows || j < 0 || j >= ncols) continue;
        trips[cnt].row = i; trips[cnt].col = j; trips[cnt].val = v; cnt++;
        if (is_sym && i != j) {
            trips[cnt].row = j; trips[cnt].col = i; trips[cnt].val = v; cnt++;
        }
    }
    fclose(f);

    /* 去重（symmetric 可能有上下三角重叠）*/
    qsort(trips, (size_t)cnt, sizeof(Triplet), triplet_cmp);
    long long ucnt = 0;
    for (long long k = 0; k < cnt; ++k) {
        if (ucnt > 0 && trips[ucnt-1].row == trips[k].row && trips[ucnt-1].col == trips[k].col)
            trips[ucnt-1].val += trips[k].val;   /* 累加重复项 */
        else
            trips[ucnt++] = trips[k];
    }
    cnt = ucnt;

    /* ── 构建 CSR ── */
    A->nrows  = nrows; A->ncols  = ncols; A->nnz = cnt;
    A->rowptr = (int*)calloc((size_t)(nrows+1), sizeof(int));
    A->colidx = (int*)malloc((size_t)cnt * sizeof(int));
    A->vals   = (double*)malloc((size_t)cnt * sizeof(double));

    for (long long k = 0; k < cnt; ++k) A->rowptr[trips[k].row + 1]++;
    for (int i = 0; i < nrows; ++i)     A->rowptr[i+1] += A->rowptr[i];
    for (long long k = 0; k < cnt; ++k) {
        A->colidx[k] = trips[k].col;
        A->vals[k]   = trips[k].val;
    }
    free(trips);
    return 0;
}

/* ─────────────────────────────────────────────
 * 均等块划分辅助
 * ───────────────────────────────────────────── */
static void make_block_offsets(int n, int np, int *off) {
    int base = n / np, rem = n % np;
    off[0] = 0;
    for (int p = 0; p < np; ++p) off[p+1] = off[p] + base + (p < rem ? 1 : 0);
}

/* ─────────────────────────────────────────────
 * 线程级负载均衡：按 nnz 均分行边界
 * 线程边界数组需分配 nthreads+1 大小
 * ───────────────────────────────────────────── */
static void compute_thread_boundaries(
        const int *rowptr, int nrows, int nthreads, int *thread_row_start) {
    long long total_nnz = (long long)rowptr[nrows] - (long long)rowptr[0];
    thread_row_start[0] = 0;
    for (int t = 1; t < nthreads; t++) {
        long long target = (long long)rowptr[0] + (total_nnz * t) / nthreads;
        int lo = thread_row_start[t-1], hi = nrows;
        while (lo < hi) {
            int mid = lo + (hi - lo) / 2;
            if ((long long)rowptr[mid] < target) lo = mid + 1;
            else hi = mid;
        }
        thread_row_start[t] = lo;
    }
    thread_row_start[nthreads] = nrows;
}

/* ─────────────────────────────────────────────
 * 【修复问题2】balanced boundary
 *
 * 与论文 Algorithm 1 对齐：
 *   - 起点是本进程的对角块 [cstart, cend)
 *   - 逐步将 right_bound 向右移动，
 *     直到区间内非零元数 >= lower_bound_nnz，或到达矩阵右边界
 *   - 遍历本进程负责的每一行，找到 rowptr_end[row]+1 处的下一个非零列，
 *     取所有行中最小的那个列号作为新的 right_bound（scanline）
 *
 * 参数：
 *   A          : 全局 CSR（或本地行的 CSR，行号已是全局）
 *   rstart/rend: 本进程负责的行范围 [rstart, rend)
 *   cstart/cend: 本进程对角块的初始列范围
 *   lb_frac    : lower_bound 比例 (0,1]
 *   leftBound  : 输出 —— 对角块左边界（固定为 cstart）
 *   rightBound : 输出 —— 对角块右边界（扩展后）
 *
 * 注：本函数在每个进程上独立调用，使用该进程能看到的行数据。
 * ───────────────────────────────────────────── */
static void compute_balanced_boundary(
        const int *rowptr, const int *colidx,
        int rstart, int rend, int ncols,
        int cstart, int cend,
        long long lower_bound,
        int *leftBound, int *rightBound,
        int *rowptr_start_out, int *rowptr_end_out)
{
    int local_nrows = rend - rstart;
    if (local_nrows <= 0) {
        *leftBound = cstart;
        *rightBound = cend - 1;
        return;
    }
    int *rptr_start = (int*)malloc((size_t)local_nrows * sizeof(int));
    int *rptr_end   = (int*)malloc((size_t)local_nrows * sizeof(int));
    long long diag_nnz = 0;

    /* 先统计对角块内初始非零元并初始化每行的 rptr_start/rptr_end */
    for (int lr = 0; lr < local_nrows; ++lr) {
        int gr = rstart + lr;
        rptr_start[lr] = -1;
        rptr_end[lr] = cstart - 1;
        for (int k = rowptr[gr]; k < rowptr[gr+1]; ++k) {
            int c = colidx[k];
            if (c >= cstart && c < cend) {
                if (rptr_start[lr] == -1) rptr_start[lr] = c;
                if (c > rptr_end[lr]) rptr_end[lr] = c;
                diag_nnz++;
            }
        }
        if (rptr_start[lr] == -1) rptr_start[lr] = cstart;
    }

    /* 初始化 ptr：第一个 >= cend 的位置（跳过对角块内部） */
    int *ptr = (int*)malloc((size_t)local_nrows * sizeof(int));
    for (int lr = 0; lr < local_nrows; ++lr) {
        int gr = rstart + lr;
        int k = rowptr[gr];
        while (k < rowptr[gr+1] && colidx[k] < cend) k++;
        ptr[lr] = k;
    }

    long long cur_nnz = diag_nnz;
    int cur_right = cend - 1;

    /* 单调推进 ptr，按列逐步扩展 rightBound，直到达到 lower_bound */
    while (cur_nnz < lower_bound && cur_right + 1 < ncols) {
        int next_col = ncols;
        /* 找到所有行中 ptr 指向的最小列号 */
        for (int lr = 0; lr < local_nrows; ++lr) {
            int gr = rstart + lr;
            while (ptr[lr] < rowptr[gr+1] && colidx[ptr[lr]] <= cur_right) ptr[lr]++;
            if (ptr[lr] < rowptr[gr+1] && colidx[ptr[lr]] < next_col)
                next_col = colidx[ptr[lr]];
        }
        if (next_col == ncols) break;
        cur_right = next_col;

        /* 统计本列命中的行并更新 cur_nnz */
        for (int lr = 0; lr < local_nrows; ++lr) {
            int gr = rstart + lr;
            if (ptr[lr] < rowptr[gr+1] && colidx[ptr[lr]] == next_col) {
                if (next_col > rptr_end[lr]) rptr_end[lr] = next_col;
                cur_nnz++;
                /* ptr[lr] 在下一轮的 while 中推进 */
            }
        }
    }

    /* 写回全局输出数组（按全局行索引） */
    for (int lr = 0; lr < local_nrows; ++lr) {
        int gr = rstart + lr;
        rowptr_start_out[gr] = rptr_start[lr];
        rowptr_end_out[gr]   = rptr_end[lr];
    }
    *leftBound = cstart;
    *rightBound = cur_right; /* 闭区间 */

    free(ptr);
    free(rptr_start); free(rptr_end);
}

/* ═══════════════════════════════════════════════
 * main
 * ═══════════════════════════════════════════════ */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc < 2) {
        if (rank == 0)
            fprintf(stderr, "Usage: %s matrix.mtx [niter=10] [lower_bound=0.8]\n用法: %s matrix.mtx [niter=10] [lower_bound=0.8]\n", argv[0], argv[0]);
        MPI_Finalize(); return 1;
    }
    const char *filename  = argv[1];
    int    niter          = (argc >= 3 ? atoi(argv[2]) : 10);
    /* 命令行传入的 lower_bound 被解释为比例系数（0,1]，表示以最大对角块 nnz 的比例作为阈值 */
    double lower_bound_frac    = (argc >= 4 ? atof(argv[3]) : 0.8);
    /* 模式："balanced"（默认） 或 "naive"（基线：Allgatherv 全局广播 x） */
    const char *mode       = (argc >= 5 ? argv[4] : "balanced");
    /* ── 性能计时 ── */

    /* ── 读取矩阵并做初始分区/本地数据准备（参考 dist_spmv_naive.c） ── */
    CSR A; memset(&A,0,sizeof(A));
    /* 注意：当前实现为简化复现，所有进程均读取完整矩阵文件到本地 A。
     * 在真实大规模分布式实验中应由 rank 0 读取并通过 MPI_Scatterv 等
     * 将各进程所需的行分发出去，以避免每个进程持有完整矩阵导致 OOM。
     */
    if (read_matrix_market(filename, &A) != 0) {
        if (rank==0) fprintf(stderr, "Failed to read %s\n无法读取矩阵文件：%s\n", filename, filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nrows = A.nrows, ncols = A.ncols; long long global_nnz = A.nnz;
    int symmetric_perm = (nrows == ncols);

    /* ---------- STEP 1: 使用 METIS 对矩阵行/列重排（仅由 rank 0 计算，然后广播 perm） ---------- */
    int *perm = NULL; /* old_index -> new_index */
    int metis_ok = 0;
    if (rank == 0) {
        int *col_counts_tmp = (int*)calloc((size_t)ncols, sizeof(int));
        for (int i = 0; i < nrows; ++i)
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k)
                col_counts_tmp[A.colidx[k]]++;

        int *col_disp_tmp = (int*)malloc((size_t)(ncols+1) * sizeof(int));
        col_disp_tmp[0] = 0;
        for (int c = 0; c < ncols; ++c) col_disp_tmp[c+1] = col_disp_tmp[c] + col_counts_tmp[c];

        int *rows_in_col = (int*)malloc((size_t)A.nnz * sizeof(int));
        int *tmp_idx = (int*)malloc((size_t)ncols * sizeof(int));
        for (int c = 0; c < ncols; ++c) tmp_idx[c] = col_disp_tmp[c];
        for (int i = 0; i < nrows; ++i)
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                int c = A.colidx[k]; rows_in_col[tmp_idx[c]++] = i;
            }
        free(tmp_idx);

        int *mark = (int*)malloc((size_t)nrows * sizeof(int));
        for (int i = 0; i < nrows; ++i) mark[i] = -1;
        int *xadj = (int*)malloc((size_t)(nrows+1) * sizeof(int));
        xadj[0] = 0;
        for (int i = 0; i < nrows; ++i) {
            int deg = 0;
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                int c = A.colidx[k];
                for (int t = col_disp_tmp[c]; t < col_disp_tmp[c+1]; ++t) {
                    int v = rows_in_col[t];
                    if (v == i) continue;
                    if (mark[v] != i) { mark[v] = i; deg++; }
                }
            }
            xadj[i+1] = deg;
        }
        for (int i = 0; i < nrows; ++i) xadj[i+1] += xadj[i];
        int nedges = xadj[nrows];
        int *adjncy = (int*)malloc((size_t)nedges * sizeof(int));
        for (int i = 0; i < nrows; ++i) mark[i] = -1;
        for (int i = 0; i < nrows; ++i) {
            int pos = xadj[i];
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                int c = A.colidx[k];
                for (int t = col_disp_tmp[c]; t < col_disp_tmp[c+1]; ++t) {
                    int v = rows_in_col[t];
                    if (v == i) continue;
                    if (mark[v] != i) { mark[v] = i; adjncy[pos++] = v; }
                }
            }
        }

        /* 调用 METIS 分区（仅在 rank 0） */
        idx_t nvtxs = (idx_t)nrows;
        idx_t ncon = 1;
        idx_t *xadj_idx = (idx_t*)malloc((size_t)(nrows+1) * sizeof(idx_t));
        idx_t *adjncy_idx = (idx_t*)malloc((size_t)nedges * sizeof(idx_t));
        for (int i = 0; i <= nrows; ++i) xadj_idx[i] = (idx_t)xadj[i];
        for (int i = 0; i < nedges; ++i) adjncy_idx[i] = (idx_t)adjncy[i];

        idx_t nparts = (idx_t)nprocs;
        idx_t objval = 0;
        idx_t *part_idx = (idx_t*)malloc((size_t)nrows * sizeof(idx_t));

        int need_fallback = 0;
        if (nvtxs <= 0) need_fallback = 1;
        if (nparts < 1) nparts = 1;
        if (nparts > nvtxs) nparts = nvtxs;
        if (nparts <= 1) { for (int i = 0; i < nrows; ++i) part_idx[i] = 0; need_fallback = 1; }
        if (!need_fallback) {
            for (int i = 1; i <= nrows; ++i) if (xadj_idx[i] < xadj_idx[i-1]) { need_fallback = 1; break; }
        }
        if (!need_fallback && nedges > 0) {
            for (int i = 0; i < nedges; ++i) if (adjncy_idx[i] < 0 || adjncy_idx[i] >= nvtxs) { need_fallback = 1; break; }
        }
        if (!need_fallback && xadj_idx[nrows] == 0) need_fallback = 1;

        if (!need_fallback) {
            int metis_ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj_idx, adjncy_idx,
                                               NULL, NULL, NULL, &nparts, NULL, NULL,
                                               NULL, &objval, part_idx);
            if (metis_ret != METIS_OK) need_fallback = 1;
        }
        if (need_fallback) {
            for (int i = 0; i < nrows; ++i) part_idx[i] = (idx_t)((long long)i * (long long)nprocs / (long long)nrows);
        }

        /* 根据分区生成 old->new 映射（仅在 rank 0） */
        int *count_per_part = (int*)calloc((size_t)nprocs, sizeof(int));
        for (int i = 0; i < nrows; ++i) count_per_part[(int)part_idx[i]]++;
        int *offp = (int*)malloc((size_t)(nprocs+1) * sizeof(int)); offp[0] = 0;
        for (int p = 0; p < nprocs; ++p) offp[p+1] = offp[p] + count_per_part[p];
        int *pos_arr = (int*)malloc((size_t)nprocs * sizeof(int));
        for (int p = 0; p < nprocs; ++p) pos_arr[p] = offp[p];
        int *perm_new = (int*)malloc((size_t)nrows * sizeof(int)); /* new_index -> old_index */
        for (int i = 0; i < nrows; ++i) {
            int p = (int)part_idx[i]; perm_new[ pos_arr[p]++ ] = i;
        }
        perm = (int*)malloc((size_t)nrows * sizeof(int)); /* old_index -> new_index */
        for (int ni = 0; ni < nrows; ++ni) perm[ perm_new[ni] ] = ni;

        /* cleanup rank-0 temporaries (keep perm) */
        free(xadj); free(adjncy); free(xadj_idx); free(adjncy_idx); free(part_idx);
        free(count_per_part); free(offp); free(pos_arr); free(perm_new);
        free(col_counts_tmp); free(col_disp_tmp); free(rows_in_col); free(mark);

        metis_ok = 1;
        if (rank == 0) printf("METIS partition computed on rank 0\n");
    }

    /* Broadcast whether rank0 produced a perm, then broadcast perm (old->new) to all ranks */
    MPI_Bcast(&metis_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (metis_ok) {
        if (rank != 0) perm = (int*)malloc((size_t)nrows * sizeof(int));
        MPI_Bcast(perm, nrows, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("METIS permutation applied\n");
    } else {
        /* fallback: identity mapping */
        if (rank != 0) perm = (int*)malloc((size_t)nrows * sizeof(int));
        for (int i = 0; i < nrows; ++i) perm[i] = i;
    }

    /* 所有进程使用 perm 对矩阵做相同的重排（old->new） */
    {
        fprintf(stderr, "rank %d: applying perm (nrows=%d ncols=%d nnz=%lld)\n", rank, nrows, ncols, (long long)A.nnz);
        fflush(stderr);
        Triplet *trips = (Triplet*)malloc((size_t)A.nnz * sizeof(Triplet));
        long long t = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                /* 行总是重排；当矩阵为方阵时同时对列做相同的对称重排 */
                trips[t].row = perm[i];
                if (symmetric_perm) trips[t].col = perm[A.colidx[k]];
                else                trips[t].col = A.colidx[k];
                trips[t].val = A.vals[k];
                t++;
            }
        }
        qsort(trips, (size_t)t, sizeof(Triplet), triplet_cmp);
        int *new_rowptr = (int*)calloc((size_t)(nrows+1), sizeof(int));
        int *new_colidx = (int*)malloc((size_t)t * sizeof(int));
        double *new_vals = (double*)malloc((size_t)t * sizeof(double));
        for (long long k = 0; k < t; ++k) new_rowptr[trips[k].row + 1]++;
        for (int i = 0; i < nrows; ++i) new_rowptr[i+1] += new_rowptr[i];
        for (long long k = 0; k < t; ++k) { new_colidx[k] = trips[k].col; new_vals[k] = trips[k].val; }
        free(A.rowptr); free(A.colidx); free(A.vals);
        A.rowptr = new_rowptr; A.colidx = new_colidx; A.vals = new_vals; A.nnz = t;
        free(trips);
        fprintf(stderr, "rank %d: perm applied, new nnz=%lld\n", rank, (long long)A.nnz);
        fflush(stderr);
        /* perm 已应用，保留 perm 以便后续对 x 的重排（稍后释放） */
    }

    /* 初始按均等块划分行/列（作为 column ownership 基准） */
    int *row_offsets = (int*)malloc((size_t)(nprocs+1)*sizeof(int));
    int *col_offsets = (int*)malloc((size_t)(nprocs+1)*sizeof(int));
    make_block_offsets(nrows, nprocs, row_offsets);
    make_block_offsets(ncols, nprocs, col_offsets);

    fprintf(stderr, "rank %d: row_offsets[%d,%d], col_offsets[%d,%d]\n", rank, row_offsets[rank], row_offsets[rank+1], col_offsets[rank], col_offsets[rank+1]); fflush(stderr);

    /* 统计每个对角块初始非零数，并取最大值作为基准（Algorithm 1 的第一步）
     * 注意：diag_nnz 对所有进程统一 calloc 初始化以避免未初始化访问风险 */
    long long *diag_nnz = (long long*)calloc((size_t)nprocs, sizeof(long long));
    if (rank == 0) {
        for (int p = 0; p < nprocs; ++p) diag_nnz[p] = 0;
        for (int p = 0; p < nprocs; ++p) {
            for (int i = row_offsets[p]; i < row_offsets[p+1]; ++i) {
                for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                    int c = A.colidx[k];
                    if (c >= col_offsets[p] && c < col_offsets[p+1]) diag_nnz[p]++;
                }
            }
        }
    }
    MPI_Bcast(diag_nnz, nprocs, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    /* 以各对角块 nnz 最大值为基准；命令行的比例系数 lower_bound_frac 控制阈值比例 */
    long long max_diag_nnz = 0;
    for (int p = 0; p < nprocs; ++p) if (diag_nnz[p] > max_diag_nnz) max_diag_nnz = diag_nnz[p];
    long long lower_bound_global = (long long)(max_diag_nnz * lower_bound_frac);
    if (lower_bound_global < 1) lower_bound_global = 1;

    /* 为 Algorithm 1 准备（可选的）行级缓存指针；在分发后由各进程按需填写 */
    int *rowptr_start = NULL;
    int *rowptr_end   = NULL;

    /* Rank-0 将对 A(rowptr/colidx/vals) 进行分块并使用 MPI_Scatterv 将本地 CSR 发给各进程。
     * 各进程仅保留其本地行片段（local_rowptr/local_colidx/local_vals），随后在本地运行
     * compute_balanced_boundary_local（基于本地 CSR）以确定 left/right boundary。
     */
    int *send_nnz = NULL, *displ_nnz = NULL;
    int *send_rowptr_counts = NULL, *displ_rowptr = NULL;
    int total_rowptr_counts = 0;
    int *rowptr_all = NULL;
    if (rank == 0) {
        send_nnz = (int*)malloc((size_t)nprocs * sizeof(int));
        displ_nnz = (int*)malloc((size_t)nprocs * sizeof(int));
        send_rowptr_counts = (int*)malloc((size_t)nprocs * sizeof(int));
        displ_rowptr = (int*)malloc((size_t)nprocs * sizeof(int));
        total_rowptr_counts = 0;
        for (int p = 0; p < nprocs; ++p) {
            int s = row_offsets[p], e = row_offsets[p+1];
            int lnrows = e - s;
            int lnnz = A.rowptr[e] - A.rowptr[s];
            send_nnz[p] = lnnz;
            displ_nnz[p] = A.rowptr[s];
            send_rowptr_counts[p] = lnrows + 1;
            displ_rowptr[p] = total_rowptr_counts;
            total_rowptr_counts += send_rowptr_counts[p];
        }
        rowptr_all = (int*)malloc((size_t)total_rowptr_counts * sizeof(int));
        int off = 0;
        for (int p = 0; p < nprocs; ++p) {
            int s = row_offsets[p], e = row_offsets[p+1];
            int lnrows = e - s;
            int base = A.rowptr[s];
            rowptr_all[off] = 0;
            for (int i = 1; i <= lnrows; ++i) rowptr_all[off + i] = A.rowptr[s + i] - base;
            off += lnrows + 1;
        }
    } else {
        /* 非 root 进程分配接收计数数组（会由 root 广播填充） */
        send_nnz = (int*)malloc((size_t)nprocs * sizeof(int));
        send_rowptr_counts = (int*)malloc((size_t)nprocs * sizeof(int));
    }

    /* 广播每个进程将接收的 nnz 与 rowptr 长度，以及行/列偏移信息 */
    MPI_Bcast(send_nnz, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(send_rowptr_counts, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(row_offsets, nprocs+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(col_offsets, nprocs+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&lower_bound_global, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    int local_rowptr_len = send_rowptr_counts[rank];
    int local_nrows = local_rowptr_len - 1;
    int local_nnz = send_nnz[rank];
    int *local_rowptr = (int*)malloc((size_t)local_rowptr_len * sizeof(int));
    int *local_colidx = (int*)malloc((size_t)local_nnz * sizeof(int));
    double *local_vals = (double*)malloc((size_t)local_nnz * sizeof(double));

    MPI_Scatterv(rowptr_all, send_rowptr_counts, displ_rowptr, MPI_INT,
                 local_rowptr, local_rowptr_len, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(A.colidx, send_nnz, displ_nnz, MPI_INT,
                 local_colidx, local_nnz, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(A.vals, send_nnz, displ_nnz, MPI_DOUBLE,
                 local_vals, local_nnz, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Rank-0 分发完成后可以释放全局矩阵 A 以节省内存 */
    if (rank == 0) {
        csr_free(&A);
        free(rowptr_all); free(displ_nnz); free(displ_rowptr);
        free(send_nnz); free(send_rowptr_counts);
    } else {
        free(send_nnz); free(send_rowptr_counts);
    }

    double *y_local = (double*)malloc((size_t)local_nrows * sizeof(double));

    /* 在本地 CSR 上计算 left/right boundary 并收集全局边界信息 */
    int leftBound_local = col_offsets[rank];
    int rightBound_local = col_offsets[rank+1] - 1;
    compute_balanced_boundary_local(local_rowptr, local_colidx, row_offsets[rank], local_nrows, ncols,
                                    col_offsets[rank], col_offsets[rank+1], lower_bound_global,
                                    &leftBound_local, &rightBound_local,
                                    NULL, NULL);
    int *leftBounds = (int*)malloc((size_t)nprocs * sizeof(int));
    int *rightBounds = (int*)malloc((size_t)nprocs * sizeof(int));
    MPI_Allgather(&leftBound_local, 1, MPI_INT, leftBounds, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Allgather(&rightBound_local, 1, MPI_INT, rightBounds, 1, MPI_INT, MPI_COMM_WORLD);

    /* owner_of_col 基于列均等划分（col_offsets） */
    int *owner_of_col = (int*)malloc((size_t)ncols * sizeof(int));
    for (int p = 0; p < nprocs; ++p)
        for (int j = col_offsets[p]; j < col_offsets[p+1]; ++j)
            owner_of_col[j] = p;

    /* 本地/全局 x 向量缓冲 */
    int cstart = col_offsets[rank];
    int cend = col_offsets[rank+1];
    int local_ncols = cend - cstart;
    double *x_owned = (double*)malloc((size_t)local_ncols * sizeof(double));
    /* 若进行了对称重排，需要将原始 x（此处初始化为 1）按 perm 重排到新编号上 */
    if (symmetric_perm) {
        double *x_init = (double*)malloc((size_t)ncols * sizeof(double));
        double *x_perm = (double*)malloc((size_t)ncols * sizeof(double));
        for (int i = 0; i < ncols; ++i) x_init[i] = 1.0;
        for (int i = 0; i < ncols; ++i) x_perm[ perm[i] ] = x_init[i];
        for (int i = 0; i < local_ncols; ++i) x_owned[i] = x_perm[cstart + i];
        free(x_init); free(x_perm);
    } else {
        for (int i = 0; i < local_ncols; ++i) x_owned[i] = 1.0;
    }
    double *x_buf = (double*)calloc((size_t)ncols, sizeof(double));
    for (int j = cstart; j < cend; ++j) x_buf[j] = x_owned[j - cstart];

    int *col_counts_init = (int*)malloc((size_t)nprocs * sizeof(int));
    int *col_disp_init   = (int*)malloc((size_t)nprocs * sizeof(int));
    for (int p = 0; p < nprocs; ++p) { col_counts_init[p] = col_offsets[p+1] - col_offsets[p]; col_disp_init[p] = col_offsets[p]; }

    double compute_time = 0.0, comm_time = 0.0;

    /* 填充 x_buf 的初始全局值：一次性 Allgatherv 将各进程持有的 x_owned 合并到 x_buf
     * 以确保 boundary 扩展进入的列位置有正确的 x 值（避免使用 calloc 的 0 值）。
     * 注意：此初始合并为预热/初始化通信，不计入随后的迭代通信时间（t0 之后）。
     * 也注意：若 x 在每轮迭代都会更新，则需要在每轮保证扩展区间 x_buf 的最新性；
     * 当前复现中 x 固定为 1，因此一次性合并即可。
     */
    MPI_Allgatherv(x_owned, local_ncols, MPI_DOUBLE,
                   x_buf, col_counts_init, col_disp_init, MPI_DOUBLE,
                   MPI_COMM_WORLD);
    /* 释放 perm（若存在），此前已用来对 A 和 x 进行重排 */
    if (perm) { free(perm); perm = NULL; }
    free(col_counts_init); free(col_disp_init);

    /* 释放在预处理阶段产生的全局中间数组（注意 prefix/remote_nnz 已不再使用） */
    free(diag_nnz);
    free(leftBounds); free(rightBounds);

    /* 将本进程 boundary 复制为后续迭代使用的标量 */
    int leftBound = leftBound_local;
    int rightBound = rightBound_local;

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();

    /* ════════════════════════════════════════════
     * 主迭代循环：两种模式
     *  - naive: 基线实现，使用 MPI_Allgatherv 每轮聚合全局 x 再本地计算
     *  - balanced: 原 DistSpMV_Balanced 实现（如上）
     * ════════════════════════════════════════════ */
    int is_naive = (strcmp(mode, "naive") == 0);

    if (is_naive) {
        /* 准备 Allgatherv 的计数与位移（按列块） */
        int *col_counts = (int*)malloc((size_t)nprocs * sizeof(int));
        int *col_disp   = (int*)malloc((size_t)nprocs * sizeof(int));
        for (int p = 0; p < nprocs; ++p) {
            col_counts[p] = col_offsets[p+1] - col_offsets[p];
            col_disp[p]   = col_offsets[p];
        }

        for (int iter = 0; iter < niter; ++iter) {
            /* 先交换全局 x（通信时间计入 comm_time） */
            double comm0 = MPI_Wtime();
            MPI_Allgatherv(x_owned, local_ncols, MPI_DOUBLE,
                           x_buf, col_counts, col_disp, MPI_DOUBLE,
                           MPI_COMM_WORLD);
            double comm1 = MPI_Wtime(); comm_time += (comm1 - comm0);

            /* 然后本地完成所有行的 SpMV（计入 compute_time） */
            double comp0 = MPI_Wtime();
            {
                int nthreads = omp_get_max_threads();
                int *thr_bound = (int*)malloc((nthreads + 1) * sizeof(int));
                compute_thread_boundaries(local_rowptr, local_nrows, nthreads, thr_bound);
                #pragma omp parallel num_threads(nthreads)
                {
                    int tid = omp_get_thread_num();
                    for (int lr = thr_bound[tid]; lr < thr_bound[tid+1]; ++lr) {
                        double s = 0.0;
                        for (int k = local_rowptr[lr]; k < local_rowptr[lr+1]; ++k) {
                            int col = local_colidx[k];
                            s += local_vals[k] * x_buf[col];
                        }
                        y_local[lr] = s;
                    }
                }
                free(thr_bound);
            }
            double comp1 = MPI_Wtime(); compute_time += (comp1 - comp0);
        }

        free(col_counts); free(col_disp);
    } else {
        /* 原有 balanced 实现保留 */
        /* 预分配 send_lists 与容量，循环内复用，避免每轮大量 malloc/free（修复问题5） */
        int *send_counts = (int*)malloc((size_t)nprocs * sizeof(int));
        int *send_caps   = (int*)malloc((size_t)nprocs * sizeof(int));
        int **send_lists = (int**)malloc((size_t)nprocs * sizeof(int*));
        for (int p = 0; p < nprocs; ++p) {
            send_caps[p]  = 64;
            send_lists[p] = (int*)malloc((size_t)send_caps[p] * sizeof(int));
        }

        for (int iter = 0; iter < niter; ++iter) {

            /* ── Phase A: 构建远端列请求（去重后按目标进程分组）── */
            /* 每轮只重置计数器，复用 send_lists/send_caps */
            memset(send_counts, 0, (size_t)nprocs * sizeof(int));

            for (int lr = 0; lr < local_nrows; ++lr) {
                for (int k = local_rowptr[lr]; k < local_rowptr[lr+1]; ++k) {
                    int col   = local_colidx[k];
                    int owner = owner_of_col[col];
                    /* 【修复问题3/6】只有 col 不在 [leftBound,rightBound] 才需要远端请求 */
                    if (col < leftBound || col > rightBound) {
                        if (owner == rank) continue;   /* 自己负责且在 boundary 外（理论上不发生）*/
                        int cnt = send_counts[owner];
                        if (cnt >= send_caps[owner]) {
                            send_caps[owner] *= 2;
                            send_lists[owner] = (int*)realloc(send_lists[owner],
                                                (size_t)send_caps[owner] * sizeof(int));
                        }
                        send_lists[owner][cnt] = col;
                        send_counts[owner]     = cnt + 1;
                    }
                }
            }
            /* 去重 */
            for (int p = 0; p < nprocs; ++p) {
                if (send_counts[p] == 0) continue;
                qsort(send_lists[p], (size_t)send_counts[p], sizeof(int), int_cmp);
                int u = 1;
                for (int i = 1; i < send_counts[p]; ++i)
                    if (send_lists[p][i] != send_lists[p][i-1])
                        send_lists[p][u++] = send_lists[p][i];
                send_counts[p] = u;
            }

            /* ── Phase B: 交换请求数量，投递非阻塞通信 ── */
            int *recv_counts = (int*)malloc((size_t)nprocs * sizeof(int));
            MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);

            /* 接收"别人向我请求哪些列" */
            MPI_Request *req_idx_recv = (MPI_Request*)malloc((size_t)nprocs * sizeof(MPI_Request));
            int **recv_req_lists      = (int**)malloc((size_t)nprocs * sizeof(int*));
            for (int p = 0; p < nprocs; ++p) {
                if (recv_counts[p] > 0) {
                    recv_req_lists[p] = (int*)malloc((size_t)recv_counts[p] * sizeof(int));
                    /* 【修复问题7】count>0 才 Irecv */
                    MPI_Irecv(recv_req_lists[p], recv_counts[p], MPI_INT,
                              p, 100, MPI_COMM_WORLD, &req_idx_recv[p]);
                } else {
                    recv_req_lists[p] = NULL;
                    req_idx_recv[p]   = MPI_REQUEST_NULL;   /* 【修复问题7】 */
                }
            }
            /* 发送"我需要哪些列" */
            MPI_Request *req_idx_send = (MPI_Request*)malloc((size_t)nprocs * sizeof(MPI_Request));
            for (int p = 0; p < nprocs; ++p) {
                if (send_counts[p] > 0)
                    MPI_Isend(send_lists[p], send_counts[p], MPI_INT,
                              p, 100, MPI_COMM_WORLD, &req_idx_send[p]);
                else
                    req_idx_send[p] = MPI_REQUEST_NULL;
            }

            /* ── Phase C: 【修复问题5】本地计算先行，不计入 comm_time ──
             *
             * 【修复问题3】对角区间内的所有列（不论 owner 是谁）都已在 x_buf 中（初始或上轮更新），
             * 直接计算，无需等待通信。
             */
            double comp0 = MPI_Wtime();
            {
                int nthreads = omp_get_max_threads();
                int *thr_bound = (int*)malloc((nthreads + 1) * sizeof(int));
                compute_thread_boundaries(local_rowptr, local_nrows, nthreads, thr_bound);
                #pragma omp parallel num_threads(nthreads)
                {
                    int tid = omp_get_thread_num();
                    for (int lr = thr_bound[tid]; lr < thr_bound[tid+1]; ++lr) {
                        double s = 0.0;
                        for (int k = local_rowptr[lr]; k < local_rowptr[lr+1]; ++k) {
                            int col = local_colidx[k];
                            /* 【修复问题3】条件：col 在 balanced boundary 内 */
                            if (col >= leftBound && col <= rightBound)
                                s += local_vals[k] * x_buf[col];
                        }
                        y_local[lr] = s;
                    }
                }
                free(thr_bound);
            }
            double comp1 = MPI_Wtime();
            compute_time += (comp1 - comp0);

            /* ── Phase D: 等待索引交換完成，準備 x 值回复 ── */
            double comm0 = MPI_Wtime();   /* 【修复问题5】comm 计时从这里开始 */
            MPI_Waitall(nprocs, req_idx_recv, MPI_STATUSES_IGNORE);
            MPI_Waitall(nprocs, req_idx_send, MPI_STATUSES_IGNORE);

            /* 发送别人请求的 x 值 */
            MPI_Request *req_val_send = (MPI_Request*)malloc((size_t)nprocs * sizeof(MPI_Request));
            double      **reply_bufs  = (double**)malloc((size_t)nprocs * sizeof(double*));
            for (int p = 0; p < nprocs; ++p) {
                if (recv_counts[p] > 0) {
                    reply_bufs[p] = (double*)malloc((size_t)recv_counts[p] * sizeof(double));
                    for (int i = 0; i < recv_counts[p]; ++i) {
                        int c = recv_req_lists[p][i];
                        reply_bufs[p][i] = (c >= cstart && c < cend) ? x_owned[c - cstart] : 0.0;
                    }
                    MPI_Isend(reply_bufs[p], recv_counts[p], MPI_DOUBLE,
                              p, 200, MPI_COMM_WORLD, &req_val_send[p]);
                } else {
                    reply_bufs[p]   = NULL;
                    req_val_send[p] = MPI_REQUEST_NULL;
                }
            }
            /* 接收我请求的 x 值 */
            MPI_Request *req_val_recv = (MPI_Request*)malloc((size_t)nprocs * sizeof(MPI_Request));
            double      **recv_val_bufs = (double**)malloc((size_t)nprocs * sizeof(double*));
            for (int p = 0; p < nprocs; ++p) {
                if (send_counts[p] > 0) {
                    recv_val_bufs[p] = (double*)malloc((size_t)send_counts[p] * sizeof(double));
                    MPI_Irecv(recv_val_bufs[p], send_counts[p], MPI_DOUBLE,
                              p, 200, MPI_COMM_WORLD, &req_val_recv[p]);
                } else {
                    recv_val_bufs[p] = NULL;
                    req_val_recv[p]  = MPI_REQUEST_NULL;
                }
            }

            MPI_Waitall(nprocs, req_val_recv, MPI_STATUSES_IGNORE);
            MPI_Waitall(nprocs, req_val_send, MPI_STATUSES_IGNORE);
            double comm1 = MPI_Wtime();
            comm_time += (comm1 - comm0);   /* 【修复问题5】comm 计时到这里结束 */

            /* 解包远端 x 值到 x_buf */
            for (int p = 0; p < nprocs; ++p) {
                for (int i = 0; i < send_counts[p]; ++i)
                    x_buf[send_lists[p][i]] = recv_val_bufs[p][i];
            }

            /* ── Phase E: 计算远端部分（boundary 外的列）──
             * 【修复问题6】条件与 Phase C 互斥：col NOT in [leftBound, rightBound]
             */
            double comp2 = MPI_Wtime();
            {
                int nthreads = omp_get_max_threads();
                int *thr_bound = (int*)malloc((nthreads + 1) * sizeof(int));
                compute_thread_boundaries(local_rowptr, local_nrows, nthreads, thr_bound);
                #pragma omp parallel num_threads(nthreads)
                {
                    int tid = omp_get_thread_num();
                    for (int lr = thr_bound[tid]; lr < thr_bound[tid+1]; ++lr) {
                        double s = y_local[lr];
                        for (int k = local_rowptr[lr]; k < local_rowptr[lr+1]; ++k) {
                            int col = local_colidx[k];
                            /* 【修复问题6】严格与 Phase C 互斥 */
                            if (col < leftBound || col > rightBound)
                                s += local_vals[k] * x_buf[col];
                        }
                        y_local[lr] = s;
                    }
                }
                free(thr_bound);
            }
            double comp3 = MPI_Wtime();
            compute_time += (comp3 - comp2);

            /* ── 释放本轮临时缓冲（保留 send_lists/send_caps/send_counts 以供复用） ── */
            for (int p = 0; p < nprocs; ++p) {
                if (recv_req_lists[p]) free(recv_req_lists[p]);
                if (reply_bufs[p])     free(reply_bufs[p]);
                if (recv_val_bufs[p])  free(recv_val_bufs[p]);
            }
            free(recv_counts);
            free(req_idx_recv); free(req_idx_send);
            free(req_val_send); free(req_val_recv);
            free(reply_bufs);   free(recv_val_bufs);
            free(recv_req_lists);
        }
        /* 迭代结束后释放复用的 send_lists/send_caps/send_counts */
        for (int p = 0; p < nprocs; ++p) free(send_lists[p]);
        free(send_lists); free(send_caps); free(send_counts);
    }
    /* 结束计时与汇总结果 */
    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double total_time = t1 - t0;

    double max_compute, max_comm, max_total;
    MPI_Reduce(&compute_time, &max_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time,    &max_comm,    1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time,   &max_total,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double gflops = (2.0 * (double)global_nnz * (double)niter) / (max_compute * 1.0e9);
        printf("== DistSpMV_Balanced (fixed) ==\n== DistSpMV_Balanced（修复版）==\n");
        printf("Processes   : %d  |  Threads/proc: %d\n", nprocs, omp_get_max_threads());
        printf("进程数：%d  |  每进程线程数：%d\n", nprocs, omp_get_max_threads());
        printf("Matrix      : %d x %d, nnz=%lld\n", nrows, ncols, global_nnz);
        printf("矩阵：%d x %d，非零元数：%lld\n", nrows, ncols, global_nnz);
        printf("Iterations  : %d\n", niter);
        printf("迭代次数：%d\n", niter);
        printf("lower_bound (frac) : %.2f\n", lower_bound_frac);
        printf("lower_bound（阈值比例）：%.2f\n", lower_bound_frac);
        printf("Total time  : %.6f s\n", max_total);
        printf("总时间：%.6f 秒\n", max_total);
        printf("Compute time: %.6f s\n", max_compute);
        printf("计算时间：%.6f 秒\n", max_compute);
        printf("Comm time   : %.6f s\n", max_comm);
        printf("通信时间：%.6f 秒\n", max_comm);
        printf("GFlops      : %.4f\n",   gflops);
        printf("GFlops（理论）：%.4f\n", gflops);
    }

    /* 计算并输出最终 y 向量的全局范数，便于与 naive 模式对比验证 */
    double local_norm_sq = 0.0;
    for (int i = 0; i < local_nrows; ++i) local_norm_sq += y_local[i] * y_local[i];
    double global_norm_sq = 0.0;
    MPI_Reduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        double y_norm = sqrt(global_norm_sq);
        printf("Y-norm: %.12e\n", y_norm);
    }

    /* ── 清理 ── */
    free(local_rowptr); free(local_colidx); free(local_vals);
    free(x_owned); free(x_buf); free(y_local);
    free(owner_of_col); free(row_offsets); free(col_offsets);

    MPI_Finalize();
    return 0;
}
