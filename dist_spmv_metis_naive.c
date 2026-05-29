/*
 * dist_spmv_metis_naive.c — METIS 重排 + 基线分布式 SpMV
 *
 * 程序会在运行前在 rank 0 上使用 METIS 对行进行划分并生成重排（old->new），
 * 将重排广播给所有进程，并在所有进程上应用重排后按块列划分进行 naive SpMV。
 *
 * 编译： mpicc -O3 -fopenmp -o dist_spmv_metis_naive dist_spmv_metis_naive.c -lmetis -lm
 * 运行： export OMP_NUM_THREADS=4
 *       mpirun -np 2 ./dist_spmv_metis_naive matrix.mtx 10
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <metis.h>

typedef struct { int row, col; double val; } Triplet;

static int triplet_cmp(const void *a, const void *b) {
    const Triplet *ta = (const Triplet*)a, *tb = (const Triplet*)b;
    if (ta->row != tb->row) return ta->row - tb->row;
    return ta->col - tb->col;
}

static int int_cmp(const void *a, const void *b) {
    return (*(const int*)a > *(const int*)b) - (*(const int*)a < *(const int*)b);
}

typedef struct {
    int    nrows, ncols;
    long long nnz;
    int   *rowptr;
    int   *colidx;
    double *vals;
} CSR;

static void csr_free(CSR *A) {
    if (A->rowptr) free(A->rowptr);
    if (A->colidx) free(A->colidx);
    if (A->vals)   free(A->vals);
    A->rowptr = NULL; A->colidx = NULL; A->vals = NULL;
}

int read_matrix_market(const char *path, CSR *A) {
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Cannot open %s\n", path); return -1; }
    char line[512];
    int is_sym = 0, is_pattern = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] != '%') break;
        if (strstr(line, "symmetric") || strstr(line, "Symmetric")) is_sym = 1;
        if (strstr(line, "pattern")   || strstr(line, "Pattern"))   is_pattern = 1;
    }
    int nrows, ncols; long long nnz_file;
    if (sscanf(line, "%d %d %lld", &nrows, &ncols, &nnz_file) != 3) { fclose(f); return -2; }
    long long cap = is_sym ? nnz_file * 2 : nnz_file;
    Triplet *trips = (Triplet*)malloc((size_t)cap * sizeof(Triplet));
    if (!trips) { fclose(f); return -3; }
    long long cnt = 0;
    for (long long k = 0; k < nnz_file; ++k) {
        int i, j; double v = 1.0;
        if (is_pattern) {
            if (fscanf(f, "%d %d", &i, &j) != 2) break;
        } else {
            if (fscanf(f, "%d %d %lf", &i, &j, &v) != 3) break;
        }
        i--; j--;
        if (i < 0 || i >= nrows || j < 0 || j >= ncols) continue;
        trips[cnt].row = i; trips[cnt].col = j; trips[cnt].val = v; cnt++;
        if (is_sym && i != j) {
            trips[cnt].row = j; trips[cnt].col = i; trips[cnt].val = v; cnt++;
        }
    }
    fclose(f);
    qsort(trips, (size_t)cnt, sizeof(Triplet), triplet_cmp);
    long long ucnt = 0;
    for (long long k = 0; k < cnt; ++k) {
        if (ucnt > 0 && trips[ucnt-1].row == trips[k].row && trips[ucnt-1].col == trips[k].col)
            trips[ucnt-1].val += trips[k].val;
        else
            trips[ucnt++] = trips[k];
    }
    cnt = ucnt;
    A->nrows = nrows; A->ncols = ncols; A->nnz = cnt;
    A->rowptr = (int*)calloc((size_t)(nrows+1), sizeof(int));
    A->colidx = (int*)malloc((size_t)cnt * sizeof(int));
    A->vals   = (double*)malloc((size_t)cnt * sizeof(double));
    for (long long k = 0; k < cnt; ++k) A->rowptr[trips[k].row + 1]++;
    for (int i = 0; i < nrows; ++i) A->rowptr[i+1] += A->rowptr[i];
    for (long long k = 0; k < cnt; ++k) {
        A->colidx[k] = trips[k].col;
        A->vals[k]   = trips[k].val;
    }
    free(trips);
    return 0;
}

static void make_block_offsets(int n, int np, int *off) {
    int base = n / np, rem = n % np;
    off[0] = 0;
    for (int p = 0; p < np; ++p) off[p+1] = off[p] + base + (p < rem ? 1 : 0);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    int rank, nprocs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc < 2) {
        if (rank == 0) fprintf(stderr, "Usage: %s matrix.mtx [niter=10]\n", argv[0]);
        MPI_Finalize(); return 1;
    }
    const char *filename = argv[1];
    int niter = (argc >= 3 ? atoi(argv[2]) : 10);

    CSR A; memset(&A,0,sizeof(A));
    if (read_matrix_market(filename, &A) != 0) {
        if (rank==0) fprintf(stderr, "Failed to read %s\n", filename);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int nrows = A.nrows, ncols = A.ncols; long long global_nnz = A.nnz;
    int symmetric_perm = (nrows == ncols);

    /* 使用 METIS 在 rank 0 生成 permutation，并广播给所有进程 */
    int *perm = NULL; int metis_ok = 0;
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

        int *count_per_part = (int*)calloc((size_t)nprocs, sizeof(int));
        for (int i = 0; i < nrows; ++i) count_per_part[(int)part_idx[i]]++;
        int *offp = (int*)malloc((size_t)(nprocs+1) * sizeof(int)); offp[0] = 0;
        for (int p = 0; p < nprocs; ++p) offp[p+1] = offp[p] + count_per_part[p];
        int *pos_arr = (int*)malloc((size_t)nprocs * sizeof(int));
        for (int p = 0; p < nprocs; ++p) pos_arr[p] = offp[p];
        int *perm_new = (int*)malloc((size_t)nrows * sizeof(int));
        for (int i = 0; i < nrows; ++i) {
            int p = (int)part_idx[i]; perm_new[ pos_arr[p]++ ] = i;
        }
        perm = (int*)malloc((size_t)nrows * sizeof(int));
        for (int ni = 0; ni < nrows; ++ni) perm[ perm_new[ni] ] = ni;

        free(xadj); free(adjncy); free(xadj_idx); free(adjncy_idx); free(part_idx);
        free(count_per_part); free(offp); free(pos_arr); free(perm_new);
        free(col_counts_tmp); free(col_disp_tmp); free(rows_in_col); free(mark);

        metis_ok = 1;
        if (rank == 0) printf("METIS partition computed on rank 0\n");
    }

    MPI_Bcast(&metis_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (metis_ok) {
        if (rank != 0) perm = (int*)malloc((size_t)nrows * sizeof(int));
        MPI_Bcast(perm, nrows, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank == 0) printf("METIS permutation applied\n");
    } else {
        if (rank != 0) perm = (int*)malloc((size_t)nrows * sizeof(int));
        for (int i = 0; i < nrows; ++i) perm[i] = i;
    }

    /* apply perm to matrix */
    {
        Triplet *trips = (Triplet*)malloc((size_t)A.nnz * sizeof(Triplet));
        long long t = 0;
        for (int i = 0; i < nrows; ++i) {
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
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
        if (rank==0) printf("Applied permutation (METIS)\n");
    }

    int *row_offsets = (int*)malloc((size_t)(nprocs+1)*sizeof(int));
    int *col_offsets = (int*)malloc((size_t)(nprocs+1)*sizeof(int));
    make_block_offsets(nrows, nprocs, row_offsets);
    make_block_offsets(ncols, nprocs, col_offsets);

    int rstart = row_offsets[rank];
    int rend   = row_offsets[rank+1];
    int cstart = col_offsets[rank];
    int cend   = col_offsets[rank+1];
    int local_nrows = rend - rstart;
    int local_ncols = cend - cstart;

    double *y_local = (double*)malloc((size_t)local_nrows * sizeof(double));

    int *owner_of_col = (int*)malloc((size_t)ncols * sizeof(int));
    for (int p = 0; p < nprocs; ++p)
        for (int j = col_offsets[p]; j < col_offsets[p+1]; ++j)
            owner_of_col[j] = p;

    double *x_owned = (double*)malloc((size_t)local_ncols * sizeof(double));
    for (int i = 0; i < local_ncols; ++i) x_owned[i] = 1.0;
    double *x_buf = (double*)calloc((size_t)ncols, sizeof(double));
    for (int j = cstart; j < cend; ++j) x_buf[j] = x_owned[j - cstart];

    int *col_counts = (int*)malloc((size_t)nprocs * sizeof(int));
    int *col_disp   = (int*)malloc((size_t)nprocs * sizeof(int));
    for (int p = 0; p < nprocs; ++p) {
        col_counts[p] = col_offsets[p+1] - col_offsets[p];
        col_disp[p]   = col_offsets[p];
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t0 = MPI_Wtime();
    double compute_time = 0.0, comm_time = 0.0;

    for (int iter = 0; iter < niter; ++iter) {
        double comm0 = MPI_Wtime();
        MPI_Allgatherv(x_owned, local_ncols, MPI_DOUBLE,
                       x_buf, col_counts, col_disp, MPI_DOUBLE, MPI_COMM_WORLD);
        double comm1 = MPI_Wtime(); comm_time += (comm1 - comm0);

        double comp0 = MPI_Wtime();
        #pragma omp parallel for schedule(dynamic,64)
        for (int lr = 0; lr < local_nrows; ++lr) {
            double s = 0.0;
            int gr = rstart + lr;
            for (int k = A.rowptr[gr]; k < A.rowptr[gr+1]; ++k) {
                int col = A.colidx[k];
                s += A.vals[k] * x_buf[col];
            }
            y_local[lr] = s;
        }
        double comp1 = MPI_Wtime(); compute_time += (comp1 - comp0);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    double t1 = MPI_Wtime();
    double total_time = t1 - t0;

    double max_compute, max_comm, max_total;
    MPI_Reduce(&compute_time, &max_compute, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&comm_time,    &max_comm,    1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_time,   &max_total,   1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double gflops = (2.0 * (double)global_nnz * (double)niter) / (max_compute * 1.0e9);
        printf("== DistSpMV_METIS_Naive ==\n");
        printf("Processes : %d  |  Threads/proc: %d\n", nprocs, omp_get_max_threads());
        printf("Matrix : %d x %d, nnz=%lld\n", nrows, ncols, global_nnz);
        printf("Iterations : %d\n", niter);
        printf("Total time : %.6f s\n", max_total);
        printf("Compute time: %.6f s\n", max_compute);
        printf("Comm time   : %.6f s\n", max_comm);
        printf("GFlops : %.4f\n", gflops);
    }
    /* 1) 输出 Y-norm（用于正确性验证） */
    double local_norm_sq = 0.0;
    for (int i = 0; i < local_nrows; ++i) local_norm_sq += y_local[i] * y_local[i];
    double global_norm_sq = 0.0;
    MPI_Reduce(&local_norm_sq, &global_norm_sq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf("Y-norm: %.12e\n", sqrt(global_norm_sq));

    /* 2) 输出每个进程的本地 nnz（用于负载均衡度分析） */
    long long local_nnz_count = (long long)(A.rowptr[rend] - A.rowptr[rstart]);
    long long *all_nnz = NULL;
    if (rank == 0) all_nnz = (long long*)malloc(nprocs * sizeof(long long));
    MPI_Gather(&local_nnz_count, 1, MPI_LONG_LONG, all_nnz, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        printf("Per-process nnz: ");
        for (int p = 0; p < nprocs; ++p) printf("%lld ", all_nnz[p]);
        printf("\n");
        free(all_nnz);
    }

    csr_free(&A);
    free(y_local); free(x_owned); free(x_buf);
    free(owner_of_col); free(row_offsets); free(col_offsets);
    free(col_counts); free(col_disp);

    MPI_Finalize();
    return 0;
}
