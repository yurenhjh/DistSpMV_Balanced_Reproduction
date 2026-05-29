#!/usr/bin/env python3
# Generate simple slide images and a multipage PDF from preset content.
import os
import textwrap
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import image as mpimg

out_dir = 'docs'
slides_dir = os.path.join(out_dir, 'slides')
os.makedirs(slides_dir, exist_ok=True)

slides = []
slides.append({
    'title': 'DistSpMV_Balanced — 论文复现',
    'lines': [
        '作者: yurenhjh',
        '日期: 2026-05-28',
        '仓库: https://github.com/yurenhjh/DistSpMV_Balanced_Reproduction'
    ]
})
slides.append({
    'title': '问题与动机',
    'lines': [
        '分布式稀疏矩阵向量乘（SpMV）在大规模下常受远程访问 x 的通信瓶颈限制。',
        '目标：通过重排与对角块扩展，将非零尽量聚集到对角块，减少跨进程通信。'
    ]
})
slides.append({
    'title': '论文核心算法',
    'lines': [
        'Algorithm 1: boundary scan（对角块扩展）',
        'Algorithm 2/3: 两阶段 MPI 通信（索引请求 → 回复 x 值）',
        'Algorithm 4: 按 nnz 划分线程区间，降低线程不均衡' 
    ]
})
slides.append({
    'title': '实现关键点与工程修复',
    'lines': [
        'METIS 在 rank0 运行，perm 广播并在方阵时做列重排（对称矩阵重排 x）。',
        '避免 OOM：在预处理完成后释放全矩阵，rank0 可选用 Scatterv 分发。',
        '初始 x 用 MPI_Allgatherv 合并，send_lists 预分配复用减少 malloc/free。'
    ]
})
slides.append({
    'title': '关键代码位置（快速导航）',
    'lines': [
        'compute_balanced_boundary — dist_spmv_balanced.c（boundary scan）',
        'compute_thread_boundaries — dist_spmv_balanced.c（线程划分）',
        'METIS / perm 广播 — dist_spmv_balanced.c（预处理）',
        '初始 x 合并（Allgatherv）与 send_lists 复用（见文件注释）'
    ]
})
slides.append({
    'title': '自动化脚本与数据',
    'lines': [
        'scripts/collect_metrics.sh, scripts/verify_y_norm.sh, scripts/plot_results.py',
        '数据文件：data/results_perf.csv, data/cant_ablation_medians.csv',
        '已生成图：data/ablation_trend.png（在下一页显示）'
    ]
})
slides.append({
    'title': '关键实验结果（摘要）',
    'lines': [
        'ecology1 (np=4)：通信时间由 0.670629s → 0.249984s，减少 ≈ 62.7%。',
        'cant 消融（np=4 中位数）：frac=0.5/0.8/1.0 的 comm_median 如报告所示。'
    ],
    'embed': 'data/ablation_trend.png'
})
slides.append({
    'title': '现场演示命令（可选）',
    'lines': [
        'mpicc -O3 -fopenmp -o dist_spmv_balanced dist_spmv_balanced.c -lmetis -lm',
        'export OMP_NUM_THREADS=1',
        "mpirun --oversubscribe -np 4 ./dist_spmv_balanced matrices/cant.mtx 10 0.8 balanced | tee logs/demo_balanced.log",
        "mpirun --oversubscribe -np 4 ./dist_spmv_balanced matrices/cant.mtx 10 0.8 naive    | tee logs/demo_naive.log"
    ]
})
slides.append({
    'title': '限制与后续工作',
    'lines': [
        '单机 VM 上测量有噪声；需在真实多节点集群上复核结果。',
        '后续：实现 rank0 读入 + MPI_Scatterv 分发；比较 ParMETIS / 并行划分器；在集群上做稳健消融。'
    ]
})
slides.append({
    'title': '结束与问答',
    'lines': [
        '谢谢！欢迎提问。',
        '仓库与文档： https://github.com/yurenhjh/DistSpMV_Balanced_Reproduction'
    ]
})

# helper to draw

def draw_slide(title, lines, embed_path=None, fname=None):
    fig = plt.figure(figsize=(11, 8.5))
    fig.patch.set_facecolor('white')
    plt.axis('off')
    # Title
    plt.text(0.05, 0.92, title, fontsize=36, weight='bold', va='top')
    y = 0.80
    for line in lines:
        wrapped = textwrap.wrap(line, width=80)
        for w in wrapped:
            plt.text(0.06, y, u'• ' + w, fontsize=20, va='top')
            y -= 0.06
            if y < 0.05:
                break
        if y < 0.05:
            break
    # embed image on right if exists
    if embed_path and os.path.exists(embed_path):
        try:
            img = mpimg.imread(embed_path)
            ax = fig.add_axes([0.55, 0.18, 0.40, 0.55])
            ax.imshow(img)
            ax.axis('off')
        except Exception:
            pass
    if fname:
        fig.savefig(fname, dpi=150)
    return fig

pdf_path = os.path.join(out_dir, 'presentation.pdf')
png_paths = []
with PdfPages(pdf_path) as pdf:
    for i, s in enumerate(slides, start=1):
        png = os.path.join(slides_dir, f'slide_{i:02d}.png')
        fig = draw_slide(s['title'], s.get('lines', []), s.get('embed'), png)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)
        png_paths.append(png)

print('Wrote', pdf_path)
for p in png_paths:
    print('WROTE', p)
