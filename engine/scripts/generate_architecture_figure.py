import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_figure(out_png='figures/fig1_system_architecture.png', out_jpg='figures/fig1_system_architecture.jpg'):
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    fig_w, fig_h = 8.8, 9.6
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    def draw_arrow(x1, y1, x2, y2, color='#000000', lw=1.8, style='-|>'):
        ax.annotate(
            '', xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle=style,
                color=color,
                lw=lw,
                mutation_scale=13,
            ),
            zorder=6
        )

    # -------------------------------------------------------------
    # 1. Top Shape: Biological Data & Analysis Parameters (Rounded Box)
    # -------------------------------------------------------------
    y_input = 0.932
    w_input, h_input = 0.880, 0.066
    r_box = patches.FancyBboxPatch(
        (0.50 - w_input/2.0, y_input - h_input/2.0), w_input, h_input,
        boxstyle="round,pad=0.0,rounding_size=0.025",
        facecolor='#ffffff', edgecolor='#000000', linewidth=2.0, zorder=3
    )
    ax.add_patch(r_box)
    ax.text(0.50, y_input + 0.013, "Biological Data & Analysis Parameters",
            ha='center', va='center', fontsize=11.2, fontweight='bold', color='#000000', zorder=4)
    ax.text(0.50, y_input - 0.014, "(Paired FASTQ, BAM, FASTA, H5AD AnnData, Metadata Samplesheet CSV, User Filtering Cutoffs)",
            ha='center', va='center', fontsize=8.2, color='#333333', zorder=4)

    # Arrow Input -> DAG
    y_dag = 0.814
    w_dag, h_dag = 0.880, 0.066
    draw_arrow(0.50, y_input - h_input/2.0, 0.50, y_dag + h_dag/2.0, color='#000000', lw=1.8)

    # -------------------------------------------------------------
    # 2. ComfyUI DAG Engine & Scheduling Layer
    # -------------------------------------------------------------
    dag_box = patches.Rectangle(
        (0.50 - w_dag/2.0, y_dag - h_dag/2.0), w_dag, h_dag,
        facecolor='#ffffff', edgecolor='#000000', linewidth=2.0, zorder=3
    )
    ax.add_patch(dag_box)
    ax.text(0.50, y_dag + 0.013, "ComfyUI Directed Acyclic Graph (DAG) Engine & Async Queue",
            ha='center', va='center', fontsize=10.6, fontweight='bold', color='#000000', zorder=4)
    ax.text(0.50, y_dag - 0.014, "Interactive Node Canvas (LiteGraph) • Topological Dependency Graph • Async Task Scheduling",
            ha='center', va='center', fontsize=8.0, color='#333333', zorder=4)

    # -------------------------------------------------------------
    # 3. Outer Execution Container (Dashed Box with Solid Header Banner)
    # -------------------------------------------------------------
    dash_x = 0.035
    dash_y = 0.160
    dash_w = 0.930
    dash_h = 0.550
    dash_box = patches.Rectangle(
        (dash_x, dash_y), dash_w, dash_h,
        facecolor='none', edgecolor='#000000', linestyle='--', linewidth=1.5, zorder=1
    )
    ax.add_patch(dash_box)

    # Subsystem Header Banner at top of outer container
    h_sub_hdr = 0.044
    y_sub_hdr = dash_y + dash_h - h_sub_hdr/2.0  # 0.688
    sub_hdr_box = patches.Rectangle(
        (dash_x, dash_y + dash_h - h_sub_hdr), dash_w, h_sub_hdr,
        facecolor='#ffffff', edgecolor='#000000', linewidth=1.8, zorder=3
    )
    ax.add_patch(sub_hdr_box)
    ax.text(0.50, y_sub_hdr,
            "ComfyBIOWMS Custom Node Execution Architecture (98 Active Registered Nodes)",
            ha='center', va='center', fontsize=9.6, fontweight='bold', color='#000000', zorder=4)

    # Arrow DAG -> Subsystem Header
    draw_arrow(0.50, y_dag - h_dag/2.0, 0.50, dash_y + dash_h, color='#000000', lw=1.8)

    # Column x positions and width
    x_c1 = 0.275
    x_c2 = 0.725
    w_col = 0.420

    # Arrows from Subsystem Header down into Class 1 and Class 2
    y_sub_bot = dash_y + dash_h - h_sub_hdr  # 0.666
    y_c_hdr = 0.598
    h_c_hdr = 0.058
    y_c_hdr_top = y_c_hdr + h_c_hdr/2.0     # 0.627
    draw_arrow(x_c1, y_sub_bot, x_c1, y_c_hdr_top, color='#000000', lw=1.6)
    draw_arrow(x_c2, y_sub_bot, x_c2, y_c_hdr_top, color='#000000', lw=1.6)

    # -------------------------------------------------------------
    # 4. Left Column: In-Memory Processing & Visualization
    # -------------------------------------------------------------
    hdr_c1 = patches.Rectangle(
        (x_c1 - w_col/2.0, y_c_hdr - h_c_hdr/2.0), w_col, h_c_hdr,
        facecolor='#ffffff', edgecolor='#000000', linewidth=1.8, zorder=3
    )
    ax.add_patch(hdr_c1)
    ax.text(x_c1, y_c_hdr + 0.012, "Python In-Memory Processing & Visualization",
            ha='center', va='center', fontsize=9.2, fontweight='bold', color='#000000', zorder=4)
    ax.text(x_c1, y_c_hdr - 0.013, "Direct In-Memory Execution • Low-Latency Scientific Analytics",
            ha='center', va='center', fontsize=7.4, color='#333333', zorder=4)

    # Arrow Header -> Main Box
    y_main = 0.354
    h_main = 0.344
    draw_arrow(x_c1, y_c_hdr - h_c_hdr/2.0, x_c1, y_main + h_main/2.0, color='#000000', lw=1.6)

    # Main Box 1: Analytics Engine
    box_c1_main = patches.Rectangle(
        (x_c1 - w_col/2.0, y_main - h_main/2.0), w_col, h_main,
        facecolor='#ffffff', edgecolor='#000000', linewidth=1.6, zorder=3
    )
    ax.add_patch(box_c1_main)
    x_c1_text = x_c1 - w_col/2.0 + 0.016
    ax.text(x_c1, y_main + h_main/2.0 - 0.024, "In-Memory Scientific Analytics & Visualization",
            ha='center', va='center', fontsize=9.3, fontweight='bold', color='#000000', zorder=4)

    bullets_c1 = [
        "• Direct In-Memory Arrays: Zero-IPC Python execution",
        "• Biopython: Alignment stats, GC metrics, pairwise search",
        "• Single-Cell: Scanpy, AnnData inspect, HVG, PCA, Leiden",
        "• Transcriptomics: Tximport count matrix aggregation",
        "• Interactive Visualization: 300+ DPI IMAGE canvas preview",
        "• Input Validation: Samplesheet CSV & path verification",
        "• In-Memory Passing: High-speed NumPy & Pandas flow",
        "• Ultra-Low Latency: Instantaneous node evaluation",
    ]
    y_start_c1 = y_main + h_main/2.0 - 0.056
    for i, b in enumerate(bullets_c1):
        ax.text(x_c1_text, y_start_c1 - i * 0.037, b,
                ha='left', va='center', fontsize=7.4, color='#1e293b', zorder=4)

    # -------------------------------------------------------------
    # 5. Right Column: Isolated Binary CLI & BioCommandRunner
    # -------------------------------------------------------------
    hdr_c2 = patches.Rectangle(
        (x_c2 - w_col/2.0, y_c_hdr - h_c_hdr/2.0), w_col, h_c_hdr,
        facecolor='#ffffff', edgecolor='#000000', linewidth=1.8, zorder=3
    )
    ax.add_patch(hdr_c2)
    ax.text(x_c2, y_c_hdr + 0.012, "External Tool CLI & Conda Isolation",
            ha='center', va='center', fontsize=9.2, fontweight='bold', color='#000000', zorder=4)
    ax.text(x_c2, y_c_hdr - 0.013, "Dedicated Conda Environments • Subprocess Containment",
            ha='center', va='center', fontsize=7.4, color='#333333', zorder=4)

    # Arrow Header -> Main Box
    draw_arrow(x_c2, y_c_hdr - h_c_hdr/2.0, x_c2, y_main + h_main/2.0, color='#000000', lw=1.6)

    # Main Box 2: BioCommandRunner Bridge & Conda Isolation
    box_c2_main = patches.Rectangle(
        (x_c2 - w_col/2.0, y_main - h_main/2.0), w_col, h_main,
        facecolor='#ffffff', edgecolor='#000000', linewidth=1.6, zorder=3
    )
    ax.add_patch(box_c2_main)
    x_c2_text = x_c2 - w_col/2.0 + 0.016
    ax.text(x_c2, y_main + h_main/2.0 - 0.024, "BioCommandRunner & Conda Execution Bridge",
            ha='center', va='center', fontsize=9.3, fontweight='bold', color='#000000', zorder=4)

    bullets_c2 = [
        "• Safe CLI Tokenization: shlex.split (shell=False execution)",
        "• Vulnerability Prevention: Command injection immunity",
        "• Dedicated Conda Envs: Fastp, STAR, Salmon, BWA-MEM2,",
        "  SPAdes, MACS3, DESeq2 (R 4.5), Samtools, BCFtools",
        "• Subprocess Dispatch: conda run -n [env] --no-capture-output",
        "• Process Control: Exit code check & live log stream",
        "• Audit Logging: Auto-generated run_manifest.sh & json",
        "• Provenance Replay: Independent CLI reproduction",
    ]
    y_start_c2 = y_main + h_main/2.0 - 0.056
    for i, b in enumerate(bullets_c2):
        ax.text(x_c2_text, y_start_c2 - i * 0.037, b,
                ha='left', va='center', fontsize=7.4, color='#1e293b', zorder=4)

    # -------------------------------------------------------------
    # 6. Bottom Shape: Analysis Outputs & Provenance Artifacts
    # -------------------------------------------------------------
    y_output = 0.052
    w_output, h_output = 0.880, 0.066
    out_box = patches.FancyBboxPatch(
        (0.50 - w_output/2.0, y_output - h_output/2.0), w_output, h_output,
        boxstyle="round,pad=0.0,rounding_size=0.025",
        facecolor='#ffffff', edgecolor='#000000', linewidth=2.0, zorder=3
    )
    ax.add_patch(out_box)
    ax.text(0.50, y_output + 0.013, "Analysis Outputs & Provenance Artifacts",
            ha='center', va='center', fontsize=11.2, fontweight='bold', color='#000000', zorder=4)
    ax.text(0.50, y_output - 0.014, "(Verified Artifacts: BAM, VCF, Counts TSV, run_manifest.sh • In-Memory AnnData, Tables & Reports)",
            ha='center', va='center', fontsize=8.2, color='#333333', zorder=4)

    # Arrows from Class 1 and Class 2 down into Output Box
    y_main_bot = y_main - h_main/2.0  # 0.182
    y_out_top = y_output + h_output/2.0  # 0.085
    draw_arrow(x_c1, y_main_bot, x_c1, y_out_top, color='#000000', lw=1.6)
    draw_arrow(x_c2, y_main_bot, x_c2, y_out_top, color='#000000', lw=1.6)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight', pad_inches=0.04, facecolor='#ffffff')
    plt.savefig(out_jpg, dpi=dpi, bbox_inches='tight', pad_inches=0.04, facecolor='#ffffff')
    plt.close()
    print(f"Generated Figure 1 (Streamlined Black Theme) at {out_png}")

if __name__ == '__main__':
    generate_figure()
