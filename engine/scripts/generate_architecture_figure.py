import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_figure(out_png='figures/fig1_system_architecture.png', out_jpg='figures/fig1_system_architecture.jpg'):
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)

    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

    fig_w, fig_h = 8.5, 10.3
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    def draw_rounded_box(x, y, w, h, title, subtitle=None, bg='#ffffff', border='#000000', lw=2.0, radius=0.035, title_fs=12.5, sub_fs=9.5):
        box = patches.FancyBboxPatch(
            (x - w/2.0, y - h/2.0), w, h,
            boxstyle=f"round,pad=0.0,rounding_size={radius}",
            facecolor=bg, edgecolor=border, linewidth=lw, zorder=3
        )
        ax.add_patch(box)
        if subtitle:
            ax.text(x, y + h*0.14, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
            ax.text(x, y - h*0.18, subtitle, ha='center', va='center', fontsize=sub_fs, color='#334155', zorder=4)
        else:
            ax.text(x, y, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
        return box

    def draw_rect_box(x, y, w, h, title, subtitle=None, bg='#ffffff', border='#000000', lw=2.0, title_fs=10.0, sub_fs=7.8):
        box = patches.Rectangle(
            (x - w/2.0, y - h/2.0), w, h,
            facecolor=bg, edgecolor=border, linewidth=lw, zorder=3
        )
        ax.add_patch(box)
        if subtitle:
            ax.text(x, y + h*0.16, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
            ax.text(x, y - h*0.18, subtitle, ha='center', va='center', fontsize=sub_fs, color='#334155', zorder=4)
        else:
            ax.text(x, y, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
        return box

    def draw_arrow(x1, y1, x2, y2, color='#000000', lw=1.8):
        ax.annotate(
            '', xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle='-|>',
                color=color,
                lw=lw,
                mutation_scale=15,
            ),
            zorder=5
        )

    # 1. Top Shape: Input (Rounded rectangle)
    y_input = 0.940
    w_input, h_input = 0.88, 0.070
    draw_rounded_box(0.50, y_input, w_input, h_input,
                     "Biological Data / Analysis Parameters",
                     "(FASTQ, BAM, FASTA, Metadata CSV, Filtering Cutoffs)",
                     lw=2.0, title_fs=12.5, sub_fs=9.2)

    # 2. Outer Dashed Box
    dash_x, dash_y = 0.02, 0.125
    dash_w, dash_h = 0.96, 0.770
    dash_box = patches.Rectangle(
        (dash_x, dash_y), dash_w, dash_h,
        facecolor='none', edgecolor='#000000', linestyle='--', linewidth=1.5, zorder=1
    )
    ax.add_patch(dash_box)

    # 3. Top row: Class 1, Class 2, Class 3
    y_row1 = 0.785
    h_row1 = 0.100
    w_class = 0.275

    x_c1 = 0.185
    x_c2 = 0.500
    x_c3 = 0.815

    # Class 1
    draw_rect_box(x_c1, y_row1, w_class, h_row1,
                  "Class 1: Pure Python Core",
                  "Path Validation, Metadata Parsing\n& Graph Routing",
                  lw=2.0, title_fs=9.8, sub_fs=7.8)

    # Class 2
    draw_rect_box(x_c2, y_row1, w_class, h_row1,
                  "Class 2: Direct Python Libs",
                  "In-Memory Library Operations\n(Biopython, Scanpy, Matplotlib)",
                  lw=2.0, title_fs=9.8, sub_fs=7.8)

    # Class 3
    draw_rect_box(x_c3, y_row1, w_class, h_row1,
                  "Class 3: Isolated Binary CLI",
                  "External Multi-Omics Binaries\n(BWA, SPAdes, MACS3, Kraken2)",
                  lw=2.0, title_fs=9.8, sub_fs=7.8)

    # Arrows from Top Input to Class 1, 2, 3
    y_split = 0.880
    ax.plot([0.50, 0.50], [y_input - h_input/2.0, y_split], color='#000000', lw=1.8, zorder=4)
    ax.plot([x_c1, x_c3], [y_split, y_split], color='#000000', lw=1.8, zorder=4)
    draw_arrow(x_c1, y_split, x_c1, y_row1 + h_row1/2.0)
    draw_arrow(x_c2, y_split, x_c2, y_row1 + h_row1/2.0)
    draw_arrow(x_c3, y_split, x_c3, y_row1 + h_row1/2.0)

    # 4. Conda Runner (Middle box below Class 2 & 3)
    y_conda = 0.590
    w_mid = 0.590
    h_conda = 0.100
    x_mid = 0.6575
    draw_rect_box(x_mid, y_conda, w_mid, h_conda,
                  "CondaRunner Execution Bridge",
                  "Command Tokenization (shlex), Safe Subprocess Dispatch\n& Unicode Path Normalization",
                  lw=2.0, title_fs=11.0, sub_fs=8.3)

    # Merge line from Class 2 and Class 3 into Conda Runner
    y_merge = 0.690
    ax.plot([x_c2, x_c2], [y_row1 - h_row1/2.0, y_merge], color='#000000', lw=1.8, zorder=4)
    ax.plot([x_c3, x_c3], [y_row1 - h_row1/2.0, y_merge], color='#000000', lw=1.8, zorder=4)
    ax.plot([x_c2, x_c3], [y_merge, y_merge], color='#000000', lw=1.8, zorder=4)
    draw_arrow(x_mid, y_merge, x_mid, y_conda + h_conda/2.0)

    # 5. CLI Subprocess
    y_cli = 0.420
    h_cli = 0.100
    draw_rect_box(x_mid, y_cli, w_mid, h_cli,
                  "CLI Subprocess Execution",
                  "12 Isolated Conda Environments, Process Containment\n& Exit Code Verification (exit code = 0)",
                  lw=2.0, title_fs=11.0, sub_fs=8.3)

    # Arrow Conda Runner -> CLI Subprocess
    draw_arrow(x_mid, y_conda - h_conda/2.0, x_mid, y_cli + h_cli/2.0)

    # 6. Data Management & Memory Protection (Wide box at bottom of dashed area)
    y_dm = 0.210
    w_dm = 0.880
    h_dm = 0.100
    draw_rect_box(0.50, y_dm, w_dm, h_dm,
                  "Data Management & Memory Protection",
                  "Zero-Memory-Bloat STRING Path Streaming (NVMe Storage)\n& Native PyTorch IMAGE Tensor Generation (1, H, W, 3) for Real-Time Canvas Preview",
                  lw=2.0, title_fs=11.2, sub_fs=8.4)

    # Arrow from Class 1 directly down to Data Management
    draw_arrow(x_c1, y_row1 - h_row1/2.0, x_c1, y_dm + h_dm/2.0)

    # Arrow from CLI Subprocess down to Data Management
    draw_arrow(x_mid, y_cli - h_cli/2.0, x_mid, y_dm + h_dm/2.0)

    # 7. Bottom Shape: Analysis Output (Rounded rectangle)
    y_output = 0.055
    w_output, h_output = 0.88, 0.070
    draw_rounded_box(0.50, y_output, w_output, h_output,
                     "Analysis Outputs & Interactive Previews",
                     "(Header-Stripped SHA-256 Validated BAM, VCF, TSV & 300+ DPI Canvas Tensors)",
                     lw=2.0, title_fs=12.5, sub_fs=9.2)

    # Arrow from Data Management down to Analysis Output
    draw_arrow(0.50, y_dm - h_dm/2.0, 0.50, y_output + h_output/2.0)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight', pad_inches=0.04, facecolor='#ffffff')
    plt.savefig(out_jpg, dpi=dpi, bbox_inches='tight', pad_inches=0.04, facecolor='#ffffff')
    plt.close()
    print(f"Generated Figure 1 matching user's clean layout at {out_png}")

if __name__ == '__main__':
    generate_figure()
