import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def generate_figure(out_png='figures/fig1_system_architecture.png', out_jpg='figures/fig1_system_architecture.jpg'):
    os.makedirs(os.path.dirname(os.path.abspath(out_png)), exist_ok=True)

    # Use classic academic serif typography matching reference
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman', 'Times', 'DejaVu Serif']

    fig_w, fig_h = 9.2, 8.8
    dpi = 300

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    ax.set_facecolor('#ffffff')
    fig.patch.set_facecolor('#ffffff')

    def draw_rounded_box(x, y, w, h, title, subtitle=None, bg='#ffffff', border='#000000', lw=1.3, radius=0.030, title_fs=11.0, sub_fs=8.2):
        box = patches.FancyBboxPatch(
            (x - w/2.0, y - h/2.0), w, h,
            boxstyle=f"round,pad=0.0,rounding_size={radius}",
            facecolor=bg, edgecolor=border, linewidth=lw, zorder=3
        )
        ax.add_patch(box)
        if subtitle:
            ax.text(x, y + h*0.16, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
            ax.text(x, y - h*0.19, subtitle, ha='center', va='center', fontsize=sub_fs, style='normal', color='#1e293b', zorder=4)
        else:
            ax.text(x, y, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
        return box

    def draw_rect_box(x, y, w, h, title, subtitle=None, bg='#ffffff', border='#000000', lw=1.3, title_fs=10.4, sub_fs=7.8):
        box = patches.Rectangle(
            (x - w/2.0, y - h/2.0), w, h,
            facecolor=bg, edgecolor=border, linewidth=lw, zorder=3
        )
        ax.add_patch(box)
        if subtitle:
            ax.text(x, y + h*0.16, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
            ax.text(x, y - h*0.19, subtitle, ha='center', va='center', fontsize=sub_fs, style='normal', color='#1e293b', zorder=4)
        else:
            ax.text(x, y, title, ha='center', va='center', fontsize=title_fs, fontweight='bold', color='#000000', zorder=4)
        return box

    def draw_arrow(x1, y1, x2, y2, color='#000000', lw=1.2):
        ax.annotate(
            '', xy=(x2, y2), xytext=(x1, y1),
            arrowprops=dict(
                arrowstyle='-|>',
                color=color,
                lw=lw,
                mutation_scale=13,
            ),
            zorder=5
        )

    def draw_poly_arrow(points, color='#000000', lw=1.2):
        for i in range(len(points) - 2):
            ax.plot([points[i][0], points[i+1][0]], [points[i][1], points[i+1][1]], color=color, lw=lw, zorder=4)
        draw_arrow(points[-2][0], points[-2][1], points[-1][0], points[-1][1], color=color, lw=lw)

    # 1. Top box: Input (Rounded rectangle)
    y_input = 0.940
    w_input, h_input = 0.44, 0.058
    draw_rounded_box(0.50, y_input, w_input, h_input, 
                     "Input Assembly & Analysis Parameters", 
                     "(FASTQ, BAM, FASTA, Metadata CSV, Analysis Cutoffs)",
                     bg='#ffffff', border='#000000', lw=1.4, title_fs=10.8, sub_fs=8.0)

    # 2. Top Stage: ComfyUI DAG Engine (Wide rectangle)
    y_ui = 0.830
    w_ui, h_ui = 0.70, 0.070
    draw_rect_box(0.50, y_ui, w_ui, h_ui,
                  "Open-Source ComfyUI DAG Engine & Web UI (LiteGraph)",
                  "Interactive Node Graph Canvas, Parameter Widgets & Async Queue Scheduler",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=11.0, sub_fs=8.3)

    # Arrow Input -> UI
    draw_arrow(0.50, y_input - h_input/2.0, 0.50, y_ui + h_ui/2.0)

    # 3. Outer Dashed Bounding Box (ComfyBIOWMS Extension Architecture)
    dash_x, dash_y = 0.035, 0.165
    dash_w, dash_h = 0.930, 0.575
    dash_box = patches.Rectangle(
        (dash_x, dash_y), dash_w, dash_h,
        facecolor='none', edgecolor='#000000', linestyle='--', linewidth=1.1, zorder=1
    )
    ax.add_patch(dash_box)

    # 4. Three Upper Boxes inside Dashed Framework (Row 1) with gaps between them
    y_row1 = 0.635
    w_col = 0.220
    h_col = 0.095

    # Center box at 0.50 (range 0.39 to 0.61)
    # Left box at 0.25 (range 0.14 to 0.36) -> Gap = 0.03
    # Right box at 0.75 (range 0.64 to 0.86) -> Gap = 0.03
    x_left = 0.255
    x_mid = 0.500
    x_right = 0.745

    # Left: Class 1
    draw_rect_box(x_left, y_row1, w_col, h_col,
                  "Class 1: Pure Python Core",
                  "Path Validation, Metadata Parsing,\nPair Pattern Matching & Routing",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=9.6, sub_fs=7.6)

    # Center: CondaRunner Bridge
    draw_rect_box(x_mid, y_row1, w_col, h_col,
                  "CondaRunner Execution Bridge",
                  "Safe shlex Command Tokenization,\nPath Normalization & Interception",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=9.6, sub_fs=7.6)

    # Right: Class 2
    draw_rect_box(x_right, y_row1, w_col, h_col,
                  "Class 2: Direct Python Libs",
                  "In-Memory Analysis & Plotting\n(Biopython, Scanpy, Seaborn)",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=9.6, sub_fs=7.6)

    # Down arrow from ComfyUI directly into Class 1 (Left Box)
    draw_arrow(x_left, y_ui - h_ui/2.0, x_left, y_row1 + h_col/2.0)

    # Horizontal arrows row1: Left -> Center, Right -> Center across the gaps!
    draw_arrow(x_left + w_col/2.0, y_row1, x_mid - w_col/2.0, y_row1)
    draw_arrow(x_right - w_col/2.0, y_row1, x_mid + w_col/2.0, y_row1)

    # 5. Middle Box: Class 3 (Row 2 - Wide Box)
    y_row2 = 0.475
    w_row2 = 0.710
    h_row2 = 0.088
    draw_rect_box(0.50, y_row2, w_row2, h_row2,
                  "Class 3: Isolated Binary CLI Subprocesses (12 Conda Envs)",
                  "Dedicated Environment Isolation: BWA-MEM2, SPAdes, MACS3, Kraken2, DESeq2, AlphaFold2\nStrict Process Isolation, Exit Code Tracking (exit code = 0) & Non-ASCII Path Protection",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=10.2, sub_fs=7.8)

    # Center row1 -> Middle row2 down arrow
    draw_arrow(0.50, y_row1 - h_col/2.0, 0.50, y_row2 + h_row2/2.0)

    # Left Bypass Arrow: From Left edge of ComfyUI, through Left Gutter, into Class 3
    x_gutter_left = 0.075
    pt_bypass_left = [
        (0.50 - w_ui/2.0, y_ui),
        (x_gutter_left, y_ui),
        (x_gutter_left, y_row2),
        (0.50 - w_row2/2.0, y_row2)
    ]
    draw_poly_arrow(pt_bypass_left)
    ax.text(x_gutter_left + 0.012, (y_row1 + y_row2)/2.0, 
            "STRING\nFile Paths", ha='left', va='center', fontsize=7.6, fontweight='bold', color='#1e293b')

    # 6. Lower Box inside Dashed Framework (Row 3 - Wide Box)
    y_row3 = 0.315
    w_row3 = 0.710
    h_row3 = 0.088
    draw_rect_box(0.50, y_row3, w_row3, h_row3,
                  "Dual-Channel Data Management & Memory Protection",
                  "High-Throughput STRING Path Streaming (Zero Memory Bloat on NVMe Storage)\nNative PyTorch IMAGE Tensor Generation for Real-Time Canvas Preview (1, H, W, 3)",
                  bg='#ffffff', border='#000000', lw=1.3, title_fs=10.2, sub_fs=7.8)

    # Middle row2 -> Lower row3 down arrow
    draw_arrow(0.50, y_row2 - h_row2/2.0, 0.50, y_row3 + h_row3/2.0)

    # Right Bypass Arrow: From Class 2 bottom, through Right Gutter, into Dual-Channel Box
    x_gutter_right = 0.925
    pt_bypass_right = [
        (x_right, y_row1 - h_col/2.0),
        (x_right, y_row1 - h_col/2.0 - 0.035),
        (x_gutter_right, y_row1 - h_col/2.0 - 0.035),
        (x_gutter_right, y_row3),
        (0.50 + w_row3/2.0, y_row3)
    ]
    draw_poly_arrow(pt_bypass_right)
    ax.text(x_gutter_right - 0.012, (y_row2 + y_row3)/2.0, 
            "Direct Tensor\nRendering", ha='right', va='center', fontsize=7.6, fontweight='bold', color='#1e293b')

    # 7. Bottom Output Box (Rounded rectangle - Terminal)
    y_out = 0.075
    w_out = 0.640
    h_out = 0.066
    draw_rounded_box(0.50, y_out, w_out, h_out,
                     "Verified Analysis Outputs & Real-Time Canvas View",
                     "Header-Stripped SHA-256 Validated Artifacts (BAM, VCF, TSV) & 300+ DPI Interactive IMAGE Previews",
                     bg='#ffffff', border='#000000', lw=1.4, title_fs=10.5, sub_fs=7.8)

    # Exit Arrow from Lower Box to Bottom Output Box
    draw_arrow(0.50, y_row3 - h_row3/2.0, 0.50, y_out + h_out/2.0)

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')

    plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)
    plt.savefig(out_png, dpi=dpi, bbox_inches='tight', pad_inches=0.05, facecolor='#ffffff')
    plt.savefig(out_jpg, dpi=dpi, bbox_inches='tight', pad_inches=0.05, facecolor='#ffffff')
    plt.close()
    print(f"Regenerated Reference-Style Figure 1 at {out_png} and {out_jpg}")

if __name__ == '__main__':
    generate_figure()
