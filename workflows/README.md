# ComfyBIOWMS Example Workflows

Ready-to-use ComfyUI workflow JSON files for various bioinformatics and computational biology pipelines.

## Available Workflows

1. **Bulk RNA-Seq Pipeline (`01_bulk_rnaseq_deseq2.json`)**
   - **Stages**: Sample QC (`FastpQCNode`) ➔ Trimming (`FastpTrimNode`) ➔ Indexing (`SalmonIndexNode`) ➔ Quantification (`SalmonQuantNode`) ➔ Count Matrix Generation (`TximportNode`) ➔ Differential Expression (`DESeq2AnalysisNode`) ➔ Visualization (`DESeq2VisualizationNode`) ➔ Image Preview.
   - **Output**: MultiQC/Fastp HTML, Transcripts Count Matrix, DESeq2 DEG CSV, Volcano Plot & PCA Plot.

2. **DNA Variant Calling Pipeline (`02_variant_calling_bcftools.json`)**
   - **Stages**: Input Validation ➔ BWA-MEM2 Indexing ➔ Read Alignment ➔ Duplicate Marking (`MarkDuplicatesNode`) ➔ Variant Calling (`BcftoolsCallNode`) ➔ Filtering (`BcftoolsFilterNode`) ➔ Visualization & Report.
   - **Output**: Sorted/Indexed BAM, VCF / Filtered VCF, Variant distribution plots.

3. **Publication Visualizers (`03_publication_visualizers.json`)**
   - **Stages**: Standalone visualization nodes: Volcano Plot, Clustermap Heatmap, UMAP Scatter Plot.
   - **Output**: 300+ DPI publication-ready figures output directly to ComfyUI's canvas and previewers.

## How to Load in ComfyUI
1. Open ComfyUI in your browser (`http://127.0.0.1:8188`).
2. Drag and drop any `.json` file from this directory onto the ComfyUI canvas, or click **Load** in the ComfyUI side panel.
3. Verify or adjust the input file paths in the nodes' widget fields.
4. Click **Queue Prompt** to execute the pipeline!
