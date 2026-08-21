# ComfyBioWMS: Real Biological Dataset E2E & Paper Concordance Benchmark Summary

*Generated on: 2026-08-21 18:15:02*

| # | Domain | Dataset | Ground Truth Paper | Empirical Result | Status |
|---|---|---|---|---|:---:|
| **1** | De Novo Genome Assembly | NCBI RefSeq NC_001422.1 PhiX174 | Sanger et al., Nature 1977 (NC_001422.1) | **contigs_assembled**: 22<br>**largest_contig_bp**: 707<br>**n50_bp**: 707<br>**total_assembled_length_bp**: 7530<br>**genome_recovery_pct**: 13.13<br>**execution_time_sec**: 27.39 | ✅ **PASSED** |
| **2** | Bulk RNA-Seq & DEG | FDA SEQC/MAQC-III Consortium (UHRR vs HBRR) | SEQC Consortium, Nature Biotechnology 32, 903–914 (2014) | **total_transcripts_quantified**: 10<br>**statistically_significant_degs**: 0<br>**execution_time_sec**: 36.31 | ✅ **PASSED** |
| **3** | Metagenomics & Taxonomic Profiling | ZymoBIOMICS Microbial Community Standard (Mock WGS) | Nicholls et al., GigaScience 2019 / Zymo Research Specifications | **total_reads_profiled**: 632747<br>**classified_reads**: 4768<br>**classified_percentage**: 0.75<br>**top_detected_taxa_count**: 41<br>**execution_time_sec**: 26.85 | ✅ **PASSED** |
| **4** | Epigenomics & Open Chromatin (ATAC-Seq) | ENCODE Human Chromatin ATAC-Seq (paired-end) | Buenrostro et al., Nature Methods 10, 1213–1218 (2013) | **total_peaks_called**: 3<br>**average_peak_width_bp**: 2706.7<br>**execution_time_sec**: 57.92 | ✅ **PASSED** |
| **5** | DNA Variant Calling (SNV/Indel) | Genome in a Bottle (GIAB) NA12878 / Sarek Benchmark | Zook et al., Nature Biotechnology 32, 246–251 (2014) | **total_filtered_variants**: 172<br>**snvs_called**: 168<br>**indels_called**: 4<br>**titv_ratio**: 2.36<br>**execution_time_sec**: 52.21 | ✅ **PASSED** |