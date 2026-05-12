# Spatial Profiling of HBV and HDV Infection in Human Liver

This repository contains the analytical pipelines for single-cell and spatial transcriptomics profiling of HBV and HDV infection in human liver samples using NanoString CosMx SMI and scRNA-seq integration.


### Core Dependencies
* **R Environment:** v4.3.2
* **Python Environment:** v3.9+
* **Key Packages:** Seurat (v5.4), InSituType (v2.0), Scanpy (v1.10), scikit-learn (v1.3), Squidpy (v1.2.3), clusterProfiler (v4.10.1).

---

## Methods overview

* **QC & Cell Typing:** Cells underwent count/probe/area-based filtering. Supervised annotation into 28 refined cell types was performed using **InSituType** mapped against a scRNA-seq reference.
* **Neighborhood Analysis:** Local cellular niches were defined using proximity-weighted kNN graphs (k=40) and Leiden clustering (res=0.2). Niche biological interpretation was assisted by an LLM (Gemini) and manually curated.
* **DEGs & Pathways:** Differential expression was calculated via Seurat's `FindMarkers`, with pathway enrichment mapped to MSigDB via **clusterProfiler**.
* **Spatial Interactions:** Ligand-receptor mappings were modeled using **SCOTIA**. Spatial co-localization was quantified via **Squidpy** using a custom Spatial Interaction Enrichment Score (SIES) compared against 1000 random permutations.
* **Pseudotime Trajectories:** Transitions between spatial niches were modeled using **Slingshot**, and dynamically expressed genes along these paths were fitted with GAMs using **tradeSeq**.

## License

MIT license.