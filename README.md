# Metagenomic Analysis of Human Gut Microbiome in Fecal Microbiota Transplantation (FMT) for Weight Loss in Overweight and Obese Individuals
---

## 1. Introduction
The global prevalence of overweight and obesity continues to rise, posing significant challenges for public health owing to their association with chronic metabolic disorders and cardiovascular diseases. Emerging evidence implicates the gut microbiome as a key modulator of host energy metabolism, nutrient harvesting, and adiposity regulation. The intestinal microbial community influences fermentation of dietary fibres to short-chain fatty acids (SCFAs), bile acid metabolism, and lipid biosynthesis—pathways that are dysregulated in obesity [1]. For instance, obese individuals often show reduced microbial gene richness and altered functional capacity of microbiomes compared to lean counterparts [2].

Fecal microbiota transplantation (FMT) has been successfully employed in treating recurrent *Clostridioides difficile* infection and is now under investigation as a therapeutic tool for metabolic diseases. In animal models, transplantation of microbiota from lean donors into obese recipients has been shown to influence adiposity and metabolic profiles, suggesting causal roles for gut microbiota in weight regulation [3]. However, human studies of FMT in obesity and metabolic syndrome report mixed outcomes, despite modifications in gut microbial composition and bile acid profiles; consistent weight loss in recipients remains elusive. This highlights the need to better understand which microbial functions (rather than just taxonomic changes) are transferred or engrafted and how these relate to clinical outcomes such as weight loss.

In this project, a comprehensive metagenomic workflow—from raw sequencing reads through taxonomic analysis, gene prediction, and functional annotation— was applied to compare gut microbiome functional capacity among FMT donors, recipient responders, and non-responders. The aim was to elucidate the microbial mechanisms that underlie successful FMT-mediated weight loss in overweight and obese individuals. 

---
## 2. Methods

The workflow for this metagenomic analysis comprised several key steps, from data retrieval and preprocessing to taxonomic and functional characterization. All preprocessing steps—from quality control to BIOM file generation—were performed on the UseGalaxy platform, while downstream statistical and visualization analyses were conducted in R. The workflow consisted of:

1. Data Retrieval and Curation
2. Quality Assessment and Preprocessing
3. Taxonomic Profiling and Visualization
4. Functional Annotation and Pathway Analysis
5. Diversity Analysis and Data Visualization

### **Summary of Methods**

| **Step**                                          | **Purpose**                                                                                                                                                                                                                                                                                                             |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1. Data Retrieval and Curation                | Raw metagenomic datasets of stool samples from FMT donors and recipients (responders and non-responders) were retrieved from the NCBI Sequence Read Archive (SRA) and curated for downstream analysis.                                                                                                                  |
| 2. Quality Assessment and Preprocessing      | Raw reads were assessed using *FastQC* to evaluate sequence quality, followed by *Trimmomatic* for adapter removal and trimming of low-quality bases, ensuring high-quality reads for accurate downstream analyses.                                                                                                 |
| 3. Taxonomic Profiling and Visualization      | Taxonomic classification was performed using *Kraken2* to identify microbial taxa present in each sample. Results were visualized using *Pavian* and *Krona* (not shown), and outputs were converted into *BIOM format* for use in R-based analyses.                                                                        |
| 4. Functional Annotation and Pathway Analysis | Gene prediction and functional annotation were performed using *Prokka* and *eggNOG-mapper* to assign predicted genes to *COG* and *KEGG* functional categories, enabling exploration of metabolic pathways related to SCFA production, bile acid metabolism, carbohydrate utilization, and lipid biosynthesis. |
| 5. Diversity Analysis and Visualization (R)   | BIOM and annotation files were imported into R. *phyloseq* was used to compute alpha and beta diversity metrics and integrate taxonomic and functional data. Visualizations were generated with *ggplot2* to compare microbial composition and functional potential between donors, responders, and non-responders. |

---


**References**

1. Ribeiro, G., Schellekens, H., Cuesta-Marti, C., Maneschy, I., Ismael, S., Cuevas-Sierra, A., ... & Calhau, C. (2025). A menu for microbes: unraveling appetite regulation and weight dynamics through the microbiota-brain connection across the lifespan. *American Journal of Physiology-Gastrointestinal and Liver Physiology*, 328(3), G206-G228.
2. Martínez-Álvaro, M., Zubiri-Gaitán, A., Hernández, P., Greenacre, M., Ferrer, A., & Blasco, A. (2021). Comprehensive functional core microbiome comparison in genetically obese and lean hosts under the same environment. *Communications Biology*, 4(1), 1246.
3. Lee, P., Yacyshyn, B. R., & Yacyshyn, M. B. (2019). Gut microbiota and obesity: An opportunity to alter obesity through faecal microbiota transplant (FMT). *Diabetes, Obesity and Metabolism*, 21(3), 479-490.

