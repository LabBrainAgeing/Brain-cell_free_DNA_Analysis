# Brain-cell free DNA Analysis
The analysis of plasma cfDNA is an emerging approach that allows to track cell death events in a minimally invasive way, also from inaccessible areas of the body like brain. DNA methylation (DNAm) profiles can be used to map the tissue of origin of cfDNA and to identify molecules released from brain upon cell death. 


Our laboratory developed an experimental and bioinformatic pipeline to identify 

<img width="284" alt="image" src="https://github.com/user-attachments/assets/ffc7d786-98fe-4ef2-92fc-a8f0d5f133fd">

# Dataset Analysis

**cfDNA_DNAm_Differential_Analysis**

We used the Core 15-state model elaborated by the NIH Roadmap Epigenomics Con-sortium, in which 15 possible chromatin states are predicted on the basis of 5 histone modifications, and we compared the data from brain regions with those from non-brain tissues and primary cell lines and selected CpG sites with a unique chromatine conformation, specific to the brain.

# Experimenta Validation

**Target bisulfite sequencing**

**cfDNA_AmpliMeth_Analysis_GTeX**

We validated selected target regions by analysing their DNAm via targeted bisulfite sequencing in a we used DNA samples from 44 autoptic tissues isolated from 4 healthy subjects, available in the framework of the Genotype-Tissue Expression (GTEx) project. Importantly, this collection included different brain regions for the same subject, including cortex, hippocampus, cerebellum, amygdala, caudate nucleus, nucleus accumbens, putamen and substantia nigra. Sequencing data preprocessing and DNAm extraction was conducted according to AmpliMethProfiler Pipeline. Please refer to Scala et al. for further information on AmpliMethProfiler bioinformatic pipeline. Alternatively, please check on GitHub Link for a modified version of AmpliMethProfiler tool.

**cfDNA_CAR-T_Analysis**

To explore the ability of the assay to quantify bcfDNA, we applied it to plasma cfDNA samples extracted from 5 patients affected by relapsed/refractory B cell lymphoma who received CD19.CAR T cell therapy. At the time of plasma collection, 3 patients were in the acute phase of ICANS, an adverse event in CD19.CAR T therapy characterized by encephalopathy, seizures and cognitive impairment.

