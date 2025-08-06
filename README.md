# Brain-cell free DNA Analysis
The analysis of plasma cfDNA is an emerging approach that allows to track cell death events in a minimally invasive way, also from inaccessible areas of the body like brain. DNA methylation (DNAm) profiles can be used to map the tissue of origin of cfDNA and to identify molecules released from brain upon cell death. 

This project contains workflows that our laboratory developed to identify loci with brain-specific DNAm profiles that could be detected in circulating cfDNA with the aim to assess early brain cell death events in neurological clinical conditions.
Below the overall workflow is depicted:

<img width="284" alt="image" src="https://github.com/user-attachments/assets/ffc7d786-98fe-4ef2-92fc-a8f0d5f133fd">


# In silico selection

**cfDNA_DNAm_Differential_Analysis : Dataset_Analysis**

First, Infinium 450K probes were annotated using the Core 15-state model elaborated by the NIH Roadmap Epigenomics Consortium. Annotated probes are stored in annotations450k.RData.gz. This workflow compares the data from brain regions with those from non-brain tissues and primary cell lines obtained from online repositories and selects CpG sites with i) a unique chromatine conformation in the brain and ii) with a brain-specific DNA methylation profile

# Experimental Validation

**Target bisulfite sequencing**

Target-bisulfite sequencing data were analysed with AmpliMethProfiler tool. Please refer to Scala et al. for further information on AmpliMethProfiler bioinformatic pipeline. Alternatively, please check on GitHub Link for a modified version of AmpliMethProfiler tool (10.1186/s12859-016-1380-3). Please refer to repository indicated by Scala et al. at: https://sourceforge.net/projects/amplimethprofiler/ for base guideline. Alternatively, refer to https://github.com/LabBrainAgeing/AmpliMethProfiler_RepetitiveElements for a modified AmpliMethProfiler pipeline. 

**cfDNA_AmpliMeth_Analysis_GTeX**

This workflows is designed to validate brain-specific DNAm profiles of selected target regions by analysing their DNAm via targeted bisulfite sequencing in a we used DNA samples from tissues collected in the framework of the Genotype-Tissue Expression (GTEx) project.

**cfDNA_CAR-T_Analysis**

This workflow validates the ability of selected targets to detect brain-cfDNA in the context of patients with neurotoxicity syndrome following CD19.CAR T cell therapy.

