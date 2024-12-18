# CBE target finder for generation of early STOP codons in bacterial genes
This repository provides an R-based pipeline to identify suitable targets for dCas9-guided Cytosine Base Editors (CBEs) to generate gene knockouts in plant pathogenic bacteria. The pipeline was made to easily identify the most suitable targets for introduction of a premature STOP codon using the system described in the study: **"Efficient CRISPR-Cas9 based cytosine base editors  for phytopathogenic bacteria"** by Li et al., *Communications Biology* (2023) https://doi.org/10.1038/s42003-023-04451-8

## The pipeline

The R-based script identifies optimal CBE targets to generate cytosine (C) to uracil (U) conversions that will lead to a thymine (T) being substituted after DNA replication. This method can be used to generate early STOP codons in bacterial genes. The pipeline identifies:

1. The earliest STOP codon that can be introduced and the gRNA for this target
   - For multiplexing approaches, multiple targets can be generated
2. Off-targets on the genome through local BLAST
   - Off-target scores are given to gRNAs to determine viability of each target
  
## Running the pipeline

The repository contains the R script to run the analysis for *Xanthomonas campestris* pv. *campestris* 8004 (Xcc8004). This pipeline can be adapted to any bacterial genome, as long as .fasta files of all annotated coding sequences are available for the species. To run local rBLAST analysis of potential off-targets, a full genome sequence needs to be available in .fasta format.
  
## Exploring the use of dCas9-CBE to generate gene knockouts at the full genome level of plant pathogenic bacteria

The pipeline was applied to explore the viability of this method for generation of knockouts at the full genome level for several model strains of plant pathogenic bacteria:

- *Agrobacterium tumefaciens* C58 (AtC58)
- *Erwinia amylovora* CFBP1430 (EaCFBP1430)
- *Pseudomonas syringae* pv. *tomato* DC3000 (PstDC3000)
- *Ralstonia solanacearum* MGI1000 (RsMGI1000)
- *Xanthomonas campestris* pv. *campestris* 8004 (Xcc8004)
- *Xylella fastidiosa* Temecula1 (XfTemecula1)

# Figures

**Figure 1. Exploration of suitable CBE targets for effective gene knockouts in six phytopathogenic bacterial species.** **A)** Percent of genes in each genome that feature a suitable CBE target site. These sites consist of a Cas9 PAM site with a suitable codon (editable to STOP codon) in position -16 to -20 **B)** Distribution of GC-content between all genes per species. **C)** Distribution of suitable codon (that could be changed into a STOP codon by the CBE) frequency between all genes per species. **D)** Density plot showing the correlation between suitable codon frequency and fraction of protein remaining after introduction of a premature STOP codon. **E)** Density plots showing the correlation between GC-content of each gene and the suitable codon frequency  (*top panel*) and the fraction of protein remaining after introduction of a premature STOP codon (*bottom panel*).

![CBE_figure1](https://github.com/user-attachments/assets/d7d503d2-3cd4-47ec-a9de-de4b6ee7532a)

**Figure 2. Application of the dCas9-CBE system to study bacterial type III effectors (T3Es) .** **A)** Number of predicted T3E genes for each strain shown in Figure 1. *A. tumefaciens* and *X. fastidiosa* lack a functional type III secretion system (T3SS). **B)** Remaining fraction of T3E protein after earliest STOP codon introduction in each T3E gene. A considerable amount of T3E proteins can only be truncated, leaving a large portion of the resulting protein intact.

![CBE_figure2](https://github.com/user-attachments/assets/8b0060d2-cfa8-4dc2-87c8-59d83f606a64)





