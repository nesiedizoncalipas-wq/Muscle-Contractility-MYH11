# MYH11 Protein Structure and Mutation Analysis

## Project Overview

This project investigates the human **MYH11 (myosin heavy chain 11)** gene and its assigned splice-site variant using NCBI databases, R programming, and published scientific literature.

The analysis includes gene identification, retrieval of the protein sequence, amino acid composition analysis, ClinVar variant investigation, and a schematic illustrating the affected splice-donor site.

**Project title:** MYH11 Protein Structure and Mutation Analysis  
**Assignment number:** 03 (MYH11)
**Organism:** Homo sapiens  
**Assigned variant:** NM_002474.3:c.3858+1G>A

# MYH11 Protein Structure and Mutation Analysis

## A. Assignment Information

| Item | Details |
|---|---|
| Assignment number | 03 |
| Gene | MYH11 |
| Official gene name | Myosin heavy chain 11 |
| Organism | Homo sapiens |
| NCBI Gene ID | 4629 |
| Chromosome | 16 |
| Assigned transcript | NM_002474.3 |
| Assigned variant | NM_002474.3:c.3858+1G>A |
| ClinVar Variation ID | 1467492 |
| ClinVar accession | VCV001467492.9 |
| Protein accession | NP_002465.1 |

## B. Research Question

What is the amino acid composition of the MYH11 protein, and how might the assigned splice-site variant, NM_002474.3:c.3858+1G>A, affect normal RNA processing and smooth muscle function?

## C. Methods

The analysis was performed using R and RStudio. The NCBI database was used to retrieve the gene, transcript, and protein information. ClinVar was used to examine the assigned variant and its clinical classification.

### 1. Set up the project and install the packages

I created a main project folder named `Muscle_Contractility_MYH11` with three subfolders: `data`, `results`, and `figures`. I also created an R script named `analysis.R`.

I used the `rentrez` package to retrieve information from NCBI and `ggplot2` to create the figures.

For a fresh R installation, I first ran:

```r
install.packages(c("rentrez", "ggplot2"))
```

Then I loaded the packages and set the assignment information:

```r
library(rentrez)
library(ggplot2)

gene_symbol <- "MYH11"
transcript_accession <- "NM_002474.3"
clinvar_id <- "1467492"
variant_protein <- ""

for (x in c("data", "results", "figures")) {
  dir.create(x, showWarnings = FALSE)
}
```

The empty `variant_protein` value was used because ClinVar did not provide a protein-level change for this splice-site variant.

### 2. Retrieve the gene information from NCBI

I searched the NCBI Gene database for MYH11 in humans and recorded the Gene ID, official name, and chromosome.

The recorded results were saved as `results/gene_summary.txt`.

The information retrieved was:

- Gene symbol: MYH11
- NCBI Gene ID: 4629
- Official name: myosin heavy chain 11
- Chromosome: 16

### 3. Retrieve the assigned transcript and protein sequence

I used the assigned RefSeq transcript accession, NM_002474.3, to retrieve its GenBank record. I then extracted the corresponding NP protein accession and downloaded its FASTA sequence.

The following commands were used:

```r
gb <- entrez_fetch(
  db = "nuccore",
  id = transcript_accession,
  rettype = "gb",
  retmode = "text"
)

pat <- '/protein_id="NP_[0-9]+[.][0-9]+"'
fragment <- regmatches(gb, regexpr(pat, gb))

if (length(fragment) == 0) {
  stop("NP protein ID not found; inspect the transcript record")
}

protein_accession <- gsub(
  '^/protein_id="|"$',
  '',
  fragment
)

protein_fasta <- entrez_fetch(
  db = "protein",
  id = protein_accession,
  rettype = "fasta",
  retmode = "text"
)

if (!startsWith(protein_fasta, ">")) {
  stop("Invalid FASTA")
}

fasta_file <- file.path(
  "data",
  paste0(gene_symbol, "_protein.fasta")
)

writeLines(protein_fasta, fasta_file)

cat("Matched protein:", protein_accession, "\n")

writeLines(
  paste("NP accession:", protein_accession),
  "results/protein_accession.txt"
)
```

The retrieved protein accession was NP_002465.1. The FASTA file was saved as `data/MYH11_protein.fasta`, while the accession information was saved in `results/protein_accession.txt`.

### 4. Calculate the amino acid composition

I read the downloaded FASTA file, removed its header, and counted the standard amino acids in the sequence.

```r
fasta_lines <- readLines(fasta_file)

sequence <- paste(
  fasta_lines[!grepl("^>", fasta_lines)],
  collapse = ""
)

aa <- strsplit(toupper(sequence), "")[[1]]

protein_length <- length(aa)

standard <- strsplit(
  "ACDEFGHIKLMNPQRSTVWY",
  ""
)[[1]]

counts <- table(factor(aa, levels = standard))

results <- data.frame(
  AminoAcid = standard,
  Count = as.integer(counts)
)

write.csv(
  results,
  "results/amino_acid_composition.csv",
  row.names = FALSE
)

cat("Protein length:", protein_length, "amino acids\n")

print(head(results[order(-results$Count), ], 3))
```

I used `ggplot2` to create the amino acid composition chart:

```r
fig1 <- ggplot(results, aes(AminoAcid, Count)) +
  geom_col(fill = "grey45") +
  theme_bw() +
  labs(
    title = paste(gene_symbol, "amino-acid composition"),
    x = "Amino acid",
    y = "Number of residues"
  )

print(fig1)

ggsave(
  "figures/01_amino_acid_composition.png",
  fig1,
  width = 8,
  height = 5,
  dpi = 300
)
```

### 5. Examine the assigned variant in ClinVar

I opened the ClinVar record using Variation ID 1467492 and recorded its transcript, HGVS notation, variant type, genomic location, clinical classification, and submission details.

The ClinVar record was accessed at:

https://www.ncbi.nlm.nih.gov/clinvar/variation/1467492/

I also recorded the two germline submissions and their classifications.

### 6. Create the splice-site schematic

Since the assigned variant is a splice-site variant and no protein change was provided, I followed Option B of the laboratory instructions and created an exon–intron schematic instead of a protein-position plot.

The figure was created using `ggplot2` and saved as:

`figures/02_MYH11_splice_site.png`

The schematic illustrates the affected splice-donor boundary. It is a simplified conceptual drawing and does not represent the complete MYH11 exon structure or exact genomic coordinates.

## D. Results and Figures

### Gene and Protein Results

The NCBI search identified MYH11 as the assigned gene, with Gene ID 4629 and chromosome 16.

The protein sequence retrieved using transcript NM_002474.3 corresponded to accession NP_002465.1 and contained **1,972 amino acids**.

### Amino Acid Composition Results

The three most abundant amino acids were:

| Amino acid | Count |
|---|---:|
| Glutamic acid (E) | 272 |
| Leucine (L) | 222 |
| Lysine (K) | 204 |

Glutamic acid had the highest count, followed by leucine and lysine.

### Generated Figures

**Figure 1. MYH11 Amino Acid Composition**

The chart shows the number of occurrences of each of the 20 standard amino acids in the retrieved protein sequence.

File: `figures/01_amino_acid_composition.png`

**Figure 2. MYH11 Splice-Site Variant Schematic**

The schematic shows the exon–intron boundary and the location of the assigned splice-donor change.

File: `figures/02_MYH11_splice_site.png`

## E. Variant Interpretation

The assigned variant is NM_002474.3:c.3858+1G>A, with ClinVar Variation ID 1467492. It is a single-nucleotide variant affecting the +1 position of a canonical splice-donor site.

ClinVar reports an aggregate germline classification of **likely pathogenic**, supported by two out of two submissions:

1. ARUP Laboratories, Molecular Genetics and Genomics — likely pathogenic, last evaluated September 14, 2023.
2. Labcorp Genetics (formerly Invitae) — likely pathogenic, last evaluated July 19, 2022.

The +1 position is important for normal splice-site recognition. A change at this position may interfere with the removal of introns during RNA processing and may result in abnormal MYH11 messenger RNA.

A published 2023 study examined the same variant using the alternative transcript description NM_001040113.2:c.3879+1G>A. RNA analysis identified abnormal transcripts, including alternative splice products and exon skipping. This provides experimental support for an effect on RNA splicing.

However, the exact protein-level consequence for the assigned NM_002474.3 transcript was not provided in the ClinVar record. Therefore, I did not assign a specific amino acid change or claim a particular protein alteration.

## F. Limitations and Possible Validation Experiment

There are limitations to this analysis. The amino acid composition only shows how frequently each amino acid occurs in the protein sequence. It cannot independently explain the protein's three-dimensional structure, activity, or the effect of a mutation.

The splice-site schematic is also simplified and does not show the complete exon–intron organization of MYH11. In addition, the ClinVar classification is based on submitted interpretations and should not be treated as direct proof of every possible molecular consequence.

One possible validation experiment is reverse transcription polymerase chain reaction (RT-PCR), followed by sequencing of the amplified complementary DNA. RNA from patient-derived cells could be compared with RNA from healthy control cells to identify abnormal splice products, such as exon skipping or alternative splice-site usage. This would help determine how the variant affects MYH11 RNA processing.

## Required Questions and Answers

### 1. What does your assigned gene encode, and how does this protein contribute to muscle function?

The assigned gene, **MYH11**, encodes myosin heavy chain 11, a protein found in smooth muscle cells. It interacts with actin and uses energy from ATP to produce the force needed for muscle contraction. This helps regulate the contraction and relaxation of blood vessels and supports the normal movement of internal organs.

### 2. What was the exact protein accession and length retrieved with R?

The exact protein accession retrieved using the assigned transcript NM_002474.3 was **NP_002465.1**. The protein sequence contained **1,972 amino acids**.

### 3. What did the amino-acid composition chart show? What can it NOT tell you?

The amino-acid composition chart showed that glutamic acid (E) was the most abundant amino acid, with 272 occurrences, followed by leucine (L), with 222, and lysine (K), with 204.

However, the chart only shows the frequency of amino acids in the protein sequence. It cannot tell us the protein's three-dimensional structure, its actual biological activity, or how a genetic variant affects its function.

### 4. What is your fixed ClinVar variant, and what is its current classification and evidence strength?

My assigned ClinVar variant is **NM_002474.3:c.3858+1G>A**, with ClinVar Variation ID **1467492** and accession **VCV001467492.9**.

Its current aggregate germline classification is **likely pathogenic**, supported by two out of two submissions. The submitting laboratories are ARUP Laboratories and Labcorp Genetics (formerly Invitae), and both classified the variant as likely pathogenic.

There is also experimental evidence from a published 2023 study that examined the same variant using an alternative MYH11 transcript description. RNA analysis identified abnormal transcripts, supporting an effect on RNA splicing. However, the ClinVar aggregate classification and the experimental findings should be considered separately because they represent different types of evidence.

### 5. What is the scientifically supported mechanism linking this variant to muscle function, and what remains uncertain?

The variant affects the +1 position of a canonical splice-donor site. This position is important for normal pre-mRNA splicing. The G-to-A substitution may interfere with splice-site recognition, resulting in abnormal processing of MYH11 messenger RNA.

Since MYH11 encodes a protein involved in smooth muscle contraction, abnormal RNA processing may affect the production or function of this protein. A 2023 study reported abnormal RNA transcripts associated with the same variant, providing experimental support for altered splicing.

However, the exact transcript and protein consequences for the assigned NM_002474.3 transcript remain uncertain. The ClinVar record does not provide a specific protein change, so a particular amino acid substitution or truncated protein should not be assumed.

### 6. What additional experiment would help test the proposed mechanism?

One useful experiment would be **reverse transcription polymerase chain reaction (RT-PCR), followed by sequencing of the amplified complementary DNA**. RNA from patient-derived cells could be compared with RNA from healthy control cells to identify abnormal splice products, such as exon skipping or alternative splice-site usage.

This experiment would help determine whether the variant affects normal MYH11 RNA processing and provide additional evidence about its possible effect on smooth muscle function.

## G. Database Links, Accession Versions, Access Dates, and References

### Databases and Accession Versions

| Database or record | Accession or ID | Link |
|---|---|---|
| NCBI Gene | Gene ID 4629 | https://www.ncbi.nlm.nih.gov/gene/4629 |
| NCBI Nucleotide | NM_002474.3 | https://www.ncbi.nlm.nih.gov/nuccore/NM_002474.3 |
| NCBI Protein | NP_002465.1 | https://www.ncbi.nlm.nih.gov/protein/NP_002465.1 |
| ClinVar Variation | VCV001467492.9 | https://www.ncbi.nlm.nih.gov/clinvar/variation/1467492/ |
| ARUP submission | SCV004562495.1 | Available through the ClinVar record |
| Labcorp/Invitae submission | SCV002256942.5 | Available through the ClinVar record |

**Date accessed:** October 31, 2026.

### References

1. National Center for Biotechnology Information. (n.d.). MYH11: Myosin heavy chain 11. NCBI Gene. https://www.ncbi.nlm.nih.gov/gene/4629

2. National Center for Biotechnology Information. (n.d.). NM_002474.3, MYH11 transcript record. NCBI Nucleotide. https://www.ncbi.nlm.nih.gov/nuccore/NM_002474.3

3. National Center for Biotechnology Information. (n.d.). NP_002465.1, MYH11 protein record. NCBI Protein. https://www.ncbi.nlm.nih.gov/protein/NP_002465.1

4. National Center for Biotechnology Information. (2026). NM_002474.3(MYH11):c.3858+1G>A, ClinVar Variation ID 1467492. https://www.ncbi.nlm.nih.gov/clinvar/variation/1467492/

5. Arnaud, P., et al. (2023). Early-Onset Aortic Dissection: Characterization of a New Pathogenic Splicing Variation in the MYH11 Gene with Several In-Frame Abnormal Transcripts. https://pmc.ncbi.nlm.nih.gov/articles/PMC11919123/

## H. Errors and Troubleshooting

The R analysis was completed using the `rentrez` and `ggplot2` packages. The protein accession was successfully extracted from the assigned transcript, and the FASTA sequence was retrieved and analyzed.

No major analysis errors remained in the final workflow. If the packages are not installed on a fresh R setup, run the installation command in Method 1 before running the analysis.

If the protein accession cannot be extracted, check that the transcript accession is correct and that the NCBI GenBank record contains a protein ID. If the FASTA file cannot be read, check that the file exists in the `data` folder and that the R working directory is the main project folder.

---
