# MYH11 Protein Structure and Mutation Analysis

## Project Overview

This project investigates the human **MYH11 (myosin heavy chain 11)** gene and its assigned splice-site variant using NCBI databases, R programming, and published scientific literature.

The analysis includes gene identification, retrieval of the protein sequence, amino acid composition analysis, ClinVar variant investigation, and a schematic illustrating the affected splice-donor site.

## 1. Gene and Assignment Information

| Item | Information |
|---|---|
| Gene symbol | MYH11 |
| Official gene name | Myosin heavy chain 11 |
| Organism | Homo sapiens |
| NCBI Gene ID | 4629 |
| Chromosome | 16 |
| Assigned transcript | NM_002474.3 |
| Assigned variant | NM_002474.3:c.3858+1G>A |
| ClinVar Variation ID | 1467492 |
| Variant type | Single-nucleotide variant (splice donor) |
| Protein accession | NP_002465.1 |
| Protein length | 1,972 amino acids |

## 2. Biological Function of MYH11

The MYH11 gene encodes smooth muscle myosin heavy chain 11, a major contractile protein found in smooth muscle cells. It is part of a myosin complex that converts chemical energy from ATP hydrolysis into mechanical energy. This allows smooth muscle to contract and relax, supporting functions such as the regulation of blood vessel diameter and movement of internal organs.

Changes affecting MYH11 may interfere with normal smooth muscle function and are associated with conditions including familial thoracic aortic aneurysm and dissection.

Source: [NCBI MYH11 Gene](https://www.ncbi.nlm.nih.gov/datasets/gene/4629/)

## 3. Protein Sequence Retrieval and Amino Acid Composition

The protein sequence was retrieved from NCBI using the assigned RefSeq transcript NM_002474.3. The corresponding protein accession identified by the R script was **NP_002465.1**.

The retrieved protein sequence contained **1,972 amino acids**.

The three most abundant amino acids were:

| Amino acid | Count |
|---|---:|
| Glutamic acid (E) | 272 |
| Leucine (L) | 222 |
| Lysine (K) | 204 |

Glutamic acid was the most abundant amino acid, followed by leucine and lysine.

The amino acid composition describes the frequency of residues in the protein sequence. However, composition alone cannot establish the protein's three-dimensional structure, biological activity, or the functional effect of a genetic variant.

## 4. ClinVar Variant Investigation

The assigned variant is:

**NM_002474.3:c.3858+1G>A**

- ClinVar Variation ID: 1467492
- ClinVar accession: VCV001467492.9
- Genomic location: chromosome 16, GRCh38 position 15,726,847
- Molecular consequence: splice-donor variant
- Protein change: Not provided
- Germline classification: Likely pathogenic
- Submissions contributing to classification: 2 out of 2

The two contributing submitters are:

1. ARUP Laboratories, Molecular Genetics and Genomics — Likely pathogenic.
2. Labcorp Genetics (formerly Invitae) — Likely pathogenic.

The aggregate classification is based on the submitted clinical interpretations. The evidence includes the variant's canonical splice-donor location, computational predictions of altered splicing, and clinical observations reported in the literature. The classification should not be interpreted as proof of a specific protein-level consequence.

[View the ClinVar record](https://www.ncbi.nlm.nih.gov/clinvar/variation/1467492/)

## 5. Discussion: Possible Molecular Effect of the Variant

The assigned MYH11 variant, NM_002474.3:c.3858+1G>A, affects the first nucleotide of an intron at a canonical splice-donor site. The +1 position is important for recognition of the splice site during pre-mRNA processing. A nucleotide substitution at this position may interfere with normal splicing and produce an abnormal messenger RNA transcript.

ClinVar currently classifies the variant as likely pathogenic, with two out of two submissions contributing to the aggregate classification. The ARUP submission reports that the variant disrupts a canonical splice-donor site and has been described in an individual with patent ductus arteriosus, aortic dissection, and intracranial vessel stenosis. The Labcorp/Invitae submission also reports computational predictions suggesting disruption of the consensus splice site. These findings support concern about its potential effect, although the clinical significance should be interpreted with the limitations of the available evidence.

Additional experimental evidence was reported in a 2023 study examining the same genomic variant using the alternative MYH11 transcript description NM_001040113.2:c.3879+1G>A. RNA analysis from cultured patient fibroblasts identified abnormal transcripts, including alternative splice products affecting the coiled-coil region of the protein. Some abnormal transcripts retained the reading frame, demonstrating that the consequences of altered splicing may be more complex than simply producing a truncated protein.

Because the assigned transcript does not have a protein-level change provided in ClinVar, a specific amino acid substitution or protein alteration should not be assumed. Further investigation is needed to establish the precise transcript and protein consequences for the assigned reference transcript.

One suitable follow-up experiment is reverse transcription polymerase chain reaction (RT-PCR), followed by sequencing of the amplified complementary DNA. Comparing patient-derived samples with healthy controls could identify abnormal splice products and help determine how the variant affects MYH11 RNA processing.

Overall, the findings support a potential disruption of normal MYH11 RNA splicing, which may affect the function of a protein essential for smooth muscle contraction.

## 6. Proposed Additional Experiment

**Experiment:** Reverse transcription polymerase chain reaction (RT-PCR) followed by Sanger sequencing.

**Purpose:** To investigate whether the assigned splice-site variant produces abnormal MYH11 RNA transcripts.

**Expected approach:** RNA would be isolated from patient-derived cells and healthy control cells, converted into complementary DNA, and amplified using primers surrounding the affected splice region. The amplified products would then be sequenced and compared to identify possible exon skipping, partial exon deletion, or alternative splice-site usage.

## 7. Methods and Tools

The project used:

- NCBI Gene database for gene identification.
- NCBI Nucleotide database for transcript retrieval.
- NCBI Protein database for protein FASTA retrieval.
- ClinVar for variant classification and clinical submission information.
- R and RStudio for sequence processing and analysis.
- `rentrez` for retrieving NCBI records.
- `ggplot2` for visualizing amino acid composition and the splice-site schematic.

### Analysis Workflow

1. Identified the assigned MYH11 gene using NCBI.
2. Retrieved the assigned transcript and corresponding protein accession.
3. Downloaded the protein FASTA sequence.
4. Calculated amino acid counts and protein length using R.
5. Generated the amino acid composition chart.
6. Examined the assigned variant in ClinVar.
7. Created a conceptual exon–intron schematic showing the splice-donor variant.

### Date Accessed

NCBI and ClinVar records were accessed on **October 31, 2026**.

## 8. Project File Organization

The repository contains the following files:

```text
Muscle-Contractility-MYH11/
│
├── README.md
├── analysis.R
│
├── data/
│   └── MYH11_protein.fasta
│
├── results/
│   ├── gene_summary.txt
│   ├── protein_accession.txt
│   └── amino_acid_composition.csv
│
└── figures/
    ├── 01_amino_acid_composition.png
    └── 02_MYH11_splice_site.png
```

## 9. References

1. National Center for Biotechnology Information (NCBI). MYH11: Myosin heavy chain 11. [NCBI Gene ID 4629](https://www.ncbi.nlm.nih.gov/datasets/gene/4629/).
2. National Center for Biotechnology Information (NCBI). Myosin-11 isoform SM1A, Homo sapiens. [RefSeq NP_002465.1](https://www.ncbi.nlm.nih.gov/protein/NP_002465.1).
3. National Center for Biotechnology Information (NCBI). ClinVar: NM_002474.3(MYH11):c.3858+1G>A, Variation ID 1467492. [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/variation/1467492/).
4. Arnaud, P., et al. (2023). *Early-Onset Aortic Dissection: Characterization of a New Pathogenic Splicing Variation in the MYH11 Gene with Several In-Frame Abnormal Transcripts*. [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC11919123/).
5. Larson, A., et al. (2020). *Intracranial Vessel Stenosis in a Young Patient with an MYH11 Mutation: A Case Report and Review of 2 Prior Cases*. World Neurosurgery, 137, 243–246. [PMID: 32081817](https://pubmed.ncbi.nlm.nih.gov/32081817/).

---

**Project:** MYH11 Protein Structure and Mutation Analysis  
**Organism:** Homo sapiens  
**Assigned variant:** NM_002474.3:c.3858+1G>A
