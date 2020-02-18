# RWC23_elt2_regulated_genes

The purpose of this project is to define targets of C. Elegans transcription factor ELT-2. This project will have the general workflow:

- Download data
    - RNA-seq of wildtype L1 worms or elt-2(-) L1 worms
    - ELT-2 ChIP-seq in L1 worms
- Differential expression (DESeq2)
- Reformat data (R)
- Integrate RNA-seq and ChIP-seq (Cistrome BETA)

# Data sources

## Differential expression data

This data was originally published in the following paper:

    Dineen A, Osborne Nishimura E, Goszczynski B, Rothman JH, McGhee JD. Quantitating transcription factor redundancy: The relative roles of the ELT-2 and ELT-7 GATA factors in the C. elegans endoderm. Dev Biol. January 2018. doi:10.1016/j.ydbio.2017.12.023

The file `Table_S1_Raw_Read_Counts.xlsx` was downloaded and converted to a `.csv` file.

## ChIP-seq

This data was originally published in the following paper:

    Kudron MM, Victorsen A, Gevirtzman L, et al. The modern resource: genome-wide binding profiles for hundreds of Drosophila and Caenorhabditis elegans transcription factors. Genetics. 2018;208(3):937-949. doi:10.1534/genetics.117.300657

ELT-2 ChIP-seq raw reads were processed by David King, and generated the file `L1_1_L1_2.IDR_0.05.narrowPeak`. See the following reposity for more information: https://github.com/meekrob/onish_ChIP_R_Analysis
