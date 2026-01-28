## FDR Evaluation for DIA

This repository contains fasta and scripts for the evaluation of FDR in Library-Free search by DIA-NN using an *in vitro* human proteome. 

**Reference:** Evaluation of the False Discovery Rate in Library-Free Search by DIA-NN Using In Vitro Human Proteome  
https://doi.org/10.1021/acs.jproteome.5c00036

###  Data Availability
The raw MS data and DIA-NN search results generated for this study have been deposited in the **jPOST repository** under the following accession numbers:
* JPST003404/PXD056519.

###  File Descriptions 

#### 1. FASTA Files for Library-Free Search
* **`clean_Mix-x.fasta`**: These file contain the predigest protein sequences used for the library-free search in DIA-NN.

#### 2. Protein Sequences
* **`ProteinSequenceMix-x.fasta`**: Contain the protein sequences present in each mixture.

**Note:** To facilitate computational processing, proteins present in multiple mixtures are assigned a **merged header**. A concatenated name indicates that the protein sequence is identical and present across the specified mixtures.

*Example:*
```text
>Mix4_FLJ93078AAAF/Mix11_FLJ93078AAAF
MAPPSVFAEVPQAQPVLVFKLTADFREDPDPRKVNLGVGAYRTDDCHPWVLPVVKKVEQK...'''

#### 3. `InSilico_TheoreticalPeptides`
InSilico_TheoreticalPeptides files contain the theoretical peptides of proteins in each mixture, used for comparing against diann search results.
