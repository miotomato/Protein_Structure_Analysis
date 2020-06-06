First part: Protein sequence and structure database creation
  - Parsing PDB, PDBsum and PISCES original files.
  - Extract corresponding DNA sequence using TBLASTN, clean and match the results to PDB format files.
  - Tessellated protein structure in 3-dimension.
  - Create SQL databases for further analysis.

Second part: Protein sequence, secondary and tertiary structure analysis
  - Statistical analysis on tessellated results including calculate the log-likelihood ratio for each simplex and perform hypothesis tests.
  - Secondary structure classification, clustering and assignment using machine learning algorithms and deep learning.
  - Codon level protein-protein interfaces analysis, try to find relationship between mutagenesis and log-likelihood potential on structure levels.
