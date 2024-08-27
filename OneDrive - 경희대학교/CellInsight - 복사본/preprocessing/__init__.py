"""

Package name: preprocessing

Package description:
- Aim to transform data to anndata format
- Process two databases: PanglaoDB and Single Cell Portal

For PanglaoDB:
  - Using R is inevitable, so one code is written in Python and another code in R
  - By using R, data whose type is RData is transformed into a list composed of a sparse matrix and metadata
  - By using Python, the list is transformed into an anndata object

For Single Cell Portal:
  - Only Python is required
  
"""