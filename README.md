## Dataset

Public GEO dataset:

- **GSE139369**
- Example sample: **GSM4138872_scRNA_BMMC_D1T1**

### Dataset Note

The processed count matrices provided by the original study have already been filtered to remove mitochondrial and ribosomal genes.

Therefore, mitochondrial quality control cannot be calculated from these matrices. The pipeline automatically detects this situation and skips mitochondrial-based filtering.
