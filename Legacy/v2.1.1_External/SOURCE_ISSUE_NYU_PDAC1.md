# NYU_PDAC1 source identity audit

## Decision

The files are complete and usable for an independent-study technical transfer
example. They must not be described as a pathology-confirmed independent PDAC
validation case. The disease identity conflict remains unresolved.

## Evidence checked

- GEO GSM6177618 title is `NYU_PDAC1_Vis 10x Genomics Visium` and source name is
  `Spatial transcriptomics of patient primary pancreatic ductal adenocarcinoma`.
- The same record's characteristics are `cancer type: Ovarian carcinoma` and
  `organ: Ovary`. These inconsistent attributes also appear in BioSample
  SAMN28600492. The official family SOFT and extracted sample block are saved.
- Barkley et al., Nature Genetics 2022, doi:10.1038/s41588-022-01141-9,
  publisher Supplementary Tables, worksheet `TableS1`, contains a sample named
  `PDAC NYU1`. Its `Pathology` value is the literal string `NA`.
- Thus the published sample nomenclature supports the PDAC label, but the table
  does not provide specimen pathology that resolves the repository conflict.
  No author correction or specimen-level clinical record was obtained.
- The sample label in TableS1 is `PDAC NYU1`, not an accession-level crosswalk.
  The naming correspondence is supporting evidence, not independent proof of
  the physical specimen's histologic diagnosis.

Do not infer identity from target scores, cell-marker expression, a favorable
correlation, or a non-expert visual assessment of the tissue image.

## Sources

- GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6177618
- Series: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203612
- BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN28600492
- Published article: https://www.nature.com/articles/s41588-022-01141-9
- Publisher XLSX: https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-022-01141-9/MediaObjects/41588_2022_1141_MOESM6_ESM.xlsx
- Official file directory: https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM6177nnn/GSM6177618/suppl/

The publisher workbook SHA-256 is
`a899f552e0042ab67a127c126982aab7c8ccbb8733ad387c26eb6009264c1410`.
The machine-readable identity audit records the exact worksheet row.

## Suggested manuscript wording

"We additionally executed the complete workflow on the public NYU_PDAC1_Vis
dataset (GSM6177618, GSE203612) as an independent-study technical transfer
example. Although its title, source description and published sample label
indicate PDAC, the GEO and BioSample characteristics list ovarian carcinoma
and ovary; the published sample table provides no pathology annotation for
PDAC NYU1. We therefore retain the repository sample name while treating its
tumor identity as unresolved and do not use it as disease-confirmed external
PDAC validation."

This sample's expression matrix was already used in the earlier technical-only
external audit. The new spatial execution expands that case; it is not an
untouched prospective holdout or an additional independent patient cohort.

## HD backup status

The official 10x pancreatic cancer HD sample was located but is not used in this
follow-up. Its spatial archive completed; the 7.67 GB binned archive download
was explicitly stopped at the user's request after a small partial download.
The partial file is retained in the work directory, is not a complete archive,
and must not be packaged or used as an analysis input. The stop was not evidence
of a network blockade. No HD counts, score, selected tissue region or HD
parameter-tuned result was produced.
