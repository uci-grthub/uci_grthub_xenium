use NCBI GEO as the primary repository for the Xenium dataset. 

place the complete original Xenium output folders, along with any processed analysis files, in a convenient shared location from which I can download them for submission

Include:
  1. morphology image files (geo_staging/Aged_cohort/morphology_focus, geo_staging/Adult_cohort/morphology_focus)
  2. decoded transcript-coordinate files (what does decoded mean here?)
  3. cell-feature matrices (for only 'all' for the 2 cohorts? not integrated data)
  4. cell and nucleus segmentation/boundary files (geo_staging/Aged_cohort/nucleus_boundaries.csv.gz, geo_staging, Aged_cohort/cell_boundaries.csv.gz, geo_staging/Aged_cohort/cell_feature_matrix.h5)
  5. panel and run information (geo_staging/metadata/Adult_cohort_experiment.xenium, geo_staging/metadata/Aged_cohort_experiment.xenium)
  6. sample-to-ROI/group/brain-region mapping. (geo_staging/Aged_cohort/coordinates.csv, geo_staging/Adult_cohort/coordinates.csv)

Since Xenium is imaging-based, I understand that the relevant raw files are the morphology images and transcript-coordinate outputs rather than FASTQ sequencing files. 

A checksum manifest would also be helpful, if available, to verify file integrity after transfer.

from ncbi guidelines: 
10x Genomics Xenium in Situ

For Xenium studies, GEO requires raw data, processed data, and metadata.

    Raw data files
        High-resolution morphology images
        Decoded transcripts (parquet or csv format)
    Processed data files
        Transcript count (MEX or HDF5 format)
        Cell summary (parquet or csv format)
        Cell segmentation (parquet or csv format)
        Nucleus segmentation (parquet or csv format)
    Metadata
        Use the Xenium-specific template for providing study and sample metadata:



this is an example GSE243168

