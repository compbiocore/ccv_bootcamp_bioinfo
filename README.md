# ccv_bootcamp_bioinfo
    90 mins each            
**FASTQ (30 mins total)**               
    - Intro to fastq format, quality encoding (5 minutes)               
    - SRA intro (study -> sample -> experiment -> run) (5 minutes)               
    - download reads from SRA using fasterq-dump (5 minutes)               
    - download reads from ENA using wget (5 minutes)               
    - read QC tools (fastQC, fastqscreen) (5 minutes)
               
**FASTA (30 mins total)**                      
    - Intro to fasta format (less than 5 minutes)
    - mostly for downloading reference genomes, you could make a transcriptome reference, etc. but they are all .fasta format, compare where the genes are actually positioned on ensembl and ucsc so don't mix and match GTFs.
    - Lots of places to download data -- different gene models and different naming conventins (e.g. chr1 vs 1). (5-10 minutes)               
    - refseq (5 minutes)               
    - ensembl (5 minutes)               
    - ucsc (5 minutes)               

**SAM (30 mins total)**               
    - Many sequencing pipelines (RNAseq, chipseq, bisulfite seq) will eventually give you a sam file, intro to the SAM format (and BAM) (5 minutes)               
    - flags and cigars (5-10 minutes?)               
    - samtools view, samtools sort, samtools index (5 minutes)               
    - alignment QC tools (qualimap, picard) (10 minutes)     

**Annotation file types and resources (30 mins total)**               
    - structural/positional (location on chr, intron or exon) vs functional annotations (gene symbol, GO term, etc.) (5 minutes)               
    - GTF, GFF2, GFF3 (5-10 minutes)  
    - gene ID types -- refseq, ucsc, ensembl, gencode, others? (5 minutes)               
    - GO, KEGG, Reactome -- go through w clusterprofiler (10 minues)     

**R annotation packages (15-20 minutes)** 
    - TxDb/EnsDb (5 minutes)
    - orgDBs (5 minutes)               
    - biomaRt (5 minutes)               
    - annotationHub (5 minutes)               
    - GenomeInfoDB (5 minutes)               

https://bioconductor.github.io/BiocWorkshops/introduction-to-bioconductor-annotation-resources.html
https://bioconductor.org/packages/devel/workflows/vignettes/annotation/inst/doc/Annotation_Resources.html

Check out ucsc browswer page w file format definitions
all the annotation stuff will be on the hub

put together a google doc w the outlines/high level info for the R workshops and the bioinformatics workshops

Getting Started on Oscar

Working on the Command Line
