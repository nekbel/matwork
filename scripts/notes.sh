##########################
# NGS analyses
##########################

traditional NGS analyses (e.g., QC, read mapping, variant calling, assembly etc.) .

• Transcriptome gives us the complete set of transcripts
It is in the form of a GTF (gene transfer format) file

Transcriptome reference file

Mapping results are only as good as the quality of the
reference transcriptome

• If one does not exist, it can be created using coordinates
from a GTF file and the genome sequence file 

salmon_in_serial
Request 6 cores 
Request 12 hours of runtime
Request 8G of memory

##########################
# Bootstraping
##########################
Bootstraping is any test or metric that uses random sampling with replacement and falls under the broader 
class of resampling methods. It uses sampling with replacement to estimate the sampling distribution for 
the desired estimator. This approach is used to assess the reliability of sequence-based phylogeny.

Bootstrapping is a statistical procedure that resamples a single dataset to create many simulated samples. 
This process allows for the calculation of standard errors, confidence intervals, and hypothesis testing

NOTE: --numBootstraps is necessary if performing isoform-level differential expression analysis with Sleuth, 
but not for gene-level differential expression analysis. Due to the statistical procedure required to assign 
reads to gene isoforms, in addition to the random processes underlying RNA-Seq, there will be technical 
variability in the abundance estimates output from the pseudo-alignment tool [2, 3] for the isoform level 
abundance estimates (not necessary for gene-level estimates). Therefore, we would need technical replicates 
to distinguish technical variability from the biological variability for gene isoforms.

The bootstraps estimate technical variation per gene by calculating the abundance estimates for all genes 
using a different sub-sample of reads during each round of bootstrapping. The variation in the abundance 
estimates output from each round of bootstrapping is used for the estimation of the technical variance for 
each gene.

##########################
# alleles in refgen
##########################
GRCh38, contains information about alternative alleles for various locations on the genome, use the HISAT2 aligner
as it is able to utilize this information during the alignment. 
There is a version of GRCh38 that does not have alleles represented, which is the appropriate version to use with STAR. 


##########################
RNASEQ Metrics
##########################
number of raw reads
percentage of reads aligned to genome
percentage of reads associated with genes