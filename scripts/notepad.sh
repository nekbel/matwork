##########################
# PROTOTYPE SCRIPT
##########################

#!/bin/bash
# SCRIPTS USED:
#       orgUCSC
source ~/.nek_bash_commands.sh #better to source just the function 
echo "Imported source commands"

ref_build="${1:-}"
input="${2:-}"

# check number of arguments given
if [ "$#" -ne 2 ]; then
    echo "Usage: bash <name>.sh <ref_build> <input>"
    echo "      <ref_build> organism reference build name (USCS version)"
    echo "      <input> .txt list of SRR codes"
    echo ""
    echo "Options: run script with -x option"
    exit 1
fi

gen_path=`orgUCSC ${ref_build}`
echo "Path to genome used: $gen_path"

echo "Creating directories (if they don't already exist)" 
mkdir -p


# index file first before bw creation
# do not mv .bai, bam files have to be in the same directory

##########################
# default value :-
##########################
refgen="${1:-hg38}"


##########################
# wget directory, accept/reject files
##########################
wget -r -nH -nd -np -R "*_random.fa.gz, chrUn_*" https://hgdownload.soe.ucsc.edu/goldenPath/mm10/chromosomes/
# -r : recursive
# -nH : Disable generation of host-prefixed directories
# -nd : all files will get saved to the current directory
# -np : Do not ever ascend to the parent directory when retrieving recursively. 
# -R : don't download files with this files pattern
# -A : get only *.html files (for this case)


##########################
# disk used
##########################
du -hs * | sort -h


##########################
# find string contained in file 
##########################
grep -rnwl ./ -e 'SRR11785646'


##########################
# change suffix of files
##########################
for f in ./*.txt
do
  echo $f | sed 's/txt/fq/g'
done


##########################
# GIT
# Matwork Repository
##########################
git add <files>
git commit -m "commit message"
git push -u origin main
token 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

git log 

git checkout <commit-hash> (chech past commit)
git chekout master (return to last commit)

git branch (list branches)
git branch <new-branch-name>
git checkout <branch-name> (check branch)

git merge <branch-name> # merge branch to master


##########################
# kati 8a kanei kai auto mwre
##########################
cat <file.fa> | grep ">" | perl -ne '$_ =~ s/\>//; print $_' | sort | uniq > tidlist.txt


##########################
# GeneCounts STAR
##########################
STAR \
--runThreadN 6  \
--genomeDir /home/myrto/work/refs/Rattus_norvegicus/UCSC/rn6/Sequence/STARIndex \
--sjdbGTFfile /home/myrto/work/refs/Rattus_norvegicus/UCSC/rn6/Annotation/Genes/genes.gtf \
--sjdbGTFtagExonParentTranscript transcript_id \
--readFilesIn sub* \
--quantMode GeneCounts


##########################
# GET FASTQ from SRR
##########################
curl -X GET \
"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR12080417&fields=fastq_ftp&result=read_run" |
cut -f 2 |
sed -n 2p |
tr \; "\n" |
wget -P ./fqs -i - &&
gunzip *.gz


##########################
# CHIP SEQ ALIGN BOWTIE2 (SINGLE-END)
##########################
trim_galore -o ./trimmed_fqs -q 20 ./fqs/SRR12080416; trim_galore -o ./trimmed_fqs -q 20 ./fqs/SRR12080417

bowtie2 --local -X 3000 -p 25 -x /home/myrto/work/refs/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U ./trimmed_fqs/SRR12080416_trimmed.fq | samtools sort -O BAM - > ./bowtie/SRR12080416.sorted.bam
bowtie2 --local -X 3000 -p 25 -x /home/myrto/work/refs/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome -U ./trimmed_fqs/SRR12080417_trimmed.fq | samtools sort -O BAM - > ./bowtie/SRR12080417.sorted.bam


##########################
# grep
##########################
grep
-B print n lines before match
-A print n lines after match
--no-group-separator dont print -- separator between matches
-n print line number of match
-v return everything that DOES NOT match
-c return the number of matches rather than the actual matching lines
-H report the filename


##########################
# sort
##########################
sort 
-k which column to sort by
-u only return unique values


##########################
# scp and rsync
##########################
scp [OPTION] [user@]SRC_HOST:]file1 [user@]DEST_HOST:]file2

rsync works similarly to scp.
-v 	verbosity (it prints on the screen what is being copied)

-e 	ssh is for encryption (use the ssh protocol 
	for encryption of the file transfer)


##########################
# basename
##########################
samplename=`basename <path.to.fq> .fq`
# retrieve filename without leading directories (prefix). 
# Also remove specified suffix


##########################
# .bashrc and .profile
##########################
.bashrc and .bash_profile located on ~
# You only log in once, and that's when ~/.bash_profile or ~/.profile is read and executed. 
# environment variables (LESS, PATH, MANPATH, LC_*, ...) e.g. 
# On the next runned shells, your environment variables will be passed along, so your non-login 
# shells don't need to load them anymore
# Non-login shells only execute ~/.bashrc, not /.profile or ~/.bash_profile, for this exact reason, 
# define everything that only applies to bash (functions, aliases, bash-only variables like HISTSIZE 
# (this is not an environment variable, shell options with set and shopt
# You'll need to add in my ~/.profile too: source ~/.bashrc.
.bashrc vs .bash_profile
# .bash_profile is executed for login shells, while .bashrc is executed for interactive non-login shells. 
# When you login (type username and password) to O2 the .bash_profile is executed. So if you want the alias 
# available only when you login, you will want to put it in your .bash_profile.


##########################
# wget files & dirs
##########################
wget
-r/--recursive 
-np/--no-parent & -nd no parent directories 
-R "index.html*" reject this file


##########################
# curl individual file (not dirs)
##########################
curl
-O will save the file in the current working directory with the same file name as remote. 
-o lets you specify a different file name or location
-l list contents (directory -> directories and files)
				 (file -> print content in stdout)


##########################
# curl hack for ensembl
# returs list of files and 
# directories contained
##########################
curl -l <ftp.link> | grep href= |
head -n -1 | sed 's/<a href="//g' | egrep -o '^.+"' | sed 's/.$//' | tail -n +2

##########################
# salmon
##########################
salmon quant -i /n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/indicies/salmon_index \
 -l A \
 -r $fq \
 -o ${samplename}.salmon \
 -p 6 \
 --numBootstraps 30 \
 --seqBias \
 --useVBOpt \
 --validateMappings

# qualimap before so u can add 
--gcBias to learn and correct for fragment-level GC biases in the input data
--posBias will enable modeling of a position-specific fragment start distribution

##########################
# qualimap
##########################
qualimap rnaseq \
-outdir results/qualimap/Mov10_oe_1 \
-a proportional \
-bam results/STAR/Mov10_oe_1_Aligned.sortedByCoord.out.bam \
-p strand-specific-reverse \
-gtf /n/holylfs05/LABS/hsph_bioinfo/Everyone/Workshops/Intro_to_rnaseq/reference_data/Homo_sapiens.GRCh38.92.gtf \
--java-mem-size=8G


##########################
# multiqc
##########################
multiqc -n multiqc_report_rnaseq \
full_dataset_results/fastqc/*zip \
full_dataset_results/STAR/*Log.final.out \
full_dataset_results/qualimap/* \
full_dataset_results/salmon/*salmon


##########################
# input arguments bash
##########################
sh  run_rnaseq.sh  input.fq  input.gtf  12

$0 => run_rnaseq.sh
$1 => input.fq
$2 => input.gtf
$3 => 12

wrappedProgram "$@"
# ^^^ this is correct and will hand over all arguments in the way
#     we received them, i. e. as several arguments, each of them
#     containing all the spaces and other uglinesses they have.
wrappedProgram "$*"
# ^^^ this will hand over exactly one argument, containing all
#     original arguments, separated by single spaces.
wrappedProgram $*
# ^^^ this will join all arguments by single spaces as well and
#     will then split the string as the shell does on the command
#     line, thus it will split an argument containing spaces into
#     several arguments.