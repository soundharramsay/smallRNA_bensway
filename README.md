# smallRNA_bensway

sg kleavelandlab 'srun --pty -n1 --mem=8G -p scu-cpu /bin/bash -i'

#loading individual modules/applications--can only do this within an interactive session or bash script 

#fastqc
=====
module load fastqc-0.11.7-gcc-8.2.0-kcqpnwl

#cutadapt
=======
module load miniconda3-4.8.2-gcc-8.2.0-azglaow
conda activate /software/apps/cutadapt3.4

#python
======
module load python/3.9.0

#bowtie
======
module load bowtie-1.3.1-gcc-8.2.0-rzjujdj

#samtools
========
module load samtools-1.10-gcc-8.2.0-nfo6kn5


#### 
#!/bin/bash

# List of compressed FASTQ files with full relative paths
fastq_files=(
    "trim_E10_1_Z8KO_k562_sample4_S4_L002_R1_001.fastq.gz"
    "trim_E13_3_z8ko_sample7_S7_L002_R1_001.fastq.gz"
    "trim_nt3_k562_sample3_S3_L002_R1_001.fastq.gz"
    "trim_E13_1_z8ko_sample5_S5_L002_R1_001.fastq.gz"
    "trim_nt1_k562_sample1_S1_L002_R1_001.fastq.gz"
    "trim_E13_2_z8ko_sample6_S6_L002_R1_001.fastq.gz"
    "trim_nt2_k562_sample2_S2_L002_R1_001.fastq.gz"
    # Add more file names here
)

# Loop over each FASTQ file
for fastq_file in "${fastq_files[@]}"; do
    # Extract the file name without the path and extension
    filename=$(basename "$fastq_file")

    # Your processing commands for each FASTQ file
    zcat "$fastq_file" | bowtie -S -n 0 -l 20 -m 20 -a --best --strata /athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/mirBASE/bowtie_1_3_1_spack_load_index_mirBase/mature_mirBASE -q - | 
    samtools sort -o - - | 
    samtools view -F 4 - > "${filename%.fastq.gz}_n0_mapped.sam" &
done
