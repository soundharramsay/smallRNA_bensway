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

###### sam to bam 

#!/bin/bash
specify the directory where your SAM files are located
sam_directory="./"

# Specify the path to the reference FASTA file
reference_fasta="/athena/kleavelandlab/store/sor4003/2_star_genome_index_nexflow/small_RNA_genomes/mirBASE/bowtie_1_3_1_spack_load_index_mirBase/mature_mirBASE_all.fa"

# Iterate over each SAM file in the directory
for sam_file in "$sam_directory"/*.sam; do
    # Generate the corresponding BAM filename
    bam_file="${sam_file%.sam}.bam"

    # Create a new header using samtools reheader
    samtools view -H "$sam_file" | sed '/^@SQ/s/$/ UR:unknown/' | samtools reheader - "$sam_file" > "$bam_file"

    echo "Converted $sam_file to $bam_file"
done


######## merged text file 

To perform a column-wise merge of column 1 and 7 from the provided text files and only include common lines based on column 1, you can use the join command along with some additional Unix utilities. Here's a step-by-step process:

Extract column 1 and 7 from each file.
Sort each file based on column 1.
Use the join command to merge the sorted files based on column 1.
Here's a command that accomplishes this:

cut -f 1,7 -d $'\t' E10_1_Z8KO_k562_sample4_S4_L002_R1_001.featureCounts.txt > col1_7_1.txt
cut -f 1,7 -d $'\t' E13_3_z8ko_sample7_S7_L002_R1_001.featureCounts.txt > col1_7_2.txt
cut -f 1,7 -d $'\t' nt3_k562_sample3_S3_L002_R1_001.featureCounts.txt > col1_7_3.txt
cut -f 1,7 -d $'\t' E13_1_z8ko_sample5_S5_L002_R1_001.featureCounts.txt > col1_7_4.txt
cut -f 1,7 -d $'\t' nt1_k562_sample1_S1_L002_R1_001.featureCounts.txt > col1_7_5.txt
cut -f 1,7 -d $'\t' E13_2_z8ko_sample6_S6_L002_R1_001.featureCounts.txt > col1_7_6.txt
cut -f 1,7 -d $'\t' nt2_k562_sample2_S2_L002_R1_001.featureCounts.txt > col1_7_7.txt

sort -t $'\t' -k1,1 col1_7_1.txt -o col1_7_1.txt
sort -t $'\t' -k1,1 col1_7_2.txt -o col1_7_2.txt
sort -t $'\t' -k1,1 col1_7_3.txt -o col1_7_3.txt
sort -t $'\t' -k1,1 col1_7_4.txt -o col1_7_4.txt
sort -t $'\t' -k1,1 col1_7_5.txt -o col1_7_5.txt
sort -t $'\t' -k1,1 col1_7_6.txt -o col1_7_6.txt
sort -t $'\t' -k1,1 col1_7_7.txt -o col1_7_7.txt

join -t $'\t' -1 1 -2 1 col1_7_1.txt col1_7_2.txt | \
join -t $'\t' -1 1 -2 1 - col1_7_3.txt | \
join -t $'\t' -1 1 -2 1 - col1_7_4.txt | \
join -t $'\t' -1 1 -2 1 - col1_7_5.txt | \
join -t $'\t' -1 1 -2 1 - col1_7_6.txt | \
join -t $'\t' -1 1 -2 1 - col1_7_7.txt > merged_output.txt

rm col1_7_*.txt

############  grep the microRNA from fastq files 
# Create separate .txt files for each pattern match
for file in trim_*.fastq.gz; do
  pattern="TGGAAGACTAGTGATTTTGTTGTT"
  zgrep "$pattern" "$file" > "${file}_mir7_result.txt"
done


############# counting the line 
# Create a result.txt file to store the counts
echo -e "Filename\tNumber_of_Lines\tAverage_Line_Length" > result.txt

# Iterate over each text file
for file in *.txt; do
  # Get the counts
  line_count=$(wc -l < "$file")
  avg_line_length=$(awk '{ total += length($0) } END { print total/NR }' "$file")

  # Append the results to result.txt
  echo -e "$file\t$line_count\t$avg_line_length" >> result.txt
done
