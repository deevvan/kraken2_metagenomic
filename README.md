# kraken2 metagenomic detection README
Python scripts that download, trim, align SRA files to host (human or mouse) genome, extract unaligned reads and detect pathogens using Kraken2

## 1. Creating conda environment for RSEM:
    conda create -n my_environment bioconda::RSEM
    conda activate RSEM


## 2. Creating RSEM transcriptome reference using STAR:
    
    human genome:
    rsem-prepare-reference --gtf /path/to/reference_genome_directory/star/UCSChg38_human.gtf --/path/to/reference_genome_directory/star/UCSChg38_human.fa /path/to/reference_genome_directory/star/UCSChg38_human_idx
    
    mouse genome:
    rsem-prepare-reference --gtf /path/to/reference_genome_directory/star/GRCm39_mouse.gtf --/path/to/reference_genome_directory/star/GRCm39_mouse.fna /path/to/reference_genome_directory/star/GRCm39_mouse_idx

## 3. Creating STAR genome index: 
    human star idx:
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /path/to/reference_genome_directory/star/ --genomeFastaFiles /path/to/reference_genome_directory/star/UCSChg38_human.fa --sjdbGTFfile /path/to/reference_genome_directory/star/UCSChg38_human.gtf --sjdbOverhang 100 --outFileNamePrefix /path/to/reference_genome_directory/star/UCSChg38_human_idx

    Creating mouse genome STAR index:
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /path/to/reference_genome_directory/star_mouse --genomeFastaFiles /path/to/reference_genome_directory/star_mouse/GRCm39_mouse.fna --sjdbGTFfile /path/to/reference_genome_directory/star_mouseGRCm39_mouse.gtf --sjdbOverhang 100 --outFileNamePrefix /path/to/reference_genome_directory/star_mouse/GRCm39_mouse_idx


## 4. Run rsem_expression.py
    Expression calculation:
    rsem-calculate-expression --paired-end -p 8 --output-genome-bam --keep-intermediate-files --star --append-names --estimate-rspd --time /path/to/trimmed_SRP_directory/{srr_id}_1.fastq /path/to/trimmed_SRP_directory/{srr_id}_2.fastq /path/to/reference_genome_directory/star/UCSChg38_human_idx /path/to/RSEM_SRP_directory/{srr_id}
    
Output files: In addition to the aligned BAM file (genome level and transcriptome level), this will generate the unaligned (unmapped) fastq files named {srr_id}_un_1.fq and {srr_id}_un_2.fq 
These fq files comprise reads that did not align to the host reference genome. Script will also generate a {srr_id}.genes.results which contains read counts, TPM and FPKM for host genes.


## 5. Run read_count_creator.py
Creating counts matrix of all host genes aligned RNAseq reads 


## 6. Creating SARS-CoV-2 genome STAR index:
    
    STAR --runThreadN 1 --runMode genomeGenerate --genomeDir /path/to/reference_genome_directory/star_star_cov --genomeFastaFiles /path/to/reference_genome_directory/star_star_cov/SARS-CoV-2-MN908947.3.fa --sjdbGTFfile /path/to/reference_genome_directory/star_star_cov/Sars_cov_2.ASM985889v3.101.gtf --sjdbOverhang 100 --genomeSAindexNbases 6 --outFileNamePrefix /path/to/reference_genome_directory/star_star_cov/cov_MN908947_idx
    
    STAR --runThreadN 1 \
     --runMode genomeGenerate \
     --genomeDir /path/to/reference_genome_directory/star_star_cov \
     --genomeFastaFiles /path/to/reference_genome_directory/star_star_cov/SARS-CoV-2-MN908947.3.fa \
     --sjdbGTFfile /path/to/reference_genome_directory/star_star_cov/Sars_cov_2.ASM985889v3.101.gtf \
     --sjdbOverhang 100 \
     --genomeSAindexNbases 6 \
     --outFileNamePrefix /path/to/reference_genome_directory/star_star_cov/cov_MN908947_idx


## 7. Run viral_star_align.py
Mapping {srr_ID}_un_1 and {srr_ID}_un_2 fasta to viral STAR genome index:
    
    STAR --genomeDir /mmfs1/projects/path/to/reference_genome_directory/star_star_cov/cov_MN908947_idx --readFilesIn {input_dir}/{srr_id}/{srr_id}_un_1.fq {input_dir}/{srr_id}/{srr_id}_un_2.fq --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --outFileNamePrefix /path/to/RSEM_viral_directory/{srr_id}/{srr_id}

    STAR --genomeDir /path/to/reference_genome_directory/star_star_cov/cov_MN908947_idx \
     --readFilesIn {srr}_un_1.fastq {srr}_un_2.fastq \
     --outSAMtype BAM SortedByCoordinate \
     --quantMode GeneCounts \
     --outFileNamePrefix {srr} \
     --limitBAMsortRAM 7999715745

A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. \

The counts coincide with those produced by htseq-count with default parameters. \

This option requires annotations (GTF or GFF with –sjdbGTFfile option) used at the genome generation step, or at the mapping step. \

STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options. 

    
With --quantMode GeneCounts option STAR will count number reads per gene while mapping. \

A read is counted if it overlaps (1nt or more) one and only one gene. Both ends of the pairedend read are checked for overlaps. \

The counts coincide with those produced by htseq-count with default parameters. \

This option requires annotations (GTF or GFF with –sjdbGTFfile option) used at the genome generation step, or at the mapping step. 


STAR outputs read counts per gene into ReadsPerGene.out.tab file with 4 columns which correspond to different strandedness options: \
  column 1: gene ID \
  column 2: counts for unstranded RNA-seq \
  column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes) \
  column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse) 
            
Unstranded RNA-seq: \
  The directionality of the RNA is not preserved, so you cannot tell which DNA strand the RNA came from. \
  The counts for unstranded RNA-seq represent the total number of reads that align to a gene, regardless of the strand. \
        
Stranded RNA-seq (htseq-count option -s yes): \
  The counts represent the number of reads that align to the same strand as the RNA. \
  This is typically the second (coding) strand of DNA, because RNA is synthesized from the template strand. 
        
Stranded RNA-seq (htseq-count option -s reverse): \
  The counts represent the number of reads that align to the opposite strand of the RNA. \
  This is typically the first (template) strand of DNA. 


## 8. Run reads_count_viral_creator.py
Creating count matrix of all viral gene aligned RNAseq reads 
    
OR 
    
  a. Run stringtie_bam.py \
      Convert viral STAR aligned bam files into stringtie gtf files 
        
  b. And then run srr_list_creator.py \
      Create a sample_list.txt file with a list of gtf files and their paths 
        
  c. Run python3 prepDE.py -i sample_list.txt \
      Convert all gtf files into count matrix        


## 9. Run sam_indexing.py
Convert {srr_id}.sortedByCoord.out.bam into {srr_id}.sortedByCoord.out.bam.bai and retrieve {srr_id}_idxstats.txt \
    samtools idxstats /path/to/RSEM_viral_directory/{srr_id}/{srr_id}.sortedByCoord.out.bam > /path/to/RSEM_viral_directory/{srr_id}/{srr_id}_idxstats.txt
    
samtools idxstats retrieve and print stats in the index file corresponding to the input file. \

Before calling idxstats, the input BAM file should be indexed by samtools index. \

The output is TAB-delimited with each line consisting of: \
    column 1: reference sequence name, \
    column 2: sequence length, \
    column 3: mapped read-segments and \
    column 4: unmapped read-segments. \
    
Note this may count reads multiple times if they are mapped more than once or in multiple fragments.


## 10. Run mRNA_read_lengths.py
Takes in {srr_id}_1.fastq and {srr_id}_2.fastq to calculate average readLength for each srr_id front and reverse file
    
    
## 11. Alternative to STAR: Extract unmapped reads to host genome using bowtie2:
Use bowtie2 to extract unmapped reads to host genome:
    bowtie2 -x /path/to/reference_genome_directory/bowtie2_GRCh38/GRCh38 \
    -p 8 -q -1 /path/to/trimmed_SRP_directory/{srr_id}_1.fastq \
    -2 /path/to/trimmed_SRP_directory/{srr_id}_2.fastq --un-conc {srr_id}_nonhost.fa
  

## 12. Create kraken2 database using only covid, RSV and Influenza genomes:
    conda activate kraken_env
  
    for file in *.fa
        do
    
    kraken2-build --add-to-library $file --db /path/to/krakendb_viral/
        done
    
  Then,
  
    kraken2-build --build --db /path/to/krakendb_viral/

## 13. Perform kraken2 alignment on unmapped reads to detect viral coinfections

    kraken2 --db {kraken_db} --threads 8 --minimum-hit-groups 3 \ --report-minimizer-data \
    --report {kraken_output_path}/{srr_id}.report \ --paired {RSEM_path}/{srr_id}/{srr_id}_nonhost.1.fa \
    {RSEM_path}/{srr_id}/{srr_id}_nonhost.1.fa > {kraken_output_path}/{srr_id}.kraken2

    
    example:
        
        kraken2 --db /path/to/krakendb_viral --threads 8 --minimum-hit-groups 3 --report-minimizer-data --report SRR22841133.report --paired SRR22841133_nonhuman.1.fa SRR22841133_nonhuman.2.fa > SRR22841133.kraken2
        
        






