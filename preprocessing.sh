############################## Align ONT data to Cassava genome with Minimap2###############################################

# -uf because it is Direct RNA-seq
# default limitation on intron length -G = 200Kb
gendir=/home/gvw266/groupdirs/SCIENCE-PLEN-Marquardt_lab/Swetha/Genomes/Cassava/v6.54
minimap2 -t 8 -ax splice -uf -k14 -L --cs --secondary=no ${gendir}/Manihot_esculenta_v6.54.fa ../FASTQ/TM12_part2.fastq.gz > TM12_part2_ont.sam
minimap2 -t 8 -ax splice -uf -k14 -L --cs --secondary=no ${gendir}/Manihot_esculenta_v6.54.fa ../FASTQ/Cold3h.fastq.gz > Cold3h_ont.sam


# Convert SAM to sorted BAM:
samtools view -hu TM12_part2_ont.sam | samtools sort - -o TM12_part2_ont.bam
samtools view -hu Cold3h_ont.sam | samtools sort - -o Cold3h_ont.bam

# Filter out unmapped reads and reads with low MAPQ:
samtools view -hb -q 10 -F 4 TM12_part2_ont.bam > TM12_part2_ont_mapq.bam
samtools view -hb -q 10 -F 4 Cold3h_ont.bam > Cold3h_ont_mapq.bam

# Count MAPQ filtered reads:
echo $(samtools flagstat TM12_part2_ont_mapq.bam | sed -n '1p' | awk '{print $1}')
echo $(samtools flagstat Cold3h_ont_mapq.bam | sed -n '1p' | awk '{print $1}')



############################## Align STRIPE-Seq data to Cassava genome with STAR ###############################################

# R1 reads are expected to begin with 8nt UMI followed by TATAGGG 

# Merge FASTQ files sequenced on different lanes:
for sample in TM12 Cold3h; do 
  for rep in 1 2 3 4; do 
    echo $sample $rep && 
    zcat ${sample}-rep${rep}*fastq.gz | 
    gzip > ../02-merged/${sample}_rep${rep}_R1.fastq.gz; 
  done; 
done


# Count raw reads:
for file in *fastq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 )); 
done

# Process UMI:
for file in *R1.fastq.gz; do 
  echo $file &&
  umi_tools extract --stdin=$file --bc-pattern=NNNNNNNN --stdout=../UMI/${file/R1.fastq.gz/UMI.fq.gz};
done

# Run FastQC:
for file in *fastq.gz; do fastqc $file; done


# Trim TATAGGG from 5' end and Illumina adapters from 3' end:
for file in *UMI.fastq; do 
  echo $file
  /home/gvw266/Documents/Softwares/TrimGalore-0.6.6/trim_galore -j 4 --illumina --clip_R1 7 --length 15 --gzip --no_report_file $file
done

# Count trimmed reads:
for file in *trimmed.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 )); 
done


# Align to Cassava genome using STAR:
genomedir=/home/gvw266/groupdirs/SCIENCE-PLEN-Marquardt_lab/Swetha/Genomes/Cassava/v6.54

files=`find ../FASTQ/ -maxdepth 1 -type f -name "*trimmed.fq" | sort`
for file in $files; do
  echo $file
  STAR \
    --genomeDir ${genomedir}/STAR_v6.54 \
    --readFilesIn $file \
    --runThreadN 4 \
    --outFileNamePrefix ${file/trimmed.fq.gz/} \
    --outSAMmultNmax 1 \
    --alignEndsType Extend5pOfRead1 \
    --outSAMtype BAM Unsorted \
    --sjdbGTFfile ${genomedir}/Manihot_esculenta_v6.54_woScaffold.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --outSAMstrandField intronMotif; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Count aligned reads:
for file in *bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' \
    | awk '{print $1}'); 
done

# Sort BAM files and remove low MAPQ reads:
for file in *bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Count filtered reads:
for file in *mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Generate indexes:
for file in *bam; do samtools index $file; done

# Deduplicate on UMI:
for file in *mapq.bam; do
  echo $file &&
  umi_tools dedup --stdin=${file} --stdout=${file/.bam/_dedup.bam};
done

# Count deduplicated reads:
for file in *dedup.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Merge replicates:
samtools merge Cold3h_merged_UMI_mapq_dedup.bam Cold3h_rep*bam
samtools merge TM12_merged_UMI_mapq_dedup.bam TM12_rep*bam

# Make stranded Bedgraph files (no strand switch):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="fw" || n="rev"; 
  for file in *dedup_sorted.bam; do 
    sample=${file/_UMI_mapq_dedup_sorted.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
for f1 in *fw.bg; do 
  f2=${f1/fw/rev} && 
  outf=${f1/fw.bg/fw_rev.bedgraph.gz} && 
  echo $f1 "+" $f2 "=" $outf && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | 
  cat $f1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outf; 
done

# Normalize to 1M tags:
for file in *bedgraph.gz; do 
  norm=$( zcat $file | sed '/^[#t]/d' | \
    awk 'BEGIN{SUM=0}{SUM+=sqrt($4^2)*($3-$2)}\
      END{print SUM / 1000000}' ) && 
  echo $file $norm && zcat $file | 
  awk -v norm=$norm 'BEGIN{OFS="\t"}{if ($0~/^[#t]/) print $0; \
    else print $1, $2, $3, $4 / norm}' | 
  gzip > ${file/.bedgraph.gz/_norm1M.bedgraph.gz}; 
done

############################## Align Quant-Seq data to Cassava genome with STAR ###############################################

genomedir=/home/gvw266/groupdirs/SCIENCE-PLEN-Marquardt_lab/Swetha/Genomes/Cassava/v6.54

# Count raw reads:
for file in *fastq.gz; do echo $file $(( $(zcat $file | wc -l | awk '{print $1}') / 4 )); done

# Align to Cassava genome in transcriptome-guided mode:
files=`find ../FASTQ/ -maxdepth 1 -type f -name "*.fastq" | sort`
for file in $files; do
  echo $file && 
  STAR \
    --genomeDir ${genomedir}/STAR_v6.54 \
    --readFilesIn $file \
    --runThreadN 4 \
    --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 \
    --alignEndsType Extend5pOfRead1 \
    --clip3pAdapterSeq AGATCGGAAGAGC \
    --outSAMtype BAM Unsorted \
    --sjdbGTFfile ${genomedir}/Manihot_esculenta_v6.54_woScaffold.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --outSAMstrandField intronMotif; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab

# Count aligned reads:
for file in *bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Sort BAM files and remove low MAPQ reads:
for file in *bam; do 
  echo $file && 
  samtools view -hu -q 10 $file | 
  samtools sort - -o ${file/.bam/_mapq.bam}; 
done

# Count filtered reads:
for file in *mapq.bam; do 
  echo $file $(samtools flagstat $file | sed -n '1p' | \
    awk '{print $1}'); 
done

# Merge replicates:
samtools merge Cold3h_merged_mapq.bam Cold3h_rep*mapq.bam
samtools merge TM12_remerged_mapq.bam TM12_rep*mapq.bam

# Make stranded Bedgraph files (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "-" ] && n="fw" || n="rev"; 
  for file in *mapq.bam; do 
    sample=${file/_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files for the same sample:
for f1 in *fw.bg; do 
  f2=${f1/fw/rev} && 
  outf=${f1/fw.bg/fw_rev.bedgraph.gz} && 
  echo $f1 "+" $f2 "=" $outf && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | 
  cat $f1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | gzip > $outf; 
done

# Normalize to 1M tags:
for file in *_fw_rev.bedgraph.gz; do 
  norm=$( zcat $file | sed '/^[#t]/d' | \
    awk 'BEGIN{SUM=0}{SUM+=sqrt($4^2)*($3-$2)}\
      END{print SUM / 1000000}' ) && 
  echo $file $norm && zcat $file | 
  awk -v norm=$norm 'BEGIN{OFS="\t"}{if ($0~/^[#t]/) print $0; \
    else print $1, $2, $3, $4 / norm}' | 
  gzip > ${file/.bedgraph.gz/_norm1M.bedgraph.gz}; 
done


############################## Align RNA-Seq data to Cassava genome with STAR ###############################################

genomedir=/home/gvw266/groupdirs/SCIENCE-PLEN-Marquardt_lab/Swetha/Genomes/Cassava/v6.54

# Count raw reads:
for file in *_1.fq.gz; do 
  echo $file $(( $(zcat $file | wc -l | 
    awk '{print $1}') / 4 ));
done


# Align to Cassava genome in transcriptome-guided mode:
#files=`find ./ -maxdepth 1 -type f -name "Cold*.sam"`
for f1 in ../FASTQ/sR6*_1.fq; do 
  f2=${f1/_1/_2} && 
  echo $f1 $f2 && 
  STAR \
    --genomeDir ${genomedir}/STAR_v6.54 \
    --readFilesIn $f1 $f2 \
    --runThreadN 4 \
    --outFileNamePrefix ${f1/R1.fq.gz/} \
    --outSAMmultNmax 1 \
    --alignEndsType Local \
    --clip3pAdapterSeq AGATCGGAAGAGC AGATCGGAAGAGC \
	--clip3pAdapterMMp 0.1 0.1 \
    --outSAMtype BAM Unsorted \
    --sjdbGTFfile ${genomedir}/Manihot_esculenta_v6.54_woScaffold.gff3 \
    --sjdbGTFtagExonParentTranscript Parent \
    --outSAMstrandField intronMotif; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rm -r *STARgenome *STARtmp *out *tab


# Count aligned read pairs:
for file in sR{1..8}*.bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{printf("%i", $1/2)}'); 
done

# Filter for reads with high MAPQ in proper pairs:
for file in *.bam; do 
  echo $file && 
  samtools view -hb -f 2 -q 10 $file -o ${file/.bam/_mapq.bam}; 
done

# Count MAPQ-filtered read pairs:
for file in *mapq.bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{printf("%i", $1/2)}'); 
done

# Run samtools fixmate (required by samtools markdup):
for file in *mapq.bam; do 
  echo $file && 
  samtools fixmate -m $file ${file/.bam/_fixmate.bam}; 
done


# Sort by coordinates:
for file in *fixmate.bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Deduplicate on start coordinates of paired reads:
for file in *sorted.bam; do 
  echo $file && 
  samtools markdup -r $file ${file/.bam/_dedup.bam}; 
done

# Count deduplicated read pairs:
for file in *dedup.bam; do 
  echo $file $(samtools flagstat $file | \
    sed -n '9p' | awk '{printf("%i", $1/2)}');
done

# Load deduplicated BAM files into R session and produce stranded normalized Bedgraph files using convert_RNAseq_BAM_to_normalized_bedGraph.R;




