RAW_DIR="../raw_fastq"
INDEX="../genome_index/GRCh38_index"

for fq in $RAW_DIR/*.fastq.gz
do
    base=$(basename "$fq" .fastq.gz)
    echo " Aligning: $base"

    # HISAT2 alignment
    hisat2 -p 8 -x $INDEX -U $fq -S ${base}.sam 2> ${base}_hisat2.log

    # Convert SAM → BAM
    samtools view -bS ${base}.sam > ${base}.bam

    # Sort BAM
    samtools sort ${base}.bam -o ${base}.sorted.bam

    # Index BAM
    samtools index ${base}.sorted.bam

    # Remove .sam to save space
    rm ${base}.sam

    echo "✨ Done: $base"
done

echo " ALL ALIGNMENTS COMPLETE!"
