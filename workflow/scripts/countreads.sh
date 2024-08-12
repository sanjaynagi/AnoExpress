#!/bin/bash

# Set reference files
GFF_FILE="resources/reference/Anopheles-gambiae-PEST_BASEFEATURES_AgamP4.12.gff3"

# Sample names
SAMPLES=("Kisumu1" "Kisumu2" "Kisumu3" "Kisumu4" "G28-BusiaSurvivors1" "G28-BusiaSurvivors2" "G28-BusiaSurvivors3" "G28-BusiaSurvivors4")

# Loop through each sample and run the tools
for SAMPLE in "${SAMPLES[@]}"
do
    BAM_FILE="results/alignments/${SAMPLE}.hisat2.bam"
    R1_FASTQ="resources/reads/${SAMPLE}_1.fastq"
    R2_FASTQ="resources/reads/${SAMPLE}_2.fastq"
    
    FEATURECOUNTS_OUTPUT="${SAMPLE}_featureCounts_output.txt"
    HTSEQCOUNT_OUTPUT="${SAMPLE}_htseq_counts_output.txt"
    KALLISTO_OUTPUT="${SAMPLE}_kallisto_output"

    # Run featureCounts
    echo "Running featureCounts for ${SAMPLE}..."
    featureCounts -a $GFF_FILE -o $FEATURECOUNTS_OUTPUT -p -B -t exon -g Parent -C $BAM_FILE
    if [ $? -ne 0 ]; then
        echo "featureCounts failed for ${SAMPLE}"
        exit 1
    fi

    # Run htseq-count
    echo "Running htseq-count for ${SAMPLE}..."
    htseq-count -f bam -r pos -s no -t exon -i Parent $BAM_FILE $GFF_FILE > $HTSEQCOUNT_OUTPUT
    if [ $? -ne 0 ]; then
        echo "htseq-count failed for ${SAMPLE}"
        exit 1
    fi

done

echo "All tools have completed successfully for all samples."
