#!/bin/bash
for SAMPLE in KO1_Post KO1_Pre KO2_Post KO2_Pre KO3_Post KO3_Pre WT1_Post WT1_Pre WT2_Post WT2_Pre WT3_Post WT3_Pre
do kallisto quant -i transcripts.idx -o ${SAMPLE} --single -l 50 -s 5 ${SAMPLE}_TL_S*_L002_R1_001.fastq.gz
done

