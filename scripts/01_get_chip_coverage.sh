!/usr/bin/env bash
# Use all bams in input directory to generate coverage (bigwig) files
# normalized by their respective controls (input) and library size
# Requires deeptools and a design.csv as used in nf-core pipelines
# cmdoret, 20210327

BAM_DIR='data/out/chip/bam'
OUT_DIR='data/out/chip/cov'
DESIGN='data/input/chip/design.csv'
mkdir -p "${OUT_DIR}"/norm
mkdir -p "${OUT_DIR}"/inputs


# Each line contains sample and input names
while read -a line; do
  # Compute log ratio of sample/input, normalize by library size
  bamCompare -b1 ${BAM_DIR}/${line[0]}*bam \
             -b2 ${BAM_DIR}/${line[1]}*bam \
             -o "${OUT_DIR}/norm/${line[0]}.bw" \
             -of "bigwig" \
             --binSize 1 \
             --smoothLength 15 \
             --extendReads 100 \
             --centerReads \
             --numberOfProcessors "max/2"

done < <(tail -n+2 "$DESIGN" | awk -vFS=',' '$NF != "" {print $1"_R"$2,$NF}')

# Compute raw coverage for each input file
while read line; do
 bamCoverage --bam ${BAM_DIR}/${line}*bam \
             -o "${OUT_DIR}/inputs/${line}.bw" \
             --binSize 1 \
             -p10 \
             --smoothLength 15 \
             --extendReads 100 \
             --centerReads \
             --numberOfProcessors "max/2"

done < <(tail -n+2 "$DESIGN" | awk -vFS=',' '$NF == "" {print $1"_R"$2}')
