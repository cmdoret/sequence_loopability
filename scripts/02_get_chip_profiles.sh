# List of commands used to plot profiles
# Bigwig files must have been generated
# Requires deeptools.
# cmdoret, 20210327

# Should contain norm and inputs subfolders
# norm contains input-normalized bigwigs, inputs contains inputs bigwigs
IN_DIR='data/out/chip/cov'
OUT_DIR='data/out/chip/profiles'
GENES='data/input/genome/saccer3_sgd_genes.bed'
mkdir -p ${OUT_DIR}

# Normalized histone mark coverage scaled by gene body (nucleosomes not aligned)
computeMatrix scale-regions -p 10 \
    --smartLabels \
    -R ${GENES} \
    --beforeRegionStartLength 3000 \
    --afterRegionStartLength 3000 \
    --regionBodyLength 5000 \
    --outFileName ${OUT_DIR}/norm_scaled_allmarks.mat.gz \
    -S ${IN_DIR}/norm/*bw

plotProfile -m ${OUT_DIR}/norm_scaled_allmarks.mat.gz \
    -out ${OUT_DIR}/norm_scaled_allmarks.png \
    --numPlotsPerRow 4 \
    --plotTitle "All marks, input-normalized, scaled, whole genome"


# Nucleosome occupancy (input signal, not scaled to gene body)
computeMatrix reference-point -p 10 \
    --smartLabels \
    -R ${GENES} \
    --binSize 1 \
    --outFileName ${OUT_DIR}inputs_refpoint.mat.gz \
    -S ${IN_DIR}/inputs/*bw

plotProfile -m ${OUT_DIR}/inputs_refpoint.mat.gz \
    -out ${OUT_DIR}/inputs_refpoint.png \
    --numPlotsPerRow 4 \
    --plotTitle "Inputs, unscaled, whole genome"
