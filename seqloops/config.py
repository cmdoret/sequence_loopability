"""
This acts as a config file centralizing path definitions for files that will be
accessed throughout the pipeline. This allows to change them easily at a single
place when experimenting.
"""

from pathlib import Path

# General directories
data_dir = Path(__file__).parents[1] / "data"
out_dir = Path(__file__).parents[1] / "outputs"

# Commonly used files
genome = data_dir / "genome" / "saccer3_sgd_2mu.fa"
loops_bed = data_dir / "hic" / "loops_tes_ctl.bed"
# loops_bed = data_dir / "control" / "GAL4_vs_neg.bed"
tokenizer = out_dir / "embeddings" / "yeast_bpe_tokenizer.json"
embedding = out_dir / "embeddings" / "yeast_w2v_embedding.vec"
loops_extracted = out_dir / "loops_extracted.pkl"
features = out_dir / "loops_features.pkl"
