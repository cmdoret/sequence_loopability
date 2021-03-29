"""
This acts as a config file centralizing path definitions for files that will be
accessed throughout the pipeline. This allows to change them easily at a single
place when experimenting.
"""

from pathlib import Path

# General directories
data_dir = Path(__file__).parents[1] / "data"
in_dir = data_dir / "input"
out_dir = data_dir / "out"

# Commonly used files
genome = in_dir / "genome" / "saccer3_sgd_2mu.fa"
tokenizer = in_dir / "embeddings" / "yeast_bpe_tokenizer.json"
embedding = in_dir / "embeddings" / "yeast_w2v_embedding.vec"
loops_bed = in_dir / "hic" / "loops.bed"
loops_extracted = out_dir / "loops_extracted.pkl"
features = out_dir / "loops_features.pkl"