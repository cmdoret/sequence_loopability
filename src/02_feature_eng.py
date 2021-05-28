"""This script is used to generate new features from the extracted data. This
includes summary statistics, embedded sequences, etc...
"""

import pickle
from typing import List
import numpy as np
from scipy.stats import entropy
import tokenizers as tok
from gensim.models import KeyedVectors
import seqloops.config as cfg

# embedding + tokenizer + loops -> features
# Inputs:
loops_seq_infile = cfg.loops_extracted
tokenizer_infile = cfg.tokenizer
embedding_infile = cfg.embedding
# Outputs:
features_outfile = cfg.features

# Load embedding trained on the yeast genome
w2v = KeyedVectors.load_word2vec_format(embedding_infile, binary=False)

# Tokenize input loops and compute their embedding vectors as the sum of
# their embedding tokens
tokenizer = tok.Tokenizer(tok.models.BPE())
tokenizer = tokenizer.from_file(str(tokenizer_infile))

# Load and encode input sequences with precomputed embedding
loops = pickle.load(open(loops_seq_infile, "rb"))
# Each sequence becomes a list of string tokens
encoded = loops.seq.apply(lambda s: tokenizer.encode(s).tokens)


def average_embedding(seq: List[str]) -> np.ndarray:
    """Average embedding vectors of all words in the sequence"""
    vec = np.zeros(w2v.vector_size)
    valid_words = 0
    # Ignore words if they are absent
    for word in seq:
        if w2v.has_index_for(word):
            vec += w2v.get_vector(word)
            valid_words += 1
    vec /= valid_words
    return vec


# Convert tokens to their embedding vectors
loops["embedding"] = [average_embedding(seq) for seq in encoded]

# Transforming and standardizing features to analyze later on.

# Basic sequence stats, GC, entropy
loops["GC"] = loops.seq.apply(lambda x: (x.count("G") + x.count("C")) / len(x))
loops["GC_bias"] = loops.seq.apply(
    lambda x: (x.count("G") - x.count("C")) / (x.count("G") + x.count("C"))
)
loops["entropy"] = loops.seq.apply(
    lambda x: entropy([x.count(base) / len(x) for base in ["A", "C", "T", "G"]])
)
with open(features_outfile, "wb") as fd:
    pickle.dump(loops, fd)
