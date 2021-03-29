"""This script is used to train the word2vec embedding vectors on the input genome."""

import sys
import pyfastx
import tokenizers as tok
from gensim.models import Word2Vec
import pyfastx
from seqloops.config import data_dir
from seqloops.data_loaders import SequenceIngester

GENOME = sys.argv[1]
EMBEDDING_FILE = sys.argv[2]

# Important parameters for tokenization and embedding
# Upper number of words in tokenization
vocab_size = 30000
# Vector dimension for output embedding
embedding_dim = 200
# Word neighbourhood to use for training
context = 4
# Minimum frequency of a word to be included in embedding
min_word_count = 10

# Generate tokens from genome. Most frequent adjacent nucleotides
# sequences (bytepairs) are merged until target vocab size is reached.
tokenizer = tok.Tokenizer(tok.models.BPE())
bpe_trainer = tok.trainers.BpeTrainer(vocab_size=vocab_size, min_frequency=10)
tokenizer.train_from_iterator(
    SequenceIngester(GENOME), trainer=bpe_trainer,
)

# Save tokenizer model for potential future use
tokenizer.save(str(data_dir / "out/bpe_tokenizer.json"), pretty=True)

# This ingester yields the tokenized chunks
token_ingester = SequenceIngester(
    GENOME, func=lambda x: tokenizer.encode(x).tokens
)

# Learn embedding vectors on genome
# Note: Like in dna2vec, we train in  skipgram (sg) mode
# Like in dnavec, we will feed "sentences" of 5kb. Instead of
# k-mers, we use byte-pair encoded tokens as "words".
model = Word2Vec(
    sentences=token_ingester,
    size=embedding_dim,
    window=context,
    min_count=min_word_count,
    workers=4,
    sg=1,
    iter=1,
)

# Save embedding as a text file to use later
model.wv.save_word2vec_format(EMBEDDING_FILE)
