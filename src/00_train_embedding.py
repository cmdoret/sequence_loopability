"""This script is used to train the word2vec embedding vectors on the input genome.
1. The input genome is split into 5kb chunks (for performance). These will
represent "sentences".

2. Bytepair encoding is performed on the sentences, iteratively merging the most
frequent adjacent byte pairs (oligonucleotides==words) until a target vocab. size
is reached.

3. The tokenized genome sequence (each word is replaced by a numeric ID) is fed
to the word2vec model in skipgram mode. The model learns to predict the
context (neighbouring words) of each word. The weight of the model are then saved
and represent embedding vectors for each word.
"""

import os
import tokenizers as tok
from gensim.models import Word2Vec
import seqloops.config as cfg
from seqloops.data_loaders import SequenceIngester

# Inputs
genome_infile = cfg.genome
# Outputs
tokenizer_outfile = cfg.tokenizer
embedding_outfile = cfg.embedding

# Important parameters for tokenization and embedding
# Upper number of words in tokenization (min count will be limiting before this)
vocab_size = 100000
# Vector dimension for output embedding
embedding_dim = 200
# Word neighbourhood to use for training
context = 4
# Minimum frequency of a word to be included in embedding
min_word_count = 10

# Generate tokens from genome. Most frequent adjacent nucleotides
# sequences (bytepairs) are merged until target vocab size is reached.
tokenizer = tok.Tokenizer(tok.models.BPE())
bpe_trainer = tok.trainers.BpeTrainer(
    vocab_size=vocab_size, min_frequency=min_word_count,
)
tokenizer.train_from_iterator(
    SequenceIngester(str(genome_infile)), trainer=bpe_trainer,
)

# Save tokenizer model for potential future use
tokenizer.save(str(tokenizer_outfile), pretty=True)

# This ingester yields the tokenized chunks
token_ingester = SequenceIngester(
    str(genome_infile), func=lambda x: tokenizer.encode(x).tokens
)

# Learn embedding vectors on genome
# Note: Like in dna2vec, we train in  skipgram (sg) mode
# Like in dnavec, we will feed "sentences" of 5kb. Instead of
# k-mers, we use byte-pair encoded tokens as "words".
os.environ["TOKENIZERS_PARALLELISM"] = "false"
model = Word2Vec(
    sentences=token_ingester,
    vector_size=embedding_dim,
    window=context,
    min_count=min_word_count,
    workers=4,
    sg=1,
    epochs=5,
)

# Save embedding as a text file to use later
model.wv.save_word2vec_format(embedding_outfile)
