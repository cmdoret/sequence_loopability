from typing import Dict, List, Callable, Optional, Any
import dvc.api
import numpy as np
import pyfastx


def get_dna2vec() -> Dict[str, List[float]]:
    """Stream pretrained dna2vec embedding from repo"""

    def get_coefs(word, *arr):
        return word, np.asarray(arr, dtype="float32")

    with dvc.api.open(
        "pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v",
        repo="https://github.com/pnpnpn/dna2vec",
    ) as fd:
        embedding = dict(get_coefs(*o.rstrip().rsplit(" ")) for o in fd)
    return embedding


class SequenceIngester:
    """Iterable feeding chunks of sequence data from a fasta
    
    Parameters
    ----------
    fasta : str
        Path to input genome.
    chunksize : Optional[int]
        Size of chunks to yield.
    func : Optional[Callable[[str], Any]]
        Function to apply to each sequence chunk.
        By default, the chunk sequence is returned.

    """

    def __init__(
        self,
        fasta: str,
        func: Optional[Callable[[str], Any]] = None,
        chunksize: int = 5000,
    ):
        self.fa = pyfastx.Fasta(fasta, uppercase=True)
        self.chunksize = chunksize
        if func is None:
            self.func = lambda x: x
        else:
            self.func = func

    def __iter__(self):
        for chrom in self.fa:
            chrom_len = len(chrom.seq)
            for start in range(0, chrom_len, self.chunksize):
                chunk = self.fa.fetch(
                    chrom.name, (start, start + self.chunksize)
                )
                yield self.func(chunk)

