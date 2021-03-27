import dvc.api
from typing import Dict, List
import numpy as np


def get_dna2vec() -> Dict[str, List[float]]:
    """Stream pretrained dna2vec embedding from repo"""

    def get_coefs(word, *arr):
        return word, np.asarray(arr, dtype="float32")

    with dvc.api.open(
        "pretrained/dna2vec-20161219-0153-k3to8-100d-10c-29320Mbp-sliding-Xat.w2v",
        repo='https://github.com/pnpnpn/dna2vec',
    ) as fd:
        embedding = dict(
            get_coefs(*o.rstrip().rsplit(" ")) for o in fd
        )
    return embedding
