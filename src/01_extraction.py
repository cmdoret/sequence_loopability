# In this file we extract histone and sequence data corresponding
# to loop regions
import numpy as np
import pandas as pd
import pickle
import pyBigWig as pb
import pyfastx
from seqloops.config import data_dir

# Inputs:
normcov_dir = data_dir / "out/chip/cov/norm"
loops_bed = data_dir / "input/hic/loops.bed"
genome_file = data_dir / "input/genome/saccer3_sgd_2mu.fa"
# Outputs:
out_file = data_dir / "out/loops_merged.pkl"

# Open all bigwigs and store filehandles in a dict
bigwigs = {}
for bw in normcov_dir.glob("*bw"):
    mark = bw.with_suffix("").name.split("_")[1]
    bigwigs[mark] = pb.open(str(bw))

# Loop on all loop anchors, and extract coverage from
# each bigwig in the regions.
loops_df = pd.read_csv(
    loops_bed, sep="\t", names=["chrom", "start", "end", "status"]
)


def load_bw_region(
    bw: "pb.BigWigFile", chrom: str, start: int, end: int
) -> np.ndarray:
    """Return the value array from all intervals in the requested region"""
    intervals = bw.intervals(chrom, start, end)
    coverage = [i[2] for i in intervals]
    return coverage


# Save the track of each each histone mark into a separate column
for m, bw in bigwigs.items():
    loops_df[m] = loops_df.apply(
        lambda row: load_bw_region(bw, row.chrom, row.start, row.end), axis=1
    )


# Retrieve DNA sequence in the interval of each loop anchor
fa = pyfastx.Fasta(str(genome_file))
loops_df["seq"] = loops_df.apply(
    lambda r: fa.fetch(r.chrom, (r.start, r.end)), axis=1
)

# Save table as a pickle (binary) file
with open(out_file, "wb") as fd:
    pickle.dump(loops_df, fd)