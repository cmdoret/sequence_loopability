# In this file we extract histone and sequence data corresponding
# to loop regions
import numpy as np
import pandas as pd
import pickle
import pyBigWig as pb
import pyfastx
from seqloops.config import data_dir
import seqloops.config as cfg

# Inputs:
genome_infile = cfg.genome
loops_bed_infile = cfg.loops_bed
normcov_dir = cfg.out_dir / "chip" / "cov" / "norm"
# Outputs:
loops_seq_outfile = cfg.loops_extracted

# Open all bigwigs and store filehandles in a dict
bigwigs = {}
for bw in normcov_dir.glob("*bw"):
    mark = bw.with_suffix("").name.split("_")[1]
    bigwigs[mark] = pb.open(str(bw))

# Loop on all loop anchors, and extract coverage from
# each bigwig in the regions.
loops_df = pd.read_csv(
    loops_bed_infile, sep="\t", names=["chrom", "start", "end", "status"]
)


# Save the track of each each histone mark into a separate column
for m, bw in bigwigs.items():
    loops_df[m] = loops_df.apply(
        lambda row: bw.values(row.chrom, row.start, row.end), axis=1
    )


# Retrieve DNA sequence in the interval of each loop anchor
# NOTE: Unlike bigwig, pyfastx uses 1-based coords, hence the +1 below.
fa = pyfastx.Fasta(str(genome_infile))
loops_df["seq"] = loops_df.apply(
    lambda r: fa.fetch(r.chrom, (r.start + 1, r.end + 1)), axis=1
)

# Save table as a pickle (binary) file
with open(loops_seq_outfile, "wb") as fd:
    pickle.dump(loops_df, fd)
