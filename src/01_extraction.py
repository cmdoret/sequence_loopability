# In this file we extract histone and sequence data corresponding
# to loop regions
import pathlib
import numpy as np
import pyBigWig as pb

# Inputs:
normcov_dir = pathlib.Path('data/out/chip/cov/norm')
loops_bed = 'data/input/hic/loops.bed'
# Outputs:
out_file = open('data/out/chip/cov/merged.tsv', 'w')

# Open all bigwigs and store filehandles in a dict
bigwigs = {}
for bw in norm_dir.glob('*bw'):
    mark = bw.with_suffix('').name.split('_')[1]
    bigwigs[mark] = pb.open(bw):

# Loop on all loop regions, and extract coverage from
# each bigwig in the regions.
loops_df = pd.DataFrame(
    loops_bed,
    sep='\t',
    names=['chrom', 'start', 'end', 'status']
)


def load_bw_region(bw: pb.BigWigFile, chrom: str, start: int, end: int) -> np.ndarray:
    """Return the value array from all intervals in the requested region"""

for m, bw in bigwigs.items():
    loops_df[m] = loops.apply(lambda row: , axis=1)
while True:
    try:
        line = {m: next(fd) for m, fd in bigwigs.items()}
        out.write()
    except StopIteration:
        out_file.close()

