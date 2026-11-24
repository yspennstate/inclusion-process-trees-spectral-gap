# Caterpillars and the spectral gap of the symmetric inclusion process

Yitzchak Shmalo

Kim and Sau showed that the spectral gap of the symmetric inclusion process equals that of one particle when all site
weights are at least one, and that this fails on tori of dimension at least two with small weights. We ask on which
trees it can fail. On every tree that is not a caterpillar we give conductances and site weights for which two
particles relax strictly more slowly than one. On a caterpillar two walkers killed when they come within distance one
die strictly faster than one walker relaxes, because they cannot exchange places. Combined with the small-weight limit
of Kim and Sau, this shows that a tree is a caterpillar if and only if, for all conductances, all weight proportions and
every number of particles, the gap equals that of one particle once the weights are small enough. We conjecture that on a caterpillar two particles never relax more slowly than one, at any weights; on stars
with equal conductances this follows from work of Feng, Miclo and Wang and of Kim and Sau.

`main.tex` and `main.pdf` are the paper; archived at Zenodo, doi:10.5281/zenodo.24510421.

## The spider check

`code/check_spider.py` verifies Lemma 9 on the spider S(2,2,2). It computes all eigenvalues of the symmetrised
generator for one and two particles at 50 significant digits and counts the eigenvalues below 1/200000 exactly by an
LDL^T inertia count in rational arithmetic. It needs `mpmath` and `sympy`:

    cd code
    python check_spider.py

Its output, recorded in `code/check_spider.out`, ends with:

    gap_1 = 5.9938549174820907956e-6
    gap_2 = 4.5581896672418914781e-6
    gap_2/gap_1 = 0.76047714367379180375
    exact counts below t=1/200000: k=1 1  k=2 2

Tested with mpmath 1.3.0 and sympy 1.14.0.

## License

MIT for the code; the paper is CC BY 4.0.
