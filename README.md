iBPP: integration of genes and traits for Bayesian Phylogenetics and Phylogeography
-------

- Claudia Solís-Lemus, L. Lacey Knowles and Cécile Ané (2014). 
Bayesian species delimitation combining multiple genes and traits in a unified framework. 
*Evolution* 69(2):492-507.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12927.svg)](https://doi.org/10.5281/zenodo.12927)

- Ziheng Yang and Bruce Rannala (2010). 
Bayesian species delimitation using multilocus sequence data. 
*PNAS* 107:9264–9269. [software](http://abacus.gene.ucl.ac.uk/software/).


See the [version history](versionHistory.txt)
for potentially new features.

For a **Windows** executable and instructions, see [here](man/winexe.md).

## Installation

*Warning*: the code compiles but does *not* work with the more recent gcc compilers.
- works well with, for example: `gcc` v5 (Linux),
  or `clang` v7.0.2, `clang` v12.0.0 (OSX)
- compiles (an executable is produced) but doesn't work with
  `gcc` v6.1.0 or `gcc` v7.5.0 (Unbuntu 18.04)

To install and compile on Linux or Mac (with Xcode installed):

- download and unzip the package, or clone it using git:
  `git clone https://github.com/cecileane/iBPP.git`
- navigate to the source directory: `cd iBPP/src/`
- compile with
  `gcc -o ibpp -O3 bpp.c tools.c -lm`.
  On Mac OS X, do
  `clang -o ibpp -O3 bpp.c tools.c -lm`
  to request the Apple compiler `clang` explicitly.
  With `gcc` 4.4.7 to 5.4.0, the option `-std=c89` is necessary:

  `gcc -o ibpp -O3 -std=c89 bpp.c tools.c -lm`

- move the executable to a directory in your PATH, typically `~/bin/`,
  to run iBPP from anywhere: `mv ibpp ~/bin/`
- check that the executable runs on a very short example:
  * navigate to the example folder: `cd ../examples/`
  * then run ibpp: `ibpp 5s.analysis.ctl`

A simulation tools is also available as part of the package.
To get it, compile with the option

`gcc -o ibpp-simul -DSIMULATION -O3 bpp.c tools.c -lm`

possibly adding the option `-std=c89` depending on your compiler.
