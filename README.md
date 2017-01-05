iBPP: integration of genes and traits for Bayesian Phylogenetics and Phylogeography
-------

- Claudia Solís-Lemus, L. Lacey Knowles and Cécile Ané (2014). 
Bayesian species delimitation combining multiple genes and traits in a unified framework. 
*Evolution* 69(2):492-507.
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12927.svg)](https://doi.org/10.5281/zenodo.12927)

- Ziheng Yang and Bruce Rannala (2010). 
Bayesian species delimitation using multilocus sequence data. 
*PNAS* 107:9264–9269. [software](http://abacus.gene.ucl.ac.uk/software/).


For a **Windows** executable and instructions, see [here](man/winexe.md).

To install and compile on Linux or Mac:

- download and unzip the package, or clone it using git:
  `git clone https://github.com/cecileane/iBPP.git`
- navigate to the source directory: `cd iBPP/src/`
- compile with `gcc -o ibpp -O3 bpp.c tools.c -lm`
- move the executable to a directory in your PATH, typically `~/bin/`,
  to run iBPP from anywhere: `mv ibpp ~/bin/`
- check that the executable runs on a very short example:
  * navigate to the example folder: `cd ../examples/`
  * then run ibpp: `ibpp 5s.analysis.ctl`

See the [version history](versionHistory.txt)
for recently added features.