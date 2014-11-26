* example file to simulate data: genes and traits        *
* run at the command line: ibpp-simul 5s.trait.simul.ctl *

          seed = -1     * change this to positive integer to fix the seed

*      seqfile = 5s.seqsimul.txt      * uncomment to get sequences
*     treefile = 5s.treesimul.txt     * uncomment to get the trees
      Imapfile = 5s.Imapsimul.txt
     traitfile = 5s.morphsimul.txt

    species&tree = 3  A  B  C
                     10 10 10
  ((A #0.1, B #0.2) : 0.5 #.12, C #0.3) : 1.0 #0.4;

*   species&tree = 1  C   * these 3 lines are to simulate a single population
*                    20   * 20 individuals
*   (C #0.3) : 1.0 #0.4;

*  species&tree = 5  A  B  C  D  E  * just another example tree
*                    3  1  2  1  1
* (((A #0.002, B) : 0.005 #.002, C #0.002) : 0.01 #.002, (D, E) :.015 #.002) : 0.02 #0.002;


 loci&length = 0 1000 * number of loci & number of sites at each locus
   locusrate = 1.0    * alpha for gamma distribution of locus rates

     ntraits = 10     * number of trait variables
       nindT = 21     * total # individuals with trait data
      lambda = 1: 0.7 * between-to-within species trait variance (Pagel's lambda)
                      * 0        - lambda ~ iid uniform(0,1) across traits --default
                      * 1: value - same lambda=value for all traits
  phyloModel = 1:1.0  * phylogenetic model for species means.
                      * 0        - Brownian motion (BM) --default
                      * 1: alpha - Ornstein-Uhlenbeck (OU)
  OUmuSD = 0.2        * SD of optimal values for the OU model, across branches and traits.
                      * relative to trait scale
  traitflow = 0 0.95  * probabilities at internal nodes: that sister species have equal evolution.
                      * analog of gene flow between sister species.
                      * can also simulate phenotypic plasticity, using an artificially larger tree. 
                      * here: true species tree is (AB, C) but 5% of traits would
                      * show plasticity with a different mean in A and in B,
                      * e.g. due to different geographic locations of A and B.
        
