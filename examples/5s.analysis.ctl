* example file to analyze data: genes and traits   *
* run at the command line:  ibpp 5s.trait.bpp.ctl  *

          seed = 1234     * change this to positive integer to fix the seed

       seqfile = 5s.sequences.txt
      Imapfile = 5s.Imap.txt
     traitfile = 5s.morph.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 5    * speciesdelimitation algorithm0 and finetune(e)
* speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)
  speciesdelimitation = 1 1 2 1 0 1
* speciesdelimitation=1, algorithm=1, finetune (a m) = 2,1, diagnosis=0 (no), 
* startwithroot=1, i.e. the starting tree has to have speciation at the root

  uniformrootedtrees = 1         * 0 means uniform labeled histories
  species&tree = 5  A  B  C  D  E
                    3  1  2  1  1
                  (((A, B), C), (D, E));
* compare to tree for simulation, where we need to specify 
* population sizes (like #0.002) and branch lengths (like :0.005):
* (((A #0.002, B) : 0.005 #.002, C #0.002) : 0.01 #.002, (D, E) :.015 #.002) : 0.02 #0.002;

*      usedata = 0    * 0: no data (prior); 1:seq & trait like
    useseqdata = 1    * 0: no seq data;     1:seq like
  usetraitdata = 1    * 0: no trait data;   1:trait like
         nloci = 5    * number of data sets in seqfile
       ntraits = 2    * number of trait variables
         nindT = 9    * total # individuals for which trait data is available
** required only if it differs from the number of individuals with genetic data

     cleandata = 0    * remove sites with ambiguity data? (1:yes, 0:no)

    thetaprior = 2 2000    # gamma(a, b) for theta
      tauprior = 2 20000 1  # gamma(a, b) for root tau & Dirichlet(a) for other tau's
           nu0 = 0         # parameters for prior on traits
        kappa0 = 0         # nu0=0 and kappa0=0 for non-informative prior

*      heredity = 0 4 4   # (0: No variation, 1: estimate, 2: from file) & a_gamma b_gamma (if 1)
*     locusrate = 0 2.0   # (0: No variation, 1: estimate, 2: from file) & a_Dirichlet (if 1)
* sequenceerror = 0 0 0 0 0 : 0.05 1   # sequencing errors: gamma(a, b) prior

      finetune = 0: 2 0.0002 0.001 0.03 0.1 0.005 1.0 .1
      # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr, traitHsq
*     finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01 .01

         print = 1
        burnin = 5 # 5000
      sampfreq = 2
       nsample = 50 # 50000


*** Note: Make your window wider (144 columns) before running the program.
