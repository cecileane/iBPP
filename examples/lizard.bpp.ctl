          seed =  -1

       seqfile = lizard.txt
      Imapfile = lizard.Imap.txt
       outfile = out.txt
      mcmcfile = mcmc.txt

* speciesdelimitation = 0 * fixed species tree
 speciesdelimitation = 1 0 5    * speciesdelimitation algorithm0 and finetune(e)
*  speciesdelimitation = 1 1 2 1  * speciesdelimitation algorithm1 finetune (a m)

  uniformrootedtrees = 1         * 0 means uniform labeled histories

  species&tree = 5  tri  cow  con  und  woo
                     4    3    4    5    1
                 (((tri, cow), con), (und, woo));

       usedata = 1    * 0: no data (prior); 1:seq like
         nloci = 1  * 29    * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 2 1000   # gamma(a, b) for theta
      tauprior = 2 1000   # gamma(a, b) for root tau & Dirichlet(a) for other tau's

*       finetune = 0: 5 0.0005 0.002  0.0005 0.5 0.2 1.0  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

       finetune = 1: .01 .01 .01 .01 .01 .01 .01 .01  # auto (0 or 1): finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1
        burnin = 50000
      sampfreq = 5
       nsample = 100000
