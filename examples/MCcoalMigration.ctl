          seed = -12345

       seqfile = MySeq.txt   * comment out this line if you don't want seqs
      treefile = MyTree.tre   * comment out this line if you don't want trees
      Imapfile = MyImap.txt

  species&tree = 4  A  B  C  D
                    1  2  1  1
  ((A #0.1, B #0.2) : 0.5 #.12, (C #0.3, D #0.4) :0.8 #0.34) : 1.0 #.1234;

     migration = 7   * number of pops (order fixed by program)

             A     B     C     D    ABCD    AB   CD
      A      0     1.1   1.2   1.3    0     0    -1
      B      0.1   0     1.4   1.5    0     0    -1
      C      0.2   0.4   0     1.6    0     1.7   0  
      D      0.3   0.5   0.6   0      0     1.8   0
      ABCD   0     0     0     0      0     0     0
      AB     0     0     0.7   0.8    0     0     1.9
      CD    -1    -1     0     0      0     0.9   0

   loci&length = 100 1000 * number of loci & number of sites at each locus
*     locusrate = 0.5     * alpha for gamma for locus rate
