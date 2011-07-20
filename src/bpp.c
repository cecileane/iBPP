/* bpp.c
   Markov chain Monte Carlo coalescent program for population genetics and 
   phylogeographic data.
   
   Copyright by Ziheng Yang, since July 2002
   Modified by Cecile Ane, July 2011, to also use quantitative trait data

   UNIX gcc/icc:
   cc -o bpp -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm
   cc -o MCcoal -DSIMULATION -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm

   icc -o bpp -fast bpp.c tools.c -lm

   MAC OSX intel:
   cc -o bpp -march=pentium-m -O4 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm
   cc -o MCcoal -DSIMULATION -march=pentium-m -O4 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm


   Windows MSC++ 6.0/2008:
   cl -O2 bpp.c tools.c
   cl -O2 -FeMCcoal.exe -DSIMULATION bpp.c tools.c
*/

/*
#define SIMULATION
*/


#define NSPECIES      20        /* max # of species */
#define NS            500       /* max # of sequences and of individuals for trait data*/
#define NBRANCH       NS*2-2
#define MAXNSONS      2
/* #define NGENE         213000   */  /* max # of loci */
#define NGENE         100000     /* max # of loci */
#define NTRAIT          1000     /* max # of traits */
#define LSPNAME       30         /* # characters in sequence names */

#include "paml.h"

struct CommonInfo {
   char *z[2*NS-1], *spname[NS];
   char seqf[96], Imapf[96],outf[96], mcmcf[96], locusratef[96], heredityf[96],traitf[96];
   char cleandata, oldconP[NS*2-1];
   int ns, ls, lgene[1], model, clock, simulation;
   int seqtype, ngene, *pose, npatt, np, readpattern;
   int ntime, ncatG, ncode, print, fix_rgene, posG[1];
   double alpha, pi[4], piG[1][4], *rates, rgene[1];
   double *conP;  /* not used */
   int curconP;
   size_t sconP;        /* size (# bytes) for conditional probabilities */
   double *conPin[2], *fpatt, space[1000000];
   double a_locusrate;  /* a_locusrate is duplicated in com for simulation & in data. */
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree; /* branch lengths and log likelihood of a particular gene tree */
/* Used by the simulation program to get a random gene tree in 'tree' and 'nodes' (below).
   used by the inference program to generate the starting gene tree with coalescent times.
*/
struct TREEN {
   int father, nson, sons[2], ibranch, ipop;  /* ibranch not used? */
   double branch, age, *conP, label;          /* age is uptodate but not branch */
   char fix_age, *nodeStr;                    /* not used */
}  *nodes, *gnodes[NGENE], nodes_t[2*NS-1]; 
/*  nodes is a pointer, gnodes holds the gene trees, nodes_t is temp space
    gnodes is not allocated for the simulation program 
*/

/* sptree.nseqsp is working space, copied from data.nseqsp[locus], which holds 
   the info for each locus.  It seems feasible to remove data.nseqsp[locus].  Its 
   major use is calculation of ndesc[] in UpdateGB_SPR, and this can be done by 
   checking ipop for tips of gene tree.  The same applies to other uses of nseqsp.
   npop is the number of theta's.  
   sptree.pops[] holds the node numbers in the species tree of the npop populations;
   younger pops are before ancestral pops.  
   sptree.pop_pop_table[i][j] = 1 if pop i can coalesce in pop j (that is, if j is 
   ancestral to i) and = 0 otherwise.  For modern species i, table[i][i] = 1 if and 
   only if there is a theta for pop i.  For ancestral species, table[i][j] = 1 if j
   is ancestral to i.
*/
struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nseqsp[NSPECIES]; 
   int npop, pops[NSPECIES*2-1], pop_pop_table[NSPECIES*2-1][NSPECIES*2-1], Itree;
   int speciesdelimitation;
   int nLHistories;  /* =0 for uniform Labelled Histories; =1 for uniform trees */
   int migration;
   double M[NSPECIES*2-1][NSPECIES*2-1];
   struct TREESPN {
      char name[LSPNAME*2];
      int father, nson, sons[2];
      double age, theta;       /* age is tau */
   } nodes[2*NSPECIES-1];
}  sptree;


#ifdef SIMULATION

struct DATA {
   double e_seqerr[NSPECIES][4*4];
   int nseqerr, iseqerr[NSPECIES];
}  data;

#else

struct DATA { /* locus-specific data and gene tree information, used in lnpGB */
   int maxns, ngene, lgene[NGENE]; // data.maxns is the max # sequences determined by the control file.
   int ns[NGENE], nseqsp[NGENE][NSPECIES], ls[NGENE], npatt[NGENE];
   int root[NGENE+1], conP_offset[NGENE];
   char *z[NGENE][NS], cleandata[NGENE], *Imap, *Inames[NS];
   /* Imap[i] = j if the ith individual listed in the Imap file belongs to species j.
      Inames[i] = name of the ith individual listed in the Imap file.
   */
   double *fpatt[NGENE], *lnpGBi, *lnpDi, *locusrate, *heredity;
   double theta_prior[2], tau_prior[3];
   double a_locusrate, a_heredity, b_heredity;
   double a_seqerr[NSPECIES][4*4], e_seqerr[NSPECIES][4*4];
   int est_locusrate, est_heredity, nseqerr, iseqerr[NSPECIES];
}  data;

struct TRAITDATA { /* trait data */
  int ntrait;                // number of trait variables
  int nind;                  // # of individuals for which trait data is available
  char *traitName[NTRAIT];   // names of the trait variables
  int indSpeciesMap[NS];     /* index of the species that individual i belongs to
			        The ith individual listed in the trait file is not 
			        required to be the same as the ith individual listed
			        in the Imap or sequence file. */
  char cleantrait[NTRAIT];   /* 1 if trait has no missing entries */
  char ismissing[NTRAIT][NS];/* 1 if entry missing for trait i and individual j */
  int ni[NTRAIT];            // # of individuals with no missing data for trait i
  double y[NTRAIT][NS];/* trait data, supposed to be continuous at this point
			   first read from file than standardized */
  double ybar[NTRAIT]; /* strait (non-phylogenetic) mean taken over populations,
			   not over individuals. Used to standardize each trait */
  double scale[NTRAIT];/* strait SD from ybar. Used to standardize each trait */
}  traitdata;

struct MCMCPARAMETERS {
   int resetFinetune, burnin, nsample, sampfreq, usedata, saveconP;
   int print, moveinnode, RJalgorithm;
   double finetune[7], RJfinetune[2];
}  mcmc;                /* control parameters */

#endif

int GetOptions (char *ctlf);
int GetOptionsSimulation (char *ctlf);
int ReadSpeciesTree (FILE* fctl, char *currline);
int ReadMigrationMatrix(FILE *fctl, char *currline);
int DownSptreeSetPops(int inode);
void GetRandomGtree(int locus);
void SimulateData(void);
int GetMem(char ipop[]);
void FreeMem(void);
int ResetSpeciesGeneTree(int locus);
int Coalescence1Pop(int ispecies);
int CoalescentMigration(void);
int ReadSeqData(char *seqfile, char *locusratef, char *heredityf, FILE*fout, char cleandata, char ipop[]);
int ReadTraitData(char *traitfile, FILE*fout);
int StandardizeTraitData(FILE*fout);
void printTraitData(FILE *fileout);
double lnpGB_ThetaTau(int locus);
double lnpData(double lnpDi[]);
double lnpD_locus(int locus);
int UseLocus(int locus, int copytreeconP, int useData, int setSeqName);
int AcceptLocus(int locus, int copyconP);
int GetInitials(void);
int collectx(FILE* fout, double x[]);
int MCMC(FILE* fout);
double UpdateGB_InternalNode(double* lnL, double finetune);
double UpdateGB_SPR(double* lnL, double finetune);
int NodeSlider(double eps);
double UpdateTheta(double finetune, double space[]);
double UpdateTau(double *lnL, double finetune, double space[]);
double mixing(double* lnL, double finetune, double space[]);
double UpdateSplitSpecies(double *lnL, double space[], double PrSplit);
double UpdateJoinSpecies(double *lnL, double space[], double PrSplit);
double UpdateLocusrateHeredity(double* lnL, double finetune);
double UpdateSequenceErrors(double* lnL, double finetune, double space[]);
void GraftNode(int source, int target, double age, int ipop);
int MatchGTree(void);
char *printSpItree(int itree);
void printGtree(int printBlength);
void checkGtree(void);
int getab_beta(void);

extern int noisy, IncludeNodeLabel;
char timestr[32];
double OLDAGE=999;
int debug=0, diagnosis=0, testlnL=0, NPMat=0, LASTROUND = 1;

#define REALSEQUENCE
#define NODESTRUCTURE
#define BPP
#include "treesub.c"


int main (int argc, char*argv[])
{
#ifndef SIMULATION
   char ctlf[128]="bpp.ctl";
#else
   char ctlf[128]="MCcoal.ctl";
#endif
   char VStr[64]="Version 2.1, May 2011, with trait modification, July 2011\n";
   FILE *fout;
   int i;

   noisy=0;
   printf("bpp %s\n", VStr);
   starttimer();
   if(argc>1) strcpy(ctlf, argv[1]);
   com.cleandata = 0;
   com.clock = 1; com.ncode = 4; com.model = 0;
   for(i=0; i<4; i++) com.pi[i] = 0.25;

#ifdef SIMULATION
   com.simulation = 1;
   GetOptionsSimulation(ctlf);
   SimulateData(); 
#else
   com.simulation=0;
   com.ngene = -1;  /* not used */
   data.ngene = 1;
   GetOptions(ctlf); 
   fout = gfopen(com.outf,"w");
   fprintf(fout, "bp&p (%s) %s\n", VerStr, com.seqf);
   /* The size of ipop[] is total#sequences*sizeof(char) */
   ReadSeqData(com.seqf, com.locusratef, com.heredityf, fout, com.cleandata, (char*)com.space);
   SetMapAmbiguity(); // the map from the ambiguity characters to resolved characters
   ReadTraitData(com.traitf, fout);
   StandardizeTraitData(fout);
   GetMem((char*)com.space);
   MCMC(fout);
   FreeMem();
#endif
   exit(0);
}


#ifdef SIMULATION


int GetOptionsSimulation (char *ctlf)
{
   int nopt=9, lline=4096, iopt, i, j, is, ierror;
   char line[4096],*pline, opt[32], *comment="*#";
   char *optstr[] = {"seed", "seqfile", "treefile", "Imapfile", "species&tree", 
                     "migration", "sequenceerror", "loci&length", "locusrate"};
   char name[LSPNAME];
   double t=1;
   FILE  *fctl=gfopen (ctlf, "r");

   strcpy(com.Imapf, "Imap.txt");
   if (fctl) {
      if (noisy) printf ("\nReading options from %s..\n", ctlf);
      for (;;) {
         if(fgets(line, lline, fctl) == NULL) 
            break;
         if(line[0]=='/' && line[1]=='/') 
            break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (strchr(comment,line[i])) break;
         if (t==0) continue;
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "="))==NULL) 
            continue;

         for (iopt=0; iopt<nopt; iopt++) {
            if (strncmp(opt, optstr[iopt], 8)==0)  {
               if (noisy>=9)
                  printf ("\n%3d %15s | %-20s %6.2f", iopt+1,optstr[iopt],opt,t);
               switch (iopt) {
                  case ( 0): SetSeed((int)t, 1);                 break;
                  case ( 1): sscanf(pline+1, "%s", com.seqf);    break;
                  case ( 2): sscanf(pline+1, "%s", com.outf);    break;
                  case ( 3): sscanf(pline+1, "%s", com.Imapf); break;
                  case ( 4): 
                     if((sptree.nspecies=com.ns=(int)t)>NSPECIES) error2("raise NSPECIES?");
                     if(sptree.nspecies>8*sizeof(int)) error2("NSPECIES larger than size of int.");
                     ReadSpeciesTree(fctl, pline+1);
                     break;
                  case ( 5):
                     ReadMigrationMatrix(fctl, pline);
                     break;
                  
                  case (6):               /* sequencing errors */
                     data.nseqerr = 0;  
                     if(sscanf(pline+1, "%d", &data.nseqerr) != 1)
                        error2("error in the sequenceerror line of the control file.");
                     for(ierror=0; ierror<data.nseqerr; ierror++) {
                        fscanf(fctl, "%s", name);
                        for(is=0; is<sptree.nspecies; is++)
                           if( strcmp(name, sptree.nodes[is].name) == 0) break;
                        if(is==sptree.nspecies) error2("expecting a species name");
                        data.iseqerr[is] = 1;
                        for(i=0; i<16; i++)
                           fscanf(fctl, "%lf", &data.e_seqerr[is][i]);
                        for(i=0; i<4; i++) 
                           abyx (1/sum(data.e_seqerr[is]+i*4, 4), data.e_seqerr[is]+i*4, 4);               
                        printf("\nsequence errors in %s:", sptree.nodes[is].name);
                        matout (F0, data.e_seqerr[is], 4, 4);
                        for(i=0; i<4; i++) for(j=1; j<4; j++) 
                           data.e_seqerr[is][i*4+j] += data.e_seqerr[is][i*4+j-1];
                        matout (F0, data.e_seqerr[is], 4, 4);
                     }
                     break;

                  case ( 7): sscanf(pline+1, "%d%d", &com.ngene, &com.ls);  break;
                  case ( 8): sscanf(pline+1, "%lf", &com.a_locusrate);  break;  /* gamma rates for loci */
               }
               break;
            }
         }
         if (iopt==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");
   return(0);
}

#else

int GetOptions (char *ctlf)
{
   int nopt=24, lline=4096, locfields[1000], iopt, i, is, ierror;
   char line[4096],*pline, opt[32], *comment="*#", *seqerrstr="0EF";
   char *optstr[] = {"seed", "seqfile", "Imapfile", "outfile", "mcmcfile", "speciesdelimitation", 
                     "uniformrootedtrees", "species&tree", "usedata", "nloci", "cleandata", 
                     "thetaprior", "tauprior", "heredity", "locusrate", "sequenceerror", 
                     "finetune", "print", "burnin", "sampfreq", "nsample",
		     "traitfile", "ntraits","nindT"};
   char name[LSPNAME];
   double t=1, *eps=mcmc.finetune;
   FILE  *fctl=gfopen (ctlf, "r");

   sptree.nLHistories = 1;  /* uniform prior prob for rooted trees */
   for(is=0; is<NSPECIES; is++) 
      for(i=0; i<16; i++) 
         data.a_seqerr[is][i] = data.e_seqerr[is][i] = 0;
   traitdata.nind = 0;      /* initializing trait data */
   traitdata.ntrait = 0;

   if (fctl) {
      if (noisy) printf ("\nReading options from %s..\n", ctlf);
      for (;;) {
         if(fgets(line, lline, fctl) == NULL) 
            break;
         if(line[0]=='/' && line[1]=='/') // the "//" stops the reading of the file.
            break;
         for (i=0,t=0,pline=line; i<lline&&line[i]; i++)
            if (isalnum(line[i]))  { t=1; break; }
            else if (strchr(comment,line[i])) break;
         if (t==0) continue;              // lines starting with a comment symbol are ignored
         sscanf (line, "%s%*s%lf", opt, &t);
         if ((pline=strstr(line, "="))==NULL) error2 ("option file.\nExpecting '=' ");

         for (iopt=0; iopt<nopt; iopt++) {
            if (strncmp(opt, optstr[iopt], 12)==0)  {
               if (noisy>=9)
                  printf ("\n%3d %15s | %-20s %6.2f", iopt+1,optstr[iopt],opt,t);
               switch (iopt) {
                  case ( 0): SetSeed((int)t, 1);                 break;
                  case ( 1): sscanf(pline+1, "%s", com.seqf);    break;
                  case ( 2): sscanf(pline+1, "%s", com.Imapf);   break;
                  case ( 3): sscanf(pline+1, "%s", com.outf);    break;
                  case ( 4): sscanf(pline+1, "%s", com.mcmcf);   break;
                  case ( 5): 
                     sscanf(pline+1, "%d%d%lf%lf%d", &sptree.speciesdelimitation, &mcmc.RJalgorithm, &mcmc.RJfinetune[0], &mcmc.RJfinetune[1], &diagnosis);
                     if(sptree.speciesdelimitation) {
                        printf("rj algorithm %d: new theta from ", mcmc.RJalgorithm);
                        if(mcmc.RJalgorithm) 
                           printf("G(a=%.2f, m=%.2f)\n", mcmc.RJfinetune[0], mcmc.RJfinetune[1]);
                        else 
                           printf("sliding window with c = %.2f\n", mcmc.RJfinetune[0]);
                        if(mcmc.RJfinetune[0]<=0 || (mcmc.RJalgorithm && mcmc.RJfinetune[1]<=0))
                           error2("RJfinetune <= 0?");
                     }
                     break;
                  case ( 6): sscanf(pline+1, "%d", &sptree.nLHistories);   break;
                     break;
                  case ( 7): 
                     if((sptree.nspecies=com.ns=(int)t)>NSPECIES)
                        error2("raise NSPECIES");
                     if(sptree.nspecies>8*sizeof(int)) error2("NSPECIES larger than size of int.");
                     ReadSpeciesTree(fctl, pline+1);
                     break;
                  case ( 8): mcmc.usedata=(int)t;   break;
                  case ( 9): data.ngene=(int)t;     break;
                  case (10): com.cleandata=(int)t;  break;
                  case (11): sscanf(pline+1, "%lf%lf", &data.theta_prior[0], &data.theta_prior[1]);
                     break;
                  case (12): sscanf(pline+1, "%lf%lf%lf", &data.tau_prior[0], &data.tau_prior[1], &data.tau_prior[2]);
                     break;
                  case (13):
                     sscanf(pline+1, "%d%lf%lf", &data.est_heredity, &data.a_heredity, &data.b_heredity);
                     if(data.est_heredity==1) {
                        printf("\n\nheredity multiplier ~ gamma(%.2f, %.2f)\n", data.a_heredity, data.b_heredity);
                        if(data.a_heredity<.001 || data.b_heredity<.001) error2("poor prior?");
                     }
                     else if(data.est_heredity==2) {
                        splitline(pline+1, locfields);
                        sscanf(pline+1+locfields[1], "%s", com.heredityf);
                     }
                     break;
                  case (14):               /* locus rate */
                     sscanf(pline+1, "%d%lf", &data.est_locusrate, &data.a_locusrate);
                     if(data.est_locusrate==1) {
                        printf("\nRates vary among loci, with ri/L ~ Dirichlet(%.2f)\n", data.a_locusrate);
                        if(data.a_locusrate<.001) error2("alpha very small?");
                     }
                     else if(data.est_locusrate==2) {
                        splitline(pline+1, locfields);
                        sscanf(pline+1+locfields[1], "%s", com.locusratef); /* check this? */
                     }
                     break;

                  case (15):               /* sequencing errors */
                     splitline(pline+1, locfields);
                     data.nseqerr = 0;  
                     if(sscanf(pline+1, "%d", &data.nseqerr) != 1)
                        error2("error in the sequenceerror line of the control file.");
                     for(ierror=0; ierror<data.nseqerr; ierror++) {
                        fscanf(fctl, "%s", name);
                        for(is=0; is<sptree.nspecies; is++)
                           if( strcmp(name, sptree.nodes[is].name) == 0) break;
                        if(is==sptree.nspecies) error2("expecting a species name");
                        data.iseqerr[is] = 1;
                        for(i=0; i<16; i++)
                           fscanf(fctl, "%lf", &data.a_seqerr[is][i]);
                        printf("\nalpha matrix for sequence errors in %s:", sptree.nodes[is].name);
                        matout (F0, data.a_seqerr[is], 4, 4);
                     }
                     break;

                  case (16):
                     sscanf(pline+1,"%d:%lf%lf%lf%lf%lf%lf%lf", &mcmc.resetFinetune, eps,eps+1,eps+2,eps+3,eps+4,eps+5,eps+6);
                     break;
                  case (17): mcmc.print=(int)t;     break;
                  case (18): mcmc.burnin=(int)t;    break;
                  case (19): mcmc.sampfreq=(int)t;  break;
                  case (20): mcmc.nsample=(int)t;   break;
                  case (21): sscanf(pline+1, "%s", com.traitf);   break;
                  case (22): traitdata.ntrait=(int)t; break;
                  case (23): traitdata.nind=(int)t; break;
               }
               break;
            }
         }
         if (iopt==nopt)
            { printf ("\noption %s in %s\n", opt, ctlf);  exit (-1); }
      }
      fclose(fctl);
   }
   else
      if (noisy) error2("\nno ctl file..");

   return(0);
}


#endif


int ReadSpeciesTree (FILE* fctl, char *currline)
{
/* It reads the species tree and initializes sptree.
   Note that (sptree.nseqsp[j] > 1) adds a population.
*/
   char line[16000]={'\0'};
   int i,j, lline=16000, locfields[1000], maxns=0;

   nodes = nodes_t;  /* species tree is read into the temp space nodes_t */
   for(i=0; i<NS; i++) {  /* for both species & gene trees */
      if (com.spname[i]) free(com.spname[i]);
      com.spname[i] = (char*)malloc((LSPNAME+1)*sizeof(char));
      /* for(j=0; j<LSPNAME; j++) com.spname[i][j]=0; */
   }
   splitline (currline, locfields);
   for(i=0; i<sptree.nspecies; i++) {
      sscanf(currline+locfields[i+1], "%s", com.spname[i]);
      strcpy(sptree.nodes[i].name, com.spname[i]);
   }
   for(i=0,maxns=0; i<sptree.nspecies; i++) {
      fscanf(fctl, "%d", &sptree.nseqsp[i]);
      maxns += sptree.nseqsp[i];
      if(sptree.nseqsp[i]<1)  /* sptree.migration is not read yet. */
         puts("Shall we have at least 1 seq from each species?");
   }
   if(maxns>NS) error2("raise NS in source file");
   if(maxns<2) error2("<2 seqs?  Too simple to do anything about.");
   if(traitdata.nind == 0) // it means the control file did not provide a value for this
     traitdata.nind = maxns;// so using the genetic data to fill this in.
   printf("%d species: ", sptree.nspecies);
   for(i=0; i<sptree.nspecies; i++) 
      printf(" %s (%d)",com.spname[i], sptree.nseqsp[i]);
   putchar('\n');

#ifndef SIMULATION
      data.maxns = maxns;
#endif
   fgets(line, lline, fctl);
   if(sptree.nspecies==1) {
      sptree.root = sptree.nodes[0].nson = 0;  sptree.nnode = 1; sptree.nbranch = 0; 
      sptree.nodes[0].age = 0;
      sptree.npop = 1;
      sptree.pops[0] = 0;
      sptree.pop_pop_table[0][0] = 1;
#ifdef SIMULATION
      fscanf(fctl, "%lf", &sptree.nodes[0].theta);
      printf("theta = %9.4f", sptree.nodes[0].theta);
#endif
   }
   else {
      ReadTreeN (fctl, &i, &j, 0, 1);
      OutTreeN(F0,1,com.simulation);
      FPN(F0); // defined in paml.h: File Put Newline to F0=stdout

      /* copy into sptree */
      sptree.nnode = tree.nnode; 
      sptree.nbranch = tree.nbranch;
      sptree.root=tree.root;
      for(i=0; i<sptree.nnode; i++) {
         sptree.nodes[i].father = nodes[i].father;
         sptree.nodes[i].nson = nodes[i].nson;
         for(j=0;j<sptree.nodes[i].nson;j++) 
            sptree.nodes[i].sons[j] = nodes[i].sons[j];
         sptree.nodes[i].age = nodes[i].age = nodes[i].branch;
         sptree.nodes[i].theta = nodes[i].label;
      }
      for(i=0,sptree.npop=0; i<sptree.nspecies; i++)
         if(sptree.nseqsp[i]>1)
            sptree.pops[sptree.npop++] = i;
      for(; i<2*sptree.nspecies-1; i++)
         sptree.nodes[i].name[0] = '\0';

      DownSptreeSetPops(sptree.root);

      printf("\npop by pop table showing node numbers in species tree\n\n%20s", " ");
      for(i=0; i<2*sptree.nspecies-1; i++)
         printf(" %2d", i+1);
      FPN(F0);
      for(i=0; i<2*sptree.nspecies-1; i++) 
         for(j=0; j<2*sptree.nspecies-1; j++) 
            sptree.pop_pop_table[i][j] = 0;
      for(i=0; i<2*sptree.nspecies-1; i++) {
         for(j=i; ; j=sptree.nodes[j].father) {
            if(j>=sptree.nspecies || sptree.nseqsp[j]>1)
               sptree.pop_pop_table[i][j] = 1;
            if(j == sptree.root) break;
         }
      }
      for(i=0; i<2*sptree.nspecies-1; i++,FPN(F0)) {
         printf("species %2d %-8s ", i+1, sptree.nodes[i].name);
         for(j=0; j<2*sptree.nspecies-1; j++) 
            printf(" %2d", sptree.pop_pop_table[i][j]);
         if(i<sptree.nspecies && sptree.nodes[i].age) 
            printf("\n\a  <-- age>0??\n");
#ifdef SIMULATION
         printf("  tau =%7.4f theta =%7.4f  ", sptree.nodes[i].age, sptree.nodes[i].theta);
         if((i>=sptree.nspecies || sptree.pop_pop_table[i][i]) && sptree.nodes[i].theta<=0)
            printf("  this theta must be > 0!!");
#endif
      }
   }

   if(com.simulation) {
      for(i=0,com.ns=0; i<sptree.nspecies; i++)
         com.ns += sptree.nseqsp[i];
      if(com.ns>NS)
         error2("raise NS in source");
   }
   else {
      printf("\n%d theta parameters (populations) in the order:", sptree.npop);
      for(i=0; i<sptree.npop; i++)
         printf(" %2d (%s)", sptree.pops[i]+1, sptree.nodes[sptree.pops[i]].name);
      if(sptree.nspecies>1) {
         printf("\n%d species divergence times in the order:", sptree.nspecies-1);
         for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++)
            printf(" %2d (%s)", i+1, sptree.nodes[i].name);
      }
      printf("\n");
   }
   return(0);
}


int ReadMigrationMatrix(FILE *fctl, char *pline)
{
/* Rates for impossible ancestor-descendent migrations are set to 0.
   This does not check node ages to make sure that the populations live 
   in the same epoch.
*/
   int nsp=sptree.nspecies*2-1, i,j, needtheta, status=0;
   double Agei, Agej, Adadi, Adadj;
   char spname[40];

   sscanf(pline+1, "%d", &sptree.migration );  
   if(sptree.migration == 0) 
      return(0);

   if(sptree.migration != nsp) 
      error2("migration = #species*2 - 1?");
   for(i=0; i<nsp; i++) {
      fscanf(fctl, "%s", spname);
      if( strcmp(sptree.nodes[i].name, spname) )  error2("Migr: spname mismatch");
   }
   for(i=0; i<nsp; i++) {
      fscanf(fctl, "%s", spname);
      if( strcmp(sptree.nodes[i].name, spname) )
         error2("Migr: spname mismatch");
      for(j=0; j<nsp; j++)
         fscanf(fctl, "%lf", &sptree.M[i][j]);
   }

   /* set rates for impossible migrations to 0 or -1 */
   for(i=0; i<nsp; i++) {
      Agei = sptree.nodes[i].age;
      Adadi = sptree.nodes[sptree.nodes[i].father].age;
      for(j=i; j<nsp; j++) {
         Agej = sptree.nodes[j].age;
         Adadj = sptree.nodes[sptree.nodes[j].father].age;
         if(i==j || sptree.pop_pop_table[i][j] || sptree.pop_pop_table[j][i])
            sptree.M[i][j] = sptree.M[j][i] = 0;
         else if(com.simulation && (Agei>Adadj || Agej>Adadi))
            sptree.M[i][j] = sptree.M[j][i] = -1;
      }
   }

   printf("\nmigration matrix\n%10s", "");
   for(j=0; j<nsp; j++)
      printf(" %9s", sptree.nodes[j].name);
   for(i=0,FPN(F0); i<nsp; i++,FPN(F0)) {
      printf("%-10s", sptree.nodes[i].name);
      for(j=0; j<nsp; j++)
         printf(" %9.4f", sptree.M[i][j]);
   }

   for(i=0; i<nsp; i++)
      for(j=0; j<nsp; j++)
         if(sptree.M[i][j] < 0)  sptree.M[i][j] = 0;
   for(i=0; i<sptree.nspecies; i++) {
      for(j=0,needtheta=0; j<nsp; j++)
         if(sptree.M[j][i]) needtheta = 1;
      if(needtheta && sptree.nodes[i].theta<=0) {
         printf("theta for %s should be > 0\n", sptree.nodes[i].name);
         status=-1;
      }
   }
   if(status) error2("fix the control file first.");
   return(0);
}


int DownSptreeSetPops (int inode)
{
/* This traverse down the species tree to see which nodes are visited,
   and initializes npop, pops[], and sptree.nodes[].name.
   Populations where coalescent is possible (theta's) are collected into 
   sptree.pops[] in the post-order tree traversal.
*/
   int k,ison;

   if(inode < sptree.nspecies) error2("should not be here?");

   for (k=0; k<sptree.nodes[inode].nson; k++) {
      ison = sptree.nodes[inode].sons[k];
      if(sptree.nodes[ison].nson)  
         DownSptreeSetPops(ison);
      if(strlen(sptree.nodes[inode].name) + strlen(sptree.nodes[ison].name) > 2*LSPNAME-1)
         error2("we are in trouble.  Increase LSPNAME?");
      strcat(sptree.nodes[inode].name, sptree.nodes[ison].name);
   }
   sptree.pops[sptree.npop++] = inode;
   return(0);
}

/* used by ResetSpeciesGeneTree() and Coalescence1Pop(). Those globals appear
   necessary as Coalescence1Pop is recursive.
*/
static int cNodeNumber, nin_G[NSPECIES*2-1], ins_G[NSPECIES*2-1][NS];

void GetRandomGtree(int locus)
{
/* This generates a random gene tree in nodes[].  This does not initialize ipop.
   This is used by the simulation program, and also by the inference program to 
   generate the starting gene tree with coalescent times.
*/
   int i;

   ResetSpeciesGeneTree(locus);
   cNodeNumber = com.ns;

   if(0 && !sptree.migration) {
      Coalescence1Pop(sptree.root);
      /* NodeToBranch(); */
      tree.root = cNodeNumber-1;  /* cNodeNumber = com.ns*2 - 1.  root is last node */
      nodes[tree.root].branch = 0;
      nodes[tree.root].father = -1;
      tree.nnode = com.ns*2-1;
      for(i=0; i<tree.nnode; i++) 
         if(i != tree.root) 
            nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
   }
   else
      CoalescentMigration();
}

int ResetSpeciesGeneTree (int locus)
{
/* This is called by GetRandomGtree() to reset the species tree before calling Coalescence1Pop().  
   It initializes nin_G & ins_G[].
   Also resets cNodeNumber for constructing the gene tree.
   gnodes[][].ipop for tips in the gene tree have been set before this routine already.
*/
   int IS='_', is, j, tip;

   for(is=0; is<2*sptree.nspecies-1; is++) nin_G[is] = 0;

#ifdef SIMULATION
  for(is=0; is<sptree.nspecies; is++) 
      nin_G[is] = sptree.nseqsp[is];
   for(is=0,tip=0; is<sptree.nspecies; is++) {
      for(j=0; j<sptree.nseqsp[is]; j++,tip++)
         ins_G[is][j] = tip;
   }
#else
   com.ns = data.ns[locus];
   nodes = gnodes[locus];
   for(is=0; is<sptree.nspecies; is++)
      sptree.nseqsp[is] = data.nseqsp[locus][is];
   for(j=0; j<com.ns; j++) {
      is = gnodes[locus][j].ipop;
      ins_G[is][nin_G[is]++] = j;
      if(sptree.nseqsp[is]>1) 
         sprintf(com.spname[j], "%s%d", sptree.nodes[is].name, nin_G[is]);
      else
         sprintf(com.spname[j], "%s", sptree.nodes[is].name);
   }
#endif

   for(j=0; j<com.ns; j++) ClearNode(j);
   return(0);
}


int Coalescence1Pop (int ispecies)
{
/* This simulates the coalescent process in population or species ispecies.
   It generates the random genealogy tree (possibly consisting of several 
   disconnected subtrees) with waiting times.
   t: waiting time; T: node age
   This is used by GetRandomGtree().  
   nin_G[] and ins_G[] are set in ResetSpeciesGeneTree().
*/
   int j, k,k1,k2, father = sptree.nodes[ispecies].father;
   double t, T;

   for (k=0; k<sptree.nodes[ispecies].nson; k++)
      Coalescence1Pop(sptree.nodes[ispecies].sons[k]);
   if(sptree.nodes[ispecies].theta > 0) {
      T = sptree.nodes[ispecies].age;
      if(nin_G[ispecies]>1 && sptree.nodes[ispecies].theta <= 0) {
         printf("theta for pop %s is %.6f", sptree.nodes[ispecies].name, sptree.nodes[ispecies].theta);
         puts(" theta <= 0");
      }
      for (; nin_G[ispecies]>1; nin_G[ispecies]--,cNodeNumber++) {
         j = nin_G[ispecies];  /* # of lineages */
         t = rndexp(sptree.nodes[ispecies].theta/(j*(j-1.)));
         T += t;
         if(ispecies!=sptree.root && T>sptree.nodes[father].age) break;

         /* k1 & k2 are lineages to coalesce */
         k  = (int)(j*rndu());  
         k1 = ins_G[ispecies][k]; ins_G[ispecies][k] = ins_G[ispecies][j-1];
         k  = (int)((j-1)*rndu());
         k2 = ins_G[ispecies][k]; ins_G[ispecies][k] = cNodeNumber;

         nodes[cNodeNumber].nson = 2;
         nodes[cNodeNumber].ipop = ispecies;
         nodes[cNodeNumber].sons[0] = k1;
         nodes[cNodeNumber].sons[1] = k2;
         nodes[k1].father = nodes[k2].father = cNodeNumber;
         nodes[cNodeNumber].age = T;
      }
   }
   /* outgoing lineages added to father's list */
   if(ispecies==sptree.root && nin_G[ispecies]>1) 
      error2("nin_Gtree > 1 at root");
   if(ispecies!=sptree.root) {
      for(k=0; k<nin_G[ispecies]; k++) 
         ins_G[father][nin_G[father]++] = ins_G[ispecies][k];
   }
   return (0);
}

double mtMRCA, mM[(NSPECIES*2-1)*(NSPECIES*2-1)];


int CoalescentMigration (void)
{
/* This simulates a gene tree with both coalescent and migration events.
   nspecies epochs: the coalescent and migration rates are used to generate 
   exponential waiting times.
   The routine works also when there is no migration.
   The algorithm keeps track of # of populations (npop), updates the list of pops in ipop[],
   the number of individuals in each pop nin_G[], and the individuals in each pop ins_G[][].  
   Some pops may be empty, with nin_G[] = 0.  Ci[] holds coalescent rate in pop i and Mi[]
   holds migration rate to pop i, while C & M are total coalescent and migration rates.  
   ipop[] keeps the list of pops during the current epoch, and is updated when we move 
   to the next epoch.
   ipopE holds the order of the epochs.
   n is the number of lineages ancestral to the sample (current sample size).  
   Tmax marks the end of the current epoch.
*/
   int n=com.ns, is, i,j,k, k1,k2, ipopE[NSPECIES-1]; /* ipop for each epoch node */
   int npop=sptree.nspecies, ipop[NSPECIES], *sons;
   double r, T, Tmax, C,Ci[NSPECIES], M=0,Mi[NSPECIES];
   double y, ages[NSPECIES-1], tmp1[NSPECIES];

   /* sort node ages in species tree to work out each epoch. */
   for(i=0; i<sptree.nspecies-1; i++)
      tmp1[i] = sptree.nodes[sptree.nspecies+i].age;
   /* The sorting is in increasing order, with ties broken in the original order.
      This way, the algorithm runs when nodes on the species tree are collapsed. */
   indexing (tmp1, sptree.nspecies-1, ipopE, 1, (int*)ages);  

   if(debug==9) {
      printf("\nages in sptree: ");
      for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++)
         printf(" %2d: %7.5f ", i, sptree.nodes[i].age);
      printf("\nafter ordering: ");
      for(i=sptree.nspecies-1-1; i>=0; i--) {
         is = sptree.nspecies + ipopE[i];
         printf(" %2d: %7.5f ", is, sptree.nodes[is].age);
      }
      printf("\n\n");
   }

   /* Initially the tips on the species tree are in the list ipop[] */
   for(i=0; i<npop; i++)  /* populations */
      ipop[i] = i;

   for(T=0; ; npop--) {   /*  # of epochs */
      if(npop==1)  
         Tmax = -1;
      else { /* is: the species in the species tree at the end of this epoch */
         is = sptree.nspecies + ipopE[npop-1-1];
         Tmax = sptree.nodes[is].age;
      }
      for ( ; ; ) {
         if(Tmax==0) break;
         /* calculate poisson rates: Ci[i] is coalescent rate in pop i */
         for(i=0,C=0; i<npop; i++) {
            Ci[i] = 0;
            if(nin_G[i] >= 2) {
               if(sptree.nodes[ipop[i]].theta <= 0)
                  printf("theta_%s = %.6f <= 0!\n", sptree.nodes[ipop[i]].name, sptree.nodes[ipop[i]].theta);
               C += Ci[i] = nin_G[i]*(nin_G[i]-1)/2 * 2.0/sptree.nodes[ipop[i]].theta;
            }
         }
         if(sptree.migration) 
            /* Mi[i] is the migration rate to pop i */
            for(i=0,M=0; i<npop; i++) {
               for(j=0,Mi[i]=0; j<npop; j++) 
                  Mi[i] += nin_G[i] * sptree.M[ipop[j]][ipop[i]]/sptree.nodes[ipop[i]].theta * 4;
               M += Mi[i];
            }

         if(debug==9) {
            printf("S%d (%5.3f) %d (", is, Tmax, npop);
            for(i=0; i<npop; i++) printf(" %d", nin_G[i]);
            printf(" ) rate CM %6.1f %6.1f ", C, M);
         }
         if(C+M < 1e-300) {             /* move to next epoch */
            if(debug==9) FPN(F0);
            break;
         }
         T += rndexp(1/(C+M));
         if(debug==9) printf(" T %6.3f ", T);
         if(T > Tmax && Tmax != -1) { /* move to next epoch */
            if(debug==9) FPN(F0);
            break;
         }
         r = rndu()*(C+M);
         if(r < C) {  /* coalescent in pop i of lineages k1 & k2 */
            /* r ~ U(0, C) */
            for(i=0,y=0; i<npop-1; i++) 
               if(r < (y += Ci[i]))
                  break;

            /* k1 & k2 (with k1 < k2) are the two lineages to coalesce */
            k1  = (int)(nin_G[i]*rndu());
            k2  = (int)((nin_G[i]-1)*rndu());
            if(k2==k1) k2++;
            if(k1>k2) {   /*  k1 < k2 */
               j=k1; k1=k2; k2=j; 
            }

            nodes[cNodeNumber].nson = 2;
            nodes[cNodeNumber].ipop = ipop[i];
            nodes[cNodeNumber].sons[0] = ins_G[i][k1];
            nodes[cNodeNumber].sons[1] = ins_G[i][k2];
            nodes[ins_G[i][k1]].father = nodes[ins_G[i][k2]].father = cNodeNumber;
            nodes[cNodeNumber].age = T;

            /* In pop i, replace k1 by new node, remove k2 */
            ins_G[i][k1] = cNodeNumber ++;
            -- nin_G[i];
            if(k2 != nin_G[i])
               ins_G[i][k2] = ins_G[i][nin_G[i]];

            if(debug==9) 
               printf("C: %s   ", sptree.nodes[ipop[i]].name);

            if(--n == 1) 
               break;      /* last coalescent */
         }
         else {            /* migration of lineage k from pop j into pop i */
            r -= C;        /* 0 < r < M */
            for(i=0,y=0; i<npop; i++) 
               if(r < (y += Mi[i])) break;
            if(i==npop) error2("u01 = 1!");
            y -= Mi[i];

            for(j=0; j<npop-1; j++)      /* duplicated calculation! */
               if(r < (y += nin_G[i]*sptree.M[ipop[j]][ipop[i]]/sptree.nodes[ipop[i]].theta * 4))
                  break;

            mM[ipop[j] * (sptree.nspecies*2-1) + ipop[i]] ++;

            k = (int)(nin_G[i]*rndu());      /* k is migrant from pop j to i */
            /* shift up lineages in pop i */
            ins_G[j][nin_G[j] ++] = ins_G[i][k];
            if(k != --nin_G[i])
               ins_G[i][k] = ins_G[i][nin_G[i]];

            if(debug==9)
               printf("M: %s < %s ", sptree.nodes[ipop[i]].name, sptree.nodes[ipop[j]].name);
         }
         if(debug==9) {
            for(i=0; i<npop; i++) {
               printf(" %-s:", sptree.nodes[ipop[i]].name);
               for(j=0; j<nin_G[i]; j++) 
                  printf("%d ", ins_G[i][j]);
            }
            FPN(F0);
         }
      }  /* forever loop inside for(epoch) */
      T = Tmax;
      /* To move to next epoch, update ipop[] and merge lineages from sons[0] & 
         sons[1] int pop is. 
         Replace pop k1 by is, move k2 into k1, delete k2.
      */
      if(npop==1 || n == 1) break;
      sons = sptree.nodes[is].sons;
      for(i=0,k1=k2=-1; i<npop; i++) {
         if(ipop[i] == sons[0])       k1 = i;
         else if(ipop[i] == sons[1])  k2 = i;
      }

      if(k1>k2) {
         i=k1; k1=k2; k2=i;
      }
      ipop[k1] = is;  
      for(i=0; i<nin_G[k2]; i++) 
         ins_G[k1][nin_G[k1] ++] = ins_G[k2][i]; 
      if(k2 != npop-1) {
         ipop[k2] = ipop[npop-1];  
         if((nin_G[k2] = nin_G[npop-1]) > 0)
            memmove(ins_G[k2], ins_G[npop-1], nin_G[k2]*sizeof(int));
      }
   }

   tree.root = cNodeNumber-1;  /* cNodeNumber = com.ns*2 - 1.  root is last node */
   if(cNodeNumber != com.ns*2-1) error2("cNodeNumber incorrect");
   nodes[tree.root].branch = 0;
   nodes[tree.root].father = -1;
   tree.nnode = com.ns*2-1;
   for(i=0; i<tree.nnode; i++) 
      if(i != tree.root) 
         nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;

   if(debug==9) {
      FPN(F0); OutTreeN(F0,1,1); FPN(F0); FPN(F0);
   }

   mtMRCA += nodes[tree.root].age;

   return 0;
}


void printGtree (int printBlength)
{
   int i,j, ipop, ipopTrue;
   double t, tb[2];

   for(i=0; i<tree.nnode; i++) 
      if(i!=tree.root) {
         nodes[i].branch = nodes[nodes[i].father].age-nodes[i].age;
         if(nodes[i].branch<0) 
            printf("blength = %9.6f\a\n", nodes[i].branch);
      }
   printf("\nns = %d  nnode = %d", com.ns, tree.nnode);
   printf("\n%7s%7s%11s (%s) %7s%7s","father","node","time","ipop","nson:","sons");
   for(i=0; i<tree.nnode; i++) {
      t = nodes[i].age;
      ipop = nodes[i].ipop;
      tb[0] = sptree.nodes[ipop].age;
      tb[1] = OLDAGE;
      if(ipop != sptree.root) 
         tb[1] = sptree.nodes[sptree.nodes[ipop].father].age;

      ipopTrue = (t>=tb[0] && t<=tb[1]);
      printf ("\n%7d%7d %11.6f (%2d %s %11.6f):%7d  ",
         nodes[i].father, i, t, ipop, (ipopTrue ? "OK" : "??"), 
         sptree.nodes[ipop].age, nodes[i].nson);
      for(j=0; j<nodes[i].nson; j++) printf (" %2d", nodes[i].sons[j]);
   }
   FPN(F0); OutTreeN(F0,0,0); FPN(F0); OutTreeN(F0,1,0); 
   if(printBlength) {
      FPN(F0); OutTreeN(F0,1,1); FPN(F0); 
   }
}


int MatchGTree (void)
{
/* This tests the gene tree topology
*/
   return 
     (nodes[0].father==nodes[1].father
   && nodes[3].father==nodes[4].father 
   && nodes[nodes[3].father].father==tree.root
   && nodes[nodes[0].father].age<sptree.nodes[sptree.root].age
   && nodes[nodes[2].father].age<sptree.nodes[sptree.root].age
   && nodes[nodes[3].father].age<sptree.nodes[sptree.root].age
   );  /* P(gtree) = 0.237913 */
}



#ifdef SIMULATION

static int n123marks[][3]= {{1,2,3},{1,2,3},{2,3,1},{3,1,2}};

void p0124Fromb0b1 (int itree, double p[5], double b[2])
{
/* This calculates p0,p1,p2,p3,p4 for the 5 site patterns for 3 species, 
   given branch lengths b0 and b1.  b0 is the gap, and b1 is the distance 
   from the ancestor of 1 and 2 to species 1.
*/
   double e1, e2, e3;

   e1=exp(-4./3*b[1]);  e2=exp(-8./3*(b[0]+b[1]));  e3=e1*e2;  e1=e1*e1;
   p[0]      = (1. + 3*e1 +  6*e2 +  6*e3)/16;
   p[n123marks[itree][0]]  = (3. + 9*e1 -  6*e2 -  6*e3)/16;
   p[n123marks[itree][1]]=p[n123marks[itree][2]] 
      = (3. - 3*e1 +  6*e2 -  6*e3)/16;
   p[4]      = (6. - 6*e1 - 12*e2 + 12*e3)/16;
}

void PMismatch3s (void)
{
/* This calculate the tree mismatch probabilities for figure 3 of Yang
   (2002 Genetics).
*/
   int iHCG=sptree.root, iHC=3+4-iHCG; /* node numbers in sptree */
   double theta_HC=sptree.nodes[iHC].theta;
   double t_HC=sptree.nodes[iHC].age, t_HCG=sptree.nodes[iHCG].age;
   /* S: species tree; G: gene tree; E: estimated gene tree 
      G[gtree], E[gtree][etree] are counts of resolved (used) loci.
   */
   double SG, SE, GE, G[4], E[4][4], b[2], p0[5],p[5],nused, over,under; 
   int ii, nii=7, ls0[100]={200, 400, 500, 1000, 2000, 4000, 10000};
   int nr=10000000, ir, i,j, gtree=-1, n[5]={0}, etree,etree2,etree3,every=100;

   printf("\nPr{S-G mismatch} = %f from formula.\n",2./3*exp(-2*(t_HCG-t_HC)/theta_HC));
   puts("Ties in genetree are removed in the following calculation.");
   every=max2(every,nr/1000);
   for(ii=0; ii<nii; ii++) {
      com.ls=ls0[ii];
      printf("\n# sites? (Ctrl-C to break) ");
      scanf("%d", &com.ls);
      printf("%d sites, %dK replicates.\n", com.ls, nr/1000);
      FOR(i,4) G[i]=0; FOR(i,4) FOR(j,4) E[i][j]=0;
      for (ir=0,SG=SE=GE=0,nused=0; ir<nr; ir++) {
         GetRandomGtree(-1);
         if (nodes[2].father==tree.root)      gtree=(nodes[0].branch>t_HCG);/* (HC) */
         else if(nodes[0].father==tree.root)  gtree=2; /* (CG) */
         else if(nodes[1].father==tree.root)  gtree=3; /* (GH) */
         else    error2("binary tree?");
         b[1]=nodes[nodes[tree.root].sons[0]].age;
         b[0]=nodes[tree.root].age-b[1];
         p0124Fromb0b1 (gtree, p0, b);
         p[0]=p0[0]+p0[4]; p[1]=p[0]+p0[1]; p[2]=p[1]+p0[2]; p[3]=p[2]+p0[3]; p[4]=-1;
         /* matout(F0, p, 1, 4); */
         MultiNomial2 (com.ls, 4, p, n, NULL);

         for(i=2,etree=1; i<=3; i++)   if(n[i]>n[etree]) etree=i;
         if(etree==1)      { etree2=2; etree3=3; }
         else if(etree==2) { etree2=3; etree3=1; }
         else              { etree2=1; etree3=2; }

         if(n[etree]!=n[etree2] && n[etree]!=n[etree3]) { /* exclude ties */
            nused++;
            G[gtree]++;
            E[gtree][etree]++;
            if(gtree>=2)  SG++;
            if(etree>=2)  SE++;
            if((gtree<=1 && etree!=1) || (gtree>=2 && etree!=gtree)) GE++;
         }
         if((ir+1)%every==0) {
            printf("%4.1f%% (%2d %2d %2d %2d) SG %.4f SE %.4f GE %.4f (+%.4f -%.4f) tie %.4f\r",
               (ir+1.)/nr*100, n[0],n[1],n[2],n[3], 
               SG/nused,SE/nused, GE/nused, 
               (E[0][2]+E[0][3] + E[1][2]+E[1][3])/nused,  /* T -> F */
               (E[2][1]+E[3][1])/nused,                    /* T -> F */
               (ir+1-nused)/(ir+1.));
         }
      }  /* for(ir) */
      if(G[0]+G[1]+G[2]+G[3]-nused!=0) error2("gtree counts incorrect.");
      SG/=nused; SE/=nused; GE/=nused; 
      FOR(i,4) FOR(j,4) E[i][j] /= G[i];
      FOR(i,4) G[i]/=nused;
      printf("\nfrequencies of gene trees 0123 (given ties are removed): ");
      matout2 (F0, &G[0], 1, 4, 9, 5);
      printf("transition probability matrix (gene tree 0123 -> MLtree (01)23):");
      matout2 (F0, &E[0][0], 4, 4, 9, 5);

      printf("\nThe following three should be equal:\n");
      printf("  (1) SE - SG = %.5f\n", SE-SG);
      printf("  (2) f(T0) * P0 = %.5f\n", G[0]*(E[0][2]+E[0][3]));
      over=G[0]*(E[0][2]+E[0][3])+G[1]*(E[1][2]+E[1][3]);
      under=G[2]*E[2][1]+G[3]*E[3][1];
      printf("  (3) over - under = %.5f - %.5f = %.5f\n", over,under,over-under);
      printf("\nf(2)*P2/2+f(3)*P3/2 = %.5f, which is GE - (over + under) = under.\n",
         G[2]*E[2][3]+G[3]*E[3][2]);
      printf("\nTime used: %s\n", printtime(timestr));
   }
   exit(0);
}

void SimulateData (void)
{
   char *zt[NS], *concatF="concatenated.txt", timestr[32];
   FILE *fseq=NULL, *ftree=NULL, *fImap=NULL;
   double rlocus, rlocusold=1, u;
   int nr=com.ngene, ir, i,j,k, variable_ns=0, nseqsp0[NSPECIES], ispecies, h;
   FILE *fconcat = (com.seqf[0] ? gfopen(concatF,"w") : NULL); 
 
   /* PMismatch3s(); */

   if(com.seqf[0]) {
      printf("\nsequence data go into %s & %s\nmap into %s\n", com.seqf, concatF, com.Imapf);
      fseq = gfopen(com.seqf, "w");
      fImap = gfopen(com.Imapf, "w");
   }
   if(com.outf[0]) {
      printf("trees go into %s.\n", com.outf);
      ftree=gfopen(com.outf,"w");
   }

   if(fconcat) {
      for(i=0; i<com.ns; i++) {
         if((zt[i] = (char*)malloc(com.ls*nr)) == NULL)
            error2("oom zt");
      }
   }
   if(variable_ns) error2("option variable_ns does not work with Imap?");
   for(j=0,k=0; j<sptree.nspecies; j++) {
      for(i=0; i<sptree.nseqsp[j]; i++, k++) {
         if(fImap) fprintf(fImap, "%6d %s\n", k+1, sptree.nodes[j].name);
         sprintf(com.spname[k], "%s%d^%d", sptree.nodes[j].name, i+1, k+1);
         nodes[k].ipop = j;
      }
   }
   if(fImap) fclose(fImap);

   if(fseq) {
      com.z[0] = (char*)malloc((2*com.ns-1)*com.ls*sizeof(char));
      for (i=1; i<2*com.ns-1; i++) com.z[i] = com.z[i-1]+com.ls;
      if(com.alpha) com.rates = (double*)malloc(com.ls*sizeof(double));
      if(com.z[0]==NULL || (com.alpha && com.rates==NULL)) error2("oom");
      com.cleandata = 1;
      com.alpha = 0;  com.ncatG = 5;
      for(i=0; i<4; i++) com.pi[i] = 1./4;
   }
   if(variable_ns) 
      for(j=0; j<sptree.nspecies; j++) 
         nseqsp0[j] = sptree.nseqsp[j];

   for (ir=0; ir<nr; ir++) {
      if(variable_ns) {
         for( ; ; ) { /* make sure at least 2 sequences are in the data */
            for(j=0,com.ns=0; j<sptree.nspecies; j++) 
               com.ns += sptree.nseqsp[j] = (int)(nseqsp0[j]*rndu()*0.95+0.5);
            if(com.ns>1) break;
         }
      }

      if(com.a_locusrate) {
         rlocus = rndgamma(com.a_locusrate);
         printf("rate for locus %2d = %.5f\n", ir+1, rlocus);
         for(j=0; j<sptree.nspecies*2-1; j++) {
            sptree.nodes[j].theta *= rlocus/rlocusold;
            if(j>=sptree.nspecies) sptree.nodes[j].age *= rlocus/rlocusold;
         }
         if(sptree.migration) {
            for(i=0; i<sptree.nspecies*2-1; i++)
               for(j=0; j<sptree.nspecies*2-1; j++) 
                  if(sptree.M[i][j] > 0) 
                     sptree.M[i][j] /= rlocus/rlocusold;
         }
         rlocusold = rlocus;
      }

      GetRandomGtree(-1);

      if(ftree) { 
         OutTreeN(ftree,1,1);  fprintf(ftree, " [TH = %.6f]\n", nodes[tree.root].age); 
      }

      if(fseq) {
         if (com.alpha) {
            Rates4Sites (com.rates, com.alpha, com.ncatG, com.ls, 1, com.space);
            for (j=1; j<com.ls; j++) com.rates[j] += com.rates[j-1];
            abyx (1/com.rates[com.ls-1], com.rates, com.ls);
         }
         for(i=0; i<com.ls; i++) 
            com.z[tree.root][i] = (char)(rndu()*4.);
         EvolveJC (tree.root);
         for(i=0; i<com.ns; i++) {
            if(data.iseqerr[j = nodes[i].ipop]==0) 
               continue;
            for(h=0; h<com.ls; h++) {
               for(k=0,u=rndu(); k<3; k++)
                  if(u < data.e_seqerr[j][com.z[i][h]*4+k])
                     break;
               com.z[i][h] = k;
            }
         }

         printSeqs(fseq, NULL, NULL, 0);

         if(fconcat) {
            for(i=0; i<com.ns; i++)
               memcpy(zt[i]+ir*com.ls, com.z[i], com.ls*sizeof(char));
         }
      }
      if((ir+1)%10000==0)
         printf("\r%10d replicates done... mtMRCA = %.6f", ir+1, mtMRCA/(ir+1));
   }  /* for(ir) */

   if(ftree) fclose(ftree);
   if(fseq) {
      fclose(fseq); 
      free(com.z[0]);
      if(com.alpha) free(com.rates);
      com.ls *= nr;
      for(i=0; i<com.ns; i++)  com.z[i] = zt[i];
      printSeqs(fconcat, NULL, NULL, 0);
      fclose(fconcat);
      for(i=0; i<com.ns; i++)  free(zt[i]);
   }
   if(variable_ns) 
      for(j=0; j<sptree.nspecies; j++) sptree.nseqsp[j] = nseqsp0[j];

   for(i=0; i<(sptree.nspecies*2-1)*(sptree.nspecies*2-1); i++)
      mM[i] /= nr;
   if(sptree.migration) {
      printf("Counts of migrations averaged over replicates\n", mtMRCA/nr);
      matout(F0, mM, (sptree.nspecies*2-1), (sptree.nspecies*2-1));
   }
   exit(0);
}

void EvolveJC (int inode)
{
/* Special version of Evolve that works with JC69-like (poisson) model only.
   This can be used to generate amino acid sequences also.
   For each branch in the tree, this determines the number of mutations 
   and then assigns mutations to sites at random.
   When alpha>0, com.rates[] are the accumulative probabilities.
*/
   int is,j, h, nmut, imut, ison;
   double r;

   if (com.alpha && fabs(com.rates[com.ls-1]-1)>1e-4) 
      { printf ("rates c.d.f.: 1 = %.6f?\n", com.rates[com.ls-1]); exit(-1); }
   for (is=0; is<nodes[inode].nson; is++) {
      ison = nodes[inode].sons[is];
      for(h=0; h<com.ls; h++) com.z[ison][h] = com.z[inode][h];
      nmut = rndpoisson (nodes[ison].branch*com.ls);
      for (imut=0; imut<nmut; imut++) {
         if (com.alpha==0) 
            h = (int)(rndu()*com.ls);
         else 
            for (h=0,r=rndu(); h<com.ls; h++) if (r<com.rates[h]) break;
         j = (int)(rndu()*3);
         if (j >= com.z[ison][h]) j++;
         com.z[ison][h] = (char)j;
      }
      if (nodes[ison].nson) EvolveJC(ison);
   }
}


#else

int ReadSeqData(char *seqfile, char *locusratef, char *heredityf, FILE*fout, char cleandata, char ipop[])
{
/* Read sequences at each locus. This sets data.nseqsp[ig].  The ipop info for 
   tips on gene trees is returned in ipop[], and passed to GetMem(), which allocates 
   gnodes[].

   use com.cleandata=1 to delete ambiguities at all loci.
*/
   FILE *fseq=gfopen(seqfile,"r"), *frates=NULL, *fheredity=NULL, *fImap;
   char line[10000], *pline, iname[24], ID='^';
   int locus, i,j, ind, is, lname=24;
   int maxind=data.maxns*data.ngene, nind=0, maxns;
   double mr;
   char clean0=cleandata;

   printf("\nReading Individual-Species map (Imap) from %s\n", com.Imapf);
   data.Imap = (char*)malloc(maxind*(1+lname)*sizeof(char));
   if(data.Imap == NULL) error2("oom Imap");
   memset(data.Imap, 0, maxind);
   /* initialized to Imap[i] = 0 for the first i's.
      The first maxind bytes of the memory block are used for the actual 
      mapping individual i -> species index.
      The next maxind*lname bytes of the memory block are used for the individual 
      names read from the file, which data.Inames[i] are going to point to. */ 

   if(sptree.nspecies>1) {
      fImap = gfopen(com.Imapf, "r");
      for(i=0; i<maxind; i++) {
         data.Inames[i] = data.Imap + maxind + i*lname;
         if(fscanf(fImap, "%s%s", data.Inames[i], line) != 2) break;
         if(strstr(data.Inames[i], "//")) break;
         for(j=0; j<sptree.nspecies; j++) 
            if(strcmp(line, sptree.nodes[j].name) == 0) {
               data.Imap[i] = j;
               break;
            }
         if(j==sptree.nspecies) {
            printf("\nspecies %s in map file not found in the control file.\n", line);
            exit(-1);
         }
      }
      printf("Individual -> Species map: ");
      nind=i;
      for(i=0; i<nind; i++)
         printf(" %d", data.Imap[i]+1);
      fputc('\n', F0);
   }

   printf("\nReading sequence data..  %d loci\n", data.ngene);
   if(data.ngene>NGENE) error2("raise NGENE?");
   if(data.est_locusrate && data.ngene==1)
      error2("ngene = 1 & locus rates");

   /* allocate space for heredity scalars and locus rates.  
      Some space is wasted if one of the two is used only. 
   */
   if(data.est_heredity || data.est_locusrate) {
      data.heredity = (double*)malloc(2*data.ngene*sizeof(double));
      if(data.heredity == NULL) error2("oom heredity");
      data.locusrate = data.heredity+data.ngene;
      for(locus=0; locus<data.ngene; locus++) 
         data.heredity[locus] = data.locusrate[locus]=1;
   }
   /* Read locus rates from the file locusratef.  The file name is read from 
      the control file but data.ngene was unknown there. */
   if(data.est_locusrate == 2) {
      frates = (FILE*) fopen(locusratef, "r");
      if(frates==NULL) 
         error2("locus rate file does not exist");
      for(locus=0,mr=0; locus<data.ngene; locus++) {
         if(fscanf(frates, "%lf", &data.locusrate[locus]) != 1) {
            printf("\nEOF when reading rates from %s\n", locusratef);
            error2("");
         }
         mr = (mr*locus + data.locusrate[locus])/(locus + 1.0);
      }
      fclose(frates);
      for(locus=0; locus<data.ngene; locus++) 
         data.locusrate[locus] /= mr;
      printf("using fixed rates for loci from file %s\n", locusratef);
      if(data.ngene<=200)
         matout2(F0, data.locusrate, 1, data.ngene, 8, 4);
   }

   /* Read heredity scalars for loci from heredityf.  The file name is read from 
      the control file but data.ngene was unknown there. */
   if(data.est_heredity == 2) {
      if((fheredity=(FILE*)fopen(heredityf, "r")) == NULL) 
         error2("heredity scalar file does not exist");
      for(locus=0; locus<data.ngene; locus++) {
         if(fscanf(fheredity, "%lf", &data.heredity[locus]) != 1) {
            printf("\nEOF when reading rates from %s\n", heredityf);
            error2("");
         }
      }
      fclose(fheredity);
   }

   for(locus=0,maxns=0; locus<data.ngene; ipop+=data.ns[locus++]) {
      fprintf(fout, "\n\n*** Locus %d ***\n", locus+1);
      com.cleandata = clean0;
      ReadSeq(fout, fseq, clean0); // this function updates com.ns to # of sequences read
      data.cleandata[locus] = com.cleandata;
      if(com.ns<=1) error2("one seq not useful in sequence file");
      if(data.nseqerr==0) PatternWeightJC69like (fout);
      if(com.ns>data.maxns) {
         printf("\n%d seqs at locus %d.  More than allowed by the control file.", com.ns, locus+1);
         exit(-1);
      }
      maxns = max2(maxns,com.ns);  // maxns is the actual maximum # sequences over all loci
      data.ns[locus] = com.ns;  
      data.ls[locus] = com.ls;
      data.npatt[locus] = com.npatt;
      data.fpatt[locus] = com.fpatt; 
      com.fpatt = NULL;
      for(i=0; i<com.ns; i++) {
         data.z[locus][i] = com.z[i]; com.z[i] = NULL;
      }

      if(sptree.nspecies == 1) {
         data.nseqsp[locus][0] = com.ns;
      }
      else {
         for(i=0; i<sptree.nspecies; i++) data.nseqsp[locus][i] = 0;
         for(i=0; i<com.ns; i++) {
            pline = strchr(com.spname[i], ID) + 1;
            sscanf(pline, "%s", iname);
            for(j=0; j<nind; j++)
               if(strcmp(iname, data.Inames[j])==0) break;
            if(j==nind)
               printf("individual label %s not recognised.", iname);
            else 
               is = data.Imap[j];  /* sequence i is individual ind and species is. */      
            ipop[i] = is;
            data.nseqsp[locus][is] ++;
            if(noisy>=9)
               printf("seq %2d %-20s is for indiv %d from species %s\r", i+1, com.spname[i], j+1, sptree.nodes[is].name);
         }
         for(i=0,j=0; i<sptree.nspecies; i++) {
            j += data.nseqsp[locus][i];
            if(data.nseqsp[locus][i]>1 && sptree.nseqsp[i]<=1) {
               printf("\nAt locus %d, species %d has %d seqs.", locus+1, i+1, data.nseqsp[locus][i]);
               error2("while control file says no more than one.");
            }
            fprintf(fout, "%s (%d) ", sptree.nodes[i].name, data.nseqsp[locus][i]);
         }
         if(j != data.ns[locus]) error2("ns not correct.");
      }
      printf("locus %2d: %2d sequences (", locus+1, data.ns[locus]);
      for(i=0; i<sptree.nspecies; i++) printf(" %d", data.nseqsp[locus][i]);
      printf(" ) ");

      printf("%5d patterns, %s\r", com.npatt,(com.cleandata?"clean":"messy"));
      if(data.est_heredity==2) 
         printf(" (heredity scalar = %4.2f ", data.heredity[locus]);
      printf((data.ngene>500 ? "\r" : "\n"));
   }

   fclose(fseq);
   for(i=0,com.cleandata=1; i<data.ngene; i++) 
      if(!data.cleandata[i]) com.cleandata=0;

   fflush(fout);
   if((data.maxns = maxns) > NS) error2("raise NS in source file?");

   data.lnpGBi = (double*)malloc(data.ngene*sizeof(double));
   data.lnpDi  = (double*)malloc(data.ngene*sizeof(double));
   if(data.lnpGBi==NULL || data.lnpDi==NULL) error2("oom");

   return(0);
}

void printTraitData(FILE *fileout)
{
  int i,j;
  fprintf(fileout, "Sp. ");
  for(i=0; i<traitdata.ntrait; i++)
    fprintf(fileout, " %8s", traitdata.traitName[i]);
  for (j=0; j<traitdata.nind; j++){
    fprintf(fileout, "\n%-4d", traitdata.indSpeciesMap[j]+1);
    for (i=0; i<traitdata.ntrait; i++)
      if (!traitdata.ismissing[i][j])
	fprintf(fileout, " %8.3g", traitdata.y[i][j]);
      else 
	fprintf(fileout, "       NA");
  }
  fprintf(fileout, "\n");fflush(fileout);
}

int ReadTraitData(char *traitfile, FILE*fout)
{
/* Read trait data, set traitdata.  
   There can be one or more traits, missing data allowed.
   traitdata.nind was previously initialized either in GetOptions or ReadSpeciesTree.
*/
   FILE *ftrait=gfopen(traitfile,"r");
   int lname = 24;  // max length for each trait variable name
   int lline=10000; // max line length
   char line[10000], name[24], *comment="*#[";
   int i, j, jsp;

   printf("\nReading Trait data from %s... ", traitfile);
   fflush(stdout);
   if(traitdata.ntrait>NTRAIT) {
     printf("%d traits are too many for me. ", traitdata.ntrait);
     error2("Raise NTRAIT?");
   }

   for (i=0; i<traitdata.ntrait; i++){
     traitdata.traitName[i] = (char*)malloc( (lname+1)*sizeof(char) );
     if(traitdata.traitName[i] == NULL) error2("Error while allocating memory for trait names.");
   }

   /* reading header with trait variable names */
   for (i=0; i<2;){ /* reading the first 2 names. First should be 'Individual',
		       second should be "Population", "clade", or "species" */
     fscanf(ftrait, "%s", name);
     if (strchr(comment,name[0])){ // ignoring the rest of the line
       fgets(line, lline, ftrait); continue;}
     i++;
   }
   for (i=0; i<traitdata.ntrait;){
     if( fscanf (ftrait, "%s", traitdata.traitName[i]) == 0 ){
       printf("Unable to read the name of %dth trait, in trait data file.\n",i+1);
       exit(-1);
     }
     if (strchr(comment, traitdata.traitName[i][0])){ // ignoring the rest of the line
       fgets(line, lline, ftrait); continue;}
     traitdata.cleantrait[i] = 1;
     i++;
   }
   fgets(line, lline, ftrait); // throwing away the rest of the line, if only the newline character.

   /* for each individual = line ... */
   for (j=0; j<traitdata.nind;){
     fscanf (ftrait, "%s", name);
     if(name[0]=='/' && name[1]=='/') break; // stop reading the file if "//"
     if (strchr(comment,name[0])){           // ignoring the rest of the line
       fgets(line, lline, ftrait); continue;}
     // otherwise we just read the individual ID, which we don't care about.
     /* ... read the species name */
     if( fscanf (ftrait, "%s", name) == 0 ){
       fprintf(stderr, "Unable to read the species name of individual %d, in trait data file.\n",j+1);
       exit(-1);
     }
     for (jsp=0; jsp<sptree.nspecies; jsp++) {
       if(strcmp(name, sptree.nodes[jsp].name) == 0) {
	 traitdata.indSpeciesMap[j] = jsp; 
	 break;
       }
     }
     if(jsp==sptree.nspecies) {
       fprintf(stderr, "\nspecies '%s' for the %dth individual in trait file not found in the control file.\n",name,j+1);
       exit(-1);
     }
     /* ... read all the numerical (or NA=missing) values */
     for (i=0; i<traitdata.ntrait;i++){
       if(fscanf(ftrait, "%s", name) ==0  || (sscanf(name, "%lf", &traitdata.y[i][j])==0
					      &&     strcmp(name, "NA"))) {
	 fprintf(stderr, "Unable to read trait value for individual %d, trait %d. Read '%s'\n", j+1,i+1,name);
	 exit(-1);
       } // else successfully read either a number or 'NA'
       if (!strcmp(name, "NA")) {
	 traitdata.cleantrait[i] = 0;
	 traitdata.ismissing[i][j]=1;
	 traitdata.y[i][j]=DBL_MAX; 
       }
     }
     fgets(line, lline, ftrait); // throwing away the rest of the line, if only the \n.
     j++;
   }
   traitdata.nind = j;
   fclose(ftrait);

   printf("%d individuals and %d variables.\n", traitdata.nind, traitdata.ntrait);
   if(noisy) printTraitData(stdout);
   fprintf(fout,"\n\n\nPrinting out trait data\n\n");
   printTraitData(fout);

   return(0);
}

int StandardizeTraitData(FILE*fout)
{
  double var;
  double * sumYbySp ;
  int i, j, isp, nsp, * nibySp;

  sumYbySp = (double*)malloc( sptree.nspecies * sizeof(double));
  nibySp   = (int*)   malloc( sptree.nspecies * sizeof(int));

  for (i=0; i<traitdata.ntrait; i++){
    /* initialize then calculate species-specific sums and sample sizes */
    for (isp=0; isp<sptree.nspecies; isp++){
      sumYbySp[isp] = 0.0;
      nibySp[isp] = 0;
    }
    for (j=0; j<traitdata.nind; j++){
      if (!traitdata.ismissing[i][j]){
	nibySp[  traitdata.indSpeciesMap[j] ]++;
	sumYbySp[traitdata.indSpeciesMap[j]] += traitdata.y[i][j];
      }
    }
    /* calculate weighted average */
    traitdata.ybar[i]= 0.0;
    traitdata.ni[i] = 0;
    nsp=0; // # species for which data is available
    for (isp=0; isp<sptree.nspecies; isp++){
      if (nibySp[isp]){ // data available for at least 1 individual from that species
	traitdata.ybar[i] += sumYbySp[isp] / nibySp[isp];
	traitdata.ni[i]   += nibySp[isp];
	nsp++;
      }
    }
    if (traitdata.ni[i] <= 1) { /* quit if only 0 or 1 observation for the trait */
      fprintf(stderr, "Error while standardizing trait data for %dth variable (%s): only %d observation.\n",
	     i+1, traitdata.traitName[i], traitdata.ni[i]);
      exit(-1);
    }
    traitdata.ybar[i] /= nsp;
    /* calculate standard deviation from ybar */
    var = 0.0;
    for (j=0; j<traitdata.nind; j++){
      if (!traitdata.ismissing[i][j]){
	traitdata.y[i][j] -= traitdata.ybar[i]; /* centering now */
	var += traitdata.y[i][j] * traitdata.y[i][j];
      }
    }
    traitdata.scale[i]= sqrt( var / (traitdata.ni[i]-1) );
    /* standardize the y values */
    for (j=0; j<traitdata.nind; j++){
      if (!traitdata.ismissing[i][j]){
	traitdata.y[i][j] /= traitdata.scale[i]; /* re-scaling now */
      }
    }
  }

  free(sumYbySp); free(nibySp);
  /* print means & sds to screen and out file */
  printf(      "Trait sample size, mean (across populations) and standard deviation from that mean: \n    ");
  fprintf(fout,"\nTrait sample size, mean (across populations) and standard deviation from that mean: \n    ");
  for(i=0; i<traitdata.ntrait; i++){
    printf(      " %8s", traitdata.traitName[i]);
    fprintf(fout," %8s", traitdata.traitName[i]);
  }
  printf("\nn   "); fprintf(fout,"\nn   ");
  for(i=0; i<traitdata.ntrait; i++){
    printf(      " %8d", traitdata.ni[i]);
    fprintf(fout," %8d", traitdata.ni[i]);
  }
  printf("\nmean"); fprintf(fout,"\nmean");
  for(i=0; i<traitdata.ntrait; i++){
    printf(      " %8.3g", traitdata.ybar[i]);
    fprintf(fout," %8.3g", traitdata.ybar[i]);
  }
  printf("\nsd  "); fprintf(fout,"\nsd  ");
  for(i=0; i<traitdata.ntrait; i++){
    printf(      " %8.3g", traitdata.scale[i]);
    fprintf(fout," %8.3g", traitdata.scale[i]);
  }
  printf(      "\n"); fflush(stdout);
  fprintf(fout,"\n"); fflush(fout);
  /* print standardized values to out file */
  if(noisy) printf("Standardized trait data:\n");
  if(noisy) printTraitData(stdout);
  fprintf(fout,"\nStandardized trait data\n\n");
  printTraitData(fout);

  return(0);
}

int GetMem (char ipop[])
{
/* This routine is called after the species tree and sequence data are read
   to allocate memory for gene trees (gnodes) and conditional probabilities 
   (conP).  The size of gnodes[locus] is determined using data.ns[locus], 
   no larger than nodes_t (static memory allocated using NS).  

   This routine also initializes ipop at tips in gene trees: gnodes[].ipop.
   It is not used by the simulation program.

   The conditional probabilities for internal nodes are com.conPin[2],
   allocated according to data.ns[locus] and data.npatt[locus].
   Two copies of the space are allocated, hence the [2].  The copy used for the 
   current gene trees is conPin[com.curconP] while the copy for proposed gene
   trees is the alternative copy conPin[!com.curconP].  Even in the alternative 
   copy, the conP space for different loci do not overlap.  

   UpdateGB_InternalNode:
   UpdateGB_SPR:
     Those two proposal steps use UseLocus() and AcceptLocus() to copy between the 
     two copies of conPin.  They do not change com.curconP.

   UpdateTau:
   mixing:
     Those two steps use SwitchconPin to accept a locus, instead of memcpy to copy 
     the whole conPin space.

   UpdateTheta: This does not change the lnL.

   Space for heredity scalars and locus rates is allocated in ReadSeqData().
   It would be fine to do it here.
*/
   int locus,i,j,k, s1, stree; // stree = number of bytes for gene trees
   double *conP;

   /* get mem for gnodes */
   for(locus=0,stree=0; locus<data.ngene; locus++)
      stree += (data.ns[locus]*2-1)*sizeof(struct TREEN);
   if((gnodes[0]=(struct TREEN*)malloc(stree))==NULL) 
      error2("oom gtree");
   for(locus=0,nodes=gnodes[0]; locus<data.ngene-1; locus++)
      gnodes[locus+1] = (nodes+=data.ns[locus]*2-1);
   for(locus=0; locus<data.ngene; ipop+=data.ns[locus++]) {
      for(j=0; j<2*data.ns[locus]-1; j++) {
         gnodes[locus][j].nson = 0; 
         gnodes[locus][j].age = 0; 
         gnodes[locus][j].ibranch = -1;
      }

#ifdef SIMULATION
      for(i=0,k=0; i<sptree.nspecies; i++) /* ipop for tips in gene tree */
         for(j=0; j<data.nseqsp[locus][i]; j++)
            gnodes[locus][k++].ipop = i;
#else
      for(i=0; i<data.ns[locus]; i++)      /* ipop for tips in gene tree */
         gnodes[locus][i].ipop = ipop[i];
#endif
   }

   /* get mem for conP (in nodes) */
   if(mcmc.usedata) {
      data.conP_offset[0] = 0;
      for(locus=0,com.sconP=0; locus<data.ngene; locus++) {
         s1 = com.ncode * data.npatt[locus];
         com.sconP += s1*(data.ns[locus]-1)*sizeof(double);
         if(locus<data.ngene-1) 
            data.conP_offset[locus+1] = data.conP_offset[locus]+(data.ns[locus]-1)*s1;
      }
      if((com.conPin[0]=(double*)malloc(2*com.sconP))==NULL) 
         error2("oom conP");
      com.conPin[1]        = com.conPin[0]+com.sconP/sizeof(double);
      printf("%d bytes for conP, %d bytes for gene trees.  \n\n", 2*com.sconP, stree);

      /* set gnodes[locus][].conP for tips and internal nodes */
      com.curconP = 0; conP = com.conPin[0];
      for(locus=0; locus<data.ngene; locus++) {
         if(!data.cleandata[locus]) /* Is this still necessary */
            UseLocus(locus, 0, 1, 0);
         for(j=data.ns[locus]; j<data.ns[locus]*2-1; j++,conP+=com.ncode*data.npatt[locus])
            gnodes[locus][j].conP = conP;
      }
   }
   return(0);
}

void FreeMem(void)
{
   int locus, j;

   free(gnodes[0]);
   if(mcmc.usedata) free(com.conPin[0]);
   for(locus=0; locus<data.ngene; locus++) {
      free(data.fpatt[locus]);
      for(j=0;j<data.ns[locus]; j++)
         free(data.z[locus][j]);
   }
   for (j=0; j<traitdata.ntrait; j++){
     free(traitdata.traitName[j]);
   }
}

int UseLocus (int locus, int copytreeconP, int useData, int setSeqName)
{
/* This point *nodes to the gene tree at locus gnodes[locus] and set com.z[] 
   etc. for likelihood calculation for the locus.
   It also updates sptree.nseqsp[], which is used in UpdateGB_SPR() to 
   calculate ndesc[].
   If(copytreeconP), the genetree gnodes[locus] is copied to nodes_t.  
   If(useData), com. and nodes[].conP are repositioned.
   If (copytreeconP && useData), the conP for internal nodes point to 
   the alternative space com.conPin[!com.curconP]+data.conP_offset[locus].  
   conP for different loci do not overlap, so that this routine is used by 
   all proposal steps, some of which change one locus only while others 
   change several or all loci.
*/
   int i,j,k, s1=com.ncode*data.npatt[locus], is,nseqsp[NSPECIES]={0};
   double *conPa=com.conPin[!com.curconP]+data.conP_offset[locus];

   com.ns=data.ns[locus]; tree.root=data.root[locus]; tree.nnode=2*com.ns-1;
   for(i=0; i<sptree.nspecies; i++) 
      sptree.nseqsp[i] = data.nseqsp[locus][i];

   if(copytreeconP) { 
      /* Old genetree is copied into nodes_t.  This is needed if the proposal 
      changes genetree topology, but may be unneccesary otherwise.
      Old conP is copied into alternative space.  Think whether this is needed. 
      */
      nodes = nodes_t;
      memcpy(nodes_t, gnodes[locus], (com.ns*2-1)*sizeof(struct TREEN));
      if(useData) {
         for(i=com.ns; i<com.ns*2-1; i++)
            nodes_t[i].conP = conPa + (i-com.ns)*s1;
         memcpy(conPa, gnodes[locus][com.ns].conP, s1*(com.ns-1)*sizeof(double));
      }
   }
   else           /* this works on the old gene tree */
      nodes = gnodes[locus];
   
   if(useData) {
      com.cleandata = data.cleandata[locus];
      com.npatt = data.npatt[locus];
      com.fpatt = data.fpatt[locus];
      for(i=0; i<com.ns; i++) 
         com.z[i] = data.z[locus][i];
   }
   
   if(setSeqName) {
      for(i=0; i<tree.nnode; i++) {
         if(i!=tree.root) 
            nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
      }
      for(i=0; i<com.ns; i++) {
         is = gnodes[locus][i].ipop;
         nseqsp[is]++;
         if(sptree.nseqsp[is]>1) 
            sprintf(com.spname[i], "%s%d", sptree.nodes[is].name, nseqsp[is]);
         else
            sprintf(com.spname[i], "%s", sptree.nodes[is].name);
      }

   }
   return(0);
}

int AcceptLocus (int locus, int copyconP)
{
/* If (copytree && usedata), this copies nodes_t[].conP into gnodes[].conP for inner 
   nodes 
*/
   int i,ns=data.ns[locus], s1=com.ncode*data.npatt[locus];
   double *conPt;

   /* copies node_t into gnodes[locus] */
   data.root[locus] = tree.root;
   memcpy(gnodes[locus], nodes_t, (ns*2-1)*sizeof(struct TREEN));
   if(mcmc.usedata) {  /* reposition gnodes[][].conP */
      conPt = com.conPin[com.curconP]+data.conP_offset[locus];
      for(i=ns; i<ns*2-1; i++)
         gnodes[locus][i].conP = conPt + (i-ns)*s1;
      if(copyconP) {/* used in UpdateGB_InternalNode() and UpdateGBTip() only */
         conPt = com.conPin[!com.curconP]+data.conP_offset[locus];
         memcpy(gnodes[locus][ns].conP, conPt, (ns-1)*s1*sizeof(double));
      }
   }
   return(0);
}

void SwitchconPin(void)
{
/* This resets gnodes[locus].conP to the alternative com.conPin, to avoid 
   recalculation of conP, when a proposal is accepted in UpdateTau() or mixing().
*/
   int i,locus;
   double *conP;

   com.curconP=!com.curconP;
   conP=com.conPin[com.curconP];
   for(locus=0; locus<data.ngene; locus++)
      for(i=data.ns[locus]; i<data.ns[locus]*2-1; i++,conP+=com.ncode*data.npatt[locus])
         gnodes[locus][i].conP = conP;
}

void CopyconPin (int locus, int FromCurrent)
{
/* This copies com.conPin[] from the current place to the alternative place or 
   vice versa.
*/
   int size = (data.ns[locus]-1)*com.ncode*data.npatt[locus]*sizeof(double);
   double *from = com.conPin[ com.curconP] + data.conP_offset[locus];
   double *to   = com.conPin[!com.curconP] + data.conP_offset[locus], *y;

   if(FromCurrent==0) {
      y=from; from=to; to=y;
   }
   memcpy(to, from, size);
}



static double PMat[4*4];

int ConditonalPNode (int inode)
{
   int n=com.ncode, i, j, k, h, ison, ispecies;
   double y;

   for(i=0; i<nodes[inode].nson; i++) {
      ison=nodes[inode].sons[i];
      if (nodes[ison].nson>0  && (!mcmc.saveconP || !com.oldconP[ison]))
         ConditonalPNode (nodes[inode].sons[i]);
   }
   for(h=0; h<com.npatt*n; h++)
      nodes[inode].conP[h] = 1;

   for(i=0; i<nodes[inode].nson; i++) {
      ison = nodes[inode].sons[i];
      ispecies = nodes[ison].ipop;
      NPMat ++;

      /*
      PMatTN93 (PMat, nodes[ison].branch/3, nodes[ison].branch/3, nodes[ison].branch/3, com.pi);
      */
      pijJC69 (PMat, nodes[ison].branch);

      if(nodes[ison].nson>0) {                          /* internal node */
         for(h=0; h<com.npatt; h++) 
            for(j=0; j<n; j++) {
               for(k=0,y=0; k<n; k++)
                  y += PMat[j != k] * nodes[ison].conP[h*n+k];
               nodes[inode].conP[h*n+j] *= y;
            }
      }
      else if (com.cleandata) {                         /* tip && clean */
         if(data.iseqerr[ispecies] == 0) { 
            for(h=0; h<com.npatt; h++) 
               for(j=0; j<n; j++)
                  nodes[inode].conP[h*n+j] *= PMat[ j != com.z[ison][h] ];
         }
         else {                                         /* has errors */
            for(h=0; h<com.npatt; h++) 
               for(j=0; j<n; j++) {
                  for(k=0,y=0; k<n; k++)
                     y += PMat[ j != k ] * data.e_seqerr[ispecies][k*4 + com.z[ison][h]];
                  nodes[inode].conP[h*n+j] *= y;
               }
         }
      }
      else if (!com.cleandata) {                        /* tip & unclean */
         for(h=0; h<com.npatt; h++) {
            if(!data.iseqerr[ispecies] || nChara[com.z[ison][h]]>1) { /* no errors || ambiguity */
               for(j=0; j<n; j++) {
                  for(k=0,y=0; k<nChara[com.z[ison][h]]; k++)
                     y += PMat[ j != CharaMap[com.z[ison][h]][k] ];
                  nodes[inode].conP[h*n+j] *= y;
               }
            }
            else {                                      /* errors & good base (TCAG) */
               for(j=0; j<n; j++) {
                  for(k=0,y=0; k<n; k++)
                     y += PMat[ j != k ] * data.e_seqerr[ispecies][k*4 + CharaMap[com.z[ison][h]][0] ];
                  nodes[inode].conP[h*n+j] *= y;
               }
            }
         }
      }
   }   /* for(ison) */
   return (0);
}

double lnpData (double lnpDi[])
{
/* This calculates the log likelihood, the log of the probability of the data 
   given gtree[nlocus] and blength[nlocus][2] for each locus.
   This updates gnodes[locus][].conP.
*/
   int j,locus;
   double lnL=0,y;

   if(mcmc.saveconP) 
      for(j=0; j<data.maxns*2-1; j++)  com.oldconP[j]=0;
   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus,0,1, 0);
      y = lnpD_locus(locus);
      if(testlnL && fabs(lnpDi[locus]-y)>1e-10) 
         printf("\tlnLi %.6f != %.6f at locus %d\n", lnpDi[locus], y, locus+1);
      lnpDi[locus]=y;
      lnL += y;
   }
   return(lnL);
}

double lnpD_locus (int locus)
{
/* This calculates ln{Di|Gi, Bi} using nodes[].age and tree.root.
   Note that nodes[].branch may not be kept up to date in the program.
*/
   int  h, i, is, haserr=0;
   double lnL=0, fh;

   if(!mcmc.usedata) return(0);
   for(i=0; i<tree.nnode; i++) {
      if(i==tree.root) continue;   
      nodes[i].branch = nodes[nodes[i].father].age - nodes[i].age;
      if(data.est_locusrate)
         nodes[i].branch *= data.locusrate[locus];
      if(nodes[i].branch < 0) {
         printf("branch length = %.6f < 0\nReport to Ziheng", nodes[i].branch);
      }
   }

   /*
   for(is=0; is<sptree.nspecies; is++)  
      if(data.iseqerr[is]) { haserr=1; break; }
   if(haserr) {
      for(is=0; is<sptree.nspecies; is++) {
         if(data.iseqerr[is] == 0) continue;
         for(i=0; i<com.ns; i++) 
            if(nodes[i].ipop == is) 
               nodes[i].branch += data.seqerr[is];
      }
   }
   */

   ConditonalPNode (tree.root);

   for (h=0; h<com.npatt; h++) {
      for (i=0,fh=0; i<com.ncode; i++) 
         fh += com.pi[i]*nodes[tree.root].conP[h*com.ncode+i];
      lnL += log(fh)*com.fpatt[h];
   }
   return (lnL);
}


double lnpGB_ThetaTau (int locus)
{
/* this calculates the prior of the gene tree and coalescent times (using tree 
   and nodes[]), when theta and tau are fixed from sptree.
   It initially collects information (coal) about the coalescent in each 
   population, such as nin, nout, coal, tj, by looping through all nodes[] 
   in the gene tree.  Then lnpGB is calculated by summing over the populations.
   Note that this does not use sptree.nseqsp[].
*/
   int i,j,k,inode, ip, ipop, isonpop, nt;
   double lnp=0, starttime,endtime, *t,y, H=(data.est_heredity==0?1:data.heredity[locus]);
   double small=0;
   /* double small=1e-12; */
   int nin[NSPECIES*2-1], ncoal[NSPECIES*2-1], nout[NSPECIES*2-1];
   double tj[NSPECIES*2-1][NS-1];

   if(debug==1) 
      puts("\nIn lnpGB_ThetaTau ");

   /* construct info (in, coal, out) for each population (the coal table) */
   for(i=0; i<sptree.nnode; i++) 
      nin[i] = ncoal[i] = nout[i] = 0;
   for(inode=0,nout[sptree.root]=1; inode<tree.nnode; inode++) {
      ipop = nodes[inode].ipop;
      if(nodes[inode].nson==0) /* tip */
         nin[ipop]++;
      else {                   /* internal node */
         tj[ipop][ncoal[ipop]++] = nodes[inode].age;
         for(i=0; i<nodes[inode].nson; i++) {
            isonpop = k=nodes[nodes[inode].sons[i]].ipop;
            if(isonpop == ipop) continue;
            for(; ; k=sptree.nodes[k].father) {
               if(k==isonpop) nout[k]++;
               else if(k==ipop)  { nin[k]++; break; }
               else { nin[k]++; nout[k]++; }
            }
         }
      }
      if(inode==tree.root && ipop!=sptree.root) {
         for(k=ipop; ; k=sptree.nodes[k].father) {
            if(k!=ipop)        nin[k]++;
            if(k!=sptree.root) nout[k]++;
            if(k==sptree.root) break;
         }
      }
   }
   for(ip=0; ip<sptree.npop; ip++) { /* sort tj for each pop */
      ipop = sptree.pops[ip]; 
      if(nin[ipop] - ncoal[ipop] != nout[ipop]) {
         printf("\nin out coal %d %d %d in ipop %d\n", 
            nin[ipop], nout[ipop], ncoal[ipop], ipop);
         printGtree(1);
         error2("in out ");
      }
      nt = ncoal[ipop]; 
      t = tj[ipop];
      if(nt>30)
         qsort(t, (size_t)nt, sizeof(double), comparedouble);
      else 
         for(i=0; i<nt-1; i++)  for (j=i+1; j<nt; j++)
            if (t[j]<t[i])  { y=t[i]; t[i]=t[j]; t[j]=y; } 
   }

   if(debug==1) {
      printf("\ncoalescence info in populations:\nspecies  in coal  out   tau     tj\n");
      for(i=0; i<sptree.nnode; i++,FPN(F0)) {
         printf("%4d: %5d%5d%5d %9.5f", i, nin[i], ncoal[i], nout[i],sptree.nodes[i].age);
         for(j=0; j<ncoal[i]; j++) printf(" %9.5f", tj[i][j]);
      }
      FPN(F0);
   }

   for(ip=0; ip<sptree.npop; ip++) {
      ipop = sptree.pops[ip]; 
      if(sptree.nodes[ipop].theta <= 0) continue;
      t = tj[ipop];
      starttime = sptree.nodes[ipop].age;
      if(ipop!=sptree.root) endtime = sptree.nodes[sptree.nodes[ipop].father].age;
      else                  endtime = -1;

      if(debug==1) printf("species %d: time: %8.4f --> %8.4f ", ipop,starttime,endtime);

      if(nin[ipop]>nout[ipop]) {
         lnp += (nin[ipop] - nout[ipop]) * log(2/(sptree.nodes[ipop].theta*H));
         for(j=nin[ipop],k=0; j>nout[ipop]; j--,k++) {
            y = t[k] - (k==0?starttime:t[k-1]);
            lnp -= j*(j-1)/(sptree.nodes[ipop].theta*H)*y;

            if(debug==1) printf(" j=%d: %9.5f", j,y);

            if(y<-small && noisy) {
               printGtree(1);
               printf("\n       tj = %.20f < 0 in ipop %d!\n", y, ipop);
               printf(" node age = %.20f\n", t[k]);
               printf("starttime = %.20f\n", starttime);
               error2("negative time in lnpGB_ThetaTau ");
            }
         }
      }
      if(nout[ipop]>1) {
         y = (nout[ipop]==nin[ipop] ? endtime-starttime : endtime-t[k-1]);
         lnp -= nout[ipop]*(nout[ipop]-1)/(sptree.nodes[ipop].theta*H)*y;

         if(debug==1) printf(" remained %.5f", y);
      }
      if(debug==1) FPN(F0);
   }
   if(debug==1) printf("\npGB=%.6g  lnpGB = %.6f\n", exp(lnp), lnp);
   return(lnp);
}

double UpdateGB_InternalNode (double* lnL, double finetune)
{
/* This slides a node in the gene tree without changing the gene tree
   topology.  However, the new coalescent time tnew might be in different 
   population from the current population, so nodes[inode].ipop might change.
   The algorithm loops through all the internal nodes in the gene tree at each 
   locus.  The lower bound tb[0] is determined by age of the two sons 
   and also the age of the species that is a common ancestor to both sons.  
   The upper bound tb[1] is the father node's age.
*/
   int  accepted=0, locus, i,inode, sons[2],pops[2], ninodes;
   int  ipopsource,ipoptarget, copytree;
   double lnacceptance, lnLd, lnpGBinew, lnpDinew, t,tnew;
   double tb[2];

   if(debug==2) puts("\nUpdateGB_InternalNode");
   if(finetune<=0) error2("steplength = 0 in UpdateGB_InternalNode");
   for(i=0,ninodes=0; i<data.ngene; i++) ninodes += data.ns[i]-1;

   for(locus=0; locus<data.ngene; locus++) {

      if(debug==2) printf("\nlocus %d (ns = %d)\n", locus+1,data.ns[locus]);
      for(inode=data.ns[locus],copytree=1; inode<data.ns[locus]*2-1; inode++) {

         if(copytree) UseLocus(locus, 1, mcmc.usedata, 0); /* work on a copy */
         t = nodes[inode].age;
         ipopsource = nodes[inode].ipop;
         for(i=0;i<2;i++) {
            sons[i] = nodes[inode].sons[i]; 
            pops[i] = nodes[sons[i]].ipop; 
         }
         tb[0] = max2(nodes[sons[0]].age, nodes[sons[1]].age);
         tb[1] = (inode==tree.root ? OLDAGE : nodes[nodes[inode].father].age);

#if(0)   /* no change to ipop, remain in the same population */
         tb[0] = max2(tb[0],sptree.nodes[ipopsource].age);
         if(ipopsource != sptree.root)
            tb[1] = min2(tb[1], sptree.nodes[sptree.nodes[ipopsource].father].age);
#else
         /* Possible change to ipop.  Find the population ancestral to both sons */
         if(pops[0] != pops[1]) { /* sptree.pops[] are from young to old */
            for(i=0; i<sptree.npop; i++) 
               if(sptree.pop_pop_table[pops[0]][sptree.pops[i]]
               && sptree.pop_pop_table[pops[1]][sptree.pops[i]])
                  break;
            tb[0] = max2(tb[0],sptree.nodes[sptree.pops[i]].age);
         }
#endif
         tnew = t + finetune*rnd2NormalSym(m2NormalSym);
         tnew = reflect(tnew, tb[0], tb[1]);

         /* determine ipoptarget by by going up sptree */
         for(i=ipoptarget=pops[0]; ; i=sptree.nodes[i].father) {
            if(sptree.nodes[i].age>tnew) break;
            ipoptarget = i;
            if(i==sptree.root) break;
         }

         /* two things are updated in the loop: node age and ipop */
         nodes[inode].age = tnew;
         nodes[inode].ipop = ipoptarget;

         if(mcmc.saveconP) {
            for(i=com.ns; i<com.ns*2-1; i++) com.oldconP[i]=1;
            for(i=inode; ; i=nodes[i].father)
               { com.oldconP[i]=0; if(i==tree.root) break; }
         }


         lnpGBinew = lnpGB_ThetaTau(locus);
         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew - data.lnpDi[locus];
         lnacceptance = lnpGBinew - data.lnpGBi[locus] + lnLd;

         if(debug==2) {
            printf("node %2d tb (%8.5f %8.5f) t %7.5f -> %7.5f lnLd %9.5f ",inode,tb[0],tb[1],t,tnew,lnLd);
            for(i=com.ns; i<com.ns*2-1; i++) printf("%d",com.oldconP[i]);
         }

         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
            accepted++;
            *lnL += lnLd;
            data.lnpDi[locus] = lnpDinew;
            data.lnpGBi[locus] = lnpGBinew;
            AcceptLocus(locus, 1);
            copytree=0;

            if(debug==2) printf(" Y\n");

         }
         else {
            nodes[inode].age = t;
            nodes[inode].ipop = ipopsource;
            copytree = 1;

            if(debug==2) printf(" N\n");
         }
      }
   }
   return((double)accepted/ninodes);
}


void GraftNode(int source, int target, double age, int ipop)
{
/* This performs surgery on the gene tree (nodes) to cut branch source 
   (branch ancestral to node source) and regrafts it to branch target.
   It assumes binary rooted tree.
   age & ipop are for the new node (father) to be created on the target branch.
   After the operation, the following branches are generated:
     grandpa->sib (at the source)
     father->source, father->target, targetfather->father (at the target)

   IMPORTANT: This keeps the node number for tree.root unchanged, as required 
   by UpdateGBTip.
*/
   int i,k, father, grandpa=-1, sib, targetfather, oldroot=tree.root;
   double y;

   father = nodes[source].father;
   nodes[father].age = age;  
   nodes[father].ipop = ipop;
   if(father!=tree.root) grandpa = nodes[father].father;
   sib = nodes[father].sons[0] + nodes[father].sons[1] - source;
   targetfather = nodes[target].father;

   if(target!=sib && target!=father) { /* change to tree topology */
      /* cut at source */
      if (father==tree.root) {         /* old root is lost */
         tree.root = sib;  nodes[sib].father = -1;
      }
      else {
         if(nodes[grandpa].sons[i=0] != father)
            i = 1;
         nodes[grandpa].sons[i] = sib; 
         nodes[sib].father = grandpa;
      }

      /* regraft */
      nodes[father].father = nodes[target].father;
      nodes[target].father = father;
      nodes[father].sons[0] = target;
      nodes[father].sons[1] = source;
      if (target == tree.root) {  /* new root is born */
         tree.root = father;
      }
      else
         for(i=0; i<2; i++) {
            if(nodes[targetfather].sons[i]==target) {
               nodes[targetfather].sons[i]=father; 
               break; 
            }
         }
   }
   if(tree.root != oldroot) { /* swap everything except .conP */
      /* if(debug==3)  printGtree(1); */
      swap2(nodes[tree.root].ipop, nodes[oldroot].ipop, k);
      swap2(nodes[tree.root].age, nodes[oldroot].age, y);
      for(i=0; i<2; i++) {
         nodes[nodes[tree.root].sons[i]].father = oldroot;
         nodes[nodes[oldroot].sons[i]].father = tree.root;
      }
      for(i=0; i<2; i++)
         swap2(nodes[tree.root].sons[i], nodes[oldroot].sons[i], k);
      for(i=0; i<2; i++)
         if(nodes[oldroot].sons[i]==oldroot)
            nodes[oldroot].sons[i] = tree.root;
      k = nodes[oldroot].father;
      for(i=0; i<2; i++)
         if(nodes[k].sons[i]==oldroot)
            nodes[k].sons[i] = tree.root;
      if((nodes[tree.root].father=k) == tree.root)
         nodes[tree.root].father = oldroot;
      nodes[oldroot].father = -1;

      father = nodes[source].father;
      if(grandpa == tree.root)    grandpa = oldroot;
      else if(grandpa == oldroot) grandpa = tree.root;
      tree.root = oldroot;
   }
   if(mcmc.usedata && mcmc.saveconP) {
      for(i=com.ns; i<com.ns*2-1; i++) com.oldconP[i]=1;
      for(i=father; ; i=nodes[i].father)
         { com.oldconP[i]=0; if(i==tree.root) break; }
      for(i=grandpa; i!=-1; i=nodes[i].father)
         { com.oldconP[i]=0; if(i==tree.root) break; }
   }
}

double UpdateGB_SPR (double* lnL, double finetune)
{
/* This removes a node (source) in the gene tree and regrafts it to a feasible 
   branch.  There is no change to the node times inside the subtree.  
   The step cycles through all loci and all nodes except the root.
   ndesc[] holds the number of descendent tips at the locus of each population.
   For example, ndesc[sptree.root]==com.ns.  This is used to identify feasible 
   populations (species) to determine range[0] for tnew.  It is used to deal 
   with loci at which some species are missing.  ndesc[] is calculated from 
   sptree.nseqsp[species], which is updated for the locus in UseLocus().  

   Note that inode is the root node for the subtree being moved, species is 
   the population inode is in, and t and tnew are the time at which the subtree
   coalesce with the rest of the gtree, that is, the age of inode's father.  
   To determine the lower bound tb[0] for tnew, the following requirements 
   are considered:
   (a) tnew must be older than the age of the subtree: nodes[inode].age.
   (b) The population of the target node must be ancestral to or the same as 
       the population of the subtree: (sptree.pop_pop_table[species][ipop]==1).
   (c) There must be other lineages outside the subtree passing that population
       to which the subtree can join.  For tips (inode<com.ns), this requirement 
       is satisfied if (ndesc[ipop]>1).  For internal nodes, (ndesc[ipop]>1) is 
       necessary but not sufficient, and the algorithm loops through all branches 
       (nodes) in the gene tree, and identifies those outside subtree whose 
       father node is older than tb[0].  For each of such branches, we 
       record the youngest age at which the subtree can join.
*/
   int accepted=0,nproposal, nsource, ntarget,targets[NS*2-1],target, copytree;
   int locus, i,j, inode, species, ipop, ipopsource, ipoptarget, father, sib;
   int ndesc[2*NSPECIES-1]; /* #descendents for ancestral species */
   double lnacceptance=0, lnLd, lnpGBinew,lnpDinew;
   double y,tb[2], t,tnew;

   if(debug==3) puts("\nUpdateGB_SPR ");
   if(finetune<=0) error2("steplength = 0 in UpdateGB_InternalNode");

   for(i=0,nproposal=0; i<data.ngene; i++) 
      nproposal += (mcmc.moveinnode ? data.ns[i]*2-1 : data.ns[i]-1);

   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 0, 0, 0); /* this resets sptree.nseqsp[] */

      if(debug==3) {
         printf("\nlocus %d (ns = %d)\n", locus+1,data.ns[locus]);
         for(i=sptree.nspecies; i<2*sptree.nspecies-1; i++) 
            printf("sp.time[%d] = %.6f\n", i, sptree.nodes[i].age);
         printGtree(1);
      }
      for(i=0; i<2*sptree.nspecies-1; i++) ndesc[i]=0;
      for(species=0; species<sptree.nspecies; species++) { /* update ndesc[] for each locus */
         ndesc[species] = sptree.nseqsp[species];
         for(i=sptree.nspecies; i<2*sptree.nspecies-1; i++)
            if(sptree.pop_pop_table[species][i]) 
               ndesc[i] += sptree.nseqsp[species];
      }
      if(ndesc[sptree.root] != com.ns) error2("ndesc");

      for(inode=0,copytree=1; inode<(mcmc.moveinnode?com.ns*2-1:com.ns); inode++) {
         /* NOTE: tree.root should not be changed in GraftNode() */
         if(inode==tree.root) continue; 
         if(copytree) UseLocus(locus, 1, mcmc.usedata, 0);
         species = nodes[inode].ipop;
         father = nodes[inode].father;
         sib = nodes[father].sons[0] + nodes[father].sons[1] - inode;
         ipopsource = nodes[father].ipop;
         t = nodes[father].age;

         /* locate lower bound tb[0] for tnew.  See notes above */
         tb[0] = nodes[inode].age;
         tb[1] = OLDAGE;
         for(i=0; i<sptree.npop; i++) {
            ipop = sptree.pops[i];
            if(sptree.nodes[ipop].theta>0 && sptree.pop_pop_table[species][ipop] && ndesc[ipop]>1)
               { tb[0] = max2(tb[0], sptree.nodes[ipop].age); break; }
         }
         /* For moving internal nodes, find min feasible speciation time (y) 
            in remaining subtree.  See note (c) above. 
         */
         if(inode>=com.ns && sptree.nspecies>1 && sptree.nodes[sptree.root].age>0) {
            for(i=0,y=sptree.nodes[sptree.root].age; i<tree.nnode; i++) {
               if(i==tree.root || nodes[i].father==inode || nodes[nodes[i].father].age<=tb[0])
                  continue;
               /* Is node i a descendent of inode? */
               for(j=i; j!=tree.root&&j!=inode; j=nodes[j].father) ;
               if(j==inode) continue;  /* skip node i if it is a descendent */
               for(ipop=nodes[i].ipop; ipop!=sptree.root; ipop=sptree.nodes[ipop].father)
                  if(sptree.nodes[ipop].theta>0 && sptree.pop_pop_table[species][ipop])
                     { y=min2(y, sptree.nodes[ipop].age); break; }
               if(y<tb[0]) break;  /* no point for checking further */
            }
            tb[0] = max2(tb[0],y);
         }

         tnew = t + finetune*rnd2NormalSym(m2NormalSym);
         tnew = reflect(tnew, tb[0], tb[1]);

         /* identify the target pop in which tnew is, by going up sptree */
         for(ipop=ipoptarget=species; ; ipop=sptree.nodes[ipop].father) {
            if(sptree.nodes[ipop].age > tnew) break;
            ipoptarget = ipop;
            if(ipop == sptree.root) break;
         }

         if(debug==3) printf("inode %2d father %2d tb (%6.4f %7.4f) time %7.5f (%2d) ->%7.5f (%2d) ",
                              inode,father,tb[0],tb[1],t,ipopsource,tnew,ipoptarget);

         /* count and identify feasible target branches. Watch out for root */
         ntarget=0;  nsource=1;
         if(tnew >= nodes[tree.root].age)
            targets[ntarget++] = tree.root;
         else
            for(i=0; i<tree.nnode; i++) { /* Is node i a possible target? */
               if(i!=inode && i!=tree.root && nodes[i].age<=tnew && nodes[nodes[i].father].age>tnew
                  && sptree.pop_pop_table[nodes[i].ipop][ipoptarget]) 
                  targets[ntarget++] = (i==father ? sib : i);
            }
         if(father!=tree.root) 
            for(i=0; i<tree.nnode; i++) { /* Is node i a possible source? */
               if(i!=inode && i!=tree.root && i!=sib && i!=father
                  && nodes[i].age<=t && nodes[nodes[i].father].age>t
                  && sptree.pop_pop_table[nodes[i].ipop][ipopsource])
                  nsource++;
            }

         if(nsource<1 || ntarget<1) {
            printf("\nsource target %d %d", nsource,ntarget); 
            printf("\nlocus %d node %d\n", locus, inode); 
            printGtree(1);
            exit(-1); 
         }
         target = targets[(int)(ntarget*rndu())];
         GraftNode(inode, target, tnew, ipoptarget);

         if(debug==3) {
            printf("line %d > %2d (", nsource,ntarget);
            for(i=0; i<ntarget; i++) printf(" %2d", targets[i]);
            printf(") > %2d ", target);
            printGtree(1);
         }
         lnacceptance = log((double)ntarget/nsource);
         lnpGBinew = lnpGB_ThetaTau(locus);
         lnacceptance += lnpGBinew - data.lnpGBi[locus];
         lnpDinew = lnpD_locus(locus);
         lnLd = lnpDinew-data.lnpDi[locus];
         lnacceptance += lnLd;

         if(debug==3) { 
            printf(" lnLd %9.5f ", lnLd);
            for(i=com.ns; i<com.ns*2-1; i++) printf("%d",com.oldconP[i]); 
         }

         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
            accepted++;
            *lnL += lnLd;
            data.lnpDi[locus] = lnpDinew;
            data.lnpGBi[locus] = lnpGBinew;
            AcceptLocus(locus, 1);
            copytree = 0;
            if(debug==3) printf(" Y (%4d)\n", NPMat);
         }
         else {
            copytree = 1;
            if(debug==3) printf(" N (%4d)\n", NPMat);
         }
      }
   }
   return((double)accepted/nproposal);
}


void ReScaleSubTree(int inode, double factor)
{
   int i;
   nodes[inode].age *= factor;
   for(i=0; i<nodes[inode].nson; i++)
      if(nodes[inode].sons[i] >= com.ns)
         ReScaleSubTree(nodes[inode].sons[i], factor);
}

int NodeSlider(double eps)
{
/* This algorithm perturbs the tree by sliding a node.  Right now it is written 
   for rooted trees only, and assumes that the tree is binary,  Some changes (not 
   done) are necessary to deal with unrooted trees.
   Right now this does not consider the constraints placed on the gene tree by the 
   species tree.
   It does not seem a good idea to cycle through the nodes on the tree, because the 
   node numbers are changed during each cycle.  Instead it may be better to pick 
   up a node at random for sliding.
*/
   int i,j, inode, dad0, dad, target, iround;
   double offset=0, tnew, factor;

   debug=11;
   if(nodes[tree.root].nson!=2)
      error2("NodeSlider: write code for unrooted trees");
   for(inode=0; inode<2*com.ns-1; inode++) {
      if(inode == tree.root) continue;
      offset = eps * rnd2NormalSym(m2NormalSym);

      if(debug) printf("\n\nnode %d offset %9.5f\n", inode, offset);
      
      dad0 = dad = nodes[inode].father;
      target = nodes[dad0].sons[0] + nodes[dad0].sons[1] - inode;  /* sib */
      offset += nodes[dad0].age - nodes[target].age;
      for(iround=0; ; iround++) {  /* forever reflecting */
         if(offset>0) {
            dad = nodes[target].father;
            if(dad==dad0) 
               dad = nodes[dad].father;
            if(target==tree.root)
               break;
            tnew = nodes[target].age + offset;
            if(dad==-1 || tnew<nodes[dad].age)
               break;
            else if(rndu()<0.5) {  /* up to dad */
               offset -= nodes[dad].age - nodes[target].age;
               target = dad;
            }
            else {                 /* down to sib */
               j = nodes[dad].sons[0];
               if(j == target || j == nodes[target].father)
                  j = nodes[dad].sons[1];
               target = j;
               offset = 2*nodes[dad].age - nodes[target].age - tnew;
            }
         }
         else {
            if(target<com.ns) /* tip */
               offset *= -1;
            else {
               target = nodes[target].sons[rndu()<0.5];
               offset += nodes[nodes[target].father].age - nodes[target].age;
            }
         }

         if(debug) printf("round %d target %2d offset %9.5f\n", iround+1, target, offset);

      }  /* forever loop */

      factor = (nodes[target].age+offset)/nodes[dad0].age;
      nodes[dad0].age = nodes[target].age + offset;
      /* rescale the branch lengths */
      if(inode>=com.ns) {
         ReScaleSubTree(inode, factor);
         if(debug) printf("factor = %9.5f\n", factor);
      }
      if(nodes[target].father!=dad0) {
         /* ipop is not updated and the algorithm does not work for the model anyway. */
         GraftNode(inode, target, nodes[target].age+offset, 0);
      }

      if(debug) printGtree(1);
   }  /*  for(inode) */
   return(0);
}


double UpdateTheta (double finetune, double space[])
{
/* This updates theta's one by one, using proportional expansion or shrinkage.
   Perhaps change to sliding window to avoid getting tracked at 0.
   This step does not require calculation of the likelihood.
*/
   int i, ipop, locus, accepted=0, ntheta=0;
   double thetaold, thetanew, lnacceptance;
   double *lnpGBinew=space, a=data.theta_prior[0], b=data.theta_prior[1];

   if(finetune<=0) error2("steplength = 0 in UpdateGB_InternalNode");
   for(i=0; i<sptree.npop; i++) {
      /* prior and proposal ratios */
      ipop = sptree.pops[i];
      if(sptree.nodes[ipop].theta < 0) continue;
      ntheta ++;
      thetaold = sptree.nodes[ipop].theta;
 
      thetanew = thetaold + finetune*rnd2NormalSym(m2NormalSym);
      if(thetanew<0) thetanew = -thetanew;
      sptree.nodes[ipop].theta = thetanew;

      lnacceptance = (a-1)*log(thetanew/thetaold) - b*(thetanew-thetaold);

      for(locus=0; locus<data.ngene; locus++) {
         UseLocus(locus, 1, 0, 0);
         lnpGBinew[locus] = lnpGB_ThetaTau(locus);
         lnacceptance += lnpGBinew[locus] - data.lnpGBi[locus];
      }
      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         for(locus=0; locus<data.ngene; locus++)
            data.lnpGBi[locus] = lnpGBinew[locus];
         accepted++;
      }
      else
         sptree.nodes[ipop].theta = thetaold;
   }
   return((double)accepted/ntheta);
}


double UpdateTau (double *lnL, double finetune, double space[])
{
/* This updates speciation times tau using the rubber-band algorithm.
   ntj[] are counts of nodes below and above tau (m and n in the paper).
*/
   int k, is, inode, locus, accepted=0, ntau=0;
   int ntj[2], ntj_locus[2];
   double tauold, taunew, taub[2], t, taufactor[2];
   double lnacceptance=0, lnLd;
   double *lnpGBinew=space, *lnpDinew=space+data.ngene;
   double a=data.tau_prior[0], b=data.tau_prior[1];

   if(debug==5) puts("\nUpdateTau ");
   if(finetune<=0) error2("steplength = 0 in UpdateTimes");
   for(is=sptree.nspecies; is<sptree.nnode; is++)
      if(sptree.nodes[is].age > 0) ntau++;

   for(is=sptree.nspecies; is<sptree.nnode; is++) {
      if(sptree.nodes[is].age == 0) continue;
      tauold = sptree.nodes[is].age;
      lnLd = 0;
      taub[0] = 0;
      taub[1] = OLDAGE;
      if(sptree.nodes[is].nson)
         taub[0] = max2(sptree.nodes[sptree.nodes[is].sons[0]].age,
                      sptree.nodes[sptree.nodes[is].sons[1]].age);
      if(is != sptree.root)
         taub[1] = sptree.nodes[sptree.nodes[is].father].age;

      taunew = tauold + finetune*rnd2NormalSym(m2NormalSym);
      taunew = sptree.nodes[is].age = reflect(taunew,taub[0],taub[1]);
      for(k=0;k<2;k++)
         taufactor[k] = (taunew-taub[k])/(tauold-taub[k]);
      lnacceptance = 0;
      if(is==sptree.root)
         lnacceptance = (a-1 - ntau+1)*log(taunew/tauold) - b*(taunew-tauold);

      if(debug==5) printf("species %d taub: %8.5f %8.5f tau: %8.5f%8.5f", is,taub[0],taub[1],tauold,taunew); 
	  
      ntj[0]=ntj[1]=0;
      for(locus=0; locus<data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & conP */
         ntj_locus[0] = ntj_locus[1] = 0;
         if(mcmc.saveconP) 
            for(k=com.ns; k<2*com.ns-1; k++) 
               com.oldconP[k] = 1;

         for(inode=com.ns; inode<tree.nnode; inode++) {
            t = nodes[inode].age;
            if (t>=taub[0] && t<taub[1]
            && (nodes[inode].ipop==is || sptree.nodes[nodes[inode].ipop].father==is)) {
               k = (t>=tauold && is!=sptree.root); /* k=0: below; 1: above */
               nodes[inode].age = t = taub[k]+taufactor[k]*(t-taub[k]);
               ntj_locus[k]++;
               if (t<taub[0] || (t<sptree.nodes[sptree.root].age && t>taub[1]))
                  error2("tj out of tb");

               if(mcmc.saveconP)
                  for(k=inode; ; k=nodes[k].father) {
                     com.oldconP[k] = 0; 
                     if(k == tree.root) break; 
                  }
            }
         }

         ntj[0] += ntj_locus[0];
         ntj[1] += ntj_locus[1];
         lnpGBinew[locus] = lnpGB_ThetaTau(locus);
         if(ntj_locus[0]+ntj_locus[1]) {
            lnpDinew[locus] = lnpD_locus(locus);
            lnLd += lnpDinew[locus] - data.lnpDi[locus];
         }
         else
            lnpDinew[locus] = data.lnpDi[locus];
         lnacceptance += lnpGBinew[locus] - data.lnpGBi[locus];

         if(debug==5) printf(" (%d %d) ",ntj_locus[0],ntj_locus[1]);

      }
      lnacceptance += lnLd + ntj[0]*log(taufactor[0]) + ntj[1]*log(taufactor[1]);

      if(debug==5) {
         printf(" ntj: %d %d", ntj[0],ntj[1]);
         printf(" lnLd %9.5f ", lnLd);
         for(k=com.ns; k<com.ns*2-1; k++)  printf("%d",com.oldconP[k]);
      }

      if(lnacceptance>=0 || rndu()<exp(lnacceptance)) {
         accepted++;
         for(locus=0; locus<data.ngene; locus++) {
            data.lnpDi[locus] = lnpDinew[locus];  
            data.lnpGBi[locus] = lnpGBinew[locus];
            UseLocus(locus, 0, 0, 0);

            for(inode=com.ns; inode<tree.nnode; inode++) { /* redo changes */
               t=nodes[inode].age;
               if (t>=taub[0] && t<taub[1]
               && (nodes[inode].ipop==is || sptree.nodes[nodes[inode].ipop].father==is)) {
                  k = (t>=tauold && is!=sptree.root);
                  nodes[inode].age = taub[k]+taufactor[k]*(t-taub[k]);
               }
            }
         }
         *lnL += lnLd;
         if(mcmc.usedata) SwitchconPin();
         if(debug==5) printf(" Y (%4d)\n", NPMat);
      }
      else {
         sptree.nodes[is].age = tauold;
         if(debug==5) printf(" N (%4d)\n", NPMat);
      }
   }
   return((double)accepted/ntau);
}


double CalculatetU4Splitting (int ispecies, char *flags)
{
/* This calculates tU, the upper bound for the new species age when node ispecies 
   in the species tree is split.  The routine is used by the SplitSpecies and 
   JoinSpecies moves.
   The allgorithm visits all gene trees in turn.  The two bits in flags for each node 
   in each gene tree are set if the node is ancestral to sons[0] and  sons[1].  
   If a node is anstral to both sons, its age is used to set the upper bound tU.  
   There is no need to flag the tips of a gene tree.
*/
   int locus, i,j, ninodes, *sons=sptree.nodes[ispecies].sons, dad=sptree.nodes[ispecies].father;
   char *f=flags;
   double tU = (dad==-1 ? OLDAGE : sptree.nodes[dad].age);

   for(i=0,ninodes=0; i<data.ngene; i++) ninodes += data.ns[i]-1;
   memset(f, '\0', ninodes*sizeof(char));

   for(locus=0; locus<data.ngene; f += data.ns[locus++]-1) {
      UseLocus (locus, 0, 0, 0);
      for(i=0; i<com.ns; i++) {
         /* If sons[0] is a tip in species tree, it may not have a theta and 
            pop_pop_table may be 0, so it is necessary below to check if(ipop==sons[0]).
         */
         if(nodes[i].ipop==sons[0] || sptree.pop_pop_table[nodes[i].ipop][sons[0]]) {
            for(j=nodes[i].father; f[j-com.ns]==0; j=nodes[j].father) {
               f[j-com.ns] = 1;  /* 1st bit */
               if(j==tree.root) break;
            }
         }
      }
      for(i=0; i<com.ns; i++) {
         if(nodes[i].ipop==sons[1] || sptree.pop_pop_table[nodes[i].ipop][sons[1]])
            for(j=nodes[i].father; f[j-com.ns]>>1==0; j=nodes[j].father) {
               f[j-com.ns] += 2;       /* 2nd bit */
               if(f[j-com.ns] == 3) {  /* if both bits are set */
                  if(tU > nodes[j].age) 
                     tU = nodes[j].age;
                  break;
               }
               if(j == tree.root) break;
            }
      }
   }
   return(tU);
}


int CountLHistories2 (void)
{
/* This counts the number of labeled histories for the rooted species tree.
   This is modified from CountLHistories().
*/
   int ns=sptree.nspecies, i,k, nLH, nLR[NSPECIES-1][2], change, *sons, j;
   double y=0;
   int debug=0;

   for(i=ns; i<sptree.nnode; i++) 
      nLR[i-ns][0] = nLR[i-ns][1] = -1;
   for(k=0; k<ns; k++) {
      if(debug) {
         for(i=0; i<sptree.nnode-ns; i++) 
            printf("\nnode %2d %d (%2d %2d): %2d %2d ", i+ns, sptree.nodes[i+ns].age>0, sptree.nodes[i+ns].sons[0], sptree.nodes[i+ns].sons[1], nLR[i][0], nLR[i][1]);
         FPN(F0);
      }
      for(i=ns,change=0; i<sptree.nnode; i++) {
         if(sptree.nodes[i].age == 0) continue;
         sons = sptree.nodes[i].sons;
         for(j=0; j<2; j++) {
            if(nLR[i-ns][j] != -1) continue;
            if(sons[j]<ns || sptree.nodes[sons[j]].age == 0) { /* sons[j] is tip */
               nLR[i-ns][j] = 0;
               change = 1;
            }
            else if(nLR[sons[j]-ns][0] != -1 && nLR[sons[j]-ns][1] != -1) {
               nLR[i-ns][j] = nLR[sons[j]-ns][0] + nLR[sons[j]-ns][1] + 1;
               change = 1;
            }
         }
      }
      if(!change) break;
   }
   for(i=0,nLH=1; i<sptree.nnode-ns; i++) {
      if(sptree.nodes[i+ns].age == 0) continue;
      if(debug)
         printf("\nnode %2d %d (%2d %2d): %2d %2d ", i+ns, sptree.nodes[i+ns].age>0, sptree.nodes[i+ns].sons[0], sptree.nodes[i+ns].sons[1], nLR[i][0], nLR[i][1]);
      
      if(nLR[i][0]==-1 || nLR[i][1]==-1)
         error2("nLR = -1");
      if(nLR[i][0] && nLR[i][1]) {
         nLH *= (int)Binomial((double)(nLR[i][0]+nLR[i][1]), nLR[i][0], &y);
         if(y) error2("y!=0 not expected");
      }
   }

   if(debug) {
      printf("\nnLH = %6d", nLH);
      exit(0);
   }
   return(nLH);   
}


#define lnPDFgamma(x, a, b)  ( (a)*log(b) - LnGamma(a) + ((a)-1)*log(x) - (b)*(x) )

double UpdateSplitSpecies (double *lnL, double space[], double PrSplit)
{
/* This is a rjMCMC move that splits a node in the guide species tree.  A node 
   is feasible if it is joined but its mother node is not joined.
   This move may change nodes[].ipop in gene trees but does not change the gene tree 
   topology or branch lengths or the likelihood.
   The old ipop is kept in ipop0[], and flags[] are used to determine the upper bound tU.   
   The number of elements in ipop0 and in flags is sum of (ns-1) over loci, the total 
   number of internal nodes (ninodes).
*/
   int nfeasible[2]={0}, feasibles[NSPECIES]={0}, locus, is, i,j, *sons,accepted=0;
   int i1,i2, ninodes, ns, ntau, nLHnew=1, ithetajk[2]={0};
   /* make sure e is the same in the split and join moves */
   double atheta=data.theta_prior[0], btheta=data.theta_prior[1];
   double atau=data.tau_prior[0], btau=data.tau_prior[1];
   double lnacceptance = log((1-PrSplit)/PrSplit), tnew, tU=1e300, *lnpGBinew=space, lnGd;
   double thetafactor=1, y, thetai, aRJ=mcmc.RJfinetune[0], mRJ=mcmc.RJfinetune[1];
   char *ipop0=(char*)(lnpGBinew+data.ngene), *ip, ipopchanged=0;  /* old ipop */
   char *flags, *f;

   if(debug==6) {
      printf("\nUpdateSplitSpecies\nSpecies tree\n");
      for(i=0; i<sptree.nnode; i++)
         printf("node %2d theta = %9.6f tau = %9.6f\n", i, sptree.nodes[i].theta, sptree.nodes[i].age);
   }

   for(i=0,ninodes=0; i<data.ngene; i++)
      ninodes += data.ns[i]-1;
   flags = ipop0 + ninodes;

   /* count nfeasible for splitting at source species tree */
   for(i=sptree.nspecies, ntau=0; i<sptree.nnode; i++) {
      if(sptree.nodes[i].age > 0)
         ntau++;
      else if(i==sptree.root || sptree.nodes[sptree.nodes[i].father].age>0)
         feasibles[nfeasible[0] ++] = i;
   }
   if(nfeasible[0]==0)
      return(0);

   /* split species is, generate theta_j & theta_k, if they can exist, and 
      calculate their prior.
   */
   is = feasibles[(int)(nfeasible[0]*rndu())];
   sons = sptree.nodes[is].sons;
   tU = CalculatetU4Splitting(is, flags);

   if(debug==6) {
      printf("Split species node %d (sons: %d %d)\n", is,sons[0],sons[1]);
      printf("tU = %.6f\nLeft & right flags to determine tU\n", tU);
      for(locus=0,f=flags; locus<data.ngene; f+=ns-1,locus++) {
         printf("locus %d\n", locus+1);
         UseLocus(locus, 0, 0, 1);
         printGtree(1);
         for(i=ns=data.ns[locus]; i<ns*2-1; i++)
            printf("node %2d: %d = %d %d\n", i, f[i-ns], f[i-ns]&1, f[i-ns]>>1);
      }
   }

   sptree.nodes[is].age = tnew = tU*pow(rndu(), 1/3.0);
   /* new theta_j & theta_k */
   thetai = sptree.nodes[is].theta;
   for(j=0; j<2; j++) 
      for(i=0; i<sptree.npop; i++) {
         if(sptree.pops[i] == sons[j]) {
            ithetajk[j] = 1;

            if(mcmc.RJalgorithm==0) {
               sptree.nodes[sons[j]].theta = thetai * exp(aRJ*(rndu()-0.5));
               thetafactor *= aRJ * sptree.nodes[sons[j]].theta;
            }
            else {
               sptree.nodes[sons[j]].theta = rndgamma(aRJ)/(aRJ/(mRJ*thetai));
               thetafactor /= exp( lnPDFgamma(sptree.nodes[sons[j]].theta, aRJ, aRJ/(mRJ*thetai)) );
            }
            /* prior on theta */
            lnacceptance += lnPDFgamma(sptree.nodes[sons[j]].theta, atheta, btheta);
            
            if(sptree.nodes[sons[j]].theta < 1e-10)
               printf("\nnew proposed theta = %20.9g\n", sptree.nodes[sons[j]].theta);
            break;
         }
      }

   /* prior on tau */
   if(is==sptree.root) { /* ntau changed from 0 to 1.  */
      y = lnPDFgamma(tnew, atau, btau);
      lnacceptance += y;
   }
   else {
      if(sptree.nLHistories)
         lnacceptance += log(sptree.nLHistories/(double)(nLHnew = CountLHistories2()));
      lnacceptance += log(ntau/sptree.nodes[sptree.root].age);
   }

   /* count nfeasible for joining at target.  This cannot be calculated before we know 
      which node to split (is).  A node is feasible for joining if both sons are either 
      tip or already joined.  
   */
   for(i=sptree.nspecies; i<sptree.nnode; i++) {
      if(sptree.nodes[i].age == 0) continue;
      i1 = sptree.nodes[i].sons[0];
      i2 = sptree.nodes[i].sons[1];
      if((i1<sptree.nspecies || sptree.nodes[i1].age==0)
      && (i2<sptree.nspecies || sptree.nodes[i2].age==0))
         nfeasible[1] ++;
   }

   y = tU*tU*tU/(3*square(sptree.nodes[is].age));
   lnacceptance += log(nfeasible[0]*y/nfeasible[1] * thetafactor);

   /* make a copy of ipop in gene trees. 
      change some node ipop from is to sons[0] or sons[1], using flags[]. */
   ip=ipop0; f=flags;
   for(locus=0,lnGd=0; locus<data.ngene; ip+=ns-1,f+=ns-1,locus++) {
      UseLocus (locus, 0, 0, 0);
      for(i=ns=data.ns[locus]; i<tree.nnode; i++) {
         ip[i-ns] = nodes[i].ipop;
         if(nodes[i].ipop!=is || nodes[i].age>tnew) continue;
         if(f[i-ns] == 3)   /* impossible for both bits to be set. */
            error2("strange here");
         if(f[i-ns] == 1) { nodes[i].ipop = sons[0]; ipopchanged=1; }
         if(f[i-ns] == 2) { nodes[i].ipop = sons[1]; ipopchanged=1; }
      }
      lnpGBinew[locus] = lnpGB_ThetaTau(locus);
      lnGd += lnpGBinew[locus] - data.lnpGBi[locus];
   }
   lnacceptance += lnGd;

   if(debug==6)
      printf("tnew = %.6f lnacceptance = %.6f\n", tnew, lnacceptance);

   if(diagnosis) {
      int i;
      double m, thetajk[2], bestm=-1, bestlnGd=-1e99;

      thetajk[0] = sptree.nodes[sons[0]].theta;
      thetajk[1] = sptree.nodes[sons[1]].theta;
      for(i=0,m=0.1; i<10; i++,m*=2) {
         if(ithetajk[0]) sptree.nodes[sons[0]].theta = thetai*m;
         if(ithetajk[1]) sptree.nodes[sons[1]].theta = thetai*m;
         for(locus=0,lnGd=0; locus<data.ngene; locus++) {
            UseLocus (locus, 0, 0, 0);
            lnGd += lnpGB_ThetaTau(locus) - data.lnpGBi[locus];
         }
         if(lnGd>bestlnGd) { bestm=m; bestlnGd=lnGd; }
      }
      sptree.nodes[sons[0]].theta = thetajk[0];
      sptree.nodes[sons[1]].theta = thetajk[1];      

      printf("S: ->qi ti qj qk: %.6f %.6f (%.6f)", sptree.nodes[is].theta, tnew, tU);
      for(i=0; i<2; i++) if(ithetajk[i]) printf(" %.6f", thetajk[i]);
      printf(" lnR %6.1f lnGd %6.1f (%6.1f mRJ %5.1f)\n", lnacceptance, lnGd, lnacceptance+bestlnGd-lnGd, bestm);
   }

   if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
      for(locus=0; locus<data.ngene; locus++) {
         data.lnpGBi[locus] = lnpGBinew[locus];
      }
      sptree.Itree |= 1<<(sptree.nspecies*2 - 2 - is);
      com.np += 1 + (sptree.nodes[sons[0]].theta>0) + (sptree.nodes[sons[1]].theta>0);
      sptree.nLHistories = nLHnew;
      accepted = 1;
   }
   else  {
      sptree.nodes[is].age = 0;
      sptree.nodes[sons[0]].theta = -1;
      sptree.nodes[sons[1]].theta = -1;

      if(ipopchanged)
         for(locus=0,ip=ipop0; locus<data.ngene; ip+=ns-1,locus++) {
            for(i=ns=data.ns[locus]; i<ns*2-1; i++)
               gnodes[locus][i].ipop = ip[i-ns];
         }
   }

   return(accepted);
}


double UpdateJoinSpecies (double *lnL, double space[], double PrSplit)
{
/* This is a rjMCMC move that joins a node in the guide species tree.  A node 
   is feasible for joining if its two daughter nodes are tips or collapsed nodes.
   This move does not change the likelihood.
   See notes in UpdateSplitSpecies().
*/
   int nfeasible[2]={0}, feasibles[NSPECIES]={0}, locus, is, i, *sons,accepted=0;
   int ninodes, ns, ntau, nLHnew=1;
   /* make sure e is the same in the split and join moves */
   double atheta=data.theta_prior[0], btheta=data.theta_prior[1];
   double atau=data.tau_prior[0], btau=data.tau_prior[1];
   double lnacceptance = log(PrSplit/(1-PrSplit)), y;
   double told, tU=1e300, *lnpGBinew=space;
   double thetafactor=1, thetai, thetajk0[2], aRJ=mcmc.RJfinetune[0], mRJ=mcmc.RJfinetune[1];
   char *ipop0=(char*)(lnpGBinew+data.ngene), *ip, ipopchanged=0;  /* old ipop */
   char *flags, *f;

   if(debug==7) puts("\nUpdateJoinSpecies ");
   for(i=0,ninodes=0; i<data.ngene; i++) 
      ninodes += data.ns[i]-1;
   flags = ipop0 + ninodes;

   /* count nfeasible for joining at source species tree.  A node is feasible 
      for joining if both sons are either tip or already joined. 
   */
   for(i=sptree.nspecies,ntau=0; i<sptree.nnode; i++) {
      if(sptree.nodes[i].age == 0) continue;
      ntau ++;
      sons = sptree.nodes[i].sons;
      if((sons[0]<sptree.nspecies || sptree.nodes[sons[0]].age==0)
      && (sons[1]<sptree.nspecies || sptree.nodes[sons[1]].age==0))
         feasibles[nfeasible[0] ++] = i;
   }
   if(nfeasible[0]==0) return(0);

   /* join species is */
   is = feasibles[(int)(nfeasible[0]*rndu())];
   sons = sptree.nodes[is].sons;

   if(debug==7) printf("Joining species node %d (sons: %d %d)\n", is,sons[0],sons[1]);

   /* get tU, for coming back though a split move.  flags are not useful later. */
   tU = CalculatetU4Splitting(is, flags);
   if(debug==7) {
      printf("tU = %.6f\nLeft & right flags to determine tU\n", tU);
      for(locus=0,f=flags; locus<data.ngene; f+=ns-1,locus++) {
         printf("locus %d\n", locus+1);
          
         UseLocus(locus, 0, 0, 1);
         printGtree(1); 

         for(i=ns=data.ns[locus]; i<ns*2-1; i++)
            printf("node %2d: %d = %d %d\n", i, f[i-ns], f[i-ns]&1, f[i-ns]>>1);
      }
   }

   /* tau_i */
   told = sptree.nodes[is].age;
   sptree.nodes[is].age = 0;
   /* theta_j & theta_k */
   thetai = sptree.nodes[is].theta;
   for(i=0; i<2; i++) {
      thetajk0[i] = sptree.nodes[sons[i]].theta;
      sptree.nodes[sons[i]].theta = -1; 
      if(thetajk0[i] > 0) {
         if(mcmc.RJalgorithm==0)
            thetafactor /= aRJ * thetajk0[i];
         else 
            thetafactor *= exp( lnPDFgamma(thetajk0[i], aRJ, aRJ/(mRJ*thetai)) );

         /* prior on old theta */         
         lnacceptance -= lnPDFgamma(thetajk0[i], atheta, btheta);
         
      }
   }
   if(is==sptree.root)  /* ntau changed from 1 to 0.  */
      lnacceptance -= lnPDFgamma(told, atau, btau);
   else {
      if(sptree.nLHistories)
         lnacceptance += log(sptree.nLHistories/(double)(nLHnew = CountLHistories2()));
      lnacceptance -= log((ntau-1)/sptree.nodes[sptree.root].age);
   }

   /* make a copy of ipop in gene trees. Change node ipop from sons[0] or sons[1] to is. */
   for(locus=0,ip=ipop0; locus<data.ngene; ip+=ns-1,locus++) {
      UseLocus (locus, 0, 0, 0);
      for(i=ns=com.ns; i<tree.nnode; i++) {
         ip[i-ns] = nodes[i].ipop;
         if(nodes[i].ipop==sons[0] || nodes[i].ipop==sons[1]) {
            nodes[i].ipop = is;
            ipopchanged = 1;
         }
      }
      lnpGBinew[locus] = lnpGB_ThetaTau(locus);
      lnacceptance += lnpGBinew[locus] - data.lnpGBi[locus];
   }

   /* count nfeasible for splitting at target.  A node is feasible if its father 
      is already split.  This is calculated after the species to join (is) is known.  
   */
   for(i=sptree.nspecies; i<sptree.nnode; i++) {
      if(sptree.nodes[i].age ==0 &&
         (i==sptree.root || sptree.nodes[sptree.nodes[i].father].age>0))
            nfeasible[1] ++;
   }

   y = tU*tU*tU/(3*told*told);
   lnacceptance += log(nfeasible[0]/(nfeasible[1]*y) * thetafactor);

   if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
      for(locus=0; locus<data.ngene; locus++) {
         data.lnpGBi[locus] = lnpGBinew[locus];
      }
      sptree.Itree &= ~(1<<(sptree.nspecies*2 - 2 - is));
      com.np -= 1 + (thetajk0[0]>0) + (thetajk0[1]>0);
      sptree.nLHistories = nLHnew;
      accepted = 1;
   }
   else  {
      sptree.nodes[is].age = told;
      for(i=0; i<2; i++)
         sptree.nodes[sons[i]].theta = thetajk0[i];
      if(ipopchanged)
         for(locus=0,ip=ipop0; locus<data.ngene; ip+=ns-1,locus++) {
            for(i=ns=data.ns[locus]; i<ns*2-1; i++)
               gnodes[locus][i].ipop = ip[i-ns];
         }
   }
   return(accepted);
}

double mixing (double* lnL, double finetune, double space[])
{
/* This multiplies all theta, tau, and branch lengths in all gene trees by c.
   This move can bring a node age in a gene tree to be older than OldAge, 
   which may not be nice.
*/
   int accepted=0, locus, i, k, ntheta=sptree.npop, ntau=sptree.nspecies-1;
   double xold,xnew, *lnpDinew=space, *lnpGBinew=lnpDinew+data.ngene;
   double c, lnc, lnacceptance=0, lnLd;
   double atheta=data.theta_prior[0], btheta=data.theta_prior[1];
   double atau=data.tau_prior[0], btau=data.tau_prior[1];

   if(finetune<=0) error2("steplength = 0 in mixing");
   if(sptree.speciesdelimitation) {
      for(i=sptree.nspecies,ntau=0; i<2*sptree.nspecies-1; i++)
         ntau += (sptree.nodes[i].age>0);
      for(i=0,ntheta=0; i<sptree.npop; i++) {
         if(sptree.nodes[sptree.pops[i]].theta > 0) 
            ntheta++;
      }
   }
   lnc = finetune*rnd2NormalSym(m2NormalSym);
   c = exp(lnc);
   for(locus=0,k=0; locus<data.ngene; locus++)
      k += data.ns[locus] - 1;
   lnacceptance = (ntheta+ntau+k)*lnc;

   for(i=0; i<sptree.npop; i++) {
      if(sptree.nodes[sptree.pops[i]].theta <= 0) continue;
      xold = sptree.nodes[sptree.pops[i]].theta;
      sptree.nodes[sptree.pops[i]].theta = xnew = xold*c;
      lnacceptance += (atheta-1)*lnc - btheta*(xnew-xold);
   }
   for(i=0; i<sptree.nspecies-1; i++) {
      if(sptree.nodes[sptree.nspecies+i].age == 0) continue;
      xold = sptree.nodes[sptree.nspecies+i].age;
      sptree.nodes[sptree.nspecies+i].age = xnew = xold*c;
      if(i==0)  /* root in species tree */
         lnacceptance += (atau-1 - ntau+1)*log(xnew/xold) - btau*(xnew-xold);
   }

   if(mcmc.saveconP) 
      for(i=0; i<data.maxns*2-1; i++) com.oldconP[i] = 0;
   for(locus=0,lnLd=0; locus<data.ngene; locus++) {
      UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & conP */
      for(i=com.ns; i<tree.nnode; i++)  nodes[i].age *= c;
      lnpGBinew[locus] = lnpGB_ThetaTau(locus);
      lnacceptance += lnpGBinew[locus] - data.lnpGBi[locus];

      lnpDinew[locus] = lnpD_locus(locus);
      lnLd += lnpDinew[locus] - data.lnpDi[locus];
      lnacceptance += lnpDinew[locus] - data.lnpDi[locus];
   }

   if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
      for(locus=0; locus<data.ngene; locus++) {
         UseLocus(locus, 0, 1, 0);
         for(i=com.ns; i<tree.nnode; i++)  
            nodes[i].age *= c;
         data.lnpDi[locus]  = lnpDinew[locus];
         data.lnpGBi[locus] = lnpGBinew[locus];
      }
      if(mcmc.usedata) SwitchconPin();
      *lnL += lnLd;
      accepted = 1;
   }
   else  {
      for(i=0; i<sptree.npop; i++)
         if(sptree.nodes[sptree.pops[i]].theta > 0) 
            sptree.nodes[sptree.pops[i]].theta /= c;
      for(i=sptree.nspecies; i<sptree.nspecies*2-1; i++)
         if(sptree.nodes[i].age > 0)
            sptree.nodes[i].age /= c;
   }
   return(accepted);
}


double UpdateLocusrateHeredity (double* lnL, double finetune)
{
/* This updates locus-specific rates and heredity multipliers.
   There is probably no need to update both when both are estimated for the 
   same locus, as the posterior for rates and heredity multipliers should be 
   quite flat.
*/
   int  accepted=0, locus, j, locusref=0;
   double lnacceptance, lnLd, lnpDinew, lnpGBinew;
   double h,hnew,r,rnew,rref,rrefnew, lnpDrefnew;

   if(debug==2) puts("\nUpdateLocusrateHeredity ");
   if(finetune<=0) error2("steplength = 0 in UpdateLocusrateHeredity");

   if(mcmc.saveconP) FOR(j,data.maxns*2-1) com.oldconP[j]=0;

   if(data.est_locusrate == 1) {
      for(j=1,locusref=0; j<data.ngene; j++)
         if(data.npatt[j] > data.npatt[locusref]) locusref = j;

      for(locus=0; locus<data.ngene; locus++) {
         if(locus == locusref) continue;
         r = data.locusrate[locus];
         rnew = r + finetune*rnd2NormalSym(m2NormalSym);
         rnew = data.locusrate[locus] = reflect(rnew, 0, r+data.locusrate[locusref]);
         rref = data.locusrate[locusref];
         rrefnew = data.locusrate[locusref] -= rnew-r;
   
         lnacceptance = (data.a_locusrate-1)*log((rnew*rrefnew)/(r*rref));
         /* lnacceptance = (data.a_locusrate-1)*log(rnew/r); */

         UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree & conP */
         lnpDinew = lnpD_locus(locus);
         UseLocus(locusref, 1, mcmc.usedata, 0);
         lnpDrefnew = lnpD_locus(locusref);

         lnLd = lnpDinew-data.lnpDi[locus] + lnpDrefnew-data.lnpDi[locusref];
         lnacceptance += lnLd;
         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
            accepted ++;
            if(mcmc.usedata) {
               *lnL += lnLd;
               data.lnpDi[locus] = lnpDinew;
               data.lnpDi[locusref] = lnpDrefnew;
               CopyconPin(locus, 0);
               CopyconPin(locusref, 0);
            }
         }
         else {
            data.locusrate[locus] = r;
            data.locusrate[locusref] = rref;
         }
      }
   }

   if(data.est_heredity == 1) { /* this step does not change likelihood */
      for(locus=0; locus<data.ngene; locus++) {
         UseLocus(locus, 1, mcmc.usedata, 0);  /* copy gtree but not conP */
         h = data.heredity[locus];
         hnew = h + finetune*rnd2NormalSym(m2NormalSym);
         if(hnew<0) hnew *= -1;
         data.heredity[locus] = hnew;
         lnacceptance = (data.a_heredity-1)*log(hnew/h) - data.b_heredity*(hnew-h);

         lnpGBinew = lnpGB_ThetaTau (locus);
         lnacceptance += lnpGBinew - data.lnpGBi[locus];
      
         if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
            accepted ++;
            data.lnpGBi[locus] = lnpGBinew;
         }
         else 
            data.heredity[locus] = h;
      }
   }
   return((double)accepted/(data.ngene*2-1.0));
}

double UpdateSequenceErrors (double* lnL, double finetune, double space[])
{
/* This updates the parameters for sequencing errors
*/
   int  accepted=0, locus, is, i, j;
   double lnacceptance, lnLd, *lnpDinew=space, eold, enew, eDold, eDnew, c, lnc, *a;

   if(finetune<=0) error2("steplength = 0 in UpdateSequenceErrors");
   for(i=com.ns; i<com.ns*2-1; i++) 
      com.oldconP[i] = 0;

   for(is=0; is<sptree.nspecies; is++) {
      if(data.iseqerr[is] == 0) continue;
      for(i=0; i<4; i++) {
         a = data.a_seqerr[is] + i*4;
         for(j=0; j<4; j++) {
            if(i==j) continue;
            eold  = data.e_seqerr[is][i*4+j];
            eDold = data.e_seqerr[is][i*4+i];
            enew  = eold + finetune*rnd2NormalSym(m2NormalSym);
            enew  = data.e_seqerr[is][i*4+j] = reflect(enew, 1e-20, eDold + eold);
            eDnew = data.e_seqerr[is][i*4+i] = eDold + eold - enew;

            lnacceptance = (a[j]-1)*log(enew/eold) + (a[i]-1)*log(eDnew/eDold);

            for(locus=0,lnLd=0; locus<data.ngene; locus++) {
               UseLocus(locus, 1, mcmc.usedata, 0);   /* copy gtree & conP */
               lnpDinew[locus] = lnpD_locus(locus);
               lnLd += lnpDinew[locus] - data.lnpDi[locus];
               lnacceptance += lnpDinew[locus] - data.lnpDi[locus];
            }

            if(lnacceptance>=0 || rndu()<exp(lnacceptance)) { /* accept */
               accepted ++;
               for(locus=0; locus<data.ngene; locus++)
                  data.lnpDi[locus] = lnpDinew[locus];
               if(mcmc.usedata) SwitchconPin();
               *lnL += lnLd;
            }
            else {
               data.e_seqerr[is][i*4+j] = eold;
               data.e_seqerr[is][i*4+i] = eDold;
            }
         }
      }
   }
   return((double)accepted/(data.nseqerr*12.0));
}


void checkGtree (void)
{
   int locus, inode,j, ipop, father;
   double t, tb[2];
   /* double small=1e-12, verysmall=1e-17; */
   double small=0, verysmall=0;

   for(locus=0; locus<data.ngene; locus++) {
      UseLocus(locus, 0, 0, 1);
      for(inode=com.ns; inode<tree.nnode; inode++) {
         t = nodes[inode].age;
         father = nodes[inode].father;
         ipop = nodes[inode].ipop;
         tb[0] = sptree.nodes[ipop].age;  
         tb[1] = OLDAGE;
         for(j=0; j<nodes[inode].nson; j++) 
            tb[0] = max2(tb[0], nodes[nodes[inode].sons[j]].age);
         if(ipop != sptree.root) 
            tb[1] = sptree.nodes[sptree.nodes[ipop].father].age;
         if(inode != tree.root) 
            tb[1] = min2(tb[1], nodes[nodes[inode].father].age);

         if(tb[1]<tb[0] || t<tb[0]-small || t>tb[1]+small) {
            printf("\nlocus %d node %d t=%9.6f tb: (%9.6f %9.6f)", locus,inode,t,tb[0],tb[1]);
            printf("\nspecies tree and gene tree:\n");
            for(j=0; j<2*sptree.nspecies-1; j++,FPN(F0)) {
               printf("species %d %-12s ", j+1,sptree.nodes[j].name);
               printf("age %10.6g  theta %10.6g  ", sptree.nodes[j].age, sptree.nodes[j].theta);
            }
            printGtree(1);
            printf("\ninconsistent gene tree at locus %d node %d\ttime %.5g (%.5g %.5g)\n",
               locus,inode,t,tb[0],tb[1]);
            puts("Enter to continue. ");
            getchar();
         }
         else if(t<tb[0] + verysmall)
            nodes[inode].age = tb[0] + verysmall;
         else if(t>tb[1] - verysmall)
            nodes[inode].age = tb[1] - verysmall;
      }
   }
}

int DownSptreeSetTime (int inode)
{
/* This moves down species tree to specify node ages, starting from root age.
   It sets the theta to -1 if the tip species is collapsed.
*/
   int j, ison;

   for (j=0; j<sptree.nodes[inode].nson; j++) {
      ison = sptree.nodes[inode].sons[j];
      if(sptree.nodes[inode].age == 0) {   /* inode collapsed */
         sptree.nodes[ison].theta = -1;
      }
      if(sptree.nodes[ison].nson) {        /* ison is an interior node */
         if(sptree.nodes[inode].age == 0)
            sptree.nodes[ison].age = 0;  

         if(sptree.nodes[ison].age > 0)
            sptree.nodes[ison].age = sptree.nodes[inode].age *(0.1+0.9*rndu());

         DownSptreeSetTime(ison);
      }
   }
   return(0);
}

int GetInitials (void)
{
/* This sets the initial values for starting the MCMC, and returns np, the 
   number of parameters in the MCMC, to be collected in collectx().

   If (data.est_locusrate == 1), the rate at locus 1 is printed out.  
   See also collectx().
   This is right now not correct, as ntau and ntheta are not correct.
*/
   int np, is, i,j,k, ntheta=sptree.npop, ntau=sptree.nspecies-1;
   double atau=data.tau_prior[0], btau=data.tau_prior[1];
   double atheta=data.theta_prior[0], btheta=data.theta_prior[1], y;

   /* initial theta's.  Some may be set to -1 later. */
   for(i=0,k=0; i<sptree.npop; i++) {
      sptree.nodes[sptree.pops[i]].theta = atheta/btheta*(0.6+0.8*rndu());
   }
   /* initial tau's.  Some theta's are reset to -1. */
   if(sptree.nspecies > 1) {
      /* initially age is used to indicate the presence of the tau parameter. */
      for(i=sptree.nspecies; i<sptree.nnode; i++) 
         sptree.nodes[i].age = 1; 
      sptree.nodes[sptree.nspecies].age =  atau/btau*(0.8+0.4*rndu());
      if(sptree.speciesdelimitation) {
         is = sptree.nspecies + (int)(sptree.nspecies*rndu());  /* 1/s chance for no collapse */

         /*
         is = sptree.nspecies + 1;
         is = 0;
         */

         /*
         printf("\nNode to collaps (a number between %d and %d, also %d)? ", sptree.nspecies+1, 2*sptree.nspecies-1, 2*sptree.nspecies);
         scanf("%d", &is);
         is--;
         */
         /* collaps node is in initial tree, possible not to collapse. */
         if(is>=sptree.nspecies && is<2*sptree.nspecies-1)
            sptree.nodes[is].age = 0;
      }

      DownSptreeSetTime (sptree.nspecies);
   }

   if(sptree.speciesdelimitation) {
      sptree.Itree = 0;
      sptree.nLHistories = CountLHistories2();
      for(is=sptree.nspecies; is<sptree.nspecies*2-1; is++)
         if(sptree.nodes[is].age>0)
            sptree.Itree |= (1 << (sptree.nspecies*2 - 2 - is));
      FPN(F0);
      printf("\nStarting species tree = %s\n", printSpItree(sptree.Itree));
      for(i=0, ntheta=ntau=0; i<sptree.nnode; i++) {
         printf("node %2d %-20s ", i+1, sptree.nodes[i].name);
         printf("theta = %9.6f  tau = %9.6f\n", sptree.nodes[i].theta, sptree.nodes[i].age);
         ntheta += (sptree.nodes[i].theta>0);
         ntau   += (sptree.nodes[i].age>0);
      }
   }
   if(data.est_heredity == 1) {
      for(i=0; i<data.ngene; i++)
         data.heredity[i] = data.a_heredity/data.b_heredity*(0.8+0.4*rndu());
   }
   if(data.est_locusrate == 1) {
      for(i=0; i<data.ngene; i++) {
         data.locusrate[i] = 0.8+0.4*rndu();
      }
      y = sum(data.locusrate,data.ngene)/data.ngene;
      for(i=0; i<data.ngene; i++) 
         data.locusrate[i] /= y;
   }

   /* sequence errors */
   if(data.nseqerr) {
      printf("\nInitials for seqerrors");

      for(is=0; is<sptree.nspecies; is++) {
         if(data.iseqerr[is]) {
            for(i=0; i<4; i++) {
               for(j=0,y=0; j<4; j++) y += data.a_seqerr[is][i*4+j];
               for(j=0; j<4; j++)
                  data.e_seqerr[is][i*4+j] = data.a_seqerr[is][i*4+j]/y * (0.8+0.4*rndu());
               for(j=0,y=0; j<4; j++) y += data.e_seqerr[is][i*4+j];
               for(j=0; j<4; j++) data.e_seqerr[is][i*4+j] /= y;
            }
            matout(F0, data.e_seqerr[is], 4, 4);
         }
      }
   }

   np = ntheta + ntau + data.nseqerr*16
      + (data.est_locusrate == 1)
      + (data.est_heredity == 1)*min2(10, data.ngene);
   if(sptree.nspecies == 1) np += data.ngene;   /* t_MRCA */
   return(np);
}


int collectx (FILE* fout, double x[])
{
/* This collects parameters into x[] for printing and summarizing.
   It checks the number of parameters.
   If(fout), it prints the header line into fout.  Otherwise it does not print.
   if(com.np == -1) com.np is assigned a value here.
*/
   int printr1=(data.est_locusrate==1), i,j,k=0, is, ipop;

   if(fout) fprintf(fout, "Gen");
   if(fout && sptree.speciesdelimitation) fprintf(fout, "\tnp\ttree");

   for(i=0; i<sptree.npop; i++) {
      ipop = sptree.pops[i];

      /* Ziheng 27/04/2011.  Changed from "if(sptree.nodes[ipop].theta <= 0)". 
      */
      if(sptree.nodes[ipop].theta == -1)
         continue;

      if(fout) fprintf(fout, "\ttheta_%s", sptree.nodes[ipop].name);
      x[k++] = sptree.nodes[ipop].theta;
   }
   for(i=sptree.nspecies; i<sptree.nnode; i++) {
      if(sptree.nodes[i].age==0)  continue;
      if(fout) fprintf(fout, "\ttau_%s", sptree.nodes[i].name);
      x[k++] = sptree.nodes[i].age;
   }
   if(data.est_heredity == 1) {
      for(i=0; i<min2(10, data.ngene); i++) {
         if(fout) fprintf(fout, "\theredity_L%d", i+1);
         x[k++] = data.heredity[i];
      }
   }
   if(printr1) {
      if(fout) fprintf(fout, "\tr_L1");
      x[k++] = data.locusrate[0];
   }
   if(data.nseqerr) {
      for(is=0; is<sptree.nspecies; is++) 
         if(data.iseqerr[is]) {
            for(i=0; i<4; i++)  for(j=0; j<4; j++) {
               if(fout) fprintf(fout, "\tseqerr_[%s]%c%c", sptree.nodes[is].name, BASEs[i], BASEs[j]);
               x[k++] = data.e_seqerr[is][i*4+j];
            }
         }
   }

   if(sptree.nspecies == 1) /* t_MRCA */
      for(i=0; i<data.ngene; i++) {
         if(fout) fprintf(fout, "\tt_MRCA_L%d", i+1);
         x[k++] = gnodes[i][data.root[i]].age;
      }
   if(k != com.np) {
      if(com.np!=-1) 
         printf("np %d != %d,  np reset in collectx().", k, com.np);
      com.np = k;
   }
   if(fout && mcmc.usedata) fprintf(fout, "\tlnL");
   if(fout) fprintf(fout, "\n");
   return(0);
}


char *printSpItree(int itree)
{
   static char treestr[NSPECIES];
   int i=0, bit=sptree.nspecies-2;

   for(; i<sptree.nspecies-1; i++,bit--)
      treestr[i] = (itree & (1 << bit)) ? '1' : '0';
   treestr[i] = '\0';
   return(treestr);
}
 
void copySptree (void)
{
/* This copies sptree into nodes = nodes_t, for printing or editing
*/
   int i,j;

   nodes = nodes_t;
   com.ns = sptree.nspecies;   tree.root = sptree.root;
   tree.nnode = sptree.nnode;  tree.nbranch = sptree.nbranch; 
   for(i=0; i<sptree.nnode; i++) {
      if(i<com.ns) com.spname[i] = sptree.nodes[i].name;
      nodes[i].father  =sptree.nodes[i].father;
      nodes[i].nson = sptree.nodes[i].nson;
      for(j=0;j<nodes[i].nson;j++) 
         nodes[i].sons[j] = sptree.nodes[i].sons[j];
      nodes[i].age = sptree.nodes[i].age;
   }
}


int MCMC (FILE* fout)
{
   FILE *fmcmc=gfopen(com.mcmcf,"w"), *fmcmctmp;
   int nsteps = 5+(data.nseqerr>0);
   char BtreeStr[NSPECIES], line[16000], *mcmctmp="mcmc.tmp";
   int locus, j,k, ir, Btree=0, bit, lline=16000;
   double pBtree=0;
   double *x, *mx, *freqtree, lnL, PrSplit=0.5, postnodes[NSPECIES-1]={0};
   double Pmodel=0, Pjump[7]={0}, nround=0;
   double PtauThreshold[NSPECIES][2]={{0}}, tauThreshold[2]={2E-5, 2E-4};

   noisy=3;

   if(sptree.speciesdelimitation) {
      k = (1<<(sptree.nspecies-1));
      if((freqtree=(double*)malloc(k*sizeof(double))) == NULL)
         error2("oom for freqtree");
      memset(freqtree, 0, k*sizeof(double));
   }
   mcmc.moveinnode = 1;     /* moves internal nodes in the gene tree as well */
   mcmc.saveconP = 1;
   if(!mcmc.usedata) mcmc.saveconP = 0;

   printf("\n%d burnin, sampled every %d, %d samples\n", 
           mcmc.burnin, mcmc.sampfreq, mcmc.nsample);
   if(mcmc.usedata) puts("Approximating posterior, using sequence data");
   else             puts("Approximating prior, not using sequence data");

   printf("(Settings: cleandata=%d print=%d saveconP=%d moveinnode=%d)\n",
          com.cleandata, mcmc.print, mcmc.saveconP, mcmc.moveinnode);

   printf("\nStarting MCMC... ");
   com.np = GetInitials();

   printf("\nprior theta ~ G(%.3f, %.3f)\n", data.theta_prior[0], data.theta_prior[1]);
   if(sptree.nspecies>1)
      printf("prior tau   ~ G(%.3f, %.3f)\n", data.tau_prior[0], data.tau_prior[1]);

   k = max2(4 + data.ngene, sptree.nspecies*3+3);
   k += data.nseqerr*16;
   x = (double*)malloc(k*2*sizeof(double));
   mx = x + k;
   if(x==NULL) error2("oom");
   zero(mx, k);

   /* initialize likelihood and prior, for each locus */
   for(locus=0,lnL=0; locus<data.ngene; locus++) {
      GetRandomGtree(locus);
      data.root[locus] = tree.root;
      if(debug) {
         printf("\nGene tree %d", locus+1);
         UseLocus(locus, 0, 0, 1);
         printGtree(1);
      }
      data.lnpGBi[locus] = lnpGB_ThetaTau(locus);
   }

   collectx(fmcmc, x);

   printf("\nInitial parameters, np = %d (gene trees generated from the prior):\n", com.np);
   for(j=0; j<com.np; j++) printf("%9.5f",x[j]);
   for(j=0; j<data.maxns*2-1; j++) com.oldconP[j]=0;
   lnL = lnpData(data.lnpDi);
   printf("\nlnL0 =%12.3f\n", lnL);

   for(ir=-mcmc.burnin,nround=0; ir<mcmc.sampfreq*mcmc.nsample; ir++) {
      if(ir==0 || (mcmc.resetFinetune && nround>=100 && mcmc.burnin>=200
               && ir<0 && ir%(mcmc.burnin/4)==0)) {
         /* reset finetune parameters.  Do this twice. */
         if(mcmc.resetFinetune && mcmc.burnin>=200) {
            ResetFinetuneSteps(fout, Pjump, mcmc.finetune, nsteps);
         }
         nround=0;  Pmodel=0;
         zero(Pjump, nsteps);
         zero(mx, com.np); 
      }
      nround++;

      if(sptree.speciesdelimitation) {
         if(rndu()<PrSplit)
            Pmodel  += UpdateSplitSpecies(&lnL, com.space, PrSplit);
         else
            Pmodel  += UpdateJoinSpecies(&lnL, com.space, PrSplit);
      }

      Pjump[0]  = (Pjump[0]*(nround-1) + UpdateGB_InternalNode(&lnL, mcmc.finetune[0]))/nround;
      Pjump[1]  = (Pjump[1]*(nround-1) + UpdateGB_SPR(&lnL, mcmc.finetune[1]))/nround;
      Pjump[2]  = (Pjump[2]*(nround-1) + UpdateTheta(mcmc.finetune[2], com.space))/nround;
      if(sptree.nspecies>1 && sptree.nodes[sptree.root].age>0)
         Pjump[3]  = (Pjump[3]*(nround-1) + UpdateTau(&lnL, mcmc.finetune[3], com.space))/nround;
      
      Pjump[4]  = (Pjump[4]*(nround-1) + mixing(&lnL, mcmc.finetune[4], com.space))/nround;
      Pjump[5]  = (Pjump[5]*(nround-1) + UpdateLocusrateHeredity (&lnL, mcmc.finetune[5]))/nround;
      if(data.nseqerr)
         Pjump[6]  = (Pjump[6]*(nround-1) + UpdateSequenceErrors (&lnL, mcmc.finetune[6], com.space))/nround;

      checkGtree();

      collectx(NULL, x);
      for(j=0; j<(sptree.speciesdelimitation ? 1 : com.np); j++) 
         mx[j] = (mx[j]*(nround-1) + x[j])/nround;
      if(ir>=0) {
         if(sptree.speciesdelimitation)
            freqtree[sptree.Itree]++;
         else 
            for(j=sptree.nspecies; j<sptree.nspecies*2-1; j++) {
               if(sptree.nodes[j].age<tauThreshold[0])  PtauThreshold[j-sptree.nspecies][0] ++;
               if(sptree.nodes[j].age<tauThreshold[1])  PtauThreshold[j-sptree.nspecies][1] ++;
            }
      }
      if(mcmc.print && ir>=0 && (ir==0 || (ir+1)%mcmc.sampfreq==0)) {
         fprintf(fmcmc,"%d", ir+1);

         if(sptree.speciesdelimitation)
            fprintf(fmcmc, "\t%d\t %s ", com.np, printSpItree(sptree.Itree));
         for(j=0;j<com.np; j++)    fprintf(fmcmc,"\t%.6f", x[j]);
         if(mcmc.usedata) fprintf(fmcmc, "\t%.3f", lnL);
         fprintf(fmcmc, "\n");
      }
      if(/* noisy && */ (ir+1)%max2(mcmc.sampfreq, mcmc.sampfreq*mcmc.nsample/1000)==0) {
         printf("\r%3.0f%%", (ir+1.)/(mcmc.nsample*mcmc.sampfreq)*100.);
         for(j=0; j<nsteps; j++) 
            printf(" %4.2f", Pjump[j]); 
         printf(" ");

         if(sptree.speciesdelimitation)
            printf(" %2d %6.4f %s ", com.np, Pmodel/nround, printSpItree(sptree.Itree));
         for(j=0; j<(sptree.speciesdelimitation ? 1 : min2(com.np,8)); j++) {
            printf(" %7.5f", mx[j]);
			   if(j>0 && j==sptree.npop-1) printf(" ");
         }

         if(mcmc.usedata) printf(" %4.2f", lnL);

         if(mcmc.sampfreq*mcmc.nsample>=50 && (ir+1)%(mcmc.sampfreq*mcmc.nsample/20)==0) {
            printf(" %s\n", printtime(timestr));
            testlnL=1;
            if(fabs(lnL - lnpData(data.lnpDi)) > 0.01) {
               printf("lnL not right: %12.6f != %12.6f ", lnL, lnpData(data.lnpDi));
            }
            testlnL=0;
            if(mcmc.print) fflush(fmcmc);
         }

         if(diagnosis) getchar();
      }
   }  /* for(ir) */

   fclose(fmcmc);
   free(data.lnpGBi);
   free(data.lnpDi);
   free(data.Imap);

   if(data.est_heredity || data.est_locusrate) 
      free(data.heredity);
   printf("\nTime used: %s", printtime(timestr));

   if(sptree.speciesdelimitation==0) {
      FPN(F0);  FPN(fout);
      for(j=sptree.nspecies; j<sptree.nspecies*2-1; j++) {
         printf("spnode %2d: P(tau <%7.2g) = %7.5f  P(tau <%7.2g) = %7.5f (%s)\n", j+1, 
            tauThreshold[0], PtauThreshold[j-sptree.nspecies][0]/nround,
            tauThreshold[1], PtauThreshold[j-sptree.nspecies][1]/nround, sptree.nodes[j].name);
         fprintf(fout, "spnode %2d: P(tau <%7.2g) = %7.5f  P(tau <%7.2g) = %7.5f (%s)\n", j+1, 
            tauThreshold[0], PtauThreshold[j-sptree.nspecies][0]/nround,
            tauThreshold[1], PtauThreshold[j-sptree.nspecies][1]/nround, sptree.nodes[j].name);
      }
      FPN(F0);  FPN(fout);

      if(mcmc.print) {
         printf("\nSummarizing, time reset.");
         fprintf(fout,"\n\nSummary of MCMC results\n");
         DescriptiveStatisticsSimple(fout, com.mcmcf, 50, 100, 1);
      }

   }   
   else {  /* sptree.speciesdelimitation==1 */
      printf("\nPrSplit = %.6f\n", PrSplit);
      printf("rj algorithm %d: new theta from ", mcmc.RJalgorithm);
      if(mcmc.RJalgorithm)  printf("G(a=%.2f, m=%.2f)\n", mcmc.RJfinetune[0], mcmc.RJfinetune[1]);
      else                  printf("sliding window with c = %.2f\n", mcmc.RJfinetune[0]);

      printf("\nTree frequencies (Ancestral nodes in order: ");
      for(k=sptree.nspecies; k<sptree.nspecies*2-1; k++)
         printf(" %2d %s", k+1, sptree.nodes[k].name);
      printf(")\n");
      fprintf(fout, "\nTree frequencies (Ancestral nodes in order: ");
      for(k=sptree.nspecies; k<sptree.nspecies*2-1; k++)
         fprintf(fout, " %2d %s", k+1, sptree.nodes[k].name);
      fprintf(fout, ")\n");
      for(j=0; j < (1<<(sptree.nspecies-1)); j++) {
         if(freqtree[j]) {
            printf(" %4d  %s %8.0f %9.6f\n", j, printSpItree(j), freqtree[j], freqtree[j]/nround);
            fprintf(fout, " %4d  %s %8.0f %9.6f\n", j, printSpItree(j), freqtree[j], freqtree[j]/nround);
            if(freqtree[j]/nround > pBtree) {
               pBtree = freqtree[j]/nround;
               Btree = j;
            }
            for(k=sptree.nspecies-1-1; k>=0; k--)
               if(j & (1 << k)) postnodes[sptree.nspecies-1-1-k] += freqtree[j]/nround;
         }
      }
      free(freqtree);

      copySptree();
      for(k=0; k<sptree.nspecies*2-1; k++)
         nodes[k].label = (k<sptree.nspecies ? 0 : postnodes[k-sptree.nspecies]);
      printf("\nGuide tree with posterior probability for presence of nodes\n");
      OutTreeN(F0, 1, PrLabel);   FPN(F0);
      fprintf(fout, "\nGuide tree with posterior probability for presence of nodes\n");
      OutTreeN(fout, 1, PrLabel); FPN(fout);

      if(mcmc.print) {
         printf("\nSummarizing parameters for the MAP tree %s\n", printSpItree(Btree));
         fprintf(fout, "\nSummarizing the posterior of parameters under the MAP tree %s\n", printSpItree(Btree));
         fmcmctmp = gfopen(mcmctmp, "w");
         com.np = -1;  /* com.np is for the last tree not the MAP tree.  Reset in collectx(). */
         for(j=0,bit=sptree.nspecies-2; j<sptree.nspecies-1; j++,bit--)
            sptree.nodes[sptree.nspecies+j].age = (Btree & (1 << bit)) ? 1 : 0;
         for(j=0; j<sptree.nspecies*2-1; j++)
            sptree.nodes[sptree.pops[j]].theta = 1;
         DownSptreeSetTime (sptree.nspecies);   /* this resets theta to -1 for collapsed tips */
         collectx(fmcmctmp, x);
         sprintf(BtreeStr, " %s ", printSpItree(Btree));

         fmcmc = gfopen(com.mcmcf,"r");
         for(j=0; j<mcmc.nsample+2; j++) {
            fgets(line, lline, fmcmc);
            if(j && strstr(line, BtreeStr))
               fprintf(fmcmctmp, "%s", line);
         }
         fclose(fmcmc);
         fclose(fmcmctmp);
         DescriptiveStatisticsSimple(fout, mcmctmp, 50, 100, 3);
      }
   }
   return(0);
}

#endif
