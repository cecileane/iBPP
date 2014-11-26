Notes by Ziheng Yang
30 May 2010

(1) To compile, try one of the following

   UNIX gcc/icc:
      cc -o bpp -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION -m64 -march=opteron -mtune=opteron -ansi -O3 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm

      icc -o bpp -fast bpp.c tools.c -lm

   MAC OSX intel:
      cc -o bpp -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm
      cc -o MCcoal -DSIMULATION -O4 -funroll-loops -fomit-frame-pointer -finline-functions bpp.c tools.c -lm

   Windows MSC++ 6.0/2008:
      cl -O2 bpp.c tools.c
      cl -O2 -FeMCcoal.exe -DSIMULATION bpp.c tools.c

(2) To run an example analysis, try 

cd examples

../bpp yu2001.bpp.ctl

../bpp ChenLi2001.bpp.ctl

../bpp lizard.bpp.ctl


Good luck.
