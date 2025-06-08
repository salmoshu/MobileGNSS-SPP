# makefile for rnx2rtkp
BINDIR  = /usr/local/bin
RTKLIB_SRC    = ../../../rtklib_src
USER_SRC = ../../../user_src
RTKLIB_INC = ../../../rtklib_inc
USER_INC = ../../../user_inc

OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN -DNFREQ=3 -DNEXOBS=3
#OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN -DNFREQ=5 -DIERS_MODEL
#OPTS    = -DENAGLO -DENAQZS -DENAGAL -DENACMP -DNFREQ=2

# for no lapack
CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(RTKLIB_SRC) $(OPTS) -std=c99
LDLIBS  = -lm
#LDLIBS  = -lgfortran -lm -lrt -pg

#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(RTKLIB_SRC) -DLAPACK $(OPTS)
#LDLIBS  = -lm -lrt -llapack -lblas

# for gprof
#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(RTKLIB_SRC) -DLAPACK $(OPTS) -pg
#LDLIBS  = -lm -lrt -llapack -lblas -pg

# for mkl
##MKLDIR  = /opt/intel/mkl
#MKLDIR  = /proj/madoca/lib/mkl
#CFLAGS  = -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(RTKLIB_SRC) $(OPTS) -DMKL
#LDLIBS  = -L$(MKLDIR)/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lgomp -lm -lrt

all        : rnx2rtkp
rnx2rtkp   : rnx2rtkp.o rtkcmn.o rinex.o rtkpos.o postpos.o solution.o
rnx2rtkp   : lambda.o geoid.o sbas.o preceph.o pntpos.o ephemeris.o options.o
rnx2rtkp   : ppp.o ppp_ar.o rtcm.o rtcm2.o rtcm3.o rtcm3e.o ionex.o tides.o trace.o

rnx2rtkp.o : $(USER_SRC)/rnx2rtkp.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/rnx2rtkp.c
rtkcmn.o   : $(RTKLIB_SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtkcmn.c
rinex.o    : $(RTKLIB_SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rinex.c
rtkpos.o   : $(RTKLIB_SRC)/rtkpos.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtkpos.c
postpos.o  : $(RTKLIB_SRC)/postpos.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/postpos.c
solution.o : $(RTKLIB_SRC)/solution.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/solution.c
lambda.o   : $(RTKLIB_SRC)/lambda.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/lambda.c
geoid.o    : $(RTKLIB_SRC)/geoid.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/geoid.c
sbas.o     : $(RTKLIB_SRC)/sbas.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/sbas.c
preceph.o  : $(RTKLIB_SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/preceph.c
pntpos.o   : $(RTKLIB_SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/pntpos.c
ephemeris.o: $(RTKLIB_SRC)/ephemeris.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ephemeris.c
options.o  : $(RTKLIB_SRC)/options.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/options.c
ppp.o      : $(RTKLIB_SRC)/ppp.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ppp.c
ppp_ar.o   : $(RTKLIB_SRC)/ppp_ar.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ppp_ar.c
rtcm.o     : $(RTKLIB_SRC)/rtcm.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm.c
rtcm2.o    : $(RTKLIB_SRC)/rtcm2.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm2.c
rtcm3.o    : $(RTKLIB_SRC)/rtcm3.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm3.c
rtcm3e.o   : $(RTKLIB_SRC)/rtcm3e.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm3e.c
ionex.o    : $(RTKLIB_SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ionex.c
tides.o    : $(RTKLIB_SRC)/tides.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/tides.c
trace.o    : $(RTKLIB_SRC)/trace.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/trace.c

rnx2rtkp.o : $(RTKLIB_INC)/rtklib.h
rtkcmn.o   : $(RTKLIB_INC)/rtklib.h
rinex.o    : $(RTKLIB_INC)/rtklib.h
rtkpos.o   : $(RTKLIB_INC)/rtklib.h
postpos.o  : $(RTKLIB_INC)/rtklib.h
solution.o : $(RTKLIB_INC)/rtklib.h
lambda.o   : $(RTKLIB_INC)/rtklib.h
geoid.o    : $(RTKLIB_INC)/rtklib.h
sbas.o     : $(RTKLIB_INC)/rtklib.h
preceph.o  : $(RTKLIB_INC)/rtklib.h
pntpos.o   : $(RTKLIB_INC)/rtklib.h
ephemeris.o: $(RTKLIB_INC)/rtklib.h
options.o  : $(RTKLIB_INC)/rtklib.h
ppp.o      : $(RTKLIB_INC)/rtklib.h
ppp_ar.o   : $(RTKLIB_INC)/rtklib.h
rtcm.o     : $(RTKLIB_INC)/rtklib.h
rtcm2.o    : $(RTKLIB_INC)/rtklib.h
rtcm3.o    : $(RTKLIB_INC)/rtklib.h
rtcm3e.o   : $(RTKLIB_INC)/rtklib.h
ionex.o    : $(RTKLIB_INC)/rtklib.h
tides.o    : $(RTKLIB_INC)/rtklib.h
trace.o    : $(RTKLIB_INC)/rtklib.h

CMD1    = ./rnx2rtkp
INPUT11 = ../../../test/data/rinex/07590920.05o ../../../test/data/rinex/30400920.05n
INPUT12 = ../../../test/data/rinex/30400920.05o
OPTS1   = -r -3978241.958 3382840.234 3649900.853

test : test1 test2 test3 test4 test5 test6 test7 test8 test9 test10
test : test11 test12 test13 test14 test15 test16 test17 test18 test19 test20
test : test21 test22 test23 test24

test1 :
	$(CMD1) $(INPUT11) -x 5 -o test1.pos
test2 :
	$(CMD1) -t -e $(OPTS1) $(INPUT11) > test2.pos
test3 :
	$(CMD1) -t -p 1 -e $(OPTS1) $(INPUT11) $(INPUT12) > test3.pos
test4 :
	$(CMD1) -t -p 3 -e $(OPTS1) $(INPUT11) $(INPUT12) > test4.pos
test5 :
	$(CMD1) -t -m 15 -e $(OPTS1) $(INPUT11) $(INPUT12) > test5.pos
test6 :
	$(CMD1) -t -f 1 -e $(OPTS1) $(INPUT11) $(INPUT12) > test6.pos
test7 :
	$(CMD1) -t -v 5 -e $(OPTS1) $(INPUT11) $(INPUT12) > test7.pos
test8 :
	$(CMD1) -t -i -e $(OPTS1) $(INPUT11) $(INPUT12) > test8.pos
test9 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) > test9.pos
test10 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) -o test10.pos
test11 :
	$(CMD1) -t -p 0 -n $(OPTS1) $(INPUT11) > test11.pos
test12 :
	$(CMD1) -t -p 0 -g $(OPTS1) $(INPUT11) > test12.pos
test13 :
	$(CMD1) -t -p 0 $(OPTS1) $(INPUT11) > test13.pos
test14 :
	$(CMD1) -t -p 0 -u $(OPTS1) $(INPUT11) > test14.pos
test15 :
	$(CMD1) -t -p 0 -d 9 $(OPTS1) $(INPUT11) > test15.pos
test16 :
	$(CMD1) -t -p 0 -s , $(OPTS1) $(INPUT11) > test16.pos
test17 :
	$(CMD1) -t -b -e $(OPTS1) $(INPUT11) $(INPUT12) > test17.pos
test18 :
	$(CMD1) -t -c -e $(OPTS1) $(INPUT11) $(INPUT12) > test18.pos
test19 :
	$(CMD1) -t -h -e $(OPTS1) $(INPUT11) $(INPUT12) > test19.pos
test20 :
	$(CMD1) -t -p 4 -a $(OPTS1) $(INPUT11) $(INPUT12) > test20.pos
test21 :
	$(CMD1) $(INPUT11) $(INPUT12) > test21.pos
test22 :
	$(CMD1) -k opts1.conf $(INPUT11) $(INPUT12) > test22.pos
test23 :
	$(CMD1) -k opts2.conf $(INPUT11) $(INPUT12) > test23.pos
test24 :
	$(CMD1) -k opts3.conf $(INPUT11) $(INPUT12) -y 2 -o test24.pos
test25 :
	$(CMD1) -k opts4.conf $(INPUT11) $(INPUT12) -y 2 -o test25.pos

clean :
	rm -f rnx2rtkp rnx2rtkp.exe *.o *.pos *.trace

prof :
	gprof rnx2rtkp > prof.txt

install :
	cp rnx2rtkp $(BINDIR)

