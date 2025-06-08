# makefile for rtkrcv
CC = gcc
BINDIR = /usr/local/bin
RTKLIB_SRC    = ../../../rtklib_src
USER_SRC = ../../../user_src
RTKLIB_INC = ../../../rtklib_inc
USER_INC = ../../../user_inc
CTARGET= -DTRACE -DENAGLO -DENAQZS -DENACMP -DENAGAL -DENAIRN -DNFREQ=3 -DNEXOBS=3 -DSVR_REUSEADDR

CFLAGS = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(RTKLIB_SRC) -I$(USER_SRC) -DTRACE $(CTARGET) -g -std=c99
LDLIBS  = -lm -lrt -lpthread

all         : rtkrcv
rtkrcv     : rtkrcv.o vt.o rtkcmn.o rtksvr.o rtkpos.o geoid.o solution.o lambda.o
rtkrcv     : sbas.o stream.o rcvraw.o rtcm.o preceph.o options.o pntpos.o ppp.o ppp_ar.o
rtkrcv     : novatel.o ublox.o crescent.o skytraq.o javad.o nvs.o binex.o
rtkrcv     : rt17.o ephemeris.o rinex.o ionex.o rtcm2.o rtcm3.o rtcm3e.o
rtkrcv     : tides.o trace.o septentrio.o swiftnav.o unicore.o
rtkrcv     : rover_station.o base_station.o at_command.o
rtkrcv     : lora_WB26.o

rover_station.o	: $(USER_SRC)/rover_station.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/rover_station.c
base_station.o	: $(USER_SRC)/base_station.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/base_station.c
at_command.o	: $(USER_SRC)/at_command.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/at_command.c
lora_WB26.o	: $(USER_SRC)/lora_WB26.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/lora_WB26.c
rtkrcv.o   : $(USER_SRC)/rtkrcv.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/rtkrcv.c
vt.o       : $(USER_SRC)/vt.c
	$(CC) -c $(CFLAGS) $(USER_SRC)/vt.c
rtkcmn.o   : $(RTKLIB_SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtkcmn.c
rtksvr.o   : $(RTKLIB_SRC)/rtksvr.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtksvr.c
rtkpos.o   : $(RTKLIB_SRC)/rtkpos.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtkpos.c
geoid.o    : $(RTKLIB_SRC)/geoid.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/geoid.c
solution.o : $(RTKLIB_SRC)/solution.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/solution.c
lambda.o   : $(RTKLIB_SRC)/lambda.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/lambda.c
sbas.o     : $(RTKLIB_SRC)/sbas.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/sbas.c
stream.o   : $(RTKLIB_SRC)/stream.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/stream.c
rcvraw.o : $(RTKLIB_SRC)/rcvraw.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcvraw.c
rtcm.o     : $(RTKLIB_SRC)/rtcm.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm.c
rtcm2.o    : $(RTKLIB_SRC)/rtcm2.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm2.c
rtcm3.o    : $(RTKLIB_SRC)/rtcm3.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm3.c
rtcm3e.o   : $(RTKLIB_SRC)/rtcm3e.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rtcm3e.c
preceph.o  : $(RTKLIB_SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/preceph.c
options.o  : $(RTKLIB_SRC)/options.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/options.c
pntpos.o   : $(RTKLIB_SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/pntpos.c
ppp.o      : $(RTKLIB_SRC)/ppp.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ppp.c
ppp_ar.o   : $(RTKLIB_SRC)/ppp_ar.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ppp_ar.c
novatel.o  : $(RTKLIB_SRC)/rcv/novatel.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/novatel.c
ublox.o    : $(RTKLIB_SRC)/rcv/ublox.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/ublox.c
swiftnav.o : $(RTKLIB_SRC)/rcv/swiftnav.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/swiftnav.c
unicore.o:   $(RTKLIB_SRC)/rcv/unicore.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/unicore.c
crescent.o : $(RTKLIB_SRC)/rcv/crescent.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/crescent.c
skytraq.o  : $(RTKLIB_SRC)/rcv/skytraq.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/skytraq.c
javad.o    : $(RTKLIB_SRC)/rcv/javad.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/javad.c
nvs.o      : $(RTKLIB_SRC)/rcv/nvs.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/nvs.c
binex.o    : $(RTKLIB_SRC)/rcv/binex.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/binex.c
rt17.o     : $(RTKLIB_SRC)/rcv/rt17.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/rt17.c
ephemeris.o: $(RTKLIB_SRC)/ephemeris.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ephemeris.c
rinex.o    : $(RTKLIB_SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rinex.c
ionex.o    : $(RTKLIB_SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/ionex.c
tides.o    : $(RTKLIB_SRC)/tides.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/tides.c
trace.o    : $(RTKLIB_SRC)/trace.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/trace.c
septentrio.o: $(RTKLIB_SRC)/rcv/septentrio.c
	$(CC) -c $(CFLAGS) $(RTKLIB_SRC)/rcv/septentrio.c


rover_station.o: $(USER_INC)/user_rtklib.h
base_station.o: $(USER_INC)/user_rtklib.h
at_command.o: $(USER_INC)/user_rtklib.h $(USER_INC)/at_command.h
lora_WB26.o: $(USER_INC)/user_rtklib.h $(USER_INC)/lora_WB26.h
rtkrcv.o   : $(RTKLIB_INC)/rtklib.h $(USER_INC)/vt.h $(USER_INC)/at_command.h
rtkcmn.o   : $(RTKLIB_INC)/rtklib.h
rtksvr.o   : $(RTKLIB_INC)/rtklib.h
rtkpos.o   : $(RTKLIB_INC)/rtklib.h
geoid.o    : $(RTKLIB_INC)/rtklib.h
solution.o : $(RTKLIB_INC)/rtklib.h
lambda.o   : $(RTKLIB_INC)/rtklib.h
sbas.o     : $(RTKLIB_INC)/rtklib.h
rcvraw.o   : $(RTKLIB_INC)/rtklib.h
rtcm.o     : $(RTKLIB_INC)/rtklib.h
rtcm2.o    : $(RTKLIB_INC)/rtklib.h
rtcm3.o    : $(RTKLIB_INC)/rtklib.h
rtcm3e.o   : $(RTKLIB_INC)/rtklib.h
preceph.o  : $(RTKLIB_INC)/rtklib.h
options.o  : $(RTKLIB_INC)/rtklib.h
pntpos.o   : $(RTKLIB_INC)/rtklib.h
ppp.o      : $(RTKLIB_INC)/rtklib.h
ppp_ar.o   : $(RTKLIB_INC)/rtklib.h
novatel.o  : $(RTKLIB_INC)/rtklib.h
ublox.o    : $(RTKLIB_INC)/rtklib.h
swiftnav.o : $(RTKLIB_INC)/rtklib.h
crescent.o : $(RTKLIB_INC)/rtklib.h
skytraq.o  : $(RTKLIB_INC)/rtklib.h
javad.o    : $(RTKLIB_INC)/rtklib.h
nvs.o      : $(RTKLIB_INC)/rtklib.h
binex.o    : $(RTKLIB_INC)/rtklib.h
rt17.o     : $(RTKLIB_INC)/rtklib.h
septentrio.o: $(RTKLIB_INC)/rtklib.h
unicore.o  : $(RTKLIB_INC)/rtklib.h
ephemeris.o: $(RTKLIB_INC)/rtklib.h
rinex.o    : $(RTKLIB_INC)/rtklib.h
ionex.o    : $(RTKLIB_INC)/rtklib.h
trace.o    : $(RTKLIB_INC)/rtklib.h

install:
	cp rtkrcv $(BINDIR)

clean:
	rm -f rtkrcv rtkrcv.exe rtkrcv.nav *.o *.out *.trace *.conf

