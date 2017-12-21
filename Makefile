#LIBS = $(pkg-config --libs gsl libxrl)
LIBS = -I/usr/local/lib -lgsl -lgslcblas -lm -lxrl

CFLAGS = -O0 -g -Wall -fopenmp -I/usr/local/include/xraylib #$(pkg-config --cflags gsl libxrl)

CC = gcc
CC_SWITCHES =	${CFLAGS} 

OBJS = polycap_v2.2.o


SRCS = polycap_v2.2.c


polycap:	$(OBJS)
	${CC} ${CC_SWITCHES} $(OBJS) $(LIBS) -o polycap_v2.2


.c.o:
	$(CC) -c $(CC_SWITCHES) $<

clean:
	rm -f $(OBJS) polycap
