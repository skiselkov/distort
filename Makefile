ACFUTILS=/home/diablos/NextDev/libacfutils

CFLAGS=-O2 -std=gnu89 -g -W -Wall -Werror -I$(ACFUTILS)/src
LDFLAGS=-L$(ACFUTILS)/qmake/lin64 -lacfutils -lm
OBJS=main.o distort.o kiss_fft.o

all : distort

distort : $(OBJS)
	gcc -o $@ $^ $(LDFLAGS)

%.o : %.c
	gcc $(CFLAGS) -c -o $@ $^ $(LDFLAGS)

clean :
	rm -f distort $(wildcard *.o)
