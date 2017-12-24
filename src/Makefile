OBJS=Heegaard1.o Heegaard2.o Heegaard3.o Heegaard4.o Heegaard5.o Heegaard6.o \
     Heegaard7.o Heegaard8.o Heegaard9.o Heegaard10.o Heegaard11.o Heegaard12.o\
     Heegaard13.o Heegaard14.o Heegaard15.o Heegaard16.o Heegaard17.o Heegaard18.o\
     Heegaard19.o Heegaard20.o Heegaard21.o utils.o qksort.o

HDRS=Heegaard.h Heegaard_Dec.h

CFLAGS=-g -Wall

Heegaard: $(OBJS) $(HDRS)
	$(CC) -o ../heegaard -lm -lc $(OBJS)

all: Heegaard

clean:
	-rm *.o
