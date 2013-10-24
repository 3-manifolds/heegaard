OBJS=Heegaard1.o Heegaard2.o Heegaard3.o Heegaard4.o Heegaard5.o Heegaard6.o \
     Heegaard7.o Heegaard8.o Heegaard9.o Heegaard10.o Heegaard11.o Heegaard12.o \
     Heegaard15.o Heegaard16.o utils.o qksort.o

OBJS2=Realizable1.o Heegaard2.o Heegaard3.o Heegaard4.o Realizable5.o Heegaard6.o \
     Heegaard7.o Heegaard8.o Heegaard9.o Heegaard10.o Heegaard11.o Realizable12.o \
     Heegaard15.o Heegaard16.o utils.o qksort.o


HDRS=Heegaard.h Heegaard_Dec.h

CFLAGS = -Wall  -Wno-return-type -Wno-implicit-int -w

heegaard: $(OBJS) $(HDRS)
	cc -o heegaard -lm -lc $(OBJS)
	cp heegaard ..

heegaard_is_realizable: $(OBJS2) $(HDRS)
	cc -o heegaard_is_realizable -lm -lc $(OBJS2)
	cp heegaard_is_realizable ..

all: heegaard heegaard_is_realizable

clean:
	-rm heegaard heegaard_is_realizable  *.o
