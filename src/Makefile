OBJS=Heegaard1.o Heegaard2.o Heegaard3.o Heegaard4.o Heegaard5.o Heegaard6.o \
     Heegaard7.o Heegaard8.o Heegaard9.o Heegaard10.o Heegaard11.o Heegaard12.o \
     Heegaard15.o Heegaard16.o utils.o qksort.o

HDRS=Heegaard.h Heegaard_Dec.h

heegaard: $(OBJS) $(HDRS)
	cc -o heegaard -lm -lc $(OBJS)
	cp heegaard ..

clean:
	-rm heegaard *.o
