LIBSO:= MolDyn_NVE.o	#da aggiungere accanto compito.cpp
LISTH:= MolDyn_NVE.h	#da aggiungere accanto %.cpp %.h
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean: 
	rm -f *.o
	rm -f main.exe
	rm -rf *epot*
	rm -rf *ekin*
	rm -rf *etot*
	rm -rf *temp*
	rm -rf old.final
	rm -rf config.final
	rm -rf frames/*.xyz

cpgas: input.gas input.dat
	cp $^

cpliquid: input.liquid input.dat
	cp $^

cpsolid: input.solid input.dat
	cp $^

esegui:	main.exe
	./main.exe $(ARGS)
