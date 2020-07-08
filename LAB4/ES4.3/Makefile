LIBSO:= MolDyn_NVE.o	#da aggiungere accanto compito.cpp
LISTH:= MolDyn_NVE.h	#da aggiungere accanto %.cpp %.h
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean: 
	rm -rf *.o
	rm -rf main.exe
	rm -rf *output*
	rm -rf *ave*	
	rm -rf old.final
	rm -rf config.final
	rm -rf frames/*.xyz

cpliquid: argon.liquid input.dat
	cp $^

cpsolid: argon.solid input.dat
	cp $^

cpgas: argon.gas input.dat
	cp $^

esegui:	main.exe
	./main.exe $(ARGS)
