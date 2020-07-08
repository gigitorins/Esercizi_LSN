LIBSO:= Monte_Carlo_NVT.o random.o		#da aggiungere accanto main.cpp
LISTH:= Monte_Carlo_NVT.h random.h		#da aggiungere accanto %.cpp %.h
CC = g++ -g
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean:
	rm -f *.o
	rm -rf main.exe seed.out
	rm -rf output.*
	rm -rf frames/*.xyz

cpliquid: input.liquid input.dat
	cp $^

cpsolid: input.solid input.dat
	cp $^

cpgas: input.gas input.dat
	cp $^

esegui:	main.exe
	./main.exe $(ARGS)