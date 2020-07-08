LIBSO:= qmc1d.o		#da aggiungere accanto main.cpp
LISTH:= qmc1d.h		#da aggiungere accanto %.cpp %.h
LIBS:= `root-config --libs`
INCS:= `root-config --cflags`
CC = g++ -g

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -O3 -Wall -o $@ ${INCS} ${LIBS}

%.o: %.cpp %.h ${LISTH}
	$(CC) -Wall -c $< ${INCS}

clean:
	rm -f *.o main.exe
	rm -rf potential.dat kinetic.dat probability.dat

esegui:	main.exe
	./main.exe $(ARGS)