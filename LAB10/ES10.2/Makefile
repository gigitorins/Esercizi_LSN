LIBSO:= genetic.o		#da aggiungere accanto main.cpp
LISTH:= genetic.h		#da aggiungere accanto %.cpp %.h
LIBS:= `root-config --libs`
INCS:= `root-config --cflags`
CC = mpicxx -g

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -O3 -Wall -o $@ ${INCS} ${LIBS}

%.o: %.cpp %.h ${LISTH}
	$(CC) -Wall -c $< ${INCS}

clean:
	rm -f *.o
	rm -rf main.exe
	rm -rf *results*
	rm -rf *best_conf*
	rm -rf *walks*
	rm -rf *averages*

esegui:	main.exe
	mpirun -np 4 ./main.exe $(ARGS)