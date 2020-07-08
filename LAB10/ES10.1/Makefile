LIBSO:= genetic.o		#da aggiungere accanto main.cpp
LISTH:= genetic.h		#da aggiungere accanto %.cpp %.h
LIBS:= `root-config --libs`
INCS:= `root-config --cflags`
CC = g++ -g

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -O3 -Wall -o $@ ${INCS} ${LIBS}

%.o: %.cpp %.h ${LISTH}
	$(CC) -Wall -c $< ${INCS}

clean:
	rm -f *.o
	rm -rf main.exe
	rm -rf averages.dat solution.dat walks.dat

esegui:	main.exe
	./main.exe $(ARGS)