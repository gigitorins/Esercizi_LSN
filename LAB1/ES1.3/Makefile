LIBSO:= random.o	#da aggiungere accanto compito.cpp
LISTH:= random.h	#da aggiungere accanto %.cpp %.h
CC = g++ -g
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean: 
	rm -f *.o
	rm -rf main.exe seed.out

esegui:	main.exe
	./main.exe $(ARGS)