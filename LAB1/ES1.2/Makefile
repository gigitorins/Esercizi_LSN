LIBSO:= random.o	#da aggiungere accanto compito.cpp
LISTH:= random.h	#da aggiungere accanto %.cpp %.h
CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean: 
	rm -f *.o
	rm -rf main.exe seed.out
#	rm -rf file1.txt file2.txt
#	rm -rf file20.txt file100.txt

esegui:
	./main.exe
