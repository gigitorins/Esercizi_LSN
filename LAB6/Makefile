LIBSO:= Monte_Carlo_ISING_1D.o random.o		#da aggiungere accanto main.cpp
LISTH:= Monte_Carlo_ISING_1D.h random.h		#da aggiungere accanto %.cpp %.h
CC = g++ -g
CFLAGS = -Wall -O3 --std=c++11

main.exe: main.cpp ${LIBSO}
	$(CC) $^ -o $@

%.o: %.cpp %.h ${LISTH}
	$(CC) -c $<

clean:
	rm -f *.o
	rm -rf main.exe seed.out
	# rm -rf output.* config.final

cpinput: input0.dat input.dat		# 1: file da copiare --> 2: target
	cp $^

esegui:	main.exe
	./main.exe