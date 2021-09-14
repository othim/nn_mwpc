

CC = g++
CFLAGS = -Wall -g


main: test.o pot_nn_mwpc.o Term.o
	$(CC) $(CFLAGS) -o main test.o pot_nn_mwpc.o Term.o -lgsl

test.o: test.cpp pot_nn_mwpc.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c test.cpp

pot_nn_mwpc.o: pot_nn_mwpc.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c pot_nn_mwpc.cpp

Term.o: Term.h
	$(CC) $(CFLAGS) -I/usr/local/Cellar/gsl/2.7/include/gsl -c Term.cpp

clean:
	rm main pot_nn_mwpc.o test.o Term.o
