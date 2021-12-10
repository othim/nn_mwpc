

CC = g++

CFLAGS = -Wall -std=c++11 -fPIC -I/opt/anaconda3/envs/nn-mwpc-env/include -I/opt/anaconda3/include/python3.8 -I/opt/anaconda3/include/python3.8

run: header.o example.o
	$(CC) -shared -o example.so -undefined dynamic_lookup header.o example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o
   
header.o: header.cpp
	$(CC) -Wall -fPIC -std=c++11 -c header.cpp -o header.o	

clean:
	rm example.o header.o example.so
