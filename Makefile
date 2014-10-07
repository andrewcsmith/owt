all: owt

owt: main.o
	gcc -std=c99 -Wall -g -lgsl -lblas -lm -o owt main.c

clean: 
	rm main.o owt

