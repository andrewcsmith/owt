command="clang -std=c99 -Wall -g "

all: owt
	echo All done

owt: owt.o main.o
	clang -g -lgsl -lm -lblas -o $@ $^

%.o: %.c
	clang -std=c99 -g -c $<

owt.o: owt.c owt.h
main.o: main.c

clean: 
	rm main.o owt

