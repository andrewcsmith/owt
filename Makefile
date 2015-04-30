command="clang -std=c99 -Wall -g "

all: owt
	echo All done

owt: owt.o main.o owt_parser.o
	clang -g -lgsl -lm -lblas -o $@ $^

%.o: %.c
	clang -std=c99 -g -c $<

owt.o: owt.c owt.h
owt_parser.o: owt_parser.c owt_parser.h
main.o: main.c
key_weights_main.o: key_weights_main.c

clean: 
	rm -f main.o owt error.o owt.o

key_weights: owt.o key_weights_main.o owt_parser.o
	clang -g -lgsl -lm -lblas -o $@ $^
