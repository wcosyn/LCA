

.PHONY: all
all: integral integral.so

.PHONY: clean
clean:
	rm -f integral integral.so


integral: integral.c
	gcc -O3 $^ -std=c99 -lgsl -lm -lblas -o $@

integral.so: integral.c
	gcc -O3 $^ -fPIC -std=c99 -lgsl -lm -lblas -shared -o $@
