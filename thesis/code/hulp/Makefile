

.PHONY: all
all: hoWave hoWave.so

.PHONY: clean
clean:
	rm -f hoWave hoWave.so


hoWave: hoWave.c
	gcc -O3 $^ -std=c99 -lgsl -lm -lblas -o $@

hoWave.so: hoWave.c
	gcc -O3 $^ -fPIC -std=c99 -lgsl -lm -lblas -shared -o $@
