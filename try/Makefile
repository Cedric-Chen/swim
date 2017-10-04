CC=/usr/bin/gcc

CPPSETTINGS = -std=c++11

CCFLAGS= -O3 -fopenmp -fstrict-aliasing -fdiagnostics-show-option -Wall -Wno-reorder 
LDFLAGS= -O3

# Linking stage
test: testArma.o
	${CC} ${CPPSETTINGS} testArma.o -o $@  ${LDFLAGS}

# Compiling stage
testArma.o: testArma.cpp 
	${CC} ${CPPSETTINGS} ${CCFLAGS} -c -o $@ $<

clean:
	rm -f *.o
	rm -f *.s rm -f test

tidy: clean
	rm -f *.vti rm -f FFTW_report.txt
