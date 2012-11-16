# Standard gcc/linux CFLAGS and LDFLAGS
CFLAGS=-g -Wall -Wextra
LDFLAGS=-lm -lgsl -lgslcblas
MAIN_LDFLAGS=${LDFLAGS} -lconfig 
ALL_SRC=ercs.c util.c kdtree.c torus.c
CC=cc

all_test: test_torus test_kdtree test_ercs

check: all_test
	./test_torus
	./test_kdtree
	./test_ercs

ercs: main.c ${ALL_SRC} *.h
	${CC} ${CFLAGS} -o ercs main.c ${ALL_SRC} ${MAIN_LDFLAGS} 

test_ercs: test_ercs.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_ercs test_ercs.c ${ALL_SRC} ${LDFLAGS} 

test_torus: test_torus.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_torus test_torus.c ${ALL_SRC} ${LDFLAGS} 

test_kdtree: test_kdtree.c ${ALL_SRC}
	${CC} ${CFLAGS} -o test_kdtree test_kdtree.c ${ALL_SRC} ${LDFLAGS} 

tags:
	etags *.c *.h

clean:
	rm -f ercs test_torus test_kdtree test_ercs
