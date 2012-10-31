# simple makefile for development.

SRC=_ercsmodule.c lib/*.[c.h]

ext3: ${SRC}
	rm -f _ercs.so
	python3 setup.py build_ext --inplace
ext2: ${SRC}
	rm -f _ercs.so
	python setup.py build_ext --inplace

figs:
	cd docs/asy && make 

docs: ext2 figs 
	cd docs && make clean && make html

	
