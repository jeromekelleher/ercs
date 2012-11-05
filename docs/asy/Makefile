
all: ../images/oriented-tree.png ../images/oriented-forest.png


../images/oriented-tree.png: oriented-tree.asy
	asy oriented-tree.asy -o tmp.eps
	convert -density 300 -units PixelsPerInch tmp.eps ../images/oriented-tree.png
	rm -f tmp.eps

../images/oriented-forest.png: oriented-forest.asy
	asy oriented-forest.asy -o tmp.eps
	convert -density 300 -units PixelsPerInch tmp.eps ../images/oriented-forest.png
	rm -f tmp.eps
