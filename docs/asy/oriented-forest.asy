
real w = 150mm; 

real x = w / 10;
real[] y = {0, x, 2x , 3x};
pair[] nodes;
int[] labels = {1, 3, 7, 9, 10, 2, 4, 5, 6, 8, 2};

int i;
for (i = 0; i < 10; ++i) {
    nodes[labels[i]] = (i * x, y[0]);

}
nodes[11] = (nodes[7].x + x / 2, y[1]);
nodes[12] = (nodes[4].x + x / 2, y[1]);
nodes[13] = (nodes[6].x + x / 2, y[1]);
nodes[14] = (nodes[12].x + x, y[2]);
nodes[15] = (nodes[12].x, y[3]);
nodes[16] = (nodes[9].x, y[2]);


for (i = 1; i < nodes.length; ++i) {
    dot(format("%d", i), nodes[i]);
}

draw(nodes[7] -- nodes[11] -- nodes[9]);
draw(nodes[4] -- nodes[12] -- nodes[5]);
draw(nodes[6] -- nodes[13] -- nodes[8]);
draw(nodes[12] -- nodes[14] -- nodes[13]);
draw(nodes[11] -- nodes[16] -- nodes[10]);
draw(nodes[2] -- nodes[15] -- nodes[14]);

