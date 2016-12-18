#autor Jakub Kołakowksi
include("zad3.jl")


h(x) = abs(x);
a = -1.0;
b = 1.0;
n = 5;
#rysujNnfx(h, a, b, n);

n = 10;
#rysujNnfx(h, a, b, n);

n = 15;
rysujNnfx(h, a, b, n);


h(x) = 1.0 / (1.0 + x*x);
a = -5.0;
b = 5.0;
n = 5;
#rysujNnfx(h, a, b, n);

n = 10;
#rysujNnfx(h, a, b, n);

n = 15;
#rysujNnfx(h, a, b, n);
