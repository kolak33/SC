#autor Jakub Kołakowski

include("zad3.jl")

g(x) = exp(x);
a = 0.0;
b = 1.0;
n = 5;
rysujNnfx(g, a, b, n);

n = 10;
#rysujNnfx(g, a, b, n);
n = 15;
#rysujNnfx(g, a, b, n);


h(x) = x*x*sin(x);
a = -1.0;
b = 1.0;
n = 5;
#rysujNnfx(h, a, b, n);

n = 10;
#rysujNnfx(h, a, b, n);

n = 15;
#rysujNnfx(h, a, b, n);
