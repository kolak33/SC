#autor Jakub Kołakowski
include("zad3.jl")

h(x) = x * x;
a = -1.0;
b = 1.0;
n = 5;
rysujNnfx(h, a, b, n);


h(x) = cos(x);
a = -3.14;
b = 3.14;
n = 10;
#rysujNnfx(h, a, b, n);


h(x) = x;
a = -2.0;
b = 2.0;
n = 5;
#rysujNnfx(h, a, b, n);
