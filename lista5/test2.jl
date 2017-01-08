include("blocksys.jl")
using blocksys
n = 10^4;
l = 8;
ck = 100.0;
file = "duza.txt";
blockmat(n, l, ck, file);
A, n, l = wczytaj_A("duza.txt");


x = Array(Float64, 1, n);
b = Array(Float64, 1, n);
x[1 : end] = 1.0;
x = x';
tic()
b = A*x;
x1 = A \ b;
toc()
blad_wzgl = norm(x1 - x) / norm(x);
println("blad wzgledny tradycyjnie: ", blad_wzgl);




tic()
b = A*x;
x1 = A \ b;
toc()
blad_wzgl = norm(x1 - x) / norm(x);
println("blad wzgledny tradycyjnie: ", blad_wzgl);
