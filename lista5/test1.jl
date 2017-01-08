include("blocksys.jl")
using blocksys
n = 10^4;
l = 4;
ck = 100.0;
file = "duza.txt";
blockmat(n, l, ck, file);
A, n, l = wczytaj_A("duza.txt");

#=
bcopy = copy(b);

@time gauss(A, b, n, l);



println("gauss");
@time gauss(A, bcopy, n, l);
=#

#Acopy = copy(A);

println("gauss");
blad_wzgl, b, x1 = oblicz_b(A, n, l, 0)
println("blad wzgledny gauss: ", blad_wzgl);

println("gauss");
blad_wzgl, b, x1 = oblicz_b(A, n, l, 0)
println("blad wzgledny gauss: ", blad_wzgl);

println("czesciowy_gauss2");
blad_wzgl, b, x1 = oblicz_b(A, n, l, 1)
println("blad wzgledny czesciowy_gauss: ", blad_wzgl);
