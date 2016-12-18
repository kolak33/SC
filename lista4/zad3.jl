#autor Jakub Kołakowski
include("zad1.jl");
include("zad2.jl");
using PyPlot
#   funkcja ktora interpoluje zadana funkcje f(x) w przedziale [a,b] za pomoca wielomianu interpolacyjnego stopnia n w postaci Newtona,
#   a nastepnie narysuje wielomian interpolacyjny i interpolowana funkcje
#   Dane:
#   f - funkcja zadana jako funkcja anonimowa
#   a,b - przedzial interpolacji
#   n - stopien wielomianu interpolowanego
#   Wyniki:
#   funkcja rysuje wielomian interpolacyjny i interpolowana funkcje w przedziale [a,b]


function rysujNnfx(f, a::Float64, b::Float64, n::Int)
    x = Vector{Float64}(n+1); #wektor dlugosci n+1 zawierajacy wezly x0, x1, ..., xn
    fx = Vector{Float64}(n+1); #wektor dlugosci n+1 zawierajacy obliczone ilorazy roznicowe
    wart = Vector{Float64}(n+1); #kolejne wartosci f(xi), dla i = 0,1,...,n


    przesuniecie = 0.0; #przesuniecie kolejnych punktow na przedziale [a,b]
    h = (b-a) / n;
    for i=1 : n+1
        x[i] = a + przesuniecie;
        wart[i] = f(x[i]);
        przesuniecie += h;
    end

    fx = ilorazyRoznicowe(x, wart);

    real_wart = Vector{Float64}(3n+3); #wartosci funkcji w kolejnych wezlach
    real_x =  Vector{Float64}(3n+3); #kolejne wezly
    horner_wart = Vector{Float64}(3n+3); #wartosci wielomianu interpolacyjnego w kolejnych punktach x0, x1, ..., xn

    przesuniecie = 0.0; #przesuniecie kolejnych punktow na przedziale [a,b]
    h = (b-a) / (3n + 2);
    for i=1 : 3*(n+1)
        real_x[i] = a + przesuniecie;
        horner_wart[i] = warNewton(x, fx, real_x[i]);
        real_wart[i] = f(real_x[i]);
        przesuniecie += h;
    end


    plot(real_x, real_wart, real_x, horner_wart, linewidth=2.0)
    #title("funkcje")


    #p = plot(wart);
end
