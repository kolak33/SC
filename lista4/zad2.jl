#autor Jakub Kołakowski

#   funkcja wartosc wielomianu interpolacyjnego stopnia n w postaci Newtona w punkcje x = t, za pomoca uogolnionego algorytmu Hornera
#   Dane:
#   x - wektor dlugosci n+1 zawierajacy wezly x0, x1, ..., xn
#   fx - wektor dlugosci n+1 zawierajacy obliczone ilorazy roznicowe
#   t - punkt w ktorym trzeba obliczyc wartosc wielomianu
#   Wyniki:
#   nt - wartosc wielomianu w punkcie t

function warNewton(x::Vector{Float64}, fx::Vector{Float64}, t::Float64)
    nt = fx[1];
    wspol = t - x[1];
    for i=2 : length(fx)
        nt += wspol * fx[i];
        wspol = wspol * (t - x[i]);
    end
return nt;
end
