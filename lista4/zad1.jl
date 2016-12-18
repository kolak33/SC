#autor Jakub Kołakowski

#   funkcja ilorazy roznicowe
#   Dane:
#   x - wektor dlugosci n+1 zawierajacy wezly x0, x1, ..., xn
#   f - wektor dlugosci n+1 zawierajacy wartosci interpolowanej funkcji w wezlach f(x0), ..., f(xn)
#   Wyniki:
#   fx - wektor dlugosci n+1 zawierajacy obliczone ilorazy roznicowe

function ilorazyRoznicowe(x::Vector{Float64}, f::Vector{Float64})
    fx = copy(f); # kopia wektora f
    for i = 1 : length(f) - 1
        for j = length(f) : -1 : i+1
            fx[j] = (fx[j] - fx[j-1]) / (x[j] - x[j-i]);
        end
    end
    return fx;
end
