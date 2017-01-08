#autor Jakub Kołakowski

module blocksys

export wczytaj_A, wczytaj_A_oraz_b, gauss, czesciowy_gauss, oblicz_b, zapisz_do_pliku_b, zapisz_do_pliku_sam_x, zapisz_do_pliku_blad_wzgl, test_calc_Ax_b, blockmat

function wczytaj_A(Afile::String)
    f = open(Afile,"r");

    linia1 = readline(f);
    n, l = split(linia1, " ");

        #inicjalizacja
        n = parse(Int64, n);
        l = parse(Int64, l);
        A = spzeros(Float64, n, n); # sparse vector nxn
        #wczytywanie do tablicy
        for line in eachline(f)
            i, j, wart = split(line, " ");
            A[parse(Int64, i), parse(Int64, j)] = parse(Float64, wart);
        end
        close(f);

        return (A,n,l);
end

function wczytaj_A_oraz_b(Afile::String, bfile::String)
#Afile = "lista5/A.txt";
f = open(Afile,"r");

linia1 = readline(f);
n, l = split(linia1, " ");

    #inicjalizacja
    n = parse(Int64, n);
    l = parse(Int64, l);
    A = spzeros(Float64, n, n); # sparse vector nxn
    b = Vector{Float64}(n);
    #wczytywanie do tablicy
    for line in eachline(f)
        i, j, wart = split(line, " ");
        A[parse(Int64, i), parse(Int64, j)] = parse(Float64, wart);
    end
    close(f);
    #wczytywanie wektora prawych stron b
    #f = open("lista5/Dane1000_1_1000/b.txt", "r");
    #f = open("lista5/b.txt", "r");
    f = open(bfile, "r");
    line = readline(f); #wyrzucamy pierwsza linijke w ktorej jest liczba kolejnych linii bo juz i tak ją mamy
    iter = 1;
    for line in eachline(f)
        b[iter] = parse(Float64, line);
        iter += 1;
    end
    close(f);

    return (A,b,n,l);
end

function gauss(A, b::Vector{Float64}, n::Int64, l::Int64)
    x = Vector{Float64}(n);
    for k=1 : n-1
        for i=k+1 : k + l - k%l
            wspol = A[i, k] / A[k, k];
            #A[i, k] = 0.0;
            for j=k+1 : min(n, k+l)
                A[i, j] = A[i, j] - wspol * A[k, j];
            end
            b[i] = b[i] - wspol*b[k];
        end
    end


    for i=n : -1 : 1
        suma = 0;
        for j=i+1 : min(n, i+l)
            suma += A[i,j]*x[j];
        end
        x[i] = (b[i] - suma) / A[i,i];
    end

return x;
end

function czesciowy_gauss(A, b::Vector{Float64}, n::Int64, l::Int64)
    x = Vector{Float64}(n);
    #temp = Array(Float64, 1, 2*l+1);
    #tempb = 0.0;

    max = 0.0
    maxrow = 0
    for k=1 : n-1
        max = abs(A[k, k]);
        maxrow = k;
        for i=k+1 : k + l - k%l
            if(abs(A[i,k]) > max)
                max = abs(A[i,k]);
                maxrow = i;
            end
        end

        if maxrow > k
            for i = k : min(n, k+2*l)
                A[k,i], A[maxrow,i] = A[maxrow,i], A[k,i]
            end
            b[k], b[maxrow] = b[maxrow], b[k]
        end


        for i = k+1 : k+l-k%l
            Iik = A[i,k] / A[k,k]
            A[i,k] = 0.0 #Iik
            for j = k+1 : min(k+2*l, n)
                A[i,j] = A[i,j] - Iik * A[k,j]
            end
            b[i] = b[i] - Iik * b[k]
        end
    end


    #obliczanie x
        for i=n : -1 : 1
            suma = 0;
            for j=i+1 : min(n, i+l +l) # musimy isc o +l dalej
                suma += A[i,j]*x[j];
            end
            x[i] = (b[i] - suma) / A[i,i];
        end

    return x;
end



function oblicz_b(A, n::Int64, l::Int64, ktora_funkcja::Int64)
    x = Vector{Float64}(n);
    x[1 : end] = 1.0;
    b = Vector{Float64}(n);
    zeros(b);
    x1 = Vector{Float64}(n);

    #pierwszy blok
    for i=1 : l
        suma = 0;
        for j=1 : l+i
            suma += A[i,j]*x[j];
        end
        b[i] = suma;
    end

    #srodkowe bloki
    blok = 0;
    count = 1;
    for i=l+1 : n-l
        suma = 0;
        for j=l-1 + l*blok : l + l*blok +l+l # +1+l bo vector c(+l) oraz b(+1) dodatkowy
            suma += A[i,j]*x[j];
        end
        b[i] = suma;

        if(count % l == 0)
            blok += 1;
        end
        count += 1;
    end
    #ostatni blok

    for i=n-l+1 : n
        suma = 0;
        for j=n-l : n
            suma += A[i,j]*x[j];
        end
        b[i] = suma;
    end
    #c = A*x;
    #println(c);
    #println(b);
    #Acopy = copy(A);
    Acopy = copy(A);

    if(ktora_funkcja == 0)
        x1 = @time czesciowy_gauss(Acopy, b, n, l); # bo czesciowy gauss modyfikuje A
    else
        x1 = @time gauss(Acopy, b, n, l); #Acopy bo gauss modyfikuje A
    end

blad_wzgl = norm(x1 - x) / norm(x);
#println(blad_wzgl);
return (blad_wzgl, b, x1);
end

function wczytaj_A_oblicz_b(Afile::String, ktora_funkcja::Int64)
    f = open(Afile,"r");

    linia1 = readline(f);
    n, l = split(linia1, " ");

    #inicjalizacja
    n = parse(Int64, n);
    l = parse(Int64, l);
    A = spzeros(Float64, n, n); # sparse vector nxn

    #wczytywanie do tablicy
    for line in eachline(f)
        i, j, wart = split(line, " ");
        A[parse(Int64, i), parse(Int64, j)] = parse(Float64, wart);
    end
    close(f);

    blad_wzgl, b, x1 = oblicz_b(A, n, l, ktora_funkcja);

    return (blad_wzgl, b, x1, n);
end

#sprawdzanie poprawnosci obliczania
function test_calc_Ax_b(Afile::String, ktora_funkcja::Int64)
    granica_bledu = 10.0^-5;
    blad_wzgl, b, x1, n = wczytaj_A_oblicz_b(Afile, ktora_funkcja);
    if(blad_wzgl < granica_bledu) println("test OK");
    else println("test FAILED, za duzy blad wzgledy")
    end
end

#funkcja wczytuje wektory A oraz b, a nastepnie oblicza wektor x zadana metoda: 0 - gauss, 1 - czesciowy_gauss i zapisuje go do pliku
function zapisz_do_pliku_sam_x(plikA::String, plikb::String, plik_do_zapisu::String, ktora_funkcja::Int64)

    A, b, n, l = wczytaj_A_oraz_b(plikA, plikb);
    if(ktora_funkcja == 0)
        Acopy = copy(A);
        bcopy = copy(b);
        x = gauss(A, bcopy, n, l);
    else
        Acopy = copy(A);
        bcopy = copy(b);
        x = czesciowy_gauss(Acopy, bcopy, n, l);
    end

    open(plik_do_zapisu, "w") do f
		for i = 1 : n
			println(f, x[i]);
		end
	end
end

#funkcja wczytuje wektor A, oblicza wektor b, oblicza x zadana metoda: 0 - gauss, 1 - czesciowy_gauss i zapisuje blad_wzgl oraz x do pliku
function zapisz_do_pliku_blad_wzgl(plikA::String,  plik_do_zapisu::String, ktora_funkcja::Int64)

    blad_wzgl, b, x1, n = wczytaj_A_oblicz_b(plikA, ktora_funkcja);

    open(plik_do_zapisu, "w") do f
        println(f, blad_wzgl);
		for i = 1 : n
			println(f, x1[i]);
		end
	end
end

#zapisuje wektor prawych stron b do pliku
function zapisz_do_pliku_b(b::Vector{Float64}, n::Int64,  l::Int64,  plik_do_zapisu::String)


    open(plik_do_zapisu, "w") do f
        print(f, n);
        print(f, " ");
        println(f, l);
		for i = 1 : n
			println(f, b[i]);
		end
	end
end


#println(pwd());
#A,b,n,l = wczytaj_A_oraz_b("A.txt", "b.txt");
#A,b,n,l = wczytaj_A_oraz_b("Dane1000_1_1000/A.txt", "Dane1000_1_1000/b.txt");
#blad_wzgl, b = wczytaj_A_oblicz_b("A.txt", 1);
#println(b);
#x = czesciowy_gauss(A, b, n, l);
#x = gauss(A, b, n, l);
#println(x);
#blad, bb = oblicz_b(A, n, l);


#autor Paweł Zieliński
function matcond(n::Int, c::Float64)
		# Function generates a random square matrix A of size n with
		# a given condition number c.
		# Inputs:
		#	n: size of matrix A, n>1
		#	c: condition of matrix A, c>= 1.0
		#
		# Usage: matcond (10, 100.0);
		#

        if n < 2
         error("size n should be > 1")
        end
        if c< 1.0
         error("condition number  c of a matrix  should be >= 1.0")
        end
        (U,S,V)=svd(rand(n,n))
        return U*diagm(linspace(1.0,c,n))*V'
	end # matcond

#autor Paweł Zieliński
  function blockmat(n::Int, l::Int, ck::Float64, outputfile::String)
		# Function generates a random block sparse matrix A of size n with
		# a given condition number ck of inner block Ak and it save the output
		# matrix in a text file.
		# Inputs:
		#	n: size of block matrix A, n>1
		# l: size of inner matrices Ak, n mod l =0 (n is  divisible by l)
		#	ck: condition of inner matrix Ak, ck>= 1.0
		# outputfile: name of the output text file
		#
		# Usage: blockmat(100, 4 ,10.0, "A.txt")
		#
		#
		#  the output file format
	  #  n  l              <--- the size of block matrix A, the size of inner matrices Ak
		#  i1  j1   A[i1,j1] <--- a non-zero element of block matrix A
		#  i2  j2   A[i2,j2] <--- a non-zero element of block matrix A
		#  i3  j3   A[i3,j3] <--- a non-zero element of block matrix A
		#  ...
		#  ...
		#  EOF
		#

    if n < 2
     error("size n should be > 1")
    end
		if n%l!=0
			error("n is not divisible by l")
		end

		nb=div(n,l)
		Ak=Array(Float64,l,l)
		open(outputfile, "w") do f
			println(f, n," ",l)
			for k in 1:nb
				Ak=matcond(l, ck)
				for i in 1:l, j in 1:l
					println(f,(k-1)*l+i," ",(k-1)*l+j," ", Ak[i,j])
				end
				if k<nb
			  	 for i in 1:l
						println(f,(k-1)*l+i," ",k*l+i," ",0.3*rand())
				 	 end
				end
				if k>1
			   for i in 1:l
					 println(f,(k-1)*l+i," ",(k-1)*l," ",0.3*rand())
				 end
				end
			end
		end	 # do
	end # blockmat



end #blocksys
