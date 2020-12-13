using Printf, Plots
gr(size=(600,400))
using PrettyTables
using Clp
using JuMP
using GLPK
using CSV
using DataFrames
CSV.read("C:\\Users\\Victor\\Desktop\\data\\export_df.csv", DataFrame)
CSV.File("C:\\Users\\Victor\\Desktop\\data\\export_df.csv"; normalizenames=true)
p = CSV.File("C:\\Users\\Victor\\Desktop\\data\\export_df.csv")
function PCV(p,gasolina,cpl,io::IO = stdout)
    m=length(p)
    s = rand(m,m)
    for i = 1:m
        for j= 1:m
            x = 6371*acos((cos((90-p[i][2])/180*pi))*cos((90-p[j][2])/180*pi) + sin((90-p[i][2])/180*pi)*sin((90-p[j][2])*pi/180)*cos((p[i][3]-p[j][3])*pi/180))*1.15
            # formula de Haversine, utiliza as leis do cossenos considerando o modelo a curvatura da terra onde o raio é 6371 
            # poderiamos também utilizar a formula m(x)=129,3799*x -34,8839 que foi adquirida atraves de uma regresão linear
            # onde meu x erra a distabncia entre pontos( no caso distancia entre latitudes e longitudes) e o meu y era a distancia em km
            s[i,j] = x
        end
    end
    f = Model(with_optimizer(GLPK.Optimizer))
    @variable(f, x[1:m,1:m], Bin)
    @objective(f, Min, sum(x[i,j]*s[i,j] for i=1:m,j=1:m))
    for i=1:m
        @constraint(f, x[i,i] == 0)
        @constraint(f, sum(x[i,1:m]) == 1)
    end
    for j=1:m
        @constraint(f, sum(x[1:m,j]) == 1)
    end
    # for i = 1:m
    #    for j = 1:m
            #coloquei @constraint (f, sum(x[i,j]) <= (m-1)) 
    #    end
    #end 
    # mas ainda ocorreu de dar o X com sub-rota então coloquei ela no final.
    for i = 1:m, j = 1:m
        @constraint(f, x[i,j]+x[j,i] <= 1) # isso impede de que você passe pelo mesmo local, 
        #no caso impede que a matriz s não pegue mesmo valores na matriz distancia como os valores de x[i,j] e x[j,i] são os mesmo na matriz s
        # um exemplo seria que nas minhas latitudes e longitudes, consegui uma matriz x com as distancias sem atribuir essa restrição
        # essa matriz x teve valores onde x[i,j] + x[j,i] = 2 isso é o mesmo que fazer ele percorrer uma rota duas vezes
        #Criou um ciclo de 2 para 4 e 4 para 2 e o mesmo para 1->5->3->1.
    end
    optimize!(f)
    km = 0.0
    custo = 0.0
    (l,t) = size(x)
    MatrizX = JuMP.value.(x)
    Solução = Int[]
    push!(Solução, 1.0)#começamos considerando que o 1 é a empresa.
    for i=1:l
        indice = argmax(MatrizX[Solução[end],1:l])#argmax pega o indice que é exatamente do que preciso para definir de que cidade
        #ele sai(no caso i) e chega( no caso j).
        if indice == Solução[1]# aqui é só para pular o primeiro elemento que é a empresa.
            break
        else
            push!(Solução,indice)
        end
    end
    if length(Solução) < l
        @constraint(f, sum(x[Solução,Solução]) <= length(Solução)-1)#Ultima restrição, ela garante que se causar uma sub-rota
        #ele ira restringir para que faça novamente o for mas com essa restrição estabelecida, que no caso se o tamanho da Solução for 
        #menor do que a quantidade de cidades(ocorreu sub-rota), a restrição faz a soma de todos os lementos de x e confere se
        # ela é menor do que a quantidade de cidades menos 1, se ocorreu sub-rota como por exemplo
        # 1->5->3->1 temos o sum dando 3 que é menor que 4
    end
    println(io, "|Cidades percorridas em ordem|     Km rodados     |    Custo da viajem   |")
    println(io, "|----------------------------|--------------------|----------------------|")
    println(io, "|  1                         |    0.0             |     0.0              |")
    for i=2:m
        km = km + s[Solução[i],Solução[i-1]]
        custo = custo + ((gasolina/cpl)*km)
        ç = @sprintf("| %2d                         | %17.15lf|     %10.4e       |", Solução[i], km, custo)
        println(io,ç)  
    end 
end