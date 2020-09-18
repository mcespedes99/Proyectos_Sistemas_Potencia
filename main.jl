using LinearAlgebra
using LinearAlgebra
println("Inserte el número de barras")
barras = parse(Int8, readline()) #cantidad de barras del problema
println("Inserte el número de barras oscilantes")
oscilante = parse(Int8, readline()) #cantidad de barras oscilantes del problema
println("Inserte el número de barras generadoras")
generadoras = parse(Int8, readline()) #cantidad de barras generadoras del problema
println("Inserte el número de barras de carga")
cargas = parse(Int8, readline()) #cantidad de barras cargas del problema
admitancia = zeros(Int8,barras,barras) #crea una matriz de dimensión de la cantidad de barras y la llena de ceros.

for i in 1:barras
  for n in 1:barras
    if n == i
      println("Inserte la admitancia entre $i y tierra")
      admitancia[i,i] = admitancia[i,i] + parse(Int8, readline())
    else
      println("Inserte la admitancia entre $i y $n")
      admitancia[i,n] = -parse(Int8, readline())
      admitancia[i,i] = admitancia[i,i] -admitancia[i,n]
    end
  end 
end# Este código arma la matriz de admitancias.
