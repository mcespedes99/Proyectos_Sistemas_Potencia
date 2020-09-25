using LinearAlgebra
import Pkg 
Pkg.add("IterativeSolvers")
using IterativeSolvers

barras = 6 #se define la cantidad de barras
r = zeros(barras,barras) #es una matriz con las resitencias, donde r[x,y] es la resitencia entre la barra x y la y.
x = zeros(barras,barras) #es una matriz con las impedancias, donde x[x,y] es la impedancia entre la barra x y la y.
bcap = zeros(barras,barras) #es una matriz con los bcap, donde bcap[x,y] es el bcap entre la barra x y la y.
r[1,2]=0.1
x[1,2]=0.2
bcap[1,2]=0.02
r[1,4]=0.05
x[1,4]=0.2
bcap[1,4]=0.02
r[1,5]=0.08
x[1,5]=0.3
bcap[1,5]=0.03
r[2,3]=0.05
x[2,3]=0.25
bcap[2,3]=0.03
r[2,4]=0.05
x[2,4]=0.1
bcap[2,4]=0.01
r[2,5]=0.1
x[2,5]=0.3
bcap[2,5]=0.02
r[2,6]=0.07
x[2,6]=0.2
bcap[2,6]=0.025
r[3,5]=0.12
x[3,5]=0.26
bcap[3,5]=0.025
r[3,6]=0.02
x[3,6]=0.1
bcap[3,6]=0.01
r[4,5]=0.2
x[4,5]=0.4
bcap[4,5]=0.04
r[5,6]=0.1
x[5,6]=0.3
bcap[5,6]=0.03
# debido a que estas tres matrices son simétricas se añade esta simetría, sin la necesidad de otras 32 líneas de código.
for i in 1:barras
  for j in 1:barras
    if i > j
      r[i,j]= r[j,i]
      x[i,j]= x[j,i]
      bcap[i,j]= bcap[j,i]
    end
  end
end
# ahora se forma la matriz de admitancias utilizando las fórmulas conocidas, donde Y=r/(r^2+x^2)+jx/(r^2+x^2)
#Y utilizando el modelo pi entre las líneas, por lo que bcap es una admitancia que va de la línea a tierra.

matriz = zeros(Complex{Float64},barras,barras)
for i in 1:barras
  for j in 1:barras
    if j == i
      matriz[i,i] = matriz[i,i] -(bcap[i,i])im
    else
      if r[i,j] != 0
        matriz[i,j]= -r[i,j]/(r[i,j]^2+x[i,j]^2) + (x[i,j]/(r[i,j]^2+x[i,j]^2))im
        cap = (bcap[i,j])*1im
        matriz[i,i]+=cap-matriz[i,j]
      end    
    end
  end
end
#atriz[1,1]=real(matriz[1,1])-11.8179152im
#matriz[2,2]=real(matriz[2,2])-23.1954965796322im
#matriz[3,3]=real(matriz[3,3])-16.5672696214169im
#matriz[4,4]=real(matriz[4,4])-14.6358820740134im
#matriz[5,5]=real(matriz[5,5])-14.1377649139613im
#matriz[6,6]=real(matriz[6,6])-17.0047269444913im
# 1. Definición del vector X: (vectores de tensiones y ángulos de las barras)
# X = [V1 V2 V3 V4 V5 V6 cita_1 cita_2 cita_3 cita_4 cita_5 cita_6]
# Aquellas V desconocidas se asumen como 1 en la primera iteración y las cita 0.
X = [1.05 1.05 1.07 1 1 1; 0 0 0 0 0 0]
#X = [1.05 1.05 1.07 0.9896 0.98565 1.004; 0 -0.0645772 -0.0750492 -0.0733038 -0.0925025 -0.102974] 

#Potencias de las barras. Las que no se conocen, se asignan como cero; de manera que la dimensión del vector sea la correcta:
Pbarras= [0 0.5 0.6 -0.7 -0.7 -0.7]
Qbarras = [0 0 0 -0.7 -0.7 -0.7]

#Definición del vector de missmatch

#Se especifica un valor para incializar el ciclo while
global max = 1
global contador = 0



while max > 1e-6
  global vector_missmatch = zeros(8)
  #Se recorre para las barras de la 2 a la 6
  for j = 2:6
    #Si son las barras 2 y 3 (j<4) entonces solo se toma en cuenta la ec. de P:
    if j < 4
      for i in 1:barras  
        vector_missmatch[j-1] += X[1,j]*X[1,i]*(real(matriz[i,j])*cos(X[2,j]-X[2,i])+imag(matriz[i,j])*sin(X[2,j]-X[2,i]))
      end
      vector_missmatch[j-1] -= Pbarras[j]
    #Para las demás barras se usa tanto la ecuación de P como la de Q
    else
      for i in 1:barras
        #Se forma la ecuación de P para la barra j
        vector_missmatch[j-1] += X[1,i]*X[1,j]*(real(matriz[j,i])*cos(X[2,j]-X[2,i])+imag(matriz[j,i])*sin(X[2,j]-X[2,i])) 
      end
      vector_missmatch[j-1] -= Pbarras[j]
      for i in 1:barras 
        #Se forma la ecuación de Q para la barra j
        vector_missmatch[j+2] += X[1,i]*X[1,j]*(real(matriz[j,i])*sin(X[2,j]-X[2,i])-imag(matriz[j,i])*cos(X[2,j]-X[2,i])) 
      end
      vector_missmatch[j+2] -= Qbarras[j]
    end
  end



  # 3. Cálculo del Jacobiano para la iteración correspondiente:
  J = zeros(8,8)

  # Primero se agrega la submatriz de dP/d_cita que es 5x5
  for i = 2:6
    for j = 2:6
      if i==j
        for k = 1:6
          J[i-1,j-1] += -X[1,i]*X[1,k]*(real(matriz[i,k])*sin(X[2,i]-X[2,k])-imag(matriz[i,k])*cos(X[2,i]-X[2,k]))
        end
        #println(i)
        #println(J[i-1,j-1])
        J[i-1,j-1] -= X[1,i]^2 * imag(matriz[i,i])
        #println(J[i-1,j-1])
      else
        J[i-1,j-1] = X[1,i]*X[1,j]*(real(matriz[i,j])*sin(X[2,i]-X[2,j]) - imag(matriz[i,j])*cos(X[2,i]-X[2,j]))
      end
    end
  end

  # A continuación se añade la submatriz dP/dV que es de 5x3
  for i = 2:6
    for j = 4:6
      if i==j
        for k = 1:6
          J[i-1,j-1+3] += X[1,k]*(real(matriz[i,k])*cos(X[2,i]-X[2,k])+imag(matriz[i,k])*sin(X[2,i]-X[2,k]))
        end
        J[i-1,j-1+3] += X[1,i]*real(matriz[i,i])
      else
        J[i-1,j-1+3] = X[1,i]*(real(matriz[i,j])*cos(X[2,i]-X[2,j]) + imag(matriz[i,j])*sin(X[2,i]-X[2,j]))
      end
    end
  end

  # Como tercer punto, de añade la submatriz dQ/d_cita que es 3x5
  for i = 4:6
    for j = 2:6
      if i==j
        for k in 1:6
          J[i-1+3,j-1] += X[1,i]*X[1,k]*(real(matriz[i,k])*cos(X[2,i]-X[2,k])+imag(matriz[i,k])*sin(X[2,i]-X[2,k]))
        end
        J[i-1+3,j-1] += -1* X[1,i]^2 *real(matriz[i,i])
      else
        J[i-1+3,j-1] = -X[1,i]*X[1,j]*(real(matriz[i,j])*cos(X[2,i]-X[2,j]) + imag(matriz[i,j])*sin(X[2,i]-X[2,j]))
      end
    end
  end

  #Por último, se añade la submatriz dQ/dV que es 3x3
  for i = 4:6
    for j = 4:6
      if i==j
        for k = 1:6
          J[i-1+3,j-1+3] += X[1,k]*(real(matriz[i,k])*sin(X[2,i]-X[2,k])-imag(matriz[i,k])*cos(X[2,i]-X[2,k]))
        end
        J[i-1+3,j-1+3] += -1*X[1,i]*imag(matriz[i,i])
      else
        J[i-1+3,j-1+3] = X[1,i]*(real(matriz[i,j])*sin(X[2,i]-X[2,j]) - imag(matriz[i,j])*cos(X[2,i]-X[2,j]))
      end
    end
  end
  #Fin del Jacobiano

  #Cálculo de LU=J: 
  L,U,p = lu(J)


  y = lsmr(L,-vector_missmatch)

  delta_x = lsmr(U,y)  
  delta_x2 = -inv(J)*vector_missmatch
  X[2,2:6]=X[2,2:6]+delta_x[1:5]
  X[1,4:6]=X[1,4:6]+delta_x[6:8]
  global contador +=1
  vector_missmatch_absoluto = broadcast(abs,vector_missmatch)
  global max = maximum(vector_missmatch_absoluto)
  println("\n Iteración número ", contador)
  for x =1:barras
    println("Tensión de la barra ", x,": ", round(X[1,x]*230,digits=4)," kV <",round(X[2,x]*(180/pi),digits = 4),"°")
  end
end

#Impresión de potencias:
P_finales = zeros(1,6)
Q_finales = zeros(1,6)
for k = 1:6
  for i = 1:6
    P_finales[1,k] += X[1,k]*X[1,i]*(real(matriz[k,i])*cos(X[2,k]-X[2,i])+imag(matriz[k,i])*sin(X[2,k]-X[2,i]))
    Q_finales[1,k] += X[1,k]*X[1,i]*(real(matriz[k,i])*sin(X[2,k]-X[2,i])-imag(matriz[k,i])*cos(X[2,k]-X[2,i]))
  end
end
println("\n Flujos de potencias de las barras hacia el sistema:")
for i = 1:barras
  println("Barra ", i, ": ", "P = ", round(P_finales[i]*100,digits = 4), " MW  Q = ", round(Q_finales[i]*100,digits=4), " MVAr")
end