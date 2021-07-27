# Mapeos proyectivos en sistemas de varios qubits

__José Alfredo de León, estudiante de Licenciatura en 
Física Aplicada, ECFM-USAC. 
Asesorado por Dr. Carlos Pineda (IFUNAM) y M.Sc. Juan Diego Chang (ECFM-USAC)__

##Proyecto de prácticas finales 

__*Descripción*__:
Paquete con funciones implementadas para el método numérico que calcula la 
forma matricial, la matriz de Choi y evalúa la positividad de las operaciones PCE.
Además se agregó una función para mostrar una figura con la esfera de Bloch
transformada por la operación PCE.

__*Archivos*__
* quantum_JA.m: paquete de Mathematica con rutinas para calcular el superoperador
  de una operación PCE, hacer el _reshuffle_ de una matriz, evaluar si una matriz 
  es positiva semidefinida y graficar la acción de una operación de 1 qubit sobre 
  la esfera de Bloch.
* PCE_example: cuaderno de Mathematica en el que mostramos un ejemplo de cómo 
  construir las operaciones PCE de 1 qubit, evaluar su completa positividad y 
  mostrar cómo transforman a la esfera de Bloch. 

##Proyecto de tesis

__*Descripción*__: 
Se diseñaron rutinas en Wolfram para evaluar numéricamente
la completa positividad y con ellas buscar los canales cuánticos 
PCE de 2 y 3 qubits. 

__*Archivos*__:
* pce.m: Paquete de Mathematica con rutinas relacionadas a las operaciones PCE y
  canales cuánticos.
* tesis.nb: Cuaderno de Mathematica en el que mostramos cómo utilizar las
  rutinas de pce.m y otras para calcular canales PCE. Adicionalmente, 
  mostramos la manera en la que reproducimos nuestros resultados 
  sobre los canales PCE de 2 y 3 qubits. 
  
##**Instrucciones para ejecutar los .nb:**
1. Clone este repositorio en su computadora. 
2. Según esté interesado en el proyecto de prácticas o tesis abra 
   el archivo .nb correspondiente en Wolfram Mathematica y evalúe 
   el cuaderno completo. 
