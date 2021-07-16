//Além disso, ao final de cada bloco (após a matriz inversa) deve constar:
//O tempo de execução da triangularização;
//O tempo de cálculo de Y e X nas fases Ly=b e Ux=y
//O valor da norma L2 do resíduo de A.A' = I, isto é, a norma L2 do resíduo de cada coluna da matriz 
//identidade com a coluna correspondente da matriz A.A', onde A'  é o resultado do cálculo de A-1 pelo Método da Fatoração LU.

//Desenvolvolvedores:
	//Anna Caroline Bozzi - acb17 -  GRR20173532
	//João Vitor Gualberto Figueiredo - jvgf17 - GRR20177037

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matriz.h"

int main ()
{
	int tamanho;
	real_t *res;
	real_t tTotal[2];

	Matriz_t *SL;
	SL = lerMatriz();

	tamanho = SL->n;
	Matriz_t *original;
	original = alocaMatriz(tamanho);
	copiaMatriz(SL,original);

	//IMPRIMIR ORIGINAL
	


	//decomposicaoLUpivoteamento(SL,tTotal);
	decomposicaoLU(SL,tTotal);

	res = normaL2residuo(original,SL); 


	printf("%d\n",tamanho );
	prnMatriz(original); // original é a original
	printf("#\n");
	prnMatriz(SL); //SL agora é a inversa
	printf("###########\n");
	printf("# Tempo Triangularização: %.5g milissegundos. \n",tTotal[0]+tTotal[1]);
	printf("# Tempo cálculo de Y: %.5g milissegundos. \n",tTotal[0]);
	printf("# Tempo cálculo de X: %.5g milissegundos. \n", tTotal[1]);
	printf("# Norma L2 do Residuo: %g\n",res);
	printf("###########\n");
	
}
