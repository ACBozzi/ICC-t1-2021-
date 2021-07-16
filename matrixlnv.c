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
#include "utils.h"

int main ()
{
	int tamanho;
	real_t *resi;
	double *tTotal;

	Matriz_t *SL;
	SL = lerMatriz();

	tamanho = SL->n;
	Matriz_t *original;
	original = alocaMatriz(tamanho);
	copiaMatriz(SL,original);

	//IMPRIMIR ORIGINAL
	
	resi = (real_t *) malloc(tamanho*sizeof(real_t));

	tTotal = (double *) malloc(tamanho*sizeof(double));

	for(int i=0; i<tamanho;i++){
		resi[i]=0;
		tTotal[i] = 0;
	}

	
	float teste	= 1.0 / 0;
	printf("%f\n",teste );
	if (isinf(teste)&& isnan(teste)){ 
		printf("Ocorreram perações com zero",teste);
	}

	double time = timestamp();
	decomposicaoLUpivoteamento(SL,tTotal);
    double elapsed = timestamp() - time;
    tTotal[2] += elapsed;

	//double time = timestamp();
	//decomposicaoLU(SL,tTotal);
    //double elapsed = timestamp() - time;
    //tTotal[2] += elapsed;

	normaL2residuo(original,SL,resi); 


	printf("%d\n",tamanho );
	prnMatriz(original); // original é a original
	printf("#\n");
	prnMatriz(SL); //SL agora é a inversa
	printf("###########\n");
	printf("# Tempo Triangularização: %.5g milissegundos. \n",tTotal[2]);
	printf("# Tempo cálculo de Y: %.5g milissegundos. \n",tTotal[0]);
	printf("# Tempo cálculo de X: %.5g milissegundos. \n", tTotal[1]);
	printf("# Norma L2 do Residuo:");
	for(int i=0; i<tamanho;i++){
		printf("%g ",resi[i]);
	}
	printf("\n");
	printf("###########\n");
	
}
