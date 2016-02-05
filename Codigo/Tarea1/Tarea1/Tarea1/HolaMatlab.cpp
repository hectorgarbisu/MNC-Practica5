/*
Fichero prototipo de función MEX
Métodos Numéricos para la Computación
*/
// incluir todas las cabeceras necesarias
#include <cstdio>
#include <cstring>
#include <cmath>
#include <mkl.h>
#ifdef _CHAR16T
#define CHAR16_T
#endif
#include <mex.h>
void mexFunction(int nlhs, mxArray *plhs[],
	int nrhs, const mxArray *prhs[]){
	/*Declarar las variables locales*/
	mexPrintf("hios\n"); //Tarea1 termina aqui
	double *A, *B, determinante;
	int *pivot, info, Nfilas, Ncolumnas;
	/*Insertar el código */
	if (nrhs != 1){ // nº args diferente de 1
		mexErrMsgTxt("Error. myla, Debe tener un arg de entrada");
	}
	if (!mxIsNumeric(prhs[0])){
		mexErrMsgTxt("Error. El argumento de entrada debe ser una matriz");
	}
	Nfilas = mxGetM(prhs[0]);
	Ncolumnas = mxGetN(prhs[0]);
	if (Nfilas != Ncolumnas){
		mexErrMsgTxt("Error. La matriz debe ser cuadrada");
	}
	if (Nfilas == 0){
		mexErrMsgTxt("Error. La matriz debe no ser vacía");
	}
	if (nlhs > 2){
		mexErrMsgTxt("Error. Debe haber uno o dos args de salida");
	}
	// copia de las variables
	A = mxGetPr(prhs[0]);
	B = (double *)mkl_malloc(Nfilas*Ncolumnas*sizeof(double), 64);
	memcpy(B, A, Nfilas*Ncolumnas*sizeof(double));
	pivot = (int *)mkl_malloc(Nfilas*sizeof(int), 32);
	//procesos computacionales
	info = LAPACKE_dgetrf(LAPACK_COL_MAJOR, Nfilas, Ncolumnas, B, Ncolumnas, pivot);
	determinante = 1.0;
	for (int i = 0; i < Nfilas; i++){
		if (pivot[i] != (i+1)){
			determinante *= -B[i*Ncolumnas + i];
		}
		else{
			determinante *= B[i*Ncolumnas + i];
		}
	}
	// crear los resultados de salida
	plhs[0] = mxCreateDoubleScalar(determinante);
	if (nlhs == 2){
		if (fabs(determinante) < 1.0e-8){
			mexWarnMsgTxt("Matriz singular o casi singular");
		}
		LAPACKE_dgetri(LAPACK_COL_MAJOR, Nfilas, B, Ncolumnas, pivot);
		plhs[1] = mxCreateDoubleMatrix(Nfilas, Ncolumnas, mxREAL);
		double *C = mxGetPr(plhs[1]);
		memcpy(C, B, Nfilas*Ncolumnas*sizeof(double));
	}
	mkl_free(pivot);
	mkl_free(B);
}


