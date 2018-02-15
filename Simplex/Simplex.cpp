#include <iostream>
#include <vector>
#include "Simplex.h"
#include "MatrixWork.h"
#include <fstream>
using namespace std;


/* Creating N_k index array */
vector<int> createNk(vector<double> x_k)
{
	vector<int> N_k(x_k.size());
	int indexOfElement = 0;
	for (int i = 0; i < x_k.size(); i++)
	{
		if (x_k[i] > 0)
		{
			N_k[indexOfElement] = i;
			indexOfElement++;
		}
	}
	N_k.resize(indexOfElement);
	return N_k;
}

/* Creating L_k index array */
vector<int> createLk(vector<double> x_k)
{
	vector<int> L_k(x_k.size());
	int indexOfElement = 0;
	for (int i = 0; i < x_k.size(); i++)
	{
		if (x_k[i] == 0)
		{
			L_k[indexOfElement] = i;
			indexOfElement++;
		}
	}
	L_k.resize(indexOfElement);
	return L_k;
}

/* Creating vector C_T on Lk index array */
vector<double> createC_T_Lk(vector<int> L_k, vector<double> C)
{
	vector<double> result(C.size(), 0);
	for (int i = 0; i < L_k.size(); i++)
	{
		result[i] = C[L_k[i]];
	}
	return result;
}

/* Creating vector C_T on Nk index array */
vector<double> createC_T_Nk(vector<int> N_k, vector<double> C)
{
	vector<double> result(C.size(), 0);
	for (int i = 0; i < N_k.size(); i++)
	{
		result[i] = C[N_k[i]];
	}
	return result;
}

/* Creating and initializing A[M][Lk] */
double** createA_M_Lk(double** A, vector<int> L_k, int M)
{
	double **result = new double*[M];
	for (int i = 0; i < M; i++)
	{
		result[i] = new double[L_k.size()];
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < L_k.size(); j++)
		{
			result[i][j] = A[i][L_k[j]];
		}
	}
	return result;
}

/* Creating and initializing A[M][Nk] */
double** createA_M_Nk(double** A, vector<int> N_k, int M)
{
	double **result = new double*[M];
	for (int i = 0; i < N_k.size(); i++)
	{
		result[i] = new double[N_k.size()];
	}

	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N_k.size(); j++)
		{
			result[i][j] = A[i][N_k[j]];
		}
	}
	return result;
}

void printMatrix(double** Mat, int N, int M)
{
	for (int i = 0; i < M; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << " " << Mat[i][j];
		}
		cout << endl;
	}
}

vector<double> vector_dot_matrix(vector<double> vec, double** Matrix, int N, int M)
{
	vector<double> result(N);
	for (int i = 0; i < N; i++)
	{
		result[i] = 0;
		for (int j = 0; j < M; j++)
		{
			result[i] += vec[j] * Matrix[j][i];
		}
	}
	return result;
}

vector<double> matrix_dot_vector(double** Matrix, vector<double> vec, int N, int M)
{
	vector<double> result(N);
	for (int i = 0; i < M; i++)
	{
		result[i] = 0;
		for (int j = 0; j < N; j++)
		{
			result[i] += vec[j] * Matrix[i][j];
		}
	}
	return result;
}


vector<double> makeA_jk(double **A, int jk, int M)
{
	vector<double> result(M);
	for (int i = 0; i < M; i++)
	{
		result[i] = A[i][jk];
	}
	return result;
}

vector<double> createDK(vector<double> C_T_Nk, vector<double> C_T_Lk, double** B,double** A_M_Lk, vector<int> N_k, vector<int> L_k, int M)
{
	vector<double> Dk(L_k.size());
	/* Calculating temp vectors for Dk */
	vector<double> temp(N_k.size());
	temp = vector_dot_matrix(C_T_Nk, B, N_k.size(), M);
	vector<double> tempO(L_k.size());
	tempO = vector_dot_matrix(temp, A_M_Lk, L_k.size(), M);
	/* End of temp calculations */
	for (int i = 0; i < L_k.size(); i++)
	{
		Dk[i] = C_T_Lk[i] - tempO[i];
	}
	return Dk;
}

vector<double> Simplex(double **A, vector<double> b, vector<double> C, vector<double> x_k, int N, int M) 
{
	while (true)
	{
		vector<int> N_k(N); /* Index array for elements > 0 */
		vector<int> L_k(N); /* Index array for elements = 0 */
		N_k = createNk(x_k); /* Creating N_k index array */
		L_k = createLk(x_k); /* Creating L_k index array */


		vector<double> C_T_Lk(N, 0); /* Creating zero - vector */
		C_T_Lk = createC_T_Lk(L_k, C); /* Initializing C_t[Lk] */

		vector<double> C_T_Nk(N, 0); /* Creating zero - vector */
		C_T_Nk = createC_T_Nk(N_k, C); /* Initializing C_t[Nk] */

		double** A_M_Lk = createA_M_Lk(A, L_k, M); /* Creating A[M][Lk] */

		double** A_M_Nk = createA_M_Nk(A, N_k, M); /* Creating A[M][Lk] */

		double** B = new double*[N_k.size()]; /* Creating matrix for inversed matrix A */
		for (int i = 0; i < N_k.size(); i++)
		{
			B[i] = new double[N_k.size()];
		}
		inverse(A_M_Nk, B, N_k.size()); /* Inverting matrix A[M][Nk] */




		vector<double> Dk(L_k.size()); /* Creating Dk vector */
		Dk = createDK(C_T_Nk, C_T_Lk, B, A_M_Lk, N_k, L_k, M);	/* Calculating Dk vector */


		int jk = 0; /* Catching jk from L_k */
		for (int i = 0; i < L_k.size(); i++)
		{
			if (Dk[i] < 0)
			{
				jk = L_k[i];
				break;
			}
			else if (i == L_k.size() - 1)
			{
				return x_k; /* OPTIMUS */
			}
		} /* jk caught */

		vector<double> Uk(N); /* Creating Uk[N] vector*/
		vector<double> A_M_jk(M); /* Column on jk position in A */
		vector<double> temp3(M); /* temp vector for multiply
								 */
		Uk[jk] = -1; /* Setting jk pos as -1 */
		A_M_jk = makeA_jk(A, jk, M); /* Creating vector-column */
		temp3 = matrix_dot_vector(B, A_M_jk, N_k.size(), M); /* Calculations */
		for (int i = 0; i < N_k.size(); i++)
		{
			Uk[N_k[i]] = temp3[i];
		} /* Uk[N] has been built */

		/* Finding theta */
		double Theta = 1000;
		for (int i = 0; i < N_k.size(); i++)
		{
			if (Uk[N_k[i]] > 0)
			{
				double temp = x_k[N_k[i]] / Uk[N_k[i]];
				if (temp < Theta)
				{
					Theta = temp;
				}
			}
		} /* Theta found */

		/* Multiplying vector to scalar */
		for (int i = 0; i < N; i++)
		{
			x_k[i] = x_k[i] - Theta * Uk[i];
		}

	}
	return x_k;
}


void readFile(void)
{
	/* Open file */
	ifstream inputFile;

	/* TODO: READ FILE TO */
	const int M = 3;
	const int N = 8;
	vector<double> x_k(N);
	x_k[0] = 0;
	x_k[1] = 0;
	x_k[2] = 0;
	x_k[3] = 0;
	x_k[4] = 0;
	x_k[5] = 5; /* y1 */
	x_k[6] = 7; /* y2 */
	x_k[7] = 2; /* y3 */
	vector<double> b(M);
	b[0] = 2;
	b[1] = 3;
	vector<double> cT(N);
	cT[0] = 0;
	cT[1] = 0;
	cT[2] = 0;
	cT[3] = 0;
	cT[4] = 0;
	cT[5] = 1;
	cT[6] = 1;
	cT[7] = 1;

	double** A_M_N = new double*[M];
	for (int i = 0; i < M; i++)
	{
		A_M_N[i] = new double[N];
	}
	A_M_N[0][0] = 2;
	A_M_N[0][1] = 1;
	A_M_N[0][2] = 1;
	A_M_N[0][3] = 1;
	A_M_N[0][4] = 3;

	A_M_N[0][5] = 1;
	A_M_N[0][6] = 0;
	A_M_N[0][7] = 0;



	A_M_N[1][0] = 3;
	A_M_N[1][1] = 0;
	A_M_N[1][2] = 2;
	A_M_N[1][3] = -1;
	A_M_N[1][4] = 6;

	A_M_N[1][5] = 0;
	A_M_N[1][6] = 1;
	A_M_N[1][7] = 0;

	
	

	A_M_N[2][0] = 1;
	A_M_N[2][1] = 0;
	A_M_N[2][2] = -1;
	A_M_N[2][3] = 2;
	A_M_N[2][4] = 1;

	A_M_N[2][5] = 0;
	A_M_N[2][6] = 0;
	A_M_N[2][7] = 1;

	
	printMatrix(A_M_N, N, M);
	x_k = Simplex(A_M_N, b, cT, x_k, N, M);
	double** A_M_K = new double*[M];
	for (int i = 0; i < M; i++)
	{
		A_M_K[i] = new double[N - M];
	}


	A_M_K[0][0] = 2;
	A_M_K[0][1] = 1;
	A_M_K[0][2] = 1;
	A_M_K[0][3] = 1;
	A_M_K[0][4] = 3;



	A_M_K[1][0] = 3;
	A_M_K[1][1] = 0;
	A_M_K[1][2] = 2;
	A_M_K[1][3] = -1;
	A_M_K[1][4] = 6;




	A_M_K[2][0] = 1;
	A_M_K[2][1] = 0;
	A_M_K[2][2] = -1;
	A_M_K[2][3] = 2;
	A_M_K[2][4] = 1;

	vector<double> cTK(N - M);
	cTK[0] = 0;
	cTK[1] = 0;
	cTK[2] = -3;
	cTK[3] = 2;
	cTK[4] = 1;


	x_k.resize(N - M);
	x_k = Simplex(A_M_K, b, cTK, x_k, N-M, M);
	return;
}

