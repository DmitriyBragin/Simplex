#include <iostream>
#include <vector>
#include "Simplex.h"
#include "MatrixWork.h"
#include <fstream>
using namespace std;

double eps = 0.000001;

/* Creating N_k index array */
vector<int> createNk(vector<double> x_k)
{
	vector<int> N_k(x_k.size());
	int indexOfElement = 0;
	for (int i = 0; i < x_k.size(); i++)
	{
		if (x_k[i] > eps)
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
		if (x_k[i] < eps)
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
		cout << "current min:";
		double min = x_k[0] * C[0];
		for (int i = 1; i < x_k.size(); i++)
		{
			min += x_k[i] * C[i];
		}
		cout << min << endl;
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



		printMatrix(B, N_k.size(), M);
		vector<double> Dk(L_k.size()); /* Creating Dk vector */
		Dk = createDK(C_T_Nk, C_T_Lk, B, A_M_Lk, N_k, L_k, M);	/* Calculating Dk vector */


		int jk = 0; /* Catching jk from L_k */
		for (int i = 0; i < L_k.size(); i++)
		{
			if (Dk[i] < -eps)
			{

				jk = L_k[i];
				break;
			}
			else if (i == L_k.size() - 1)
			{
				vector<double> tempVec(M);
				tempVec = vector_dot_matrix(C_T_Nk, B, N_k.size(), M);
				vector<double> tempVec2(N);
				tempVec2 = vector_dot_matrix(tempVec, A, N, M);
				vector<double> optimium(N);
				for (int i = 0; i < N; i++)
				{
					optimium[i] = C[i] - tempVec2[i];
				}
				double res = 0;
				for (int i = 0; i < N; i++)
				{
					res += optimium[i] * x_k[i]; /* Checking optimum conditions */
				}

				vector<double> solutionOtherTask(M);
				solutionOtherTask = vector_dot_matrix(C_T_Nk, B, N_k.size(), M);
				double difference = 0;
				for (i = 0; i < N; i++)
				{
					difference += 0.1 * x_k[i];
				}
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
		vector<double> thetaUk(N);
		for (int i = 0; i < N; i++)
		{
			thetaUk[i] = Theta * Uk[i];
		}
		/* Multiplying vector to scalar */
		for (int i = 0; i < N; i++)
		{
			x_k[i] = x_k[i] - Theta * Uk[i];
		}
		vector<double> allowableVector(N);
		vector<double> curiousVector(N);
		allowableVector = matrix_dot_vector(A, x_k, N, M);
		curiousVector = matrix_dot_vector(A, thetaUk, N, M);

		printMatrix(A_M_Nk, N_k.size(), M);
	}
	return x_k;
}



void readFile(void)
{

	/* TODO: READ FILE TO */
	const int M = 4;
	const int N = 15;



	//parseFile(inputFile, N, M);
	vector<double> x_k(N, 0);

	vector<double> b(M);

	vector<double> cT(N);

	double** A_M_N = new double*[M];
	for (int i = 0; i < M; i++)
	{
		A_M_N[i] = new double[N];
	}

	/* Open file */
	ifstream inputFile;
	inputFile.open("InputDataHelper.txt");
	while (!inputFile.eof())
	{
		/* CT[N] */
		for (int i = 0; i < N; i++)
		{
			inputFile >> cT[i];
		}
		/* A[M][N] */
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N; j++)
			{
				inputFile >> A_M_N[i][j];
			}
		}
		/* b[M] */
		for (int i = 0; i < M; i++)
		{
			inputFile >> b[i];
		}
	}

	int k = 0;
	for (int i = N - M; i < N; i++)
	{
		x_k[i] = b[k];
		k++;
	}
	inputFile.close();
	cout << "Matrix A_M_N" << endl;
	printMatrix(A_M_N, N, M);
	x_k = Simplex(A_M_N, b, cT, x_k, N, M);





	double** A_M_K = new double*[M];
	for (int i = 0; i < M; i++)
	{
		A_M_K[i] = new double[N - M];
	}


	vector<double> cTK(N - M);


	x_k.resize(N - M);
	ifstream inputFile1;
	inputFile1.open("InputData.txt");
	while (!inputFile1.eof())
	{
		/* CT[N] */
		for (int i = 0; i < N - M; i++)
		{
			inputFile1 >> cTK[i];
		}
		/* A[M][N] */
		for (int i = 0; i < M; i++)
		{
			for (int j = 0; j < N - M; j++)
			{
				inputFile1 >> A_M_K[i][j];
			}
		}
		/* b[M] */
		for (int i = 0; i < M; i++)
		{
			inputFile1 >> b[i];
		}
	}

	cout << "Matrix A_M_K" << endl;
	printMatrix(A_M_K, N - M, M);
	x_k = Simplex(A_M_K, b, cTK, x_k, N-M, M);


	return;
}

