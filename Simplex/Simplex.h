#pragma once
#include <vector>
using namespace std;
void readFile(void);
vector<int> createNk(vector<double> x_k);
vector<int> createLk(vector<double> x_k);
vector<double> createC_T_Lk(vector<int> L_k, vector<double> C);
vector<double> createC_T_Nk(vector<int> N_k, vector<double> C);
double** createA_M_Lk(double** A, vector<int> L_k, int M);
double** createA_M_Nk(double** A, vector<int> N_k, int M);
void printMatrix(double** Mat, int N, int M);
vector<double> vector_dot_matrix(vector<double> vec, double** Matrix, int N, int M);
vector<double> makeA_jk(double **A, int jk, int M);
vector<double> Simplex(double **A, vector<double> b, vector<double> C, vector<double> x_k, int N, int M);