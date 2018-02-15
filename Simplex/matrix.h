#pragma once

class Matrix
{
public:
	Matrix(int n, int m, int *arrVec)
	{
		num = n;
		for (int i = 0; i < num; i++)
		{
			vec[i] = arrVec[i];
		}
	}

	Matrix* truncateMatrix(Matrix m);
private:
	int * vec;
	int num;
};

