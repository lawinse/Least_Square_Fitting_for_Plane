#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#define MAX 10
#define NMAX 100

double alCo(double **mat, int rank, int row, int col);
double cofactor(double **mat, int rank, int row, int col);


void inverse(double **mat1, double **mat2, int n, double d) {
	int i,j;
	for (i=0; i<n; ++i)
		mat2[i] = new double[n];
	for (i=0; i<n; ++i)
		for (j=0; j<n; ++j)
			mat2[i][j] = alCo(mat1,n,i,j)/d;
}

double determinant(double ** mat, int n) {
	double ret = 0, tmp;
	int i;
	if (n == 1) return mat[0][0];
	for (i=0; i<n; ++i) {
		tmp = alCo(mat,n,n-1,i);
		ret += (mat[n-1][i])*tmp;
	}
	return ret;
}

double alCo(double **mat, int rank, int row, int col) {
	if ((row + col) % 2 == 0) return cofactor(mat, rank, row, col);
	else return (-1.0)*cofactor(mat,rank, row, col); 
}

double cofactor(double **mat, int rank, int row, int col) {
	double ret;
	int i,j;
	double *smallmatr[MAX-1];
	for (i=0; i<rank-1; ++i) smallmatr[i] = new double[rank-1];
	for (i=0; i<row; ++i)
		for (j=0; j<col; ++j) smallmatr[i][j] = mat[i][j];
	for (i=row; i<rank-1; ++i)
		for (j=0; j<col; ++j) smallmatr[i][j] = mat[i+1][j];
	for (i=0; i<row; ++i)
		for (j=col; j<rank-1; ++j) smallmatr[i][j] = mat[i][j+1];

	for (i=row; i<rank-1; ++i)
		for (j=col; j<rank-1; ++j)
			smallmatr[i][j] = mat[i+1][j+1];
	ret = determinant(smallmatr,rank-1);
	for (i=0; i<rank-1; ++i) delete [] smallmatr[i];
	return ret;
}

void standard(double &a, double &b, double &c) {
	double mo = sqrt(a*a + b*b + c*c);
	if (mo == 0) return;
	a /= mo;
	b /= mo;
	c /= mo;
}

int main(int argc, char const *argv[])
{	
	srand(time(0));
	double arr[NMAX][3],Y[3];
	double A,B,C;
	A = B = C = 0.0;
	memset(arr, 0, sizeof(arr));
	memset(Y, 0, sizeof(Y));
	for (int i=0; i<NMAX; ++i){
		printf("(");
		for (int j=0; j<3; ++j)	arr[i][j] = (-1+2*(bool)(rand()&1))*(double)(rand()%1000)/500, printf("%.2f, ", arr[i][j]);
		printf("\b\b)\n");
	}
	
	double *Matrix[3], *IMatrix[3];
	for (int i=0; i<3; ++i) {
		Matrix[i] = new double[3];
		memset(Matrix[i], 0, sizeof(double)*3);
		IMatrix[i] = new double[3];
	}

	for (int j=0; j<3; ++j) 
		for (int i=0; i<NMAX; ++i) {
			Matrix[0][j] += arr[i][0] * arr[i][j];
			Matrix[1][j] += arr[i][1] * arr[i][j];
			Matrix[2][j] += arr[i][2] * arr[i][j];
			Y[j] -= arr[i][j];
		}

	double d = determinant(Matrix, 3);
	if (fabs(d) < 0.00001) {
		printf("\nMatrix is singular\n");
		getchar();
		return -1;
	}
	inverse(Matrix, IMatrix, 3, d);
	for (int i=0; i<3; ++i) {
		A += IMatrix[0][i] * Y[i];
		B += IMatrix[1][i] * Y[i];
		C += IMatrix[2][i] * Y[i];
	}

	standard(A,B,C);

	printf("\nA = %5.3f, B = %5.3f, C = %5.3f\n", A,B,C);

	for (int i=0; i<3; ++i) {
		delete [] Matrix[i];
		delete [] IMatrix[i];
	}

	return 0;
}