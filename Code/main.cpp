#include <iostream>
#include <cmath>
using namespace std;
double F(int k, double x, int key);
void GetMatr(double** mas, double** p, int i, int j, int m);
double Determinant(double** mas, int m);
void print(double** A, int n);

//Ввод
void Enter(double* x, int m) {
	for (int i = 0; i < m; i++)
		cin >> x[i];
}
//получение матрицы без i-той строки и j-того столбца
void GetMatr(double** mas, double** p, int i, int j, int m) {
	int ki, kj, di = 0, dj;
	for (ki = 0; ki < m - 1; ki++) { // проверка индекса строки
		if (ki == i) di = 1;
		dj = 0;
		for (kj = 0; kj < m - 1; kj++) { // проверка индекса столбца
			if (kj == j) dj = 1;
			p[ki][kj] = mas[ki + di][kj + dj];
		}
	}
}
//нахождение детерминанта
double Determinant(double** mas, int m) {
	int i, j, k, n;
	double** p, d;
	p = new double* [m];
	for (i = 0; i < m; i++)
		p[i] = new double[m];
	j = 0; d = 0;
	k = 1; //(-1) в степени i
	n = m - 1;
	if (m == 2) 
		return(mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]));
	if (m > 2) {
		for (i = 0; i < m; i++) {
			GetMatr(mas, p, i, 0, m);
			d = d + k * mas[i][0] * Determinant(p, n);
			k = -k;
		}
	}
	return(d);
}
//нахождение обратной матрицы
double** invMatr(double** a, int n) {
	double** q = new double* [n];
	double** T = new double* [n];
	double** h = new double* [n];
	double d;
	for (int i = 0; i < n; i++) {
		q[i] = new double[n];
		T[i] = new double[n];
		h[i] = new double[n];
	}

	d = 1 / Determinant(a, n);
	//Транспонирование
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			T[j][i] = a[i][j];
	}

	//Матрица алгебраических дополнений
	cout << endl;
	for (int i = 0; i < n; i++) 
		for (int j = 0; j < n; j++) {
			GetMatr(T, h, i, j, n);
			q[i][j] = d * Determinant(h, 2);
			if ((i + j + 2) % 2 == 1)
				q[i][j] = q[i][j] * (-1);
		}
	return q;
}
//Умножение матриц
double** prodMatr(double** a, double** b, int n)
{
	double** res;
	double sum = 0;
	res = new double* [n];
	int i, j;
	res = new double* [n];
	for (i = 0; i < n; i++)
		res[i] = new double[n];

	for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			res[i][j] = 0;
			for (int k = 0; k < n; k++)
				res[i][j] += a[i][k] * b[k][j];
		}
	}

	return res;
}
//Расчёт коэффициентов
double** AklF(double* x, int n, int m, int key) {
	double** A = new double* [n];
	for (int i = 0; i < n; i++) 
		A[i] = new double[n];
	for (int k = 0; k < n; k++) {// 3 - количество базисных функций
		for (int l = 0; l < n; l++) {
			A[k][l] = 0;
			for (int i = 0; i < m; i++) { // 5 - количество точек аппроксимации 
				A[k][l] += F(k, x[i], key) * F(l, x[i], key);
			}
		}
	}
	return A;
}
double** BkF(double* x, double* y, int n, int m, int key) {
	double** B = new double* [n];
	for (int i = 0; i < n; i++)
		B[i] = new double[n];
	for (int k = 0; k < n; k++) {
		for (int i = 0; i < n; i++) {
			B[k][i] = 0;
		}
	}
	for (int k = 0; k < n; k++) {// 3 - количество базисных функций
		for (int i = 0; i < m; i++) { // 5 - количество точек аппроксимации 
				B[k][0] += y[i] * F(k, x[i], key);
			}
	}
	return B;
}
double F(int k, double x, int key) {
	switch (key)
	{
	case (1):
		if (k == 0)
			return cos(x);
		if (k == 1)
			return sin(x);
		else
			return 1;
		break;
	case (2):
		if (k == 0)
			return 1;
		if (k == 1)
			return (x);
		else
			return (3*pow(x,2)-1);
		break;
	case (3):
		if (k == 0)
			return 1;
		if (k == 1)
			return exp(-x);
		else
			return x;
		break;
	}
}
void print(double** C, int n) {
	for (int i = 0; i < n; i++)
		cout << C[i][0] << " ";
	cout << endl;
}
int main() {
	setlocale(LC_ALL, "Russian");
	int n = 3, m, key;
	double J = 0, max = 0;
	cout << "Введите количество точек апроксимации"
	cin >> m;
	double** A;
	double** B;
	double** C;
	double* x = new double[m];
	double* y = new double[m];
	double* f = new double[m];
	double* delta = new double[m];
	Enter(x, m);
	Enter(y, m);
	while (true) {
		cout << "Выберите набор функций:" << endl
		<< "1)sin, cos, 1 " << endl 
		<< "2)x, 1, 3*x^2 - 1" << endl 
		<< "3)1, e^(-x), x " << endl
		<< "4)Выход" << endl;
		cin >> key;
		if ((key == 1) or (key == 2) or (key == 3)) {
			A = AklF(x, n, m, key);
			B = BkF(x, y, n, m, key);
			A = invMatr(A, n);
			C = prodMatr(A, B, n);
			print(C, n);
			cout << "Значения апросимирующей функции:  ";
			for (int i = 0; i < m; i++) {
				f[i] = C[0][0] * F(0, x[i], key) + C[1][0] * F(1, x[i], key) + C[2][0]* F(2, x[i], key);
				cout << f[i] << " ";
			}
			cout << endl << "Значения отклонения: ";
			for (int i = 0; i < m; i++) {
				delta[i] = abs(f[i] - y[i]);
				cout << delta[i] << " ";
			}
			max = delta[0];
			for (int i = 1; i < m; i++) {
				if (max < delta[i])
					max = delta[i];
			}
			cout << endl << "Максимальное по модулю отклонение:" << max;
			for (int i = 0; i < m; i++)
				J = J + pow(delta[i], 2);
			cout << endl << "Значения критерия аппроксимации: " << J << endl;
		}
		else if (key == 4)
			return 0;
		else
			cout << "Неправильные данные" << endl;
		J = 0;
		system("pause");
		system("cls");
	}
}
/*
Набор 1
5
0 0.78 1.57 2.35 3.14
0 4 5 8 9
Набор 2
 5
 -1.0 -0.6 -0.1 0.2 0.7
 0.4 0.6 1.0 1.3 1.8
 Набор 3
 5
-1.0 -0.5 0.0 0.67 1.0
 1.0 -0.25 0.0 0.25 1.0
 Набор 4
 6
 0.0 0.75 1.57 2.36 3.14 6.28
 -2.0 0.1 0.0 -2.0 -4.0 -6.0
*/
