#include "NumberSolution.h"
#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
using namespace std;


int main(void)
{
	while (true)
	{
		setlocale(LC_ALL, "rus");

		int		Nmax = 10000;
		int		S = 0;
		double	eps = 0.0000000000001;
		double	epsMax = 0;
		int		n = 0, m = 0;
		double** V = NULL;
		double** F = NULL;
		double	a, b, c, d;
		int kostyl;
		int Exit = 1;


		a = 1;
		b = 2;
		c = 2;
		d = 3;

		n = 0;
		m = 0;
		Nmax = 0;

		system("cls");
		cout << "������� ����������� �����:";
		while (n <= 0)
		{
			cout << endl << "n = ";
			cin >> n;
		}
		while (m <= 0)
		{
			cout << endl << "m = ";
			cin >> m;
		}
		while (Nmax <= 0)
		{
			cout << endl << "������� ������������ ���������� ��������: ";
			cin >> Nmax;
			cout << endl;
		}

		V = MemoryAllocator(n + 1, m + 1);
		FillStartSolution(V, n, m, a, b, c, d);
		ZeidelsMethod(V, n, m, a, b, c, d, eps, Nmax, epsMax, S);

		if (n * m < 420)
		{
			cout << "������� �� ��������� ��������" << endl;
			ShowSolution(V, n, m);
			cout << endl;
		}

		cout << "---------------------------------------------" << endl;
		cout << "������� ����������� X: [" << a << "; " << b << "]" << endl;
		cout << "������� ����������� Y: [" << c << "; " << d << "]" << endl;
		cout << "��� ����� �� ��� Ox: h = " << (b - a) / n << endl;
		cout << "��� ����� �� ��� Oy: k = " << (d - c) / m << endl;
		cout << "��������� ��������: " << eps << endl;
		cout << "����������� ��������: " << epsMax << endl;
		cout << "��������� ��������: " << S << endl;
		cout << "������� �������: " << DiscrepancyOfSolution(V, n, m, a, b, c, d) << endl;
		cout << "����� ����������� �������: " << CheckNumSolution(V, n, m, a, b, c, d) << endl;
		cout << "---------------------------------------------" << endl << endl;

		MemoryCleaner(V, n);

		cout << "���� ������ ���������� ������� 1" << endl;
		cout << endl;
		cin >> kostyl;
		if (kostyl != 1)
		{
			break;
		}
	}
}