#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;
	y = -cos(0.1 * x(0)) * exp(-pow(0.1*x(0)-2*M_PI, 2)) + 0.002 * pow(0.1 * x(0), 2);
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera wartość funkcji celu
	matrix Y0 = matrix(3, 1),								// Y0 zawiera warunki początkowe
		MT = matrix(1, new double[1] { m2d(x) * 1e-4 });	// MT zawiera pole przekroju DA w m^2
	Y0(0) = 5.0;											// VA_start = 5 m^3
	Y0(1) = 1.0;											// VB_start = 1 m^3  
	Y0(2) = 20.0;											// TB_start = 20 C
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, MT);	// rozwiązujemy równanie różniczkowe
	int n = get_len(Y[0]);									// długość rozwiązania
	double T_max = 0;										// szukamy maksymalnej temperatury w zbiorniku B
	for (int i = 0; i < n; ++i)
		if (Y[1](i, 2) > T_max)
			T_max = Y[1](i, 2);
	y = abs(T_max - 50.0);									// wartość funkcji celu (minimalizujemy różnicę od 50°C)
	Y[0].~matrix();											// usuwamy z pamięci rozwiązanie RR
	Y[1].~matrix();
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(3, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double PA = 2.0, PB = 1.0, g = 9.81;					// definiujemy parametry modelu
	double a = 0.98, b = 0.63;								// współczynniki lepkości i zwężenia strumienia
	double TA_in = 95.0, TB_in_temp = 20.0;				// temperatury dopływu
	double Fin_B = 10.0 / 1000.0;							// dopływ zewnętrzny do zbiornika B [m^3/s]
	double DB = 36.5665 * 1e-4;							// pole powierzchni odpływu ze zbiornika B [m^2]
	double DA = ud2(0);									// pole powierzchni przekroju DA [m^2]
	double VA = Y(0), VB = Y(1), TB = Y(2);				// aktualny stan: objętości i temperatura
	double hA = VA / PA, hB = VB / PB;						// wysokości słupa wody w zbiornikach
	double Fout_A = (hA > 0) ? a * b * DA * sqrt(2 * g * hA) : 0;		// odpływ ze zbiornika A (prawo Torricellego)
	double Fout_B = (hB > 0) ? a * b * DB * sqrt(2 * g * hB) : 0;		// odpływ ze zbiornika B (prawo Torricellego)
	dY(0) = -Fout_A;										// zmiana objętości w zbiorniku A
	dY(1) = Fout_A + Fin_B - Fout_B;						// zmiana objętości w zbiorniku B
	dY(2) = (abs(VB) < 1e-9) ? 0 : (Fout_A * (TA_in - TB) + Fin_B * (TB_in_temp - TB)) / VB;	// zmiana temperatury w zbiorniku B
	return dY;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego Lab 2
{
	matrix y;
	// f(x1, x2) = x1^2 + x2^2 - cos(2.5*pi*x1) - cos(2.5*pi*x2) + 2
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2.0;
	return y;
}
