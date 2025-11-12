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

matrix ff2R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego Lab 2 - ramię robota
{
	matrix y;
	// x(0) = k1, x(1) = k2 - współczynniki wzmocnienia regulatora
	// Warunki początkowe: alpha(0) = 0, dalpha(0) = 0
	matrix Y0 = matrix(2, 1);
	Y0(0) = 0.0;		// alpha(0) = 0 rad
	Y0(1) = 0.0;		// dalpha/dt(0) = 0 rad/s
	
	// Parametry: k1, k2, alpha_ref = pi rad, omega_ref = 0 rad/s
	matrix MT = matrix(2, 1);
	MT(0) = x(0);		// k1
	MT(1) = x(1);		// k2
	
	// Symulacja dla czasu t_end = 100s z krokiem dt = 0.1s
	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, ud1, MT);
	
	int n = get_len(Y[0]);
	
	// Obliczenie funkcji celu: całka z (10*(alpha_ref - alpha(t))^2 + (omega_ref - omega(t))^2 + (M(t))^2) dt
	double integral = 0.0;
	double alpha_ref = M_PI;
	double omega_ref = 0.0;
	
	// Całkowanie metodą prostokątów (lewy punkt)
	for (int i = 0; i < n - 1; ++i)
	{
		double alpha = Y[1](i, 0);
		double omega = Y[1](i, 1);
		
		// Obliczenie momentu siły M(t) dla aktualnego stanu
		double alpha_error = alpha_ref - alpha;
		double omega_error = omega_ref - omega;
		double M = x(0) * alpha_error + x(1) * omega_error;
		
		// Wartość podcałkowa
		double integrand = 10.0 * pow(alpha_error, 2) + pow(omega_error, 2) + pow(M, 2);
		
		// Krok czasowy
		double dt = Y[0](i + 1) - Y[0](i);
		
		// Metoda prostokątów (lewy punkt)
		integral += integrand * dt;
	}
	
	y = integral;
	
	Y[0].~matrix();
	Y[1].~matrix();
	
	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)		// równania różniczkowe dla ramienia robota Lab 2
{
	matrix dY(2, 1);
	
	// Parametry ramienia
	double m_c = 5.0;		// kg - masa ciężarka na platformie
	double m_r = 1.0;		// kg - masa ramienia
	double l = 2.0;			// m - długość ramienia
	double b = 0.25;		// N·m·s - współczynnik tarcia
	
	// Moment bezwładności: I = (1/3)*m_r*l^2 + m_c*l^2
	double I = (1.0 / 3.0) * m_r * pow(l, 2) + m_c * pow(l, 2);
	
	// Aktualny stan
	double alpha = Y(0);	// [rad]
	double omega = Y(1);	// [rad/s]
	
	// Wartości referencyjne
	double alpha_ref = M_PI;	// rad
	double omega_ref = 0.0;		// rad/s
	
	// Współczynniki regulatora z ud2
	double k1 = ud2(0);
	double k2 = ud2(1);
	
	// Moment siły sterujący
	double M = k1 * (alpha_ref - alpha) + k2 * (omega_ref - omega);
	
	// Równania ruchu
	dY(0) = omega;											// dalpha/dt = omega
	dY(1) = (M - b * omega) / I;							// domega/dt = (M - b*omega) / I
	
	return dY;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego Lab 3
{
	matrix y;
	
	// Parametr a z ud1
	double a = ud1(0);
	
	// Typ funkcji kary z ud2(0): 0 = zewnętrzna, 1 = wewnętrzna
	// Współczynnik kary c z ud2(1)
	int penalty_type = (int)ud2(0);
	double c = (get_len(ud2) > 1) ? ud2(1) : 1.0;
	
	double x1 = x(0);
	double x2 = x(1);
	
	// Funkcja celu: f(x1, x2) = sin(pi*sqrt((x1/pi)^2 + (x2/pi)^2)) / (pi*sqrt((x1/pi)^2 + (x2/pi)^2))
	double r_arg = sqrt(pow(x1 / M_PI, 2) + pow(x2 / M_PI, 2));
	double f = 0.0;
	
	if (r_arg < 1e-10)
	{
		// Wartość graniczna dla r → 0: lim_{r→0} sin(pi*r)/(pi*r) = 1
		f = 1.0;
	}
	else
	{
		f = sin(M_PI * r_arg) / (M_PI * r_arg);
	}
	
	// Ograniczenia:
	// g1(x1) = -x1 + 1 <= 0  =>  x1 >= 1
	// g2(x2) = -x2 + 1 <= 0  =>  x2 >= 1
	// g3(x1, x2) = sqrt(x1^2 + x2^2) - a <= 0  =>  sqrt(x1^2 + x2^2) <= a
	
	double g1 = -x1 + 1.0;
	double g2 = -x2 + 1.0;
	double g3 = sqrt(pow(x1, 2) + pow(x2, 2)) - a;
	
	// Funkcja kary S(x)
	double S = 0.0;
	
	if (penalty_type == 0)
	{
		// ZEWNĘTRZNA FUNKCJA KARY
		// Wzór: S(x) = Σ(max(0, g_i(x)))^2
		// 
		// Dla każdego ograniczenia g_i(x) <= 0:
		// - Jeśli g_i(x) <= 0: punkt wewnątrz, max(0, g_i) = 0, brak kary
		// - Jeśli g_i(x) > 0: punkt poza obszarem, max(0, g_i) = g_i, kara = (g_i)^2
		
		S = pow(max(0.0, g1), 2) + pow(max(0.0, g2), 2) + pow(max(0.0, g3), 2);
	}
	else
	{
		// WEWNĘTRZNA FUNKCJA KARY (BARIERA LOGARYTMICZNA)
		// Wzór: S(x) = -Σ(1/g_i(x)) dla g_i(x) < 0
		// 
		// MECHANIZM:
		// - Dla punktów wewnątrz obszaru: g_i < 0
		// - Gdy punkt zbliża się do granicy: g_i → 0⁻, więc 1/g_i → -∞
		// - Stąd: -1/g_i → +∞ (bariera odbijająca od granicy)
		// - Współczynnik c ZMNIEJSZA SIĘ w kolejnych iteracjach (c → 0)
		// - Im mniejsze c, tym słabsza bariera, punkt może być bliżej granicy
		
		double epsilon_barrier = 1e-10;  // minimalna odległość od granicy
		
		// WARUNEK KONIECZNY: punkt MUSI być wewnątrz obszaru (wszystkie g_i < 0)
		if (g1 >= -epsilon_barrier || g2 >= -epsilon_barrier || g3 >= -epsilon_barrier)
		{
			// Punkt jest na granicy lub poza obszarem dopuszczalnym
			// Zwracamy bardzo dużą karę, żeby metoda optymalizacji odrzuciła ten punkt
			S = 1e20;
		}
		else
		{
			// Punkt jest wewnątrz obszaru dopuszczalnego (wszystkie g_i < 0)
			// S(x) = -Σ(1/g_i)
			// Uwaga: g_i < 0, więc 1/g_i < 0, a zatem -1/g_i > 0
			S = -(1.0 / g1 + 1.0 / g2 + 1.0 / g3);
		}
	}
	
	// F(x) = f(x) + c * S(x)
	y = f + c * S;
	
	return y;
}
