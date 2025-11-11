/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udostępniony na licencji CC BY-SA 3.0
Autor: dr inż. Łukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 30.09.2025
*********************************************/

#include"opt_alg.h"
#include<vector>
#include<algorithm>
#include<ctime>
#include<cstdlib>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		lab3();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;									// dok�adno��
	int Nmax = 10000;										// maksymalna liczba wywo�a� funkcji celu
	matrix lb(2, 1, -5), ub(2, 1, 5),						// dolne oraz g�rne ograniczenie
		a(2, 1);											// dok�adne rozwi�zanie optymalne
	solution opt;											// rozwi�zanie optymalne znalezione przez algorytm
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Wahadlo
	Nmax = 1000;											// dok�adno��
	epsilon = 1e-2;											// maksymalna liczba wywo�a� funkcji celu
	lb = 0, ub = 5;											// dolne oraz g�rne ograniczenie
	double teta_opt = 1;									// maksymalne wychylenie wahad�a
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo�anie procedury optymalizacji
	cout << opt << endl << endl;							// wypisanie wyniku
	solution::clear_calls();								// wyzerowanie licznik�w

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	ofstream Sout("symulacja_lab0.csv");					// definiujemy strumie� do pliku .csv
	Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
	Sout.close();											// zamykamy strumie�
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
}

void lab1()
{
	/*
	Link do excela:
	https://docs.google.com/spreadsheets/d/1cG12sYTO_ambfhdUiorknPooRxzYnViA0UMpRMmPSZA/edit?usp=sharing
	Minima:
	x=62.7
	x=0
	*/

	solution::clear_calls();
	double x0 = 45.0;
	double d = 1.0;
	double alpha[3] = { 1.2, 1.5, 1.7 };
	int Nmax = 100;
	double epsilon = 1e-2;
	double* interval = nullptr;

	cout << "Wstępne szacowanie przedziału poszukań:\n";
	try
	{
		interval = expansion(ff1T, x0, d, alpha[0], Nmax);
		cout << "alpha = " << alpha[2] << ": [" << interval[0] << ", " << interval[1] << "]\n";
		cout << "Liczba wywolan funkcji: " << solution::f_calls << "\n\n";
		solution::clear_calls();
		interval = expansion(ff1T, x0, d, alpha[1], Nmax);
		cout << "alpha = " << alpha[1] << ": [" << interval[0] << ", " << interval[1] << "]\n";
		cout << "Liczba wywolan funkcji: " << solution::f_calls << "\n\n";
		solution::clear_calls();
		interval = expansion(ff1T, x0, d, alpha[2], Nmax);
		cout << "alpha = " << alpha[0] << ": [" << interval[0] << ", " << interval[1] << "]\n";
		cout << "Liczba wywolan funkcji: " << solution::f_calls << "\n\n";
		solution::clear_calls();
	}
	catch (string EX_INFO)
	{
		throw EX_INFO;
	}

	cout << "Poszukiwanie minimum metodą Fibonacciego:\n";
	try
	{
		solution opt = fib(ff1T, interval[0], interval[1], epsilon);
		cout << opt << endl;
		solution::clear_calls();
	}
	catch (string EX_INFO)
	{
		throw EX_INFO;
	}

	cout << "Poszukiwanie minimum metodą Lagrange'a:\n";
	double gamma = 1e-4;
	try
	{
		solution opt = lag(ff1T, interval[0], interval[1], epsilon, gamma, Nmax);
		cout << opt << endl;
		solution::clear_calls();
	}
	catch (string EX_INFO)
	{
		throw EX_INFO;
	}

	// Dodatkowa analiza z losowymi x0
	cout << "\nAnaliza z losowymi punktami startu:\n";
	
	// Otwieramy plik CSV do zapisu wyników
	// Kolumny: x0, a_expansion, b_expansion, expansion_calls, fib_x, fib_y, fib_calls, fib_minimum_type, lag_x, lag_y, lag_calls, lag_minimum_type
	ofstream csvFile("../data/wyniki_lab1_fun_test.csv");
	
	// Globalne minimum znajduje się w x=62.7
	double global_min_x = 62.7;
	double tolerance = 5.0; // tolerancja do określenia czy minimum jest globalne
	
	// Dla każdego alpha
	for (int alpha_idx = 0; alpha_idx < 3; alpha_idx++)
	{
		cout << "Analiza dla alpha = " << alpha[alpha_idx] << "\n";
		
		// Generowanie 100 losowych x0
		vector<double> x0_values;
		srand(time(nullptr) + alpha_idx); // różne seed dla każdego alpha
		
		for (int i = 0; i < 100; i++)
		{
			double random_x0 = (rand() / (double)RAND_MAX) * 200.0 - 100.0; // losowe x0 z zakresu [-100, 100]
			x0_values.push_back(random_x0);
		}
		
		// Sortowanie x0
		sort(x0_values.begin(), x0_values.end());
		
		// Dla każdego x0
		for (int i = 0; i < 100; i++)
		{
			double current_x0 = x0_values[i];
			
			try
			{
				// Ekspansja
				solution::clear_calls();
				double* current_interval = expansion(ff1T, current_x0, d, alpha[alpha_idx], Nmax);
				int expansion_calls = solution::f_calls;
				
				// Fibonacci
				solution::clear_calls();
				solution fib_opt = fib(ff1T, current_interval[0], current_interval[1], epsilon);
				int fib_calls = solution::f_calls;
				double fib_x = m2d(fib_opt.x);
				string fib_min_type = (abs(fib_x - global_min_x) < tolerance) ? "globalne" : "lokalne";
				
				// Lagrange
				solution::clear_calls();
				solution lag_opt = lag(ff1T, current_interval[0], current_interval[1], epsilon, gamma, Nmax);
				int lag_calls = solution::f_calls;
				double lag_x = m2d(lag_opt.x);
				string lag_min_type = (abs(lag_x - global_min_x) < tolerance) ? "globalne" : "lokalne";
				
				// Zapisywanie do CSV
				csvFile << current_x0 << "," 
						<< current_interval[0] << "," 
						<< current_interval[1] << "," 
						<< expansion_calls << ","
						<< fib_x << "," 
						<< m2d(fib_opt.y) << "," 
						<< fib_calls << "," 
						<< fib_min_type << ","
						<< lag_x << "," 
						<< m2d(lag_opt.y) << "," 
						<< lag_calls << "," 
						<< lag_min_type << "\n";
				
				delete[] current_interval;
			}
			catch (string EX_INFO)
			{
				cout << "Błąd dla x0 = " << current_x0 << ", alpha = " << alpha[alpha_idx] << ": " << EX_INFO << "\n";
				// W przypadku błędu zapisujemy pusty wiersz lub pomijamy
				csvFile << current_x0 << ",ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR\n";
			}
		}
	}
	
	csvFile.close();
	cout << "Wyniki zapisane do pliku wyniki_lab1_fun_test.csv\n";

	cout << "\n=== PROBLEM RZECZYWISTY - OPTYMALIZACJA ZBIORNIKOW ===\n";
	solution::clear_calls();
	
	double x0_real = 50.0;
	double d_real = 5.0;
	double alpha_real = 1.5;
	int Nmax_real = 1000;
	double epsilon_real = 1e-2;
	
	double* interval_real = nullptr;
	try
	{
		interval_real = expansion(ff1R, x0_real, d_real, alpha_real, Nmax_real);
		cout << "Przedzial poszukiwan: [" << interval_real[0] << ", " << interval_real[1] << "] cm^2\n";
		solution::clear_calls();
	}
	catch (string EX_INFO)
	{
		cout << "Blad podczas ekspansji: " << EX_INFO << "\n";
		return;
	}
	
	solution opt_fib;
	int fib_calls = 0;
	try
	{
		opt_fib = fib(ff1R, interval_real[0], interval_real[1], epsilon_real);
		fib_calls = solution::f_calls;
		solution::clear_calls();
		cout << "Fibonacci - DA: " << m2d(opt_fib.x) << " cm^2, y*: " << m2d(opt_fib.y) << ", wywolania: " << fib_calls << "\n";
	}
	catch (string EX_INFO)
	{
		cout << "Blad Fibonacci: " << EX_INFO << "\n";
		return;
	}
	
	solution opt_lag;
	int lag_calls = 0;
	double gamma_real = 1e-6;
	try
	{
		opt_lag = lag(ff1R, interval_real[0], interval_real[1], epsilon_real, gamma_real, Nmax_real);
		lag_calls = solution::f_calls;
		solution::clear_calls();
		cout << "Lagrange - DA: " << m2d(opt_lag.x) << " cm^2, y*: " << m2d(opt_lag.y) << ", wywolania: " << lag_calls << "\n";
	}
	catch (string EX_INFO)
	{
		cout << "Blad Lagrange: " << EX_INFO << "\n";
		return;
	}
	
	cout << "\n=== TABELA 3 ===\n";
	cout << "DA_fib\ty*_fib\tcalls_fib\tDA_lag\ty*_lag\tcalls_lag\n";
	cout << m2d(opt_fib.x) << "\t" << m2d(opt_fib.y) << "\t" << fib_calls << "\t"
		 << m2d(opt_lag.x) << "\t" << m2d(opt_lag.y) << "\t" << lag_calls << "\n";
	
	cout << "\n=== SYMULACJE ===\n";
	cout << "Kolumny CSV: t, VA_fib, VA_lag, VB_fib, VB_lag, TB_fib, TB_lag\n";
	
	matrix Y0_fib = matrix(3, 1);
	Y0_fib(0) = 5.0; Y0_fib(1) = 1.0; Y0_fib(2) = 20.0;
	matrix MT_fib = matrix(1, new double[1] { m2d(opt_fib.x) * 1e-4 });
	matrix* Y_fib = solve_ode(df1, 0, 1, 2000, Y0_fib, NAN, MT_fib);
	
	matrix Y0_lag = matrix(3, 1);
	Y0_lag(0) = 5.0; Y0_lag(1) = 1.0; Y0_lag(2) = 20.0;
	matrix MT_lag = matrix(1, new double[1] { m2d(opt_lag.x) * 1e-4 });
	matrix* Y_lag = solve_ode(df1, 0, 1, 2000, Y0_lag, NAN, MT_lag);
	
	ofstream csvSym("../data/symulacja_lab1_real.csv");
	int n_sim = get_len(Y_fib[0]);
	for (int i = 0; i < n_sim; ++i)
	{
		csvSym << Y_fib[0](i, 0) << ","
				<< Y_fib[1](i, 0) << ","
				<< Y_lag[1](i, 0) << ","
				<< Y_fib[1](i, 1) << ","
				<< Y_lag[1](i, 1) << ","
				<< Y_fib[1](i, 2) << ","
				<< Y_lag[1](i, 2) << "\n";
	}
	csvSym.close();
	cout << "Wyniki symulacji zapisane do ../data/symulacja_lab1_real.csv\n";
	
	Y_fib[0].~matrix(); Y_fib[1].~matrix();
	Y_lag[0].~matrix(); Y_lag[1].~matrix();
	delete[] interval_real;
	
	// Dodatkowa symulacja dla DA = 50 cm^2
	cout << "\n=== SYMULACJA DLA DA = 50 cm^2 ===\n";
	matrix Y0_test = matrix(3, 1);
	Y0_test(0) = 5.0; Y0_test(1) = 1.0; Y0_test(2) = 20.0;
	matrix MT_test = matrix(1, new double[1] { 50.0 * 1e-4 }); // 50 cm^2 -> m^2
	matrix* Y_test = solve_ode(df1, 0, 1, 2000, Y0_test, NAN, MT_test);
	
	int n_test = get_len(Y_test[0]);
	double T_max_test = 0;
	for (int i = 0; i < n_test; ++i)
	{
		if (Y_test[1](i, 2) > T_max_test)
			T_max_test = Y_test[1](i, 2);
	}
	cout << "Maksymalna temperatura w zbiorniku B dla DA=50cm^2: " << T_max_test << " C\n";
	
	Y_test[0].~matrix(); Y_test[1].~matrix();
}

void lab2()
{
	// Excel
	// https://docs.google.com/spreadsheets/d/1H7ypVKu4YyRn_saU08J23drAsOx_1Rzn/edit?usp=sharing&ouid=117458587915467310851&rtpof=true&sd=true

	/*
	Testowa funkcja celu Lab 2:
	f(x1, x2) = x1^2 + x2^2 - cos(2.5*pi*x1) - cos(2.5*pi*x2) + 2
	Punkt startowy: x1 ∈ [-1, 1], x2 ∈ [-1, 1]
	Optymalizacja: metoda Hooke'a-Jeevesa i metoda Rosenbrocka
	Minimum globalne: f(0, 0) = 0
	*/
	
	cout << "=== LAB 2: Optymalizacja wielowymiarowa ===\n\n";
	
	// Parametry optymalizacji
	double epsilon = 1e-3;
	int Nmax = 10000;
	double alpha_HJ = 0.5;		// współczynnik redukcji kroku dla Hooke-Jeeves
	double alpha_Rosen = 2.0;	// współczynnik zwiększenia kroku dla Rosenbrocka
	double beta_Rosen = 0.5;	// współczynnik zmniejszenia kroku dla Rosenbrocka
	
	// Różne długości kroku startowego
	double step_sizes[3] = { 0.01, 0.05, 0.075 };
	
	// Globalne minimum: (0, 0) z wartością funkcji = 0
	double global_min_x1 = 0.0;
	double global_min_x2 = 0.0;
	double tolerance = 0.1;		// tolerancja do określenia czy minimum jest globalne
	
	// Pliki CSV dla wyników
	ofstream csv_tabela1("../data/lab2_tabela1.csv");
	ofstream csv_tabela2("../data/lab2_tabela2.csv");
	ofstream csv_rosen("../data/rosen_results.csv");
	
	// Opis kolumn dla tabeli 1 (bez nagłówka w pliku CSV)
	cout << "TABELA 1 - Struktura kolumn (12 kolumn, 300 wierszy):\n";
	cout << "  Kol 1-2:   x1(0), x2(0) - punkt startowy\n";
	cout << "  Kol 3-7:   Hooke-Jeeves -> x1, x2, y, fcalls, is_global\n";
	cout << "  Kol 8-12:  Rosenbrock -> x1, x2, y, fcalls, is_global\n";
	cout << "  Wiersze 1-100:   krok = 0.01\n";
	cout << "  Wiersze 101-200: krok = 0.05\n";
	cout << "  Wiersze 201-300: krok = 0.075\n\n";
	
	srand(time(nullptr));
	
	// Dla każdego rozmiaru kroku
	for (int s_idx = 0; s_idx < 3; s_idx++)
	{
		double step_size = step_sizes[s_idx];
		cout << "Analiza dla kroku startowego s = " << step_size << "\n";
		
		// Statystyki dla średnich - teraz ze WSZYSTKICH optymalizacji
		int hj_global_count = 0, rosen_global_count = 0;
		double hj_fcalls_sum = 0, rosen_fcalls_sum = 0;
		double hj_f_sum = 0, rosen_f_sum = 0;
		double hj_x1_sum = 0, hj_x2_sum = 0;
		double rosen_x1_sum = 0, rosen_x2_sum = 0;
		
		// 100 optymalizacji dla każdej metody
		for (int run = 0; run < 100; run++)
		{
			// Losowy punkt startowy w przedziale [-1, 1] x [-1, 1]
			matrix x0(2, 1);
			x0(0) = (rand() / (double)RAND_MAX) * 2.0 - 1.0;	// x1 ∈ [-1, 1]
			x0(1) = (rand() / (double)RAND_MAX) * 2.0 - 1.0;	// x2 ∈ [-1, 1]
			
			// ===== METODA HOOKE'A-JEEVESA =====
			solution::clear_calls();
			solution opt_hj = HJ(ff2T, x0, step_size, alpha_HJ, epsilon, Nmax);
			int hj_fcalls = solution::f_calls;
			
			// Sprawdzenie czy znaleziono minimum globalne
			double dist_hj = sqrt(pow(opt_hj.x(0) - global_min_x1, 2) + pow(opt_hj.x(1) - global_min_x2, 2));
			bool is_global_hj = (dist_hj < tolerance);
			
			// Statystyki dla średnich - wszystkie próby
			hj_x1_sum += opt_hj.x(0);
			hj_x2_sum += opt_hj.x(1);
			hj_fcalls_sum += hj_fcalls;
			hj_f_sum += opt_hj.y(0);
			if (is_global_hj)
				hj_global_count++;
			
			// ===== METODA ROSENBROCKA =====
			solution::clear_calls();
			matrix s0_rosen(2, 1);
			s0_rosen(0) = step_size;
			s0_rosen(1) = step_size;
			solution opt_rosen = Rosen(ff2T, x0, s0_rosen, alpha_Rosen, beta_Rosen, epsilon, Nmax);
			int rosen_fcalls = solution::f_calls;
			
			// Zapisanie do rosen_results.csv
			csv_rosen << opt_rosen.x(0) << "," << opt_rosen.x(1) << "," << opt_rosen.y(0) << "," << rosen_fcalls << "\n";
			
			// Sprawdzenie czy znaleziono minimum globalne
			double dist_rosen = sqrt(pow(opt_rosen.x(0) - global_min_x1, 2) + pow(opt_rosen.x(1) - global_min_x2, 2));
			bool is_global_rosen = (dist_rosen < tolerance);
			
			// Statystyki dla średnich - wszystkie próby
			rosen_x1_sum += opt_rosen.x(0);
			rosen_x2_sum += opt_rosen.x(1);
			rosen_fcalls_sum += rosen_fcalls;
			rosen_f_sum += opt_rosen.y(0);
			if (is_global_rosen)
				rosen_global_count++;
			
			// Zapisanie do tabeli 1 (jeden wiersz z wynikami obu metod)
			// Format: x1(0), x2(0), HJ_x1, HJ_x2, HJ_y, HJ_fcalls, HJ_global, Rosen_x1, Rosen_x2, Rosen_y, Rosen_fcalls, Rosen_global
			csv_tabela1 << x0(0) << "," << x0(1) << ","
						<< opt_hj.x(0) << "," << opt_hj.x(1) << "," << opt_hj.y(0) << "," << hj_fcalls << ",,"
						<< opt_rosen.x(0) << "," << opt_rosen.x(1) << "," << opt_rosen.y(0) << "," << rosen_fcalls << ",\n";
		}
		
		cout << "  HJ: " << hj_global_count << " optymalizacji znalazło minimum globalne\n";
		cout << "  Rosenbrock: " << rosen_global_count << " optymalizacji znalazło minimum globalne\n\n";
	}
	
	csv_tabela1.close();
	csv_tabela2.close();
	csv_rosen.close();
	
	cout << "Wyniki zapisane do:\n";
	cout << "  - ../data/lab2_tabela1.csv (wszystkie wyniki)\n";
	cout << "  - ../data/lab2_tabela2.csv (średnie dla minimum globalnego)\n";
	cout << "  - ../data/rosen_results.csv (x1, x2, y, f_calls dla Rosena)\n\n";
	
	// ===== WYKRES KONTUROWY - ŚCIEŻKA OPTYMALIZACJI =====
	cout << "Generowanie ścieżki optymalizacji dla wykresu konturowego...\n";
	cout << "TABELA 3 - Struktura kolumn (5 kolumn):\n";
	cout << "  Kol 1:     nr iteracji\n";
	cout << "  Kol 2-3:   Hooke-Jeeves -> x1, x2\n";
	cout << "  Kol 4-5:   Rosenbrock -> x1, x2\n\n";
	
	// Wybrany przypadek: krok s=0.05, punkt startowy bliżej globalnego minimum
	matrix x0_example(2, 1);
	x0_example(0) = 0.3;
	x0_example(1) = 0.2;
	
	// Plik CSV dla ścieżek optymalizacji
	ofstream csv_wykres("../data/lab2_wykres.csv");
	
	// ===== HOOKE-JEEVES Z ZAPISEM HISTORII =====
	solution::clear_calls();
	vector<matrix> history_hj;
	
	// Prosta modyfikacja HJ - zapisujemy punkty bazowe
	solution XB_hj(x0_example);
	XB_hj.fit_fun(ff2T);
	history_hj.push_back(XB_hj.x);
	
	double s_wykres = 0.05;
	int iter = 0;
	int max_iter = 1000;
	
	while (iter < max_iter)
	{
		solution X_hj = HJ_trial(ff2T, XB_hj, s_wykres);
		
		if (X_hj.y < XB_hj.y)
		{
			while (true)
			{
				solution XB_old_hj = XB_hj;
				XB_hj = X_hj;
				history_hj.push_back(XB_hj.x);  // Zapisz nowy punkt bazowy
				
				matrix x_new_hj = 2.0 * XB_hj.x - XB_old_hj.x;
				X_hj.x = x_new_hj;
				X_hj.fit_fun(ff2T);
				
				X_hj = HJ_trial(ff2T, X_hj, s_wykres);
				
				if (X_hj.y >= XB_hj.y)
					break;
			}
		}
		else
		{
			s_wykres = alpha_HJ * s_wykres;
			// Zapisz punkt bazowy nawet gdy redukowaliśmy krok
			history_hj.push_back(XB_hj.x);
		}
		
		if (s_wykres < epsilon || solution::f_calls > Nmax)
			break;
			
		iter++;
	}
	
	// ===== ROSENBROCK Z ZAPISEM HISTORII =====
	solution::clear_calls();
	vector<matrix> history_rosen;
	
	int n = 2;
	matrix l_rosen(n, n);
	for (int i = 0; i < n; ++i)
		l_rosen(i, i) = 1.0;
	
	matrix p_rosen(n, 1), lambda_rosen(n, 1);
	solution X_rosen(x0_example);
	X_rosen.fit_fun(ff2T);
	history_rosen.push_back(X_rosen.x);
	
	matrix s_rosen(n, 1);
	s_rosen(0) = 0.05;
	s_rosen(1) = 0.05;
	
	iter = 0;
	while (iter < max_iter)
	{
		for (int j = 0; j < n; ++j)
		{
			p_rosen(j) = 0;
			lambda_rosen(j) = 0;
		}
		
		bool any_change = false;
		while (true)
		{
			for (int j = 0; j < n; ++j)
			{
				solution X_new_rosen = X_rosen;
				for (int i = 0; i < n; ++i)
					X_new_rosen.x(i) = X_rosen.x(i) + s_rosen(j) * l_rosen(i, j);
				X_new_rosen.fit_fun(ff2T);
				
				if (X_new_rosen.y < X_rosen.y)
				{
					X_rosen = X_new_rosen;
					history_rosen.push_back(X_rosen.x);  // Zapisz każdy poprawny krok
					p_rosen(j) = p_rosen(j) + s_rosen(j);
					lambda_rosen(j) = lambda_rosen(j) + 1;
					s_rosen(j) = alpha_Rosen * s_rosen(j);
					any_change = true;
				}
				else
				{
					s_rosen(j) = -beta_Rosen * s_rosen(j);
					lambda_rosen(j) = lambda_rosen(j) - 1;
				}
			}
			
			if (solution::f_calls > Nmax)
				break;
			
			bool any_lambda_nonzero = false;
			for (int j = 0; j < n; ++j)
			{
				if (lambda_rosen(j) != 0)
				{
					any_lambda_nonzero = true;
					break;
				}
			}
			if (!any_lambda_nonzero)
				break;
		}
		
		// Zapisz punkt po zakończeniu wewnętrznej pętli, jeśli nie było zmian
		if (!any_change)
			history_rosen.push_back(X_rosen.x);
		
		bool change_small = true;
		for (int j = 0; j < n; ++j)
		{
			if (abs(p_rosen(j)) >= epsilon)
			{
				change_small = false;
				break;
			}
		}
		
		if (change_small || solution::f_calls > Nmax)
			break;
		
		// Gram-Schmidt (uproszczony)
		for (int j = 0; j < n; ++j)
			s_rosen(j) = 0.5 * s_rosen(j);
			
		iter++;
	}
	
	// Zapisanie do CSV
	size_t max_history = max(history_hj.size(), history_rosen.size());
	
	for (size_t i = 0; i < max_history; ++i)
	{
		csv_wykres << i << ",";
		
		if (i < history_hj.size())
			csv_wykres << history_hj[i](0) << "," << history_hj[i](1) << ",";
		else
			csv_wykres << ",,";
			
		if (i < history_rosen.size())
			csv_wykres << history_rosen[i](0) << "," << history_rosen[i](1);
		else
			csv_wykres << ",";
			
		csv_wykres << "\n";
	}
	
	csv_wykres.close();
	
	cout << "Punkt startowy: (" << x0_example(0) << ", " << x0_example(1) << ")\n";
	cout << "HJ: " << history_hj.size() << " punktów bazowych\n";
	cout << "Rosenbrock: " << history_rosen.size() << " punktów bazowych\n";
	cout << "Zapisano do: data/lab2_wykres.csv\n\n";
	
	cout << "=== LAB 2 ZAKOŃCZONE ===\n";
	
	// ================================================================
	// ===== CZĘŚĆ 2: PROBLEM RZECZYWISTY - OPTYMALIZACJA RAMIENIA =====
	// ================================================================
	
	cout << "\n\n=== LAB 2 CZĘŚĆ 2: PROBLEM RZECZYWISTY - RAMIĘ ROBOTA ===\n\n";
	
	/*
	Problem: Optymalizacja współczynników wzmocnienia k1 i k2 regulatora dla ramienia robota
	Parametry:
	- m_r = 1 kg (masa ramienia)
	- m_c = 5 kg (masa ciężarka)
	- l = 2 m (długość ramienia)
	- b = 0.25 N·m·s (współczynnik tarcia)
	- I = (1/3)*m_r*l^2 + m_c*l^2 (moment bezwładności)
	
	Równanie ruchu: I * d²α/dt² + b * dα/dt = M(t)
	Moment siły: M(t) = k1 * (α_ref - α(t)) + k2 * (ω_ref - ω(t))
	
	Funkcja celu: Q(k1, k2) = ∫[0,t_end] (10*(α_ref - α(t))² + (ω_ref - ω(t))² + (M(t))²) dt
	
	Zakres poszukiwań: k1 ∈ [0, 20] Nm, k2 ∈ [0, 20] Nms
	Warunki początkowe: α(0) = 0, dα/dt(0) = 0
	Wartości referencyjne: α_ref = π rad, ω_ref = 0 rad/s
	t_end = 100s, dt = 0.1s
	*/
	
	// Parametry optymalizacji
	double epsilon_real = 1e-2;
	int Nmax_real = 10000;
	
	// WERYFIKACJA POPRAWNOŚCI IMPLEMENTACJI
	cout << "=== WERYFIKACJA POPRAWNOŚCI IMPLEMENTACJI ===\n";
	cout << "Obliczanie wartości funkcji celu dla k1 = 5 Nm, k2 = 5 Nms\n";
	
	matrix k_test(2, 1);
	k_test(0) = 5.0;		// k1 = 5 Nm
	k_test(1) = 5.0;		// k2 = 5 Nms
	
	solution::clear_calls();
	solution sol_test(k_test);
	sol_test.fit_fun(ff2R);
	double Q_test = sol_test.y(0);
	
	cout << "Wartość funkcji celu Q(5, 5) = " << Q_test << "\n";
	cout << "Oczekiwana wartość: Q(5, 5) ≈ 775.229\n";
	
	if (abs(Q_test - 775.229) < 1.0)
	{
		cout << "✓ POPRAWNA IMPLEMENTACJA - różnica: " << abs(Q_test - 775.229) << "\n";
	}
	else
	{
		cout << "✗ UWAGA - różnica: " << abs(Q_test - 775.229) << " (może wymagać weryfikacji)\n";
	}
	
	cout << "Wywołań funkcji celu: " << solution::f_calls << "\n";
	cout << "==========================================\n\n";
	
	solution::clear_calls();
	
	// Długości kroku dla różnych prób
	double step_sizes_real[3] = { 0.5, 1.0, 2.0 };
	
	// Tabela 3: porównanie metod dla problemu rzeczywistego
	cout << "TABELA 3 - Optymalizacja ramienia robota dla różnych długości kroku\n";
	cout << "Format: krok, k1_HJ, k2_HJ, Q_HJ, fcalls_HJ, k1_Rosen, k2_Rosen, Q_Rosen, fcalls_Rosen\n\n";
	
	ofstream csv_tabela3("../data/lab2_tabela3_real.csv");
	
	// Punkt startowy dla optymalizacji (środek przedziału)
	matrix x0_real(2, 1);
	x0_real(0) = 10.0;		// k1_start = 10 Nm
	x0_real(1) = 10.0;		// k2_start = 10 Nms
	
	cout << "Punkt startowy: k1 = " << x0_real(0) << " Nm, k2 = " << x0_real(1) << " Nms\n\n";
	
	// Najlepsze rozwiązanie (do późniejszych symulacji)
	solution best_hj, best_rosen;
	bool best_hj_set = false, best_rosen_set = false;
	
	for (int s_idx = 0; s_idx < 3; s_idx++)
	{
		double step_size_real = step_sizes_real[s_idx];
		cout << "Optymalizacja dla kroku s = " << step_size_real << "\n";
		
		// ===== METODA HOOKE'A-JEEVESA =====
		solution::clear_calls();
		solution opt_hj_real = HJ(ff2R, x0_real, step_size_real, 0.5, epsilon_real, Nmax_real);
		int hj_fcalls_real = solution::f_calls;
		
		cout << "  Hooke-Jeeves:\n";
		cout << "    k1 = " << opt_hj_real.x(0) << " Nm\n";
		cout << "    k2 = " << opt_hj_real.x(1) << " Nms\n";
		cout << "    Q = " << opt_hj_real.y(0) << "\n";
		cout << "    Wywołania funkcji celu: " << hj_fcalls_real << "\n";
		
		// ===== METODA ROSENBROCKA =====
		solution::clear_calls();
		matrix s0_rosen_real(2, 1);
		s0_rosen_real(0) = step_size_real;
		s0_rosen_real(1) = step_size_real;
		solution opt_rosen_real = Rosen(ff2R, x0_real, s0_rosen_real, 2.0, 0.5, epsilon_real, Nmax_real);
		int rosen_fcalls_real = solution::f_calls;
		
		cout << "  Rosenbrock:\n";
		cout << "    k1 = " << opt_rosen_real.x(0) << " Nm\n";
		cout << "    k2 = " << opt_rosen_real.x(1) << " Nms\n";
		cout << "    Q = " << opt_rosen_real.y(0) << "\n";
		cout << "    Wywołania funkcji celu: " << rosen_fcalls_real << "\n\n";
		
		// Zapisanie do CSV
		csv_tabela3 << step_size_real << ","
					<< opt_hj_real.x(0) << "," << opt_hj_real.x(1) << "," << opt_hj_real.y(0) << "," << hj_fcalls_real << ","
					<< opt_rosen_real.x(0) << "," << opt_rosen_real.x(1) << "," << opt_rosen_real.y(0) << "," << rosen_fcalls_real << "\n";
		
		// Zachowanie najlepszego rozwiązania
		if (!best_hj_set || opt_hj_real.y(0) < best_hj.y(0))
		{
			best_hj = opt_hj_real;
			best_hj_set = true;
		}
		if (!best_rosen_set || opt_rosen_real.y(0) < best_rosen.y(0))
		{
			best_rosen = opt_rosen_real;
			best_rosen_set = true;
		}
	}
	
	csv_tabela3.close();
	cout << "Wyniki zapisane do: ../data/lab2_tabela3_real.csv\n\n";
	
	// ===== SYMULACJE DLA OPTYMALNYCH WARTOŚCI =====
	cout << "=== SYMULACJE ===\n";
	
	// Symulacja dla najlepszego rozwiązania Hooke-Jeeves
	cout << "Symulacja dla optymalnych k1, k2 (Hooke-Jeeves):\n";
	cout << "  k1 = " << best_hj.x(0) << " Nm, k2 = " << best_hj.x(1) << " Nms\n";
	
	matrix Y0_sim(2, 1);
	Y0_sim(0) = 0.0;		// alpha(0) = 0
	Y0_sim(1) = 0.0;		// omega(0) = 0
	
	matrix MT_hj_sim(2, 1);
	MT_hj_sim(0) = best_hj.x(0);
	MT_hj_sim(1) = best_hj.x(1);
	
	matrix* Y_hj_sim = solve_ode(df2, 0, 0.1, 100, Y0_sim, NAN, MT_hj_sim);
	
	// Symulacja dla najlepszego rozwiązania Rosenbrock
	cout << "Symulacja dla optymalnych k1, k2 (Rosenbrock):\n";
	cout << "  k1 = " << best_rosen.x(0) << " Nm, k2 = " << best_rosen.x(1) << " Nms\n";
	
	matrix MT_rosen_sim(2, 1);
	MT_rosen_sim(0) = best_rosen.x(0);
	MT_rosen_sim(1) = best_rosen.x(1);
	
	matrix* Y_rosen_sim = solve_ode(df2, 0, 0.1, 100, Y0_sim, NAN, MT_rosen_sim);
	
	// Zapisanie wyników symulacji do CSV (tylko 5 kolumn: t, alpha_HJ, omega_HJ, alpha_Rosen, omega_Rosen)
	ofstream csv_sim("../data/lab2_symulacja_real.csv");
	
	int n_sim_real = get_len(Y_hj_sim[0]);
	for (int i = 0; i < n_sim_real; ++i)
	{
		csv_sim << Y_hj_sim[0](i) << ","
				<< Y_hj_sim[1](i, 0) << "," << Y_hj_sim[1](i, 1) << ","
				<< Y_rosen_sim[1](i, 0) << "," << Y_rosen_sim[1](i, 1) << "\n";
	}
	
	csv_sim.close();
	cout << "\nWyniki symulacji zapisane do: ../data/lab2_symulacja_real.csv\n";
	cout << "Kolumny: t, alpha_HJ, omega_HJ, alpha_Rosen, omega_Rosen\n";
	
	// Zwolnienie pamięci
	Y_hj_sim[0].~matrix(); Y_hj_sim[1].~matrix();
	Y_rosen_sim[0].~matrix(); Y_rosen_sim[1].~matrix();
	
	cout << "\n=== LAB 2 CZĘŚĆ 2 ZAKOŃCZONA ===\n";
}


void lab3()
{
	// Link do Excela:
	// https://docs.google.com/spreadsheets/d/1vXfR9t_j6LvPTwkt_JZoQdYq9V--8BDT/edit?usp=sharing&ouid=117458587915467310851&rtpof=true&sd=true
	
	/*
	Funkcja testowa Lab 3:
	f(x1, x2) = sin(pi*sqrt((x1/pi)^2 + (x2/pi)^2)) / (pi*sqrt((x1/pi)^2 + (x2/pi)^2))
	
	Ograniczenia:
	g1(x1) = -x1 + 1 <= 0
	g2(x2) = -x2 + 1 <= 0
	g3(x1, x2) = sqrt(x1^2 + x2^2) - a <= 0
	
	Wartości parametru a: 4, 4.4934, 5
	
	Zadanie:
	- 100 optymalizacji dla każdej wartości a
	- Punkt startowy losowany w obszarze dopuszczalnym
	- Epsilon = 1e-3
	- Metoda: Nelder-Mead (simpleks)
	- Uwzględnienie ograniczeń: zewnętrzna i wewnętrzna funkcja kary
	*/
	
	cout << "=== LAB 3: Optymalizacja z ograniczeniami ===\n\n";
	
	// Parametry optymalizacji
	double epsilon = 1e-3;
	int Nmax = 10000;
	
	// Wartości parametru a
	double a_values[3] = { 4.0, 4.4934, 5.0 };
	
	// Plik CSV dla wyników (Tabela 1)
	ofstream csv_tabela1("data/lab3_tabela1.csv");
	
	// Nagłówek (opcjonalnie, ale ułatwia analizę)
	// csv_tabela1 << "x1_0,x2_0,x1_ext,x2_ext,r_ext,y_ext,fcalls_ext,x1_int,x2_int,r_int,y_int,fcalls_int\n";
	
	cout << "TABELA 1 - Struktura kolumn (12 kolumn, 300 wierszy):\n";
	cout << "  Kol 1-2:   x1(0), x2(0) - punkt startowy\n";
	cout << "  Kol 3-7:   Zewnętrzna funkcja kary -> x1*, x2*, r*, y*, fcalls\n";
	cout << "  Kol 8-12:  Wewnętrzna funkcja kary -> x1*, x2*, r*, y*, fcalls\n";
	cout << "  Wiersze 1-100:   a = 4.0\n";
	cout << "  Wiersze 101-200: a = 4.4934\n";
	cout << "  Wiersze 201-300: a = 5.0\n\n";
	
	srand(time(nullptr));
	
	// Dla każdej wartości parametru a
	for (int a_idx = 0; a_idx < 3; a_idx++)
	{
		double a = a_values[a_idx];
		cout << "Optymalizacja dla a = " << a << "\n";
		
		// 100 optymalizacji dla każdej wartości a
		for (int run = 0; run < 100; run++)
		{
			// Losowanie punktu startowego w obszarze dopuszczalnym
			// Ograniczenia: x1 >= 1, x2 >= 1, sqrt(x1^2 + x2^2) <= a
			// Losujemy punkt w pierścieniu: max(sqrt(2), 1.1) <= r <= a - 0.1
			
			double r_min = max(sqrt(2.0), 1.1);	// minimalna odległość od początku
			double r_max = a - 0.1;				// maksymalna odległość (z marginesem)
			
			if (r_min >= r_max)
			{
				// Jeśli przedział jest pusty lub zbyt wąski, użyjmy wartości domyślnej
				r_min = 1.5;
				r_max = a - 0.1;
			}
			
			double r_start = r_min + (r_max - r_min) * (rand() / (double)RAND_MAX);
			double theta_start = M_PI / 4.0 + (M_PI / 4.0) * (rand() / (double)RAND_MAX);	// kąt w przedziale [45°, 90°]
			
			matrix x0(2, 1);
			x0(0) = r_start * cos(theta_start);
			x0(1) = r_start * sin(theta_start);
			
			// Zapewnienie że x1 >= 1 i x2 >= 1
			if (x0(0) < 1.0) x0(0) = 1.0 + 0.1 * (rand() / (double)RAND_MAX);
			if (x0(1) < 1.0) x0(1) = 1.0 + 0.1 * (rand() / (double)RAND_MAX);
			
			// ===== ZEWNĘTRZNA FUNKCJA KARY =====
			solution::clear_calls();
			
			// Parametry funkcji kary zewnętrznej
			double c_ext = 1.0;			// początkowa wartość współczynnika kary
			double dc_ext = 2.0;		// współczynnik zwiększania kary
			
			matrix ud1_ext(1, 1);
			ud1_ext(0) = a;				// przekazujemy parametr a
			
			matrix ud2_ext(1, 1);
			ud2_ext(0) = 0;				// 0 = zewnętrzna funkcja kary
			
			solution opt_ext = pen(ff3T, x0, c_ext, dc_ext, epsilon, Nmax, ud1_ext, ud2_ext);
			int fcalls_ext = solution::f_calls;
			
			// Obliczenie odległości od początku układu współrzędnych
			double r_ext = sqrt(pow(opt_ext.x(0), 2) + pow(opt_ext.x(1), 2));
			
			// Obliczenie prawdziwej wartości funkcji celu (bez kary)
			matrix ud2_true(2, 1);
			ud2_true(0) = 0;
			ud2_true(1) = 0;	// c = 0, więc nie ma kary
			solution opt_ext_true(opt_ext.x);
			opt_ext_true.fit_fun(ff3T, ud1_ext, ud2_true);
			
			// ===== WEWNĘTRZNA FUNKCJA KARY =====
			solution::clear_calls();
			
			// Parametry funkcji kary wewnętrznej
			double c_int = 10.0;		// początkowa wartość współczynnika kary
			double dc_int = 0.5;		// współczynnik zmniejszania kary
			
			matrix ud1_int(1, 1);
			ud1_int(0) = a;				// przekazujemy parametr a
			
			matrix ud2_int(1, 1);
			ud2_int(0) = 1;				// 1 = wewnętrzna funkcja kary
			
			solution opt_int = pen(ff3T, x0, c_int, dc_int, epsilon, Nmax, ud1_int, ud2_int);
			int fcalls_int = solution::f_calls;
			
			// Obliczenie odległości od początku układu współrzędnych
			double r_int = sqrt(pow(opt_int.x(0), 2) + pow(opt_int.x(1), 2));
			
			// Obliczenie prawdziwej wartości funkcji celu (bez kary)
			solution opt_int_true(opt_int.x);
			opt_int_true.fit_fun(ff3T, ud1_int, ud2_true);
			
			// Zapisanie do CSV
			csv_tabela1 << x0(0) << "," << x0(1) << ","
						<< opt_ext.x(0) << "," << opt_ext.x(1) << "," << r_ext << "," << opt_ext_true.y(0) << "," << fcalls_ext << ","
						<< opt_int.x(0) << "," << opt_int.x(1) << "," << r_int << "," << opt_int_true.y(0) << "," << fcalls_int << "\n";
			
			// Postęp co 25 iteracji
			if ((run + 1) % 25 == 0)
			{
				cout << "  Ukończono " << (run + 1) << "/100 optymalizacji\n";
			}
		}
		
		cout << "  Zakończono optymalizacje dla a = " << a << "\n\n";
	}
	
	csv_tabela1.close();
	
	cout << "Wyniki zapisane do: ../data/lab3_tabela1.csv\n";
	cout << "\n=== LAB 3 ZAKOŃCZONE ===\n";
}void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
