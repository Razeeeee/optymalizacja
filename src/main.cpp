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
		lab1();
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

}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
