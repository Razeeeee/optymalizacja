#include "opt_alg.h"
#include <algorithm>
#include <vector>

solution MC(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	// Zmienne wej�ciowe:
	// ff - wska�nik do funkcji celu
	// N - liczba zmiennych funkcji celu
	// lb, ub - dolne i g�rne ograniczenie
	// epslion - zak��dana dok�adno�� rozwi�zania
	// Nmax - maksymalna liczba wywo�a� funkcji celu
	// ud1, ud2 - user data
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N); // losujemy macierz Nx1 stosuj�c rozk�ad jednostajny na przedziale [0,1]
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i); // przeskalowywujemy rozwi�zanie do przedzia�u [lb, ub]
			Xopt.fit_fun(ff, ud1, ud2);							 // obliczmy warto�� funkcji celu
			if (Xopt.y < epsilon)								 // sprawdzmy 1. kryterium stopu
			{
				Xopt.flag = 1; // flaga = 1 ozancza znalezienie rozwi�zanie z zadan� dok�adno�ci�
				break;
			}
			if (solution::f_calls > Nmax) // sprawdzmy 2. kryterium stopu
			{
				Xopt.flag = 0; // flaga = 0 ozancza przekroczenie maksymalne liczby wywo�a� funkcji celu
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution MC(...):\n" + ex_info);
	}
}

double *expansion(matrix (*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0; // 1. i = 0
		solution x0_sol(x0);
		x0_sol.fit_fun(ff, ud1, ud2);
		solution x1(x0 + d); // 2. x(1) = x(0) + d
		x1.fit_fun(ff, ud1, ud2);

		if (x1.y == x0_sol.y) // 3. if f(x(1)) = f(x(0)) then
			return new double[2]{x0_sol.x(0), x1.x(0)}; // 4. return [x(0), x(1)]
		// 5. end if

		if (x1.y > x0_sol.y) // 6. if f(x(1)) > f(x(0)) then
		{
			d = -d; // 7. d = -d
			x1.x = x0 + d; // 8. x(1) = x(0) + d
			x1.fit_fun(ff, ud1, ud2);

			if (x1.y >= x0_sol.y) // 9. if f(x(1)) ≥ f(x(0)) then
			{
				return new double[2]{x1.x(0), x0 - d}; // 10. return [x(1), x(0) - d]
			} // 11. end if
		} // 12. end if

		solution x_current = x1;
		solution x_next;

		while (true) // 13. repeat
		{
			if (solution::f_calls > Nmax) // 14. if fcalls > Nmax
			{
				throw string("Przekroczono maksymalną liczbę wywołań funkcji celu."); // 15. return error
			} // 16. end if

			i = i + 1; // 17. i = i + 1
			x_next.x = x0 + pow(alpha, i) * d; // 18. x(i+1) = x(0) + α^i·d
			x_next.fit_fun(ff, ud1, ud2);

			if (x_current.y <= x_next.y) // 19. until f(x(i)) ≤ f(x(i+1))
				break;

			x_current = x_next;
		}

		double *p = new double[2];
		double x_previous = (i == 1) ? x0 : x0 + pow(alpha, i - 1) * d;

		if (d > 0) // 20. if d > 0
		{
			p[0] = x_previous; // 21. return [x(i-1), x(i+1)]
			p[1] = m2d(x_next.x(0));
		}
		else
		{
			p[0] = m2d(x_next.x(0)); // 23. return [x(i+1), x(i-1)]
			p[1] = x_previous;
		}

		return p;
	}
	catch (string ex_info)
	{
		throw("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		// 1. znajdź najmniejszą liczbę k spełniającą nierówność φ^k > (b - a) / ε
		const double phi = (1.0 + sqrt(5.0)) / 2.0; // złoty podział
		int k = 1;
		while (pow(phi, k) <= (b - a) / epsilon)
		{
			k++;
		}

		solution A(a), B(b); // 2. a(0) = a, b(0) = b
		// 3. c(0) = b(0) - φ^(k-1) / φ^k·(b(0) - a(0))
		solution C(B.x(0) - pow(phi, k - 1) / pow(phi, k) * (B.x(0) - A.x(0)));
		solution D(A.x(0) + B.x(0) - C.x(0)); // 4. d(0) = a(0) + b(0) - c(0)
		C.fit_fun(ff, ud1, ud2);
		D.fit_fun(ff, ud1, ud2);

		for (int i = 0; i <= k - 3; ++i) // 5. for i = 0 to k – 3 do
		{
			if (C.y < D.y) // 6. if f(c(i)) < f(d(i)) then
			{
				// A(i+1) = A(i) - niezmieniane // 7. a(i+1) = a(i)
				B = D; // 8. b(i+1) = d(i)
			}
			else // 9. else
			{
				// B(i+1) = B(i) - niezmieniane // 10. b(i+1) = b(i)
				A = C; // 11. a(i+1) = c(i)
			} // 12. end if

			// 13. c(i+1) = b(i+1) - φ^(k-i-2) / φ^(k-i-1)·(b(i+1) - a(i+1))
			C.x = B.x - pow(phi, k - i - 2) / pow(phi, k - i - 1) * (B.x - A.x);
			C.fit_fun(ff, ud1, ud2);
			
			D.x = A.x + B.x - C.x; // 14. d(i+1) = a(i+1) + b(i+1) - c(i+1)
			D.fit_fun(ff, ud1, ud2);
		} // 15. end for

		solution sol = C; // 16. return x* = c(i+1)
		sol.flag = 1;
		return sol;
	}
	catch (string ex_info)
	{
		throw("solution fib(...):\n" + ex_info);
	}
}

solution lag(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		int i = 0; // 1. i = 0
		solution A(a), B(b), C((a + b) / 2.0); // 2. a(0) = a, b(0) = b, c(0) = c
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);

		solution D;
		double d_previous = NAN;

		while (true) // 3. repeat
		{
			if (solution::f_calls > Nmax) // 36. if fcalls > Nmax then
			{
				D.flag = 0; // 37. return error
				return D;
			} // 38. end if

			// 4. l = f(a(i))((b(i))^2 - (c(i))^2) + f(b(i))((c(i))^2 - (a(i))^2) + f(c(i))((a(i))^2 - (b(i))^2)
			double l = A.y(0) * (pow(B.x(0), 2) - pow(C.x(0), 2)) + B.y(0) * (pow(C.x(0), 2) - pow(A.x(0), 2)) + C.y(0) * (pow(A.x(0), 2) - pow(B.x(0), 2));
			// 5. m = f(a(i))(b(i) - c(i)) + f(b(i))(c(i) - a(i)) + f(c(i))(a(i) - b(i))
			double m = A.y(0) * (B.x(0) - C.x(0)) + B.y(0) * (C.x(0) - A.x(0)) + C.y(0) * (A.x(0) - B.x(0));

			if (m <= 0) // 6. if m ≤ 0
			{
				D.flag = 0; // 7. return error
				return D;
			} // 8. end if

			double d_val = 0.5 * l / m; // 9. d(i) = 0.5 * l / m
			D.x = d_val;
			D.fit_fun(ff, ud1, ud2);

			if (d_val > A.x(0) && d_val < C.x(0)) // 10. if a(i) < d(i) < c(i) then
			{
				if (D.y < C.y) // 11. if f(d(i)) < f(c(i)) then
				{
					// A(i+1) = A(i) // 12. a(i+1) = a(i)
					B = C; // 14. b(i+1) = c(i)
					C = D; // 13. c(i+1) = d(i)
				}
				else // 15. else
				{
					A = D; // 16. a(i+1) = d(i)
					// C(i+1) = C(i) // 17. c(i+1) = c(i)
					// B(i+1) = B(i) // 18. b(i+1) = b(i)
				} // 19. end if
			}
			else if (d_val > C.x(0) && d_val < B.x(0)) // 20. else 21. if c(i) < d(i) < b(i) then
			{
				if (D.y < C.y) // 22. if f(d(i)) < f(c(i)) then
				{
					A = C; // 23. a(i+1) = c(i)
					C = D; // 24. c(i+1) = d(i)
					// B(i+1) = B(i) // 25. b(i+1) = b(i)
				}
				else // 26. else
				{
					// A(i+1) = A(i) // 27. a(i+1) = a(i)
					// C(i+1) = C(i) // 28. c(i+1) = c(i)
					B = D; // 29. b(i+1) = d(i)
				} // 30. end if
			}
			else // 31. else
			{
				D.flag = 0; // 32. return error
				return D;
			} // 33. end if // 34. end if

			i = i + 1; // 35. i = i + 1

			// 39. until b(i) - a(i) < ε or |d(i) - d(i-1)| < γ
			if (abs(B.x(0) - A.x(0)) < epsilon || (!isnan(d_previous) && abs(d_val - d_previous) < gamma))
			{
				D.flag = 1;
				return D; // 40. return x* = d(i)
			}

			d_previous = d_val;
		}
	}
	catch (string ex_info)
	{
		throw("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution XB(x0);								// 1. XB(0) = x0
		XB.fit_fun(ff, ud1, ud2);
		
		solution XB_old;
		int n = get_len(x0);							// liczba zmiennych
		
		while (true)									// 2. repeat
		{
			solution X = HJ_trial(ff, XB, s, ud1, ud2);	// 3. X = HJ_trial(XB(i), s)
			
			if (X.y < XB.y)								// 4. if f(X) < f(XB(i)) then
			{
				while (true)							// 5. repeat
				{
					XB_old = XB;						// 6. XB_old = XB(i)
					XB = X;								// 7. XB(i+1) = X
					
					// 8. X = 2*XB(i+1) - XB_old
					matrix x_new = 2.0 * XB.x - XB_old.x;
					X.x = x_new;
					X.fit_fun(ff, ud1, ud2);
					
					// 9. X = HJ_trial(X, s)
					X = HJ_trial(ff, X, s, ud1, ud2);
					
					// 10. until f(X) >= f(XB(i+1))
					if (X.y >= XB.y)
						break;
				}
			}
			else										// 11. else
			{
				s = alpha * s;							// 12. s = alpha * s
			}											// 13. end if
			
			// 14. if fcalls > Nmax then
			if (solution::f_calls > Nmax)
			{
				XB.flag = 0;							// 15. return XB(i), fcalls > Nmax
				return XB;
			}											// 16. end if
			
			// 17. until s < epsilon
			if (s < epsilon)
			{
				XB.flag = 1;							// 18. return XB(i)
				return XB;
			}
		}
	}
	catch (string ex_info)
	{
		throw("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix (*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		int n = get_len(XB.x);							// 1. liczba zmiennych
		solution X = XB;								// 2. X = XB
		
		// 3. for j = 0 to n-1 do
		for (int j = 0; j < n; ++j)
		{
			// 4. X(j) = X(j) + s
			matrix x_plus = X.x;
			x_plus(j) = x_plus(j) + s;
			solution X_plus(x_plus);
			X_plus.fit_fun(ff, ud1, ud2);
			
			// 5. if f(X) < f(XB) then
			if (X_plus.y < XB.y)
			{
				XB = X_plus;							// 6. XB = X
			}
			else										// 7. else
			{
				// 8. X(j) = XB(j) - s
				matrix x_minus = XB.x;
				x_minus(j) = x_minus(j) - s;
				solution X_minus(x_minus);
				X_minus.fit_fun(ff, ud1, ud2);
				
				// 9. if f(X) < f(XB) then
				if (X_minus.y < XB.y)
				{
					XB = X_minus;						// 10. XB = X
				}										// 11. end if
			}											// 12. end if
			
			X = XB;
		}												// 13. end for
		
		return XB;										// 14. return XB
	}
	catch (string ex_info)
	{
		throw("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix (*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		// 1: i = 0
		int i = 0;
		int n = get_len(x0);
		
		// 2: dj(0) = ej, j = 1, 2, …, n
		matrix D = ident_mat(n);  // macierz kierunków (kolumny to wektory bazowe)
		
		// 3: λj(0) = 0, j = 1, 2, …, n
		matrix lambda(n, 1);  // suma przemieszczeń w każdym kierunku
		
		// 4: pj(0) = 0, j = 1, 2, …, n
		matrix p(n, 1);  // liczba porażek w każdym kierunku
		
		// 5: xB = x(0)
		solution xB(x0);
		xB.fit_fun(ff, ud1, ud2);
		
		matrix s = s0;  // długości kroków
		
		// 6: repeat
		while (true)
		{
			// 7: for j = 1 to n do
			for (int j = 0; j < n; ++j)
			{
				// 8: if f(xB + sj(i)·dj(i)) < f(xB) then
				matrix x_new = xB.x + s(j) * D[j];
				solution x_trial(x_new);
				x_trial.fit_fun(ff, ud1, ud2);
				
				if (x_trial.y < xB.y)
				{
					// 9: xB = xB + sj(i)·dj(i)
					xB = x_trial;
					
					// 10: λj(i+1) = λj(i) + sj(i)
					lambda(j) = lambda(j) + s(j);
					
					// 11: sj(i+1) = α·sj(i)
					s(j) = alpha * s(j);
				}
				else
				{
					// 12: else
					// 13: sj(i+1) = -β·sj(i)
					s(j) = -beta * s(j);
					
					// 14: pj(i+1) = pj(i) + 1
					p(j) = p(j) + 1;
					// 15: end if
				}
			}
			// 16: end for
			
			// 17: i = i + 1
			i = i + 1;
			
			// 18: x(i) = xB  (już mamy w xB)
			
			// 19: if λj(i) ≠ 0 and pj(i) ≠ 0, j = 1, 2, …, n then
			bool all_lambda_nonzero = true;
			bool all_p_nonzero = true;
			
			for (int j = 0; j < n; ++j)
			{
				if (lambda(j) == 0)
					all_lambda_nonzero = false;
				if (p(j) == 0)
					all_p_nonzero = false;
			}
			
			if (all_lambda_nonzero && all_p_nonzero)
			{
				// 20: zmiana bazy kierunków dj(i)
				// Budujemy nową bazę przez ortogonalizację Grama-Schmidta
				// Pierwszy wektor to suma przemieszczeń
				
				matrix D_new = ident_mat(n);
				
				// Pierwszy nowy kierunek: suma ważonych kierunków
				matrix v0(n, 1);
				for (int j = 0; j < n; ++j)
				{
					for (int k = 0; k < n; ++k)
					{
						v0(k) = v0(k) + lambda(j) * D[j](k);
					}
				}
				
				// Normalizacja pierwszego wektora
				double norm0 = 0.0;
				for (int k = 0; k < n; ++k)
					norm0 += v0(k) * v0(k);
				norm0 = sqrt(norm0);
				
				if (norm0 > 1e-15)
				{
					for (int k = 0; k < n; ++k)
						D_new[0](k) = v0(k) / norm0;
				}
				
				// Ortogonalizacja pozostałych wektorów
				for (int j = 1; j < n; ++j)
				{
					matrix vj = D[j];
					
					// Ortogonalizacja względem już znormalizowanych wektorów
					for (int k = 0; k < j; ++k)
					{
						double dot = 0.0;
						for (int m = 0; m < n; ++m)
							dot += vj(m) * D_new[k](m);
						
						for (int m = 0; m < n; ++m)
							vj(m) = vj(m) - dot * D_new[k](m);
					}
					
					// Normalizacja
					double normj = 0.0;
					for (int m = 0; m < n; ++m)
						normj += vj(m) * vj(m);
					normj = sqrt(normj);
					
					if (normj > 1e-15)
					{
						for (int m = 0; m < n; ++m)
							D_new[j](m) = vj(m) / normj;
					}
				}
				
				D = D_new;
				
				// 21: λj(i) = 0, j = 1, 2, …, n
				for (int j = 0; j < n; ++j)
					lambda(j) = 0;
				
				// 22: pj(i) = 0, j = 1, 2, …, n
				for (int j = 0; j < n; ++j)
					p(j) = 0;
				
				// 23: sj(i) = sj(0), j = 1, 2, …, n
				s = s0;
				// 24: end if
			}
			
			// 25: if fcalls > Nmax then
			if (solution::f_calls > Nmax)
			{
				// 26: return error
				xB.flag = 0;
				return xB;
				// 27: end if
			}
			
			// 28: until maxj=1,…,n(|sj(i)|) < ε
			double max_s = 0.0;
			for (int j = 0; j < n; ++j)
			{
				if (abs(s(j)) > max_s)
					max_s = abs(s(j));
			}
			
			if (max_s < epsilon)
			{
				// 29: return x* = x(i)
				xB.flag = 1;
				return xB;
			}
		}
	}
	catch (string ex_info)
	{
		throw("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	/*
	METODA FUNKCJI KARY - zgodnie z pseudokodem
	
	PSEUDOKOD:
	1: i = 0
	2: repeat
	3:    i = i + 1
	4:    wyznacz F^(i)(x) = f(x) + c^(i)S(x)
	5:    wyznacz x^(i) dla F^(i) startując z x^(i-1)
	6:    c^(i+1) = α·c^(i)
	7:    if f_calls > Nmax then
	8:       return error
	9:    end if
	10: until ||x^(i) - x^(i-1)||_2 < ε
	11: return x* = x^(i)
	
	FUNKCJE KARY:
	- Zewnętrzna: S(x) = Σ(max(0, g_i(x)))^2
	- Wewnętrzna: S(x) = -Σ(1/g_i(x))
	
	PARAMETRY:
	- Zewnętrzna: c^(1) > 0 (mała), α > 1 (zwiększamy)
	- Wewnętrzna: c^(1) > 0 (duża), 0 < α < 1 (zmniejszamy)
	*/
	
	try
	{
		solution Xopt;
		
		// 1: i = 0
		int i = 0;
		
		// x^(0) = x0
		solution x_prev(x0);
		
		// Parametry dla metody Nelder-Mead (simpleks)
		double s = 0.5;			// długość krawędzi początkowego simpleksu
		double alpha_nm = 1.0;	// współczynnik odbicia
		double beta = 0.5;		// współczynnik kontrakcji
		double gamma = 2.0;		// współczynnik ekspansji
		double delta = 0.5;		// współczynnik redukcji
		
		// c^(0) = c (współczynnik kary)
		double c_current = c;
		
		// Typ funkcji kary z ud2: 0 = zewnętrzna, 1 = wewnętrzna
		int penalty_type = (int)ud2(0);
		
		// 2: repeat
		while (true)
		{
			// 3: i = i + 1
			i = i + 1;
			
			// 4: wyznacz F^(i)(x) = f(x) + c^(i)S(x)
			// (funkcja kary jest obliczana w ff z ud2_pen)
			
			// 5: wyznacz x^(i) dla F^(i) startując z x^(i-1)
			matrix ud2_pen(2, 1);
			ud2_pen(0) = penalty_type;		// typ funkcji kary
			ud2_pen(1) = c_current;			// c^(i)
			
			solution x_current = sym_NM(ff, x_prev.x, s, alpha_nm, beta, gamma, delta, epsilon, Nmax, ud1, ud2_pen);
			
			// 6: c^(i+1) = α·c^(i)
			c_current = dc * c_current;
			
			// 7: if f_calls > Nmax then
			if (solution::f_calls > Nmax)
			{
				// 8: return error
				Xopt = x_current;
				Xopt.flag = 0;
				return Xopt;
			}
			// 9: end if
			
			// 10: until ||x^(i) - x^(i-1)||_2 < ε
			matrix diff = x_current.x - x_prev.x;
			double norm_diff = 0.0;
			for (int j = 0; j < get_len(diff); ++j)
				norm_diff += pow(diff(j), 2);
			norm_diff = sqrt(norm_diff);
			
			if (norm_diff < epsilon)
			{
				// 11: return x* = x^(i)
				Xopt = x_current;
				Xopt.flag = 1;
				return Xopt;
			}
			
			// x^(i-1) = x^(i) dla następnej iteracji
			x_prev = x_current;
			
			// Dodatkowe zabezpieczenie: maksymalna liczba iteracji zewnętrznej pętli
			if (i > 100)
			{
				Xopt = x_current;
				Xopt.flag = 1;
				return Xopt;
			}
		}
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix (*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		
		// 1. n = dim(x(0))
		int n = get_len(x0);
		
		// 2. Tworzenie początkowego simpleksu
		// Simpleks będzie miał n+1 wierzchołków
		vector<solution> simplex(n + 1);
		
		// Pierwszy wierzchołek to x0
		simplex[0].x = x0;
		simplex[0].fit_fun(ff, ud1, ud2);
		
		// Pozostałe wierzchołki tworzymy przez przesunięcie wzdłuż osi współrzędnych
		for (int i = 1; i <= n; ++i)
		{
			simplex[i].x = x0;
			simplex[i].x(i - 1) = simplex[i].x(i - 1) + s;
			simplex[i].fit_fun(ff, ud1, ud2);
		}
		
		// 3. repeat
		while (true)
		{
			// 4. Sortowanie wierzchołków według wartości funkcji celu (rosnąco)
			sort(simplex.begin(), simplex.end(), [](const solution& a, const solution& b) {
				return a.y(0) < b.y(0);
			});
			
			// Najlepszy: simplex[0]
			// Drugi najgorszy: simplex[n-1]
			// Najgorszy: simplex[n]
			
			// 5. Obliczenie centroidu (bez najgorszego punktu)
			matrix centroid(n, 1);
			for (int i = 0; i < n; ++i)
			{
				for (int j = 0; j < n; ++j)
				{
					centroid(j) = centroid(j) + simplex[i].x(j);
				}
			}
			for (int j = 0; j < n; ++j)
			{
				centroid(j) = centroid(j) / n;
			}
			
			// 6. Odbicie (reflection)
			matrix x_reflected = centroid + alpha * (centroid - simplex[n].x);
			solution X_reflected(x_reflected);
			X_reflected.fit_fun(ff, ud1, ud2);
			
			// 7. if f(x_best) <= f(x_reflected) < f(x_second_worst) then
			if (simplex[0].y(0) <= X_reflected.y(0) && X_reflected.y(0) < simplex[n - 1].y(0))
			{
				// 8. Zastąp najgorszy punkt punktem odbitym
				simplex[n] = X_reflected;
			}
			// 9. else if f(x_reflected) < f(x_best) then
			else if (X_reflected.y(0) < simplex[0].y(0))
			{
				// 10. Ekspansja (expansion)
				matrix x_expanded = centroid + gamma * (x_reflected - centroid);
				solution X_expanded(x_expanded);
				X_expanded.fit_fun(ff, ud1, ud2);
				
				// 11. if f(x_expanded) < f(x_reflected) then
				if (X_expanded.y(0) < X_reflected.y(0))
				{
					// 12. Zastąp najgorszy punkt punktem rozszerzonym
					simplex[n] = X_expanded;
				}
				else
				{
					// 13. Zastąp najgorszy punkt punktem odbitym
					simplex[n] = X_reflected;
				}
				// 14. end if
			}
			// 15. else
			else
			{
				// 16. Kontrakcja (contraction)
				matrix x_contracted = centroid + beta * (simplex[n].x - centroid);
				solution X_contracted(x_contracted);
				X_contracted.fit_fun(ff, ud1, ud2);
				
				// 17. if f(x_contracted) < f(x_worst) then
				if (X_contracted.y(0) < simplex[n].y(0))
				{
					// 18. Zastąp najgorszy punkt punktem skurczonym
					simplex[n] = X_contracted;
				}
				// 19. else
				else
				{
					// 20. Redukcja (shrink) - kurczenie całego simpleksu w kierunku najlepszego punktu
					for (int i = 1; i <= n; ++i)
					{
						simplex[i].x = simplex[0].x + delta * (simplex[i].x - simplex[0].x);
						simplex[i].fit_fun(ff, ud1, ud2);
					}
				}
				// 21. end if
			}
			// 22. end if
			
			// 23. Sprawdzenie warunku zakończenia
			if (solution::f_calls > Nmax)
			{
				Xopt = simplex[0];
				Xopt.flag = 0;
				return Xopt;
			}
			
			// Obliczenie odchylenia standardowego wartości funkcji w simpleksie
			double mean_y = 0.0;
			for (int i = 0; i <= n; ++i)
			{
				mean_y += simplex[i].y(0);
			}
			mean_y /= (n + 1);
			
			double std_dev = 0.0;
			for (int i = 0; i <= n; ++i)
			{
				std_dev += pow(simplex[i].y(0) - mean_y, 2);
			}
			std_dev = sqrt(std_dev / (n + 1));
			
			// 24. until std_dev < epsilon
			if (std_dev < epsilon)
			{
				Xopt = simplex[0];
				Xopt.flag = 1;
				return Xopt;
			}
		}
		
		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix (*ff)(matrix, matrix, matrix), matrix (*gf)(matrix, matrix, matrix),
				matrix (*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix (*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution golden(...):\n" + ex_info);
	}
}

solution Powell(matrix (*ff)(matrix, matrix, matrix), matrix x0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution Powell(...):\n" + ex_info);
	}
}

solution EA(matrix (*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, int mi, int lambda, matrix sigma0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw("solution EA(...):\n" + ex_info);
	}
}
