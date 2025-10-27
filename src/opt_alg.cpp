#include "opt_alg.h"

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
		int n = get_len(x0);								// liczba zmiennych
		matrix l(n, n);										// macierz kierunków (identity matrix)
		for (int i = 0; i < n; ++i)
			l(i, i) = 1.0;
		
		matrix p(n, 1), lambda(n, 1);						// przyrosty i mnożniki
		solution X(x0);										// 1. X(0) = x0
		X.fit_fun(ff, ud1, ud2);
		solution X_old = X;
		matrix s = s0;										// długości kroków
		
		while (true)										// 2. repeat
		{
			// 3. for j = 0 to n-1 do
			for (int j = 0; j < n; ++j)
			{
				p(j) = 0;									// 4. p(j) = 0
				lambda(j) = 0;								// 5. lambda(j) = 0
			}												// 6. end for
			
			// 7. repeat
			while (true)
			{
				// 8. for j = 0 to n-1 do
				for (int j = 0; j < n; ++j)
				{
					solution X_new = X;
					// 9. X = X + s(j) * l(j)
					for (int i = 0; i < n; ++i)
						X_new.x(i) = X.x(i) + s(j) * l(i, j);
					X_new.fit_fun(ff, ud1, ud2);
					
					// 10. if f(X) < f(X(i)) then
					if (X_new.y < X.y)
					{
						X = X_new;							// 11. X(i+1) = X
						p(j) = p(j) + s(j);					// 12. p(j) = p(j) + s(j)
						lambda(j) = lambda(j) + 1;			// 13. lambda(j) = lambda(j) + 1
						s(j) = alpha * s(j);				// 14. s(j) = alpha * s(j)
					}
					else									// 15. else
					{
						s(j) = -beta * s(j);				// 16. s(j) = -beta * s(j)
						lambda(j) = lambda(j) - 1;			// 17. lambda(j) = lambda(j) - 1
					}										// 18. end if
				}											// 19. end for
				
				// 20. if fcalls > Nmax then
				if (solution::f_calls > Nmax)
				{
					X.flag = 0;								// 21. return X, fcalls > Nmax
					return X;
				}											// 22. end if
				
				// 23. until any(lambda(j) != 0) = false
				bool any_lambda_nonzero = false;
				for (int j = 0; j < n; ++j)
				{
					if (lambda(j) != 0)
					{
						any_lambda_nonzero = true;
						break;
					}
				}
				if (!any_lambda_nonzero)
					break;
			}
			
			// Sprawdzenie warunku stopu
			bool change_small = true;
			for (int j = 0; j < n; ++j)
			{
				if (abs(p(j)) >= epsilon)
				{
					change_small = false;
					break;
				}
			}
			
			// 24. if ||p|| < epsilon then
			if (change_small)
			{
				X.flag = 1;									// 25. return X
				return X;
			}												// 26. end if
			
			// 27. Gram-Schmidt orthogonalization
			// Q = [l(0), l(1), ..., l(n-1)]
			matrix Q = l;									// kopiujemy aktualną bazę
			matrix v = p;									// v(0) = p
			
			// 28. for j = 1 to n-1 do
			for (int j = 1; j < n; ++j)
			{
				matrix v_j(n, 1);
				for (int i = 0; i < n; ++i)
					v_j(i) = Q(i, j);
				
				// 29. v(j) = Q(j)
				// 30. for k = 0 to j-1 do
				for (int k = 0; k < j; ++k)
				{
					matrix v_k(n, 1);
					for (int i = 0; i < n; ++i)
						v_k(i) = (k == 0) ? v(i) : Q(i, k);
					
					// 31. v(j) = v(j) - (v(k)^T * Q(j)) / (v(k)^T * v(k)) * v(k)
					double numerator = 0, denominator = 0;
					for (int i = 0; i < n; ++i)
					{
						numerator += v_k(i) * v_j(i);
						denominator += v_k(i) * v_k(i);
					}
					
					if (abs(denominator) > 1e-10)
					{
						double coeff = numerator / denominator;
						for (int i = 0; i < n; ++i)
							v_j(i) = v_j(i) - coeff * v_k(i);
					}
				}												// 32. end for
				
				// Zapisujemy v_j
				for (int i = 0; i < n; ++i)
				{
					if (j == 1)
						v(i) = v_j(i);						// przesuwamy v
					Q(i, j - 1) = v_j(i);
				}
			}													// 33. end for
			
			// 34. for j = 0 to n-1 do
			for (int j = 0; j < n; ++j)
			{
				double norm = 0;
				for (int i = 0; i < n; ++i)
				{
					double val = (j == 0) ? v(i) : Q(i, j - 1);
					norm += val * val;
				}
				norm = sqrt(norm);
				
				// 35. l(j) = v(j) / ||v(j)||
				if (norm > 1e-10)
				{
					for (int i = 0; i < n; ++i)
					{
						if (j == 0)
							l(i, j) = v(i) / norm;
						else
							l(i, j) = Q(i, j - 1) / norm;
					}
				}
				
				// 36. s(j) = 0.5 * s(j)
				s(j) = 0.5 * s(j);
			}													// 37. end for
		}
	}
	catch (string ex_info)
	{
		throw("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix (*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		// Tu wpisz kod funkcji

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
		// Tu wpisz kod funkcji

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
