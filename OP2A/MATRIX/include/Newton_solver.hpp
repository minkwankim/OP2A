/*
 * Open-source multi-Physics Phenomena Analyzer (OPPA) ver. 0.0
 *
 * 		Copyright (c) 2014 MINKWAN KIM
 *
 * 	Initial Developed Date: May 28, 2014
 *      			Author: Minkwan Kim
 *
 * Newton_solver.hpp
 * 			-  
 *  
 */

#ifndef NEWTON_SOLVER_HPP_
#define NEWTON_SOLVER_HPP_

#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;

#include "../../UTIL/include/OP2A_utilities.hpp"
#include "math_misc.hpp"

template <class Fn_class>
class Newton_solver
{
public:
	double		X_min;
	double		X_max;
	int 		Iter_MAX;
	double		epsilon;

	Newton_solver();
	~Newton_solver();

	double solve(Fn_class fn, int type);
};


template <class Fn_class>
Newton_solver<Fn_class>::Newton_solver()
{
	X_min		= 0.0;
	X_max		= 0.0;
	Iter_MAX	= 1e6;
	epsilon		= 1.0e-5;
}

template <class Fn_class>
Newton_solver<Fn_class>::~Newton_solver()
{
}

template <class Fn_class>
double Newton_solver<Fn_class>::solve(Fn_class fn, int type)
{
	int		i_iter;
	bool	flag;
	double	X_new, X_old, dX;
	double	f_max, f_min;
	double	f_old, f_new;
	double	sign;
	double	error;


	/*
	 * Find X_min and X_max
	 */
	switch (type)
	{
	case 1:
		// Step1 : Find X_min and X_max
		// T_min
		flag	= false;
		i_iter	= 0;
		dX		= 10.0;
		X_old	= X_min;

		while (flag == false)
		{
			X_new	= X_old + dX;
			f_min	= fn.f(X_new);

			if (f_min < 0.0)
			{
				dX		= 1.03*dX;
				X_old	= X_new;
			}
			else if (f_min > 0.0)
			{
				flag = true;
			}
			else
			{
				return (X_new);
			}

			i_iter++;
			if (i_iter > Iter_MAX)
			{
				cout << "ERROR[0004]: CAUGHT IN A LOOP IN Newton Solver!! Xmin =" << X_min << endl;
				program_error_type1("PLEASE, CHECK VIBRATIONAL ENERGY AND TEMPERATURE !!! \n");
			}
		}

		X_min	= X_old;
		f_min	= fn.f(X_min);


		// T_max
		flag	= false;
		i_iter	= 0;
		dX		= 50.0;
		X_old	= X_max;

		while (flag == false)
		{
			X_new	= X_old + dX;
			f_max	= fn.f(X_new);

			if (f_max > 0.0)
			{
				flag 	= true;
				X_old	= X_new;
			}
			else if (f_max < 0.0)
			{
				dX		= 1.1*dX;
				X_old	= X_new;
			}
			else
			{
				return (X_new);
			}

			i_iter++;
			if (i_iter > Iter_MAX)
			{
				cout << "ERROR[0004]: CAUGHT IN A LOOP IN Newton Solver!! Xmax =" << X_max << endl;
				program_error_type1("PLEASE, CHECK VIBRATIONAL ENERGY AND TEMPERATURE !!! \n");
			}
		}

		X_max	= X_old;
		f_max	= fn.f(X_max);


		if (f_min * f_max > 0.0)
		{
			cout << "ERROR[0004]: PROBLEM IN THE FINDING OF Xmin AND Xmax!! Xmin = " << X_min << " Xmax =" << X_max << endl;
			program_error_type1("PLEASE, CHECK THE ALGORITHM OF FINDING MAX AND MIN!!! \n");
		}
		break;

	case 2:
		// Step1 : Find X_min and X_max
		// T_min
		flag	= false;
		i_iter	= 0;
		dX		= 10.0;
		X_old	= X_min;

		while (flag == false)
		{
			X_new	= X_old + dX;
			f_min	= fn.f(X_new);

			if (f_min > 0.0)
			{
				dX		= 1.03*dX;
				X_old	= X_new;
			}
			else if (f_min < 0.0)
			{
				flag = true;
			}
			else
			{
				return (X_new);
			}

			i_iter++;
			if (i_iter > Iter_MAX)
			{
				cout << "ERROR[0004]: CAUGHT IN A LOOP IN Newton Solver!! Xmin =" << X_min << endl;
				program_error_type1("PLEASE, CHECK VIBRATIONAL ENERGY AND TEMPERATURE !!! \n");
			}
		}

		X_min	= X_old;
		f_min	= fn.f(X_min);


		// T_max
		flag	= false;
		i_iter	= 0;
		dX		= 50.0;
		X_old	= X_max;

		while (flag == false)
		{
			X_new	= X_old + dX;
			f_max	= fn.f(X_new);

			if (f_max < 0.0)
			{
				flag 	= true;
				X_old	= X_new;
			}
			else if (f_max > 0.0)
			{
				dX		= 1.1*dX;
				X_old	= X_new;
			}
			else
			{
				return (X_new);
			}

			i_iter++;
			if (i_iter > Iter_MAX)
			{
				cout << "ERROR[0004]: CAUGHT IN A LOOP IN Newton Solver!! Xmax =" << X_max << endl;
				program_error_type1("PLEASE, CHECK VIBRATIONAL ENERGY AND TEMPERATURE !!! \n");
			}
		}

		X_max	= X_old;
		f_max	= fn.f(X_max);


		if (f_min * f_max > 0.0)
		{
			cout << "ERROR[0004]: PROBLEM IN THE FINDING OF Xmin AND Xmax!! Xmin = " << X_min << " Xmax =" << X_max << endl;
			program_error_type1("PLEASE, CHECK THE ALGORITHM OF FINDING MAX AND MIN!!! \n");
		}
		break;
	}




	/*
	 * METHOD 1: NEWTON
	 *
	 * */
	error = 1.0;
	X_old = X_min;
	X_new = X_min;

	for (i_iter = 0; i_iter <= Iter_MAX; i_iter++)
	{
		X_old = X_new;
		if (error <= epsilon)
		{
			return (X_new);
		}
		else
		{
			dX 		= fn.f(X_old) / fn.df(X_old);
			error	= fabs(dX);
			X_new 	= X_old - dX;
		}

		if (X_new > 1.5*X_max)
		{
			X_new = X_max - 2.0*fabs(dX);
		}
		else if (X_new < 0.9*X_min)
		{
			X_new = X_min + 0.1*fabs(dX);
		}
	}

	cout << ".....FAILURE METHOD 1[NEWTON MEHTOD] TO SOLVE EQUATION. WILL TRY METHOD 2!!" << endl;


	/* METHOD 2 */
	dX = 10.0;
	X_old = X_min;
	f_old = fn.f(X_old);
	sign = fsgn(f_old);

	switch (type)
	{
	case 1:
		if (sign > 0.0)
		{
			cout << "ERROR[0004]: PROBLEM IN Xmin, f(Xmin) should be smaller than 0 Xmin=" << X_min << " fmin=" << f_old << endl;
			program_error_type1("PLEASE, CHECK THE ALGORITHM OF FINDING Tmin!!! \n");
		}

		X_new = X_old + dX;
		for (i_iter = 0; i_iter <= 10*Iter_MAX; i_iter++)
		{
			X_old = X_new;
			f_new = fn.f(X_new);

			sign = fsgn(f_new) * fsgn(f_old);

			if (sign > 0.0)
			{
				dX = 1.2*dX;
				if (f_new < 0.0)
					X_new = X_old + dX;
				else
					X_new = X_old - dX;
			}
			else if (sign < 0.0)
			{
				dX = 0.7*dX;
				if (f_new < 0.0)
					X_new = X_old + dX;
				else
					X_new = X_old - dX;
			}
			else
			{
				if (f_new == 0.0)
					return (X_new);
			}

			f_old = f_new;

			if (fabs(X_new - X_old) <= epsilon)
			{
				return (X_new);
			}

			if (X_new > X_max)
			{
				X_new = X_max - 0.1*fabs(dX);
			}
			else if (X_new < X_min)
			{
				X_new = X_min + 0.1*fabs(dX);
			}
		}
		break;

	case 2:
		if (sign < 0.0)
		{
			cout << "ERROR[0004]: PROBLEM IN Xmin, f(Xmin) should be greater than 0 Xmin=" << X_min << " fmin=" << f_old << endl;
			program_error_type1("PLEASE, CHECK THE ALGORITHM OF FINDING Tmin!!! \n");
		}

		X_new = X_old + dX;
		for (i_iter = 0; i_iter <= 10*Iter_MAX; i_iter++)
		{
			X_old = X_new;
			f_new = fn.f(X_new);

			sign = fsgn(f_new) * fsgn(f_old);

			if (sign < 0.0)
			{
				dX = 1.2*dX;
				if (f_new > 0.0)
					X_new = X_old + dX;
				else
					X_new = X_old - dX;
			}
			else if (sign > 0.0)
			{
				dX = 0.7*dX;
				if (f_new > 0.0)
					X_new = X_old + dX;
				else
					X_new = X_old - dX;
			}
			else
			{
				if (f_new == 0.0)
					return (X_new);
			}

			f_old = f_new;

			if (fabs(X_new - X_old) <= epsilon)
			{
				return (X_new);
			}

			if (X_new > X_max)
			{
				X_new = X_max - 0.1*fabs(dX);
			}
			else if (X_new < X_min)
			{
				X_new = X_min + 0.1*fabs(dX);
			}
		}
		break;
	}


	cout << "ERROR[0004]: PROBLEM IN THE CALCULATION OF VIBRATIONAL TEMP" << endl;
	program_error_type1(" PLEASE, CHECK VIBRATIONAL ENERGY AND TEMPERATURE !!! \n");
	return (0);
}


#endif /* NEWTON_SOLVER_HPP_ */

