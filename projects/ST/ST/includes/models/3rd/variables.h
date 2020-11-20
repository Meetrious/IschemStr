/* This header contains classes for variables list declaration in 3rd-model ODE_system*/

#pragma once
#include <cstdint>
#include <cmath>


namespace StraightTask
{
	class RetVar{
	public:
		/* each double_t field is to be named by following pattern:
		"DependentVariableName_HoursLag" */
		RetVar(): adh_4{0} {}
		~RetVar() = default;

		double_t adh_4;
	};

//==========================================================================================================================
	class variables	{
	public:

		variables() {}
		~variables() = default;

		// list of dependent variables retarded in time
		RetVar ret;

		// independent variable of time
		float_t tj;


		// down below is a list of dependent variables

		double_t nec, acu_c, hel ;
		double_t cy, ch, adh;
		double_t mia, mii, lm, ln;


		// subbordinate values

		double_t eps_s, eps_w;
		double_t d_F, d_ini;
		double_t dp_A, dp_N;


		uint32_t spl_gap_counter;
	};

}