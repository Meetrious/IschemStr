/* This header contains classes for variables list declaration in 4th-model ODE_system*/

#pragma once
#include <cstdint>
#include <cmath>

namespace StraightTask
{
	class RetVar{
	public:
		/* each double_t field is to be named by following pattern:
		"DependentVariableName_HoursLag" */
		RetVar() {
			*this = 0;
		}

		void operator=(double_t val) {
			adh_4 = val;
		}

		~RetVar() = default;

		double_t adh_4;
	};

//==========================================================================================================================
	class variables	{
	public:

		variables() {
			*this = 0;
		}
		
		void operator=(double_t val) {
			nec = acu_c = hel = val;
			cy = ch = adh = val;
			mia = mii = lm = ln = val;
			eps_s = eps_w = val;
			d_F = d_ini = val;
			dp_A = dp_N = val;

			ret = val;

		}

		~variables() = default;

		// list of dependent variables retarded in time
		RetVar ret;

		// independent variable of time
		double_t tj;


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