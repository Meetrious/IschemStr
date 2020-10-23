#pragma once
#include <cstdint>
#include <cmath>


namespace StraightTask
{
	class RetVar{
	public:

		RetVar(): adh_4{0} {}

		~RetVar() = default;

		double_t adh_4;
	};

//==========================================================================================================================
	class variables	{
	public:

		variables() {}

		~variables() = default;

		RetVar ret;

		double_t tj;

		double_t nec, acu_c, hel ;

		double_t cy, ch, adh ;
		double_t eps_s, eps_w;

		double_t mia, mii, lm, ln;

		double_t d_F, d_ini;

		double_t dp_A, dp_N;

		uint32_t spl_gap_counter;
	};

}