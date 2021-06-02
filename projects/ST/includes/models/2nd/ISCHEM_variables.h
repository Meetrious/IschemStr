#pragma once
#include <array>
#include <iomanip> // нужен для std::setpresicion модификатора в std::fstream потоках
#include <base/DirPaths.h>
#include <base/SplRealisation.h>

namespace StraightTask
{

	class RetVar
	{
	public:
		double_t a_12 = 0;
		double_t hel_12 = 0;
		double_t adh_4 = 0, adh_12 = 0, adh_24 = 0;
		double_t d_12 = 0;

		double_t x_1 = 0;
		double_t x_pi2 = 0, y_pi2 = 0, z_pi2 = 0;
	};

//==========================================================================================================================
	class variables
	{
	public:

		variables() {}

		~variables() = default;

		RetVar ret;

		double_t tj = 0;

		double_t nec = 0, acu_c = 0, ap_s = 0 , ap_e = 0 , hel = 0 ;

		double_t cy = 0, ch = 0, adh = 0 ;
		double_t eps_s = 0, eps_w = 0;

		double_t mia = 0, mii = 0, lm = 0, ln = 0;

		double_t d_F = 0, d_ini = 0;
		double_t dp_A = 0, dp_N = 0;

		double_t psy = 0;

		double_t x = 0, y = 0, z = 0; // for tests

		uint32_t spl_gap_counter=0;
	};

}