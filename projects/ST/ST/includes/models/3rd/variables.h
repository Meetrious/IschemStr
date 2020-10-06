#pragma once
namespace StraightTask
{

	class RetVar
	{
	public:
		double_t adh_4 = 0;
	};

//==========================================================================================================================
	class variables
	{
	public:

		variables() {}

		~variables() = default;

		RetVar ret;

		double_t tj = 0;

		double_t nec = 0, acu_c = 0, hel = 0 ;

		double_t cy = 0, ch = 0, adh = 0 ;
		double_t eps_s = 0, eps_w = 0;

		double_t mia = 0, mii = 0, lm = 0, ln = 0;

		double_t d_F = 0, d_ini = 0;
		double_t dp_A = 0, dp_N = 0;

		uint32_t spl_gap_counter=0;
	};

}