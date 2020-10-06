#pragma once

namespace StraightTask
{

	class RetVar
	{
	public:
		// Value_TimeGapBackFromCurrentMoment
		double_t hel_12 = 0;
		double_t adh_12 = 0, adh_24 = 0;
		double_t d_12 = 0;
	};

//==========================================================================================================================
	class variables
	{
	public:

		variables() {}

		~variables() = default;

		RetVar ret; // storage for current ret-values

		// current time
		double_t tj = 0;

		// population distribution among neuroglia
		double_t nec = 0, ap_s = 0 , ap_e = 0 , hel = 0 ;

		// signal protein
		double_t cy = 0, ch = 0, adh = 0 ;

		double_t eps_s = 0;

		// macrophages population
		double_t mia = 0, mii = 0, lm = 0, ln = 0;

		double_t d_F = 0, d_ini = 0;

		double_t psy = 0;

		uint32_t spl_gap_counter=0;
	};

}