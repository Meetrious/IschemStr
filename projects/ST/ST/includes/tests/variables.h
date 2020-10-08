// This header contains dependent variables of the Test-tasks.

#pragma once



namespace StraightTask
{

	class RetVar
	{
	public:
		// value_TimeGapBackFromCurrentMoment
		double_t x_1 = 0;
		double_t x_pi2 = 0, y_pi2 = 0, z_pi2 = 0;
	};

//==========================================================================================================================
	class variables 
	{
	public:

		variables() {}

		~variables() = default;

		RetVar ret; // storage for ret-values

		double_t tj = 0.0;

		double_t x = 0.0, y = 0.0, z = 0.0;

		uint32_t spl_gap_counter=0;
	};

}