/* This header contains non-member functions that implements surface algorythms
 that serve to solve given straight task for various purposes:
 i.e in order to output solution, or calculate the relative error,
 or to count the practical order of the method on given task using Runge-method // */

// This header is included in every header *_pipe.h located on top of ST/include directory

#pragma once

namespace StraightTask
{
	// Solves the straight task and outputs solution
	void ODE_solver_output()
	{
		// choosing num-method from StraightTask
		Euler SYS;

		// setting an object that will output 2 or 3 dimentional phase trajectory for solution 
		// PhaseTrajOutput PTO("APPROX");

		// Setting parameters for num-method: N, \tau, full_amount_of_gaps
		SYS.ST.Set(1500, 24.0, 4);

		SYS.PrepairTheTask();

		SYS.PrepairTheOutput();

		uint16_t current_gap = 0;
		double_t Tj;

		while (current_gap < SYS.ST.full_amount_of_gaps)
		{

			std::cout << SYS.ST.full_amount_of_gaps - current_gap << " "; // outputing amount of gaps remained to process
			SYS.ST.X_prev = SYS.ST.X_init; // previous step is defined by initial one at every begining of the gap
			Tj = SYS.ST.X_pred.tj = SYS.ST.X_prev.tj; // synchronising the independent variable
			uint32_t Nj = 1; // restating the number of the current step in num-method

			// 1st approximations for multistep-methods
			if (current_gap == 0) SYS.ApplyPrepStep(Nj, Tj); // in one-step-methods it is void{return;}

			// cycle for processing current gap in <step_method>
			for (; Nj <= SYS.ST.N; Nj++)
			{
				// shifting independent variable on one step further for predicted solution
				Tj = SYS.ST.X_pred.tj += SYS.ST.H;

				// setting ret-values for Tj-time-moment in X_pred
				if (SYS.is_SYS_deflecting()) SYS.RetUpload(Nj);

				// контролируем промежуток интерпол€ции
				//ST.X_pred.CheckShiftInterpGap(STpar);

				/* pushing presolved data in X_pred in case we want to freeze the system relatively given behaviour
				 defined in SYS.<IAggregate_member_field>.PreSavedSolData;
				 matrix<double_t> StraightTask::IOs SolDataGetter is a member function that does things // */
				{
				//for (auto const& cur : SYS.SolDataGetter) { cur(current_day, Nj, SYS.ST.X_pred); }	

					//SYS.SolDataGetter[3](current_gap, Nj, SYS.ST.N, SYS.ST.X_pred);
					//SYS.SolDataGetter[3](current_gap, Nj, SYS.ST.N, SYS.ST.X_cor); // in pred_cor scheme
				}

				SYS.ApplyMethod();

				// outputting solution in current Tj - time-moment
				for (auto const& cur : SYS.SolutionOutputter) { cur(Nj, Tj, *SYS.ST.X_sol); }

				// Phase trajectory output
				//PTO.OutputPhaseTraj(Nj, Tj, { (*SYS.ST.X_sol).x, (*SYS.ST.X_sol).y, (*SYS.ST.X_sol).z });


				// outputting equation budget-bricks
				// for (auto const& cur : SYS.BudgetOutputter) { cur(Nj, Tj); }

				// updating RetArray(s) pushing X_prev.value of solution
				if (SYS.is_SYS_deflecting()) SYS.RetDataUpdate(Nj);

				// shifting X_prev next step further in one-step-methods and
				// X[1],X[2]... next step further in multistep-methods
				SYS.NodeShift();


			} // end of the cycle processing current gap in <step_method>

			current_gap++; std::cout << " ; ";

			// Checking if the any gap remained processed 
			if (current_gap == SYS.ST.full_amount_of_gaps)
				std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			else
			{
				SYS.ST.X_init = *SYS.ST.X_sol;
				Nj = 1;
			}
		}
	}
}