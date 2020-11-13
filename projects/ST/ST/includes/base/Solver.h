/* This header contains non-member functions that implements surface algorythms
 that serve to solve given straight task for various purposes:
 i.e in order to output solution, or calculate the relative error,
 or to count the practical order of the method on given task using Runge-method // */

// This header is included in every header *_pipe.h located on top of ST/include directory

#pragma once

namespace StraightTask
{
	// Solves the straight task and outputs solution
	template<typename Method>
	void ISolver<Method>::SolveAndOutput(uint32_t N, double_t GapWidth, uint16_t full_amount_of_gaps){

		// setting an object that will output 2 or 3 dimentional phase trajectory for solution 
		// PhaseTrajOutput PTO("APPROX");

		// Setting parameters for num-method: N, \tau, full_amount_of_gaps
		ST.Set(N, GapWidth, full_amount_of_gaps);

		PrepairTheTask();

		PrepairTheOutput();

		uint16_t current_gap = 0;
		double_t Tj;

		while (current_gap < full_amount_of_gaps)
		{

			std::cout << full_amount_of_gaps - current_gap << " "; // outputing amount of gaps remained to process
			ST.X_prev = ST.X_init; // previous step is defined by initial one at every begining of the gap
			Tj = ST.X_pred.tj = ST.X_prev.tj; // synchronising the independent variable
			uint32_t Nj = 1; // restating the number of the current step in num-method

			// 1st approximations for multistep-methods
			if (current_gap == 0) ApplyPrepStep(Nj, Tj); // in one-step-methods it is void{return;}

			// cycle for processing current gap in <step_method>
			for (; Nj <= N; Nj++)
			{
				// shifting independent variable on one step further for predicted solution
				Tj = ST.X_pred.tj += ST.H;

				// setting ret-values for Tj-time-moment in X_pred
				if (is_SYS_deflecting()) RetUpload(Nj);

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

				ApplyMethod();

				// outputting solution in current Tj - time-moment
				//for (auto const& cur : SolutionOutputter) { cur(Nj, Tj, (*ST.X_sol)); }
				OutputSolution(Nj, Tj, (*ST.X_sol));

				// Phase trajectory output
				//PTO.OutputPhaseTraj(Nj, Tj, { (*SYS.ST.X_sol).x, (*SYS.ST.X_sol).y, (*SYS.ST.X_sol).z });


				// outputting equation budget-bricks
				//for (auto const& cur : BudgetOutputter) { cur(Nj, Tj); }
				OutputBudgets(Nj, Tj);

				// updating RetArray(s) pushing X_prev.value of solution
				if (is_SYS_deflecting()) RetDataUpdate(Nj);

				// shifting X_prev next step further in one-step-methods and
				// X[1],X[2]... next step further in multistep-methods
				NodeShift();


			} // end of the cycle processing current gap in <step_method>

			current_gap++; std::cout << " ; ";

			// Checking if the any gap remained processed 
			if (current_gap == ST.full_amount_of_gaps)
				std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			else
			{
				ST.X_init = (*ST.X_sol);
				Nj = 1;
			}
		}
	}

	// Solves the straight task three times and calculates the order of assignated num-method using Runge-rule
	template<typename Method>
	void ISolver<Method>::SolveForRungeAnalysis(){
		return;
	}

}