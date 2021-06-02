namespace StraightTask

{
	template<typename Method>
	double_t ISolver<Method>::SolveForBGA(ReverseTask::IAggregateControls& F) {


		PrepairTheTask();

		uint16_t current_gap = 0;
		double_t Tj;



		F.ResetStates();
		F.CollectCalc2(0, 0, Mthd.N, Mthd.X_init);


		while (current_gap < Mthd.full_amount_of_gaps) {

			Mthd.X_prev = Mthd.X_init; // previous step is defined by initial one at every begining of the gap
			Tj = Mthd.X_pred.tj = Mthd.X_prev.tj; // synchronising the independent variable
			uint32_t Nj = 1; // restating the number of the current step in num-method

			// 1st approximations for multistep-methods
			if (current_gap == 0) ApplyPrepStep(Nj, Tj); // in one-step-methods it is void{return;}

			// cycle for processing current gap in <step_method>
			for (; Nj <= Mthd.N; Nj++)
			{
				// shifting independent variable on one step further for predicted solution
				Tj = Mthd.X_pred.tj += Mthd.H;

				// setting ret-values for Tj-time-moment in X_pred
				if (is_SYS_deflecting()) RetUpload(Nj);

				
				AssignSolData(current_gap, Nj, Mthd.X_pred);
				ApplyMethod();

				//F.CollectCalc(Mthd.H, Tj, (*Mthd.X_sol));
				F.CollectCalc2(Nj, current_gap, Mthd.N, *Mthd.X_sol);
				//F.CollectCalc2(Nj, current_gap, Mthd.N, Mthd.X_pred);

				// updating RetArray(s) pushing X_prev.value of solution
				if (is_SYS_deflecting()) RetDataUpdate(Nj);

				// shifting X_prev next step further in one-step-methods and
				// X[1],X[2]... next step further in multistep-methods
				NodeShift();


			} // end of the cycle processing current gap in <step_method>

			current_gap++;

			// Checking if the any gap remained processed 
			if (current_gap == Mthd.full_amount_of_gaps) {}
			else
			{
				Mthd.X_init = (*Mthd.X_sol);
				Nj = 1;
			}
		}

		F.CountFullResult();
		return F.full_result;
	}
}