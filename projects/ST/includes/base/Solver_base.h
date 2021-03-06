/* This header serves as pre-surface shell for implementation of common functions that manifest solution algorithm;
* main points of interest are ISolver public methods (void Solve*...*(-//-));
* a little bit less important, but not it the least (void Prepair*(-//-)) methods;
* the class ISolver is to get publicly inherited in class <Method> in *Pipe.h headers, 
	where presented definitions of all declared functions(relatively to each ODE system)
	required in algorithm stated in (void Solve*...*(-//-)) // */


#pragma once
#include <base/Settings_base.h>

// #define ST_COUNT_TIME_ELAPSED

#ifdef ST_COUNT_TIME_ELAPSED
#include <base/Timer.h>
#endif

// namespace and class declaration for reverse-task purposes
namespace ReverseTask {
	class IAggregateControls;
}

namespace StraightTask
{
	class PhaseTrajOutput
	{
	public:
		std::ofstream str;
		PhaseTrajOutput(const char* name) { AllocateOutputStream(name); }

		void AllocateOutputStream(const char* name)
		{
			str.open(output_dir + "ST/SOL/PTraj/" + name + ".txt");
			if (!str)
			{
				std::cout << "\n For some reason budgets output stream for <" << name << "> wasn't allocated \n";
				getchar();
				return;
			}
		}
		void OutputPhaseTraj(uint32_t Nj, double_t Tj, vector<double_t> Vals) {
			str << Nj << "\t\t\t" << std::setprecision(5) << Tj << "\t\t\t" << std::setprecision(15);
			for (auto const& cur : Vals)
				str << cur << "\t\t\t";
			str << std::endl;
		}

		~PhaseTrajOutput() { if (str.is_open()) str.close(); }
	};

	// is defined somewhere above through inclusion of this header in another one
	class IAggregate;

	template<typename Method>
	class ISolver : public IAggregate {

	public:
		Method Mthd;

		// Solves the straight task and outputs solution
		void SolveAndOutput(uint32_t N, double_t GapWidth, uint16_t full_amount_of_gaps) {

#ifdef ST_COUNT_TIME_ELAPSED
			Timer dot;
#endif

			// setting an object that will output 2 or 3 dimentional phase trajectory for solution 
			// PhaseTrajOutput PTO("APPROX");

			// Setting parameters for calc-method:
			Mthd.Set(N, GapWidth, full_amount_of_gaps);

			PrepairTheTask();

			PrepairTheOutput();

			uint16_t current_gap = 0;
			double_t Tj;

#ifdef ST_COUNT_TIME_ELAPSED
			std::cout << "\n Preparation took somewhere about: ";
			dot.CountInterval();
#endif // ST_COUNT_TIME_ELAPSED


			while (current_gap < full_amount_of_gaps)
			{

				std::cout << full_amount_of_gaps - current_gap << " "; // outputting amount of gaps remained to process
				Mthd.X_prev = Mthd.X_init; // previous step is defined by initial one at every begining of the gap
				Tj = Mthd.X_pred.tj = Mthd.X_prev.tj; // synchronising the independent variable
				uint32_t Nj = 1; // restating the number of the current step in num-method

				// 1st approximations for multistep-methods
				if (current_gap == 0) ApplyPrepStep(Nj, Tj); // for one-step-methods it is void{return;}

				// cycle for processing current gap in <step_method>
				for (; Nj <= N; Nj++)
				{
					// shifting independent variable on one step further for predicted solution
					Tj = Mthd.X_pred.tj += Mthd.H;

					// setting ret-values for Tj-time-moment in X_pred
					if (is_SYS_deflecting()) RetUpload(Nj);

					// interoplation gap monitor
					//Mthd.X_pred.CheckShiftInterpGap(STpar);

					// defining solution with presolved data
					AssignSolData(current_gap, Nj, Mthd.X_pred);

					ApplyMethod();

					// outputting solution in current Tj - time-moment
					//for (auto const& cur : SolutionOutputter) { cur(Nj, Tj, (*Mthd.X_sol)); }
					OutputSolution(Nj, Tj, (*Mthd.X_sol));

					// Phase trajectory output
					//PTO.OutputPhaseTraj(Nj, Tj, { (*SYS.Mthd.X_sol).x, (*SYS.Mthd.X_sol).y, (*SYS.Mthd.X_sol).z });


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
				if (current_gap == Mthd.full_amount_of_gaps) {
					std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";

#ifdef COUNT_TIME_ELAPSED
					std::cout << "\n Solution took somewhere about: ";
					dot.CountInterval();
#endif // COUNT_TIME_ELAPSED
				}

				else
				{
					Mthd.X_init = (*Mthd.X_sol);
					Nj = 1;
				}
			}
		}

		// Solves the straight task in order to calculate aberration value in BGA
		double_t SolveForBGA(ReverseTask::IAggregateControls& F);

		// Solves the straight task for sensibility analysis
		void SolveForSensAn(vector<variables>& RawData);

		// Solves the straight task three times and calculates the order of assignated num-method using Runge-rule
		void SolveForRungeAnalysis() {
			// some day in another life
			return;
		}


		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };


	
		void InitialiseIniData() {
			vector<double_t> T0s;

			SetIniData(T0s);
			//SetIniDataFromOutside(T0s);
			double_t mean = 0;
			for (size_t i = 0; i < T0s.size(); i++) {
				mean += T0s[0];
			}
			mean = mean / T0s.size();

			Mthd.X_init.tj = Mthd.t_0 = mean;

			bool trg = false;

			// making sure, that each initial datum is set in the same t0 
			for (size_t i = 0; i < T0s.size(); i++) {
				if (mean == T0s[i])
					continue;
				else trg = true;
			}
			if (trg)
			{
				std::cout << "\n WARNING: initial data for the ODE system is not consistent:"
					<< "\n initial time moments are not the same: \n {_";
				for (auto const& cur : T0s) std::cout << cur << '_';
				std::cout << "} \n t_0 will be defined as mean value over listed ones."
					<< "\n Do you wish to proceed? :\n  0. NO; \n 1. YES; \n ans: " << std::endl;
				std::cin >> trg;
				if (!trg)
				{
					inconsistent_initial_data const exception;
					throw(exception);
				}
				else std::cout << "Approved \n Proceeding. \n";
				
			}
			T0s.clear();
			ExpressSubValues(Mthd.X_init);

		}

		void PrepairTheTask()
		{
			try {

				// prepairing equations
				SynchronizeAllCoefs();

				InitialiseIniData();

				// collecting presolved solution to freeze the system relatively given behaviour
				CollectData(); // full
				
				// setting initial data from presolved_solution_data
				AssignSolData(0, 0, Mthd.X_init);// full

				//X_0 = Mthd.X_init;

				if (is_SYS_deflecting())
				{
					// initialise ret-value storage
					
					InitialiseRetArrays();

					// setting first ret-values required on the first step of calculation
					//for (auto const& cur : IniRetInitialiser) { cur(Mthd.X_init); }
					InitialiseIniRetValues();
				}

			}
			catch (const std::exception & exception) {
				std::cout << "\n Terminating.";
				throw( &exception );
			}

		}

		void PrepairTheOutput()
		{
			try {
				// streams for output
				AllocateOutputStreams();

				// outputting initial solution data
				OutputSolution(0, Mthd.X_init.tj, Mthd.X_init);
			}
			catch (std::exception & exception) {
				std::cout << "\n Closing all output streams ";
				DeallocateOutputStreams();
				throw(exception);
			}
		}

		/* other members listed below is to be defined or redefined in
		particular cases of ODE system somewhere down in another header */

		bool is_SYS_deflecting();

		protected:

		virtual void NodeShift() { Mthd.X_prev = *Mthd.X_sol; }

		void RetUpload(uint32_t Nj);
		void RetDataUpdate(uint32_t Nj);

		virtual void ApplyMethod() = 0;

		void CollectData();

		void SetIniData(vector<double_t>& T0s);
		void SetIniDataFromOutside(vector<double_t>& T0s);

		void AssignSolData(uint16_t day, uint32_t Nj, variables& X);
		void ExpressSubValues(variables& X);

		void InitialiseRetArrays();
		void InitialiseIniRetValues();

		void AllocateOutputStreams();
		void OutputSolution(uint32_t Nj, double_t Tj, variables const& X);
		void OutputBudgets(uint32_t Nj, double_t Tj);
		void DeallocateOutputStreams();

		void SynchronizeAllCoefs();
		

		virtual void ApplyPrepMethod() { return; }
	};
}