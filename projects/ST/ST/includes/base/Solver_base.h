/* This header serves as pre-surface shell for implementation of common functions that manifest solution algorithm;
* main points of interest are ISolver public methods (void Solve*...*(-//-));
* a little bit less important, but not it the least (void Prepair*(-//-)) methods;
* the class ISolver is to be publicly inherited in class <Method> in *Pipe.h headers, 
where presented definitions of all functions(relatively to each ODE system)
required in algorythm stated in (void Solve*...*(-//-)) // */


#pragma once
#include <base/Settings_base.h>

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

	// is defined somewhere above throught inclusion of this header in another one
	class IAggregate;

	template<typename Method>
	class ISolver : public IAggregate {

	public:
		Method ST; // ST is for Straight Task

		// Solves the straight task and outputs solution
		void SolveAndOutput(uint32_t N, double_t GapWidth, uint16_t full_amount_of_gaps) {

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
		void SolveForRungeAnalysis() {
			return;
		}


		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };


	protected:
		void InitialiseIniData() {
			vector<double_t> T0s;
			/*	for (auto const& cur : IniDataInitialiser) {
					cur(ST.X_init); T0s.emplace_back(ST.X_init.tj);
				}*/
			SetIniData(T0s);
			double_t mean = 0;
			for (size_t i = 0; i < T0s.size(); i++) {
				mean += T0s[0];
			}
			mean = mean / T0s.size();

			ST.X_init.tj = ST.t_0 = mean;

			bool trg = false;

			// making sure, that each initial datum is set in the same t0 
			for (size_t i = 0; i < T0s.size(); i++) {
				if (mean == T0s[i])
					continue;
				else trg = true;
			}
			if (trg)
			{
				std::cerr << "\n WARNING: initial data for the ODE system is not consistent:"
					<< "\n initial time moments are not the same: \n {_";
				for (auto const& cur : T0s) std::cerr << cur << '_';
				std::cerr << "} \n t_0 will be defined as mean."
					<< "\n Do you wish to proceed? :\n  0. NO; \n 1. YES;" << std::endl;
				std::cin >> trg;
				if (!trg) throw(" Initial data is inacceptable. ");
				else std::cerr << "Approved \n Proceeding. \n";
			}

			ExpressSubValues(ST.X_init);// */

			return;
		}

		void PrepairTheTask()
		{
			try {
				InitialiseIniData();

				// collecting presolved solution to freeze the system relatively given behaviour
				//CollectData() // full
				//DataCollector[3]();

				// setting initial data from presolved_solution_data
				//AssignSolData(ST.X_init);// full
				//SolDataGetter[3](0, 0, ST.N, ST.X_init);

				if (is_SYS_deflecting())
				{
					// initialise ret-value storage
					//for (auto const& cur : RetInitialiser) { cur(ST.N, ST.t_0, ST.gap_width); }
					InitialiseRetArrays();

					// setting first ret-values required on the first step of calculation
					//for (auto const& cur : IniRetInitialiser) { cur(ST.X_init); }
					InitialiseIniRetValues();
				}

			}
			catch (const char* exception) {
				std::cerr << "WARNING:" << exception << "\n Terminating.";
				throw(exception);
			}

		}

		void PrepairTheOutput()
		{
			try {
				// streams for output
				// for (auto const& cur : OutStreamAllocator) { cur(); }
				AllocateOutputStreams();

				// outputting initial solution data
				// for (auto const& cur : SolutionOutputter) { cur(0, ST.X_init.tj, ST.X_init); }
				OutputSolution(0, ST.X_init.tj, ST.X_init);
			}
			catch (const char* exception) {
				std::cerr << "WARNING:" << exception << "\n Terminating.";
				throw(exception);
			}
		}

		/* other members listed below is to be defined or redefined in
		particular cases of ODE system somewhere down in another header */

		bool is_SYS_deflecting();

		virtual void NodeShift() { ST.X_prev = *ST.X_sol; }

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

		virtual void ApplyPrepMethod() { return; }
	};
}