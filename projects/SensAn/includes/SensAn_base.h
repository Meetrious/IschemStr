#pragma once
#include<ExpPipe.h>


namespace StraightTask {
	template<typename Method>
	void ISolver<Method>::SolveForSensAn(vector<variables>& RawData) {

		PrepairTheTask();

		RawData[0] = Mthd.X_init;

		uint16_t current_gap = 0;
		double_t Tj;


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

				// collecting solution for the analysis
				RawData[current_gap * Nj + Nj] = *Mthd.X_sol;

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
	}

	class IVariable {
	public:
		IVariable(const char* name) :
			name(name), value(0.0)
		{}
		const char* name;
		double_t value;
	};
}

namespace SensAnalysisTask {

	using namespace StraightTask;

	template<typename T>
	using vector = std::vector<T>;

	template<typename T>
	using matrix = vector<vector<T>>; // */

	class Parameters {
	public:
		Parameters(
			// parameters for SensAnTask
			float_t Var_Percent = 0.1,

			// parameters for StraightTask 

			uint32_t ST_gridPow = 1500,
			double_t ST_gapWidth = 24.0,
			uint16_t ST_amountOfGaps = 1)
			:
			VarPercent(Var_Percent),
			//t_0{ t0 },
			N{ ST_gridPow },
			gap_width{ ST_gapWidth },
			full_amount_of_gaps{ ST_amountOfGaps }

		{}

		// parameters for Sensitivity analysis
		float_t VarPercent;


		// parameters for StraightTask
		uint32_t N; // amount of nods in a grid
		double_t gap_width;
		uint16_t full_amount_of_gaps; // crucial parameter in terms of dealing with delayed(retarded) arguments in equations

	};



	template <typename ST_Method>
	class SensAnTask : public Parameters {
	public:

		/* one and only constructor;
		* we don't want to have a default constructor in this class */
		SensAnTask(Parameters CurrentPars);

		/* the method speaks for itself*/
		void SolveForOutput();

		// default destructor
		~SensAnTask() = default;

		// coefficient collection
		vector<IConstant*> Coefs;

		/* vatiations correction 
		* may be omited in the algorithm,
		* but it is posible that I'd like to output them */
		vector<double_t> variation;

		/* 1st index is for time;
		* this collecion may be omited in the algorithm */
		vector <variables> DefaultSolData;

		/* collection for max-values of the default solution
		* also for keeping names of solution components */
		vector <IVariable> Max;

		/* 1st index is for coefficient, 
		* 2nd - is for time */
		matrix<variables> SolData;

		/* collection for preparatory sensitivity coefficient values;
		* they are to be normalize along t-parameter in terms of l_2 norm */
		matrix<vector<double_t>> PreS;

		/* 1st index is for coefficient 
		* 2nd - is for solution-member */
		matrix<double_t> RS;

		// system agreggator, and num_method of the straight task
		ST_Method STM; 

		// capacity of the Coefs-collection
		size_t amount_of_attributes;

	private:
		size_t ListCoefsInTouch();
		void CalculateMax();
		void ConfigureMax();
	};

	template <typename ST_Method>
	SensAnTask<ST_Method>::SensAnTask(Parameters CurrentPars) :
		Parameters(CurrentPars)
	{
		/* Prepating Sensitivity analysis task 
		*allocation required amount of memory that need 
		*realise for solution */

		amount_of_attributes = ListCoefsInTouch();

		variation.assign(amount_of_attributes, 0.0);
		variation.shrink_to_fit();

		DefaultSolData.assign(N * full_amount_of_gaps + 1, variables());

		SolData.assign(2, std::vector<variables>());
		SolData.shrink_to_fit();

		PreS.assign(2, matrix<double_t>());
		PreS.shrink_to_fit();

		RS.assign(amount_of_attributes, vector<double_t>(9, 0.0)); RS.shrink_to_fit();

		/* at this point we have to know all parameters of the straight task
		* we got them from the arguments of the constructor */
		for (auto& cur : SolData) {
			cur.assign(N * full_amount_of_gaps + 1, variables());
			cur.shrink_to_fit();
		}

		for (auto& cur : PreS) {
			cur.assign(N * full_amount_of_gaps + 1, vector<double_t>(9, 0.0));
			cur.shrink_to_fit();
		}

		// filling container with default solution
		STM.SolveForSensAn(DefaultSolData);
		
		CalculateMax();

	}

	template<typename ST_Method>
	void SensAnTask<ST_Method>::SolveForOutput() {

		uint16_t ST_counter = 0;

		for (size_t j = 0; j < amount_of_attributes; j++) {

			// saving default coef_value for later
			double_t ini_coef_value = Coefs[j]->value; 

			//saving variations for later
			variation[j] = Coefs[j]->value * VarPercent; 

			/* to figure out whether we need to add or substract the variation to current coefficient
			* in order to get max-posible sensitivity value, we do both, 
			* and calculate solution for each case	*/

			Coefs[j]->value = ini_coef_value + variation[j];
			STM.SolveForSensAn(SolData[0]);

			Coefs[j]->value = ini_coef_value - variation[j];
			STM.SolveForSensAn(SolData[1]);

			/* counting amount of times the programme would solve the StraightTask;
			* not really useful in this algorith, and is here just for curiosity purposes*/
			ST_counter += 2;

			Coefs[j]->value = ini_coef_value; // returning default coef_value


		//  defining sensitivity parameters for each solution component at each time moment
#define COUNT_PRECs(NJ, COMP_SOL_NAME, COMP_ID, COEF_ID) \
\
	PreS[0][NJ][COMP_ID] = ((SolData[0][NJ].COMP_SOL_NAME - DefaultSolData[NJ].COMP_SOL_NAME) / variation[COEF_ID]) *\
			(Coefs[COEF_ID]->value / Max[COMP_ID].value);\
\
	PreS[1][NJ][COMP_ID] = ((SolData[1][NJ].COMP_SOL_NAME - DefaultSolData[NJ].COMP_SOL_NAME) / variation[COEF_ID]) *\
			(Coefs[COEF_ID]->value / Max[COMP_ID].value)

			for (size_t n = 0; n <= N * full_amount_of_gaps; n++) {
				COUNT_PRECs(n, nec, 0, j);
				COUNT_PRECs(n, acu_c, 1, j);
				COUNT_PRECs(n, hel, 2, j);

				COUNT_PRECs(n, cy, 3, j);
				COUNT_PRECs(n, adh, 4, j);

				COUNT_PRECs(n, lm, 5, j);
				COUNT_PRECs(n, ln, 6, j);

				COUNT_PRECs(n, mia, 7, j);
				COUNT_PRECs(n, mii, 8, j);
			}


			// caclulating the l_2-norm value along time variable for sensitivity parameters
			for (size_t i = 0; i < 9; i++) {
				double_t SumCollector0 = 0.0;
				double_t SumCollector1 = 0.0;
				for (size_t n = 0; n <= N * full_amount_of_gaps; n++) {
					SumCollector0 += PreS[0][n][i] * PreS[0][n][i];
					SumCollector1 += PreS[1][n][i] * PreS[1][n][i];
				}
				
				/* figuring out, which variant of sign for variation of the parameter was better*/
				if (SumCollector0 < SumCollector1)
					RS[j][i] = std::sqrt(SumCollector1) / (N * full_amount_of_gaps);
				else
					RS[j][i] = std::sqrt(SumCollector0) / (N * full_amount_of_gaps);
			}
		}

		std::cout << ST_counter << " times I solved given StraightTask" << std::endl;
		getchar();

		// outputting data required for analysis
		std::ofstream out;
		out.open(output_dir + "/SensAn/CoefNames.txt", std::ios_base::out);
		for (auto& cur : Coefs) {
			out << cur->name << "\n";
		}
		out.close();

		out.open(output_dir + "/SensAn/SolCompsNames.txt", std::ios_base::out);
		for (auto& cur : Max) {
			out << cur.name << "\n";
		}
		out.close();

		out.open(output_dir + "/SensAn/data.txt", std::ios_base::out);
		for (size_t j = 0; j < amount_of_attributes; j++) {
			for (auto const& cur : RS[j]) {
				out << cur << "\t\t";
			}
			out << std::endl;
		}
		out.close();

	}


	//------------------------------------------------------------------------------------------------------------------------------------------------------------------------


	template<typename ST_Method>
	size_t SensAnTask<ST_Method>::ListCoefsInTouch() {

		// Defining Coeficients

		// for 1st block
		Coefs.emplace_back(&STM.NEUR_C.q_N);
		Coefs.emplace_back(&STM.NEUR_C.q_H);
		Coefs.emplace_back(&STM.NEUR_C.q_epN);
		Coefs.emplace_back(&STM.NEUR_C.q_epA);
		Coefs.emplace_back(&STM.NEUR_C.q_A_n);
		Coefs.emplace_back(&STM.NEUR_C.q_A);
		Coefs.emplace_back(&STM.NEUR_C.p_r);
		//Coefs.emplace_back(&STM.NEUR_C.l1);
		//Coefs.emplace_back(&STM.NEUR_C.l2);
		//Coefs.emplace_back(&STM.NEUR_C.l3);
		Coefs.emplace_back(&STM.NEUR_C.c_H);
		Coefs.emplace_back(&STM.NEUR_C.c_A_h);
		Coefs.emplace_back(&STM.NEUR_C.c_A);

		//for Cytokines
		Coefs.emplace_back(&STM.CYTO_C._A);
		Coefs.emplace_back(&STM.CYTO_C.p_Macy);
		Coefs.emplace_back(&STM.CYTO_C.p_Lncy);
		Coefs.emplace_back(&STM.CYTO_C.p_Lmcy);
		//Coefs.emplace_back(&STM.CYTO_C.l3);
		//Coefs.emplace_back(&STM.CYTO_C.l2);
		//Coefs.emplace_back(&STM.CYTO_C.l1);
		Coefs.emplace_back(&STM.CYTO_C.e_cy);
		Coefs.emplace_back(&STM.CYTO_C.C_Ma);
		Coefs.emplace_back(&STM.CYTO_C.C_Ln);
		Coefs.emplace_back(&STM.CYTO_C.C_Lm);

		// for Adhesion molecules
		Coefs.emplace_back(&STM.ADH_C.o_cy_2);
		Coefs.emplace_back(&STM.ADH_C.o_cy_1);
		Coefs.emplace_back(&STM.ADH_C.e_adh);

		// for leucoMacrophags
		//Coefs.emplace_back(&STM.LM_C.l1);
		Coefs.emplace_back(&STM.LM_C.K_Lm);
		Coefs.emplace_back(&STM.LM_C.d_Lm);
		Coefs.emplace_back(&STM.LM_C.c_Lm);
		Coefs.emplace_back(&STM.LM_C.c_dLm);

		// for leuNeutrophils
		//Coefs.emplace_back(&STM.LN_C.K_Ln2);
		Coefs.emplace_back(&STM.LN_C.K_Ln1);
		Coefs.emplace_back(&STM.LN_C.d_Ln);
		Coefs.emplace_back(&STM.LN_C.c_Ln);
		//Coefs.emplace_back(&STM.LN_C.c_dLn2);
		Coefs.emplace_back(&STM.LN_C.c_dLn1);

		// for Microglia
		Coefs.emplace_back(&STM.MI_C.T_M1);
		Coefs.emplace_back(&STM.MI_C.K_Mi);
		Coefs.emplace_back(&STM.MI_C.c_pro);
		Coefs.emplace_back(&STM.MI_C.c_m_N);
		Coefs.emplace_back(&STM.MI_C.c_dMi);
		Coefs.emplace_back(&STM.MI_C.c_m_A);

		// for ToxDamage subValues
		Coefs.emplace_back(&STM.D_C.c_q1);
		Coefs.emplace_back(&STM.D_C.c_q2);
		Coefs.emplace_back(&STM.D_C.c_q3);
		Coefs.emplace_back(&STM.D_C.c_q4);

		Coefs.emplace_back(&STM.D_C.p_q1);
		Coefs.emplace_back(&STM.D_C.p_q2);
		Coefs.emplace_back(&STM.D_C.p_q3);
		Coefs.emplace_back(&STM.D_C.p_q4);

		Coefs.emplace_back(&STM.D_C.D_0);

		// for Phagocytosis subValues
		Coefs.emplace_back(&STM.EPS_C.e_Mi);
		Coefs.emplace_back(&STM.EPS_C.e_Ma);
		Coefs.emplace_back(&STM.EPS_C.e_Ln);
		Coefs.emplace_back(&STM.EPS_C.e_Lm); // 54 that is, but 45 for 4th model

		Coefs.shrink_to_fit();

		return Coefs.size();
	}

	template<typename ST_Method>
	void SensAnTask<ST_Method>::ConfigureMax() {

		Max.clear();

		Max.emplace_back(IVariable("nec"));
		Max.emplace_back(IVariable("acu\\_c"));
		Max.emplace_back(IVariable("hel"));

		Max.emplace_back(IVariable("cy"));
		Max.emplace_back(IVariable("adh"));

		Max.emplace_back(IVariable("lm"));
		Max.emplace_back(IVariable("ln"));

		Max.emplace_back(IVariable("mia"));
		Max.emplace_back(IVariable("mii"));

		Max.shrink_to_fit();
	}

	template<typename ST_Method>
	void SensAnTask<ST_Method>::CalculateMax() {
		ConfigureMax();

		for (size_t n = 0; n <= N * full_amount_of_gaps; n++) {
			if (DefaultSolData[n].nec > Max[0].value) Max[0].value = DefaultSolData[n].nec;
			if (DefaultSolData[n].acu_c > Max[1].value) Max[1].value = DefaultSolData[n].acu_c;
			if (DefaultSolData[n].hel > Max[2].value) Max[2].value = DefaultSolData[n].hel;

			if (DefaultSolData[n].cy > Max[3].value) Max[3].value = DefaultSolData[n].cy;
			if (DefaultSolData[n].adh > Max[4].value) Max[4].value = DefaultSolData[n].adh;

			if (DefaultSolData[n].lm > Max[5].value) Max[5].value = DefaultSolData[n].lm;
			if (DefaultSolData[n].ln > Max[6].value) Max[6].value = DefaultSolData[n].ln;

			if (DefaultSolData[n].mia > Max[7].value) Max[7].value = DefaultSolData[n].mia;
			if (DefaultSolData[n].mii > Max[8].value) Max[8].value = DefaultSolData[n].mii;
		}


	}


}