#pragma once
#include <ExpPipe.h>

#include <base/RT_exceptions.h>

#include <base/Aberration_base.h>

namespace ReverseTask
{
	class IAggregateControls {
	public:

		IAggregateControls(const char* defaultname = "mem") :
			NEC("N"),
			AC("AC"),// */
			
			CY("CY"),
			LN("Ln")
			//LM("LM") 
		{}
		~IAggregateControls() = default;


		double_t full_result = 0;

		void GatherData() {
			AC.GatherData(input_dir + "preserved_solution/2/Ac_ch.txt");
			NEC.GatherData(input_dir + "preserved_solution/2/Necr.txt");
			CY.GatherData(input_dir + "preserved_solution/2/Cy.txt");
			LN.GatherData(input_dir + "preserved_solution/2/Ln.txt");
			//LM.GatherData(input_dir + "preserved_solution/2/Lm.txt");// */
		}

		void CollectCalc1(float_t H, float_t Tj, StraightTask::variables& X_sol) {
			//one.CollectCalc(H, Tj, X_sol.cy);

		}

		void CollectCalc2(uint32_t Nj, uint16_t gap, uint32_t N, StraightTask::variables& X_sol) {

			AC.CollectCalc(Nj, gap, N, X_sol.acu_c);
			NEC.CollectCalc(Nj, gap, N, X_sol.nec); // */

			CY.CollectCalc(Nj, gap, N, X_sol.cy);
			
			LN.CollectCalc(Nj, gap, N, X_sol.ln);
			//LM.CollectCalc(Nj, gap, N, X_sol.lm);// */
		}

		void CountFullResult() {
			full_result = 0.0;
		
#ifdef _DEBUG
			//std::cout << "\n CY_res = " << CY.CountResult() << std::endl;
			//std::cout << " LN_res = " << LN.CountResult() << std::endl;
#endif // _DEBUG

			full_result += NEC.CountResult();
			full_result += AC.CountResult(); // */

			full_result += CY.CountResult();
			full_result += LN.CountResult();
			//full_result += LM.CountResult();// */
		}


		void ResetStates() {
			AC.ResetState();
			NEC.ResetState(); // */

			CY.ResetState();
			
			LN.ResetState();
			//LM.ResetState();// */
		}


	private:

		Continuous NEC;
		Continuous AC; // */
		//Continuous HEL;

		Continuous CY;
		Continuous LN;
		//Continuous LM;// */
	};
}

#include <base/ST_SolverForBGA.h>

#include <base/BGA_Base.h>


namespace ReverseTask
{
	namespace BGA
	{
		template <typename ST_Method>
		size_t Task<ST_Method>::SetCVlist() {

#define CONF_COEFS(NUM, COEF_PACK_NAME, COEF_NAME, L, R) \
			for (auto & Worker: Workers) { \
				Worker.CoefsToVariate[NUM] = &Worker.STM.COEF_PACK_NAME.COEF_NAME; \
			}\
			bounds[0].emplace_back(L);\
 			bounds[1].emplace_back(R)

			CONF_COEFS(0, NEUR_C, q_N, 0.02, 4.0);
			CONF_COEFS(1, NEUR_C, q_A_n, 6e-2, 4.0);
			CONF_COEFS(2, NEUR_C, c_H, 2e-3, 5.0);
			CONF_COEFS(3, NEUR_C, c_A, 1e-2, 3.0);
			CONF_COEFS(4, NEUR_C, q_epN, 2e-3, 1.5);
			CONF_COEFS(5, NEUR_C, q_epA, 1e-1, 1.0);
			//CONF_COEFS(6, NEUR_C, q_A, 6e-1, 4.0);
			//CONF_COEFS(7, NEUR_C, c_A_h, 1e-3, 3.0);
			//CONF_COEFS(8, NEUR_C, p_r, 2e-2, 2.0);
			CONF_COEFS(6, NEUR_C, q_H, 7.0, 16.0);
			//CONF_COEFS(7, NEUR_C, l1, 0.5, 4.0); // */
			
			//CONF_COEFS(7, LN_C, c_Ln, 1e-3, 1.5);
			CONF_COEFS(7, LN_C, d_Ln, 10.0, 16.0);
			CONF_COEFS(8, LN_C, c_dLn1, 1.0, 6.0);
			CONF_COEFS(9, LN_C, K_Ln1, 4.5, 9.0);
			//CONF_COEFS(11, LN_C, l1, 0.5, 3.0);

			CONF_COEFS(10, CYTO_C, p_Macy, 1e-2, 1.0);
			CONF_COEFS(11, CYTO_C, C_Ma, 2.0, 8.0);
			//CONF_COEFS(7, CYTO_C, l1, 1.0, 6.0);
			CONF_COEFS(12, CYTO_C, _A, 3e-3, 0.5);
			CONF_COEFS(13, CYTO_C, p_Lmcy, 3.0, 8.0);
			CONF_COEFS(14, CYTO_C, C_Lm, 1.0, 6.0);
			//CONF_COEFS(11, CYTO_C, l2, 1.0, 5.0);
			CONF_COEFS(15, CYTO_C, p_Lncy, 6.5, 11.0);
			CONF_COEFS(16, CYTO_C, C_Ln, 1e-2, 2.0);
			//CONF_COEFS(14, CYTO_C, l3, 1.0, 5.0);
			CONF_COEFS(17, CYTO_C, e_cy, 1.0, 4.0);// */

			/*CONF_COEFS(19, LM_C, c_Lm, 0.1, 3.0);
			CONF_COEFS(20, LM_C, c_dLm, 3.0, 7.0);
			CONF_COEFS(18, LM_C, K_Lm, 5.0, 15.0);
			CONF_COEFS(19, LM_C, d_Lm, 1e-2, 1.0); // */
			//CONF_COEFS(20, LM_C, l1, 1.0, 5.0);


#undef CONF_COEFS

			/* we return the number of the last defined coefficient setter;
			* it is rather an individual process and each time is to be set by hand */
			return bounds[0].size() - 1;


		}		
	}
}

