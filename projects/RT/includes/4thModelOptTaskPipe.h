#pragma once
#include <4thModelPipe.h>

#include <base/Aberration_base.h>

namespace ReverseTask
{
	class IAggregateControls {
	public:

		IAggregateControls(const char* defaultname = "mem") :
			NEC("N"),
			AC("AC")// */
			
			//CY("CY"),
			//LN("Ln"),
			//LM("LM") 
		{}
		~IAggregateControls() = default;


		double_t full_result = 0;

		void GatherData() {
			AC.GatherData(input_dir + "preserved_solution/2/Ac_ch.txt");
			NEC.GatherData(input_dir + "preserved_solution/2/Necr.txt");
			//CY.GatherData(input_dir + "preserved_solution/2/Cy.txt");
			//LN.GatherData(input_dir + "preserved_solution/2/Ln.txt");
			//LM.GatherData(input_dir + "preserved_solution/2/Lm.txt");// */
		}

		void CollectCalc1(float_t H, float_t Tj, StraightTask::variables& X_sol) {
			//one.CollectCalc(H, Tj, X_sol.cy);

		}

		void CollectCalc2(uint32_t Nj, uint16_t gap, uint32_t N, StraightTask::variables& X_sol) {

			AC.CollectCalc(Nj, gap, N, X_sol.acu_c);
			NEC.CollectCalc(Nj, gap, N, X_sol.nec); // */

			//CY.CollectCalc(Nj, gap, N, X_sol.cy);
			
			//LN.CollectCalc(Nj, gap, N, X_sol.ln);
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

			//full_result += CY.CountResult();
			//full_result += LN.CountResult();
			//full_result += LM.CountResult();// */
		}


		void ResetStates() {
			AC.ResetState();
			NEC.ResetState(); // */

			//CY.ResetState();
			
			//LN.ResetState();
			//LM.ResetState();// */
		}


	private:

		Continuous NEC;
		Continuous AC; // */
		//Continuous HEL;

		//Continuous CY;
		//Continuous LN;
		//Continuous LM;// */
	};
}

#include <base/ST_SolverForBGA.h>

#include <base/BGA_Base.h>

namespace ReverseTask
{
	namespace BGA
	{
		size_t Parameters::SetCTVlist() {
#define CONF_COEFS(NUM, EQ_BLOCK, NAME, ID, L, R) \
			CoefsToVariate[NUM] = std::make_pair(NAME, &StraightTask::EQ_BLOCK::ID);\
			bounds[0].emplace_back(L); bounds[1].emplace_back(R)

			CONF_COEFS(0, Neurons, "c_A", c_A, 6e-3, 1.0);
			CONF_COEFS(1, Neurons, "q_{A_n}", q_A_n, 2e-1, 2);
			CONF_COEFS(2, Neurons, "q_{ep_N}", q_epN, 5e-2, 0.7);
			CONF_COEFS(3, Neurons, "l_1", l1, 0.5, 3.0);
			//CONF_COEFS(4, ToxDamage, "D_0", D_0, 1e-3, 0.1);// */

		/*	CONF_COEFS(0, Neurons, "q_N", q_N, 1e-3, 2);
			CONF_COEFS(1, Neurons, "q_{A_n}", q_A_n, 2e-1, 2);
			CONF_COEFS(2, Neurons, "c_H", c_H, 1e-3, 1);
		
			CONF_COEFS(5, Neurons, "q_{ep_A}", q_epA, 1e-1, 0.7);
			CONF_COEFS(6, Neurons, "q_A", q_A, 0.5, 1.5);
			CONF_COEFS(7, Neurons, "c_{A_h}", c_A_h, 1e-3, 0.5);
			CONF_COEFS(8, Neurons, "q_H", q_H, 7.0, 9.0);
			


			/*CONF_COEFS(0, Cytokines, "p_{Ma_{cy}}", p_Macy, 5.0, 8.0);
			CONF_COEFS(1, Cytokines, "C_{Ma}", C_Ma, 4.0, 5.5);

			CONF_COEFS(2, Cytokines, "p_{Lm_{cy}}", p_Lmcy, 6.7, 9.2);
			CONF_COEFS(3, Cytokines, "C_{Lm}", C_Lm, 80.0, 110.0);

			CONF_COEFS(4, Cytokines, "p_{Ln_{cy}}", p_Lncy, 5.5, 7.5);
			CONF_COEFS(5, Cytokines, "C_{Ln}", C_Ln, 1e-3, 1.0);

			CONF_COEFS(6, Cytokines, "e_{cy}", e_cy, 2.5, 6.0);*/

			//CONF_COEFS(7, Cytokines, "A_{btm}", _A, 1e-3, 10.0);
			
			/*CONF_COEFS(13, Cytokines, "l_1", l1, 0.5, 10.0);
			CONF_COEFS(14, Cytokines, "l_2", l2, 0.5, 10.0);
			CONF_COEFS(15, Cytokines, "l_3", l3, 0.5, 10.0);// */

			

			/*CONF_COEFS(7, Adhesion, "o_{cy_1}", o_cy_1, 0.003, 1.7);
			CONF_COEFS(8, Adhesion, "o_{cy_2}", o_cy_2, 0.003, 1.5);
			CONF_COEFS(9, Adhesion, "e_{adh}", e_adh, 0.003, 1.2);// */

			/*CONF_COEFS(0, LeuNeutrophils, "c_{Ln}", c_Ln, 0.03, 30.0);
			CONF_COEFS(1, LeuNeutrophils, "d_{Ln}", d_Ln, 0.03, 30.0);
			CONF_COEFS(2, LeuNeutrophils, "c_{d_{Lm_1}}", c_dLn1, 0.03, 30.0);
			CONF_COEFS(3, LeuNeutrophils, "K_{Ln_1}", K_Ln1, 0.03, 30.0);
			CONF_COEFS(4, LeuNeutrophils, "c_{d_{Ln_2}}", c_dLn2, 0.03, 30.0);
			CONF_COEFS(5, LeuNeutrophils, "K_{Ln_2}", K_Ln2, 0.03, 30.0); // */


			/*CONF_COEFS(0, LeuMacrophags, "c_{Lm}", c_Lm, 0.005, 10.0);
			CONF_COEFS(1, LeuMacrophags, "c_{d_{Lm}}", c_dLm, 0.5, 10.0);
			CONF_COEFS(2, LeuMacrophags, "K_{Lm}", K_Lm, 0.003, 10.0);
			CONF_COEFS(3, LeuMacrophags, "d_{Lm}", d_Lm, 0.003, 10.0);
			CONF_COEFS(4, LeuMacrophags, "l_{1}", l1, 4.0, 15.0); //*/


			

			/*CONF_COEFS(14, ToxDamage, "p_q1", p_q1, 1e-3, 2);
			CONF_COEFS(15, ToxDamage, "c_q1", c_q1, 5, 15);*/


#undef CONF_COEFS

			/* we return the number of the last defined coefficient setter;
			* it is rather an individual process and each time is to be set by hand */
			return bounds[0].size() - 1;
		}
		
	}
}

