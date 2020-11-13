#pragma once
#include <3rdModelPipe.h>

#include <base/Aberration_base.h>

namespace ReverseTask
{
	class IAggregateControls {
	public:

		IAggregateControls(const char* name_1 = "cyto") :
			one{ name_1 }, two{name_1} {}
		~IAggregateControls() = default;


		double_t full_result;

		void GatherData() {
			one.GatherData(input_dir + "exp/cyto4SPL.txt");
		}

		void CollectCalc(float_t H, float_t Tj, StraightTask::variables& X_sol) {
			one.CollectCalc(H, Tj, X_sol.cy);
		}

		void CountFullResult() {
			full_result = 0;
			full_result += one.CountResult();
		}


		void ResetStates() {
			one.ResetState();
		}


	private:
		Discrete one;
		Continuous two;
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

			CONF_COEFS(0, Cytokines, "p_{Ma_{cy}}", p_Macy, 0.0, 10.0);
			CONF_COEFS(1, Cytokines, "C_{Ma}", C_Ma, 0.03, 10.9);

			CONF_COEFS(2, Cytokines, "p_{Lm_{cy}}", p_Lmcy, 1.0, 10.0);
			CONF_COEFS(3, Cytokines, "C_{Lm}", C_Lm, 0.03, 10.9);

			CONF_COEFS(4, Cytokines, "p_{Ln_{cy}}", p_Lncy, 0.1, 10.5);
			CONF_COEFS(5, Cytokines, "C_{Ln}", C_Ln, 0.03, 10.9);

			CONF_COEFS(6, Cytokines, "e_{cy}", e_cy, 0.03, 10.0);

			CONF_COEFS(7, Cytokines, "l_1", l1, 0.5, 10.0);
			CONF_COEFS(8, Cytokines, "l_2", l2, 0.5, 10.0);
			CONF_COEFS(9, Cytokines, "l_3", l3, 0.5, 10.0);

			CONF_COEFS(10, Cytokines, "A_{btm}", _A, 0.05, 1.0);

#undef CONF_COEFS

			/* we return the number of the last defined coefficient setter;
			* it is rather an individual process and each time is to be set by hand */
			return 10;
		}
	}
}

