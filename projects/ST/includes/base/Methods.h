/* This header constains classes with functions that implement ODE_calculation methods*/

#pragma once
#include <vector>


namespace StraightTask
{
	namespace Methods
	{
		class Parameters {
		public:

			double_t t_0; // initial time moment
			double_t H; // grid step in calc-scheme
			uint32_t N; // amount of nods in a grid

			double_t gap_width;
			uint16_t full_amount_of_gaps; // crucial parameter in terms of dealing with delayed arguments in equations

			// default constructor
			Parameters() :t_0(0.5), N(1500), gap_width(24.0), full_amount_of_gaps(1) { H = gap_width / N; }

			// custom constructor
			Parameters(double_t t0, uint32_t N, uint16_t full_amount_of_gaps)
				:t_0(t0), N(N), full_amount_of_gaps(full_amount_of_gaps) {
				H = gap_width / this->N; 
			}

			void Sync(Parameters& method) { method = *this; }

			void Set(uint32_t N, double_t GapWidth, uint16_t full_amount_of_gaps){
				this->N = N; this->gap_width = GapWidth; this->full_amount_of_gaps = full_amount_of_gaps;
				this->H = gap_width / N;
			}


		};
		class IMethod: public Parameters {
		public:
			IMethod() { X_init.tj = t_0; X_sol = nullptr; }
			~IMethod() = default;

			variables X_init; // initial data storage 
			variables X_prev; // solution on previous step
			variables X_pred; // predicted solution on current step
			variables* X_sol; // a pointer to variables obj with final solution 
		};

		//1-st order numerical Eulers method 
		class Euler : public IMethod {
		public:

			Euler() { X_sol = &X_pred; }
			~Euler() = default;

			double_t Predictor // value returns to in-X_pred-double_t-member field
			(double_t const prev_U, IEqMember<variables const &> & RP) {
				return prev_U + H * RP.Expression(X_prev); 
			}

		};

		// 2-nd order modified numerical Euler method 
		class ModEuler :public Euler {
		public:
			ModEuler() { X_sol = &X_cor; }
			~ModEuler() = default;

			variables X_cor; // corrected solution in Predictor-Corrector scheme

			double_t Corrector // value returns to in-X_cor-double_t-member field
			(double_t const prev_U, IEqMember<variables const &> & RP) {
				return prev_U + H * (RP.Expression(X_prev) + RP.Expression(X_pred)) * 0.5; }
		};

		// 4-th order numerical Runge-Kutta method 
		class RunKut4 : public IMethod {
		public:
			RunKut4() { X_sol = &X_pred; }
			~RunKut4() = default;

			variables X_sub; // additional solution storage for scheme purposes
			double_t k1 = 0.0, k2 = 0.0, k3 = 0.0, k4 = 0.0;

			double_t Predictor // value returns to in-X_pred-double_t-member field
			(double_t & sub_U, 
				double_t const prev_U,
				IEqMember <variables const&> & RP){

				//sub_U is ref on a value of current component from X_sub obj

				k1 = H * RP.Expression(X_prev);

				X_sub.tj = X_prev.tj + H * 0.5F;		sub_U = prev_U + k1 * 0.5;
				k2 = H * RP.Expression(X_sub);

				/*X_sub.tj = X_prev.tj + H * 0.5;*/		sub_U = prev_U + k2 * 0.5;
				k3 = H * RP.Expression(X_sub);

				X_sub.tj = X_prev.tj + H;		sub_U = prev_U + k3;
				k4 = H * RP.Expression(X_sub);

				return prev_U + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
			}
		};

		// 4-th order Gear Method
		class RKGear : public RunKut4 {
		public:

			RKGear() { X_sol = &X_cor; }
			~RKGear() = default;

			std::array<variables,4> X;

			variables X_cor; // corrected solution in Predictor-Corrector scheme

			double_t GPred // value returns to in-X_pred-double_t-member field
			(double_t const U_1,
				double_t const U_2,
				double_t const U_3,
				double_t const prev_U, IEqMember<variables const&>& RP) {
				//return prev_U + H * RP.Expression(X_prev);
				return 4 * H * RP.Expression(X_prev) + (U_1 - 10.0 * prev_U)/3.0 + 6.0 * U_3 - 2.0 * U_2;
			}

			double_t Corrector // value returns to in-X_cor-double_t-member
			(double_t const U_1,
				double_t const U_2,
				double_t const U_3,
				double_t const prev_U,
				IEqMember<variables const &> & RP) {
				return 0.04 * (12.0 * H * RP.Expression(X_pred) + 48.0 * prev_U - 36.0 * U_3 + 16.0 * U_2 - 3.0 * U_1);
			}

		};

		// 4-th order multistep Adams method
		class Adams : public RunKut4 {
		public:
			Adams() { X_sol = &X_pred; }
			~Adams() = default;

			std::array<variables, 4> X;

			// solution 1 time-step back from current is kept in X_prev

			double_t A_Predictor // value returns to in X_pred double_t-member
			(double_t const prev_U,
				IEqMember<variables const &> & RP) {
				return prev_U + H * (55.0 * RP.Expression(X_prev) - 59.0 * RP.Expression(X[3]) + 37.0 * RP.Expression(X[2]) - 9.0 * RP.Expression(X[1])) / 24.0;
			}

		};

		// 4-th order multistep Adams-Bashforth-Moulton method
		class ABM : public Adams {
		public:
			ABM() { X_sol = &X_cor; }
			~ABM() = default;

			variables X_cor; // corrected solution in Predictor-Corrector scheme

			double_t A_Corrector // value returns to in-X_cor-double_t-member
			(double_t const prev_U,
				IEqMember<variables const&>& RP) {
				return prev_U + H * (9 * RP.Expression(X_pred) + 19.0 * RP.Expression(X_prev) - 5.0 * RP.Expression(X[3]) + RP.Expression(X[2])) / 24.0;
			}
		
		};
	}
}