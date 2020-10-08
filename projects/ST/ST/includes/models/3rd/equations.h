#pragma once

#include <models/3rd/variables.h>
#include <base/Eq_base.h>
#include <cmath>


double_t pi = 3.14159265358979323;

namespace StraightTask
{
	template<typename T>
	using vector = std::vector<T>;

	template<typename T>
	using matrix = vector<vector<T>>;


	[[nodiscard]] inline double_t Hill (double_t C, double_t ed, double_t K, int32_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	[[nodiscard]] inline double_t CHill (double_t C, double_t ed, double_t K, double_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	[[nodiscard]] inline double_t GetPositive (double_t val) {
		if (val > 0.0) return val; else return 0;
	}

	[[nodiscard]] inline double_t indicator(double_t val, double_t border){
		if (val > border) return 1.0; else return 0;
	}


	template<typename argType> class IEquation: public IRightPart<argType>{};


	namespace Neurons
	{
		static double_t
			q_N = 1.5, q_A_n = 0.406, c_H = 0.006, c_A = 1.2e-2, q_epN = 0.15, q_epA = 0.001,
			q_A = 1.0, q_A_h = 0.006, p_r = 0.1, q_H = 8.13;


		namespace NecroticCells
		{
			class Equation : public IEquation<variables const&>{
			public:
				Equation() { B.insert(B.end(), 4, 0.0); B.shrink_to_fit(); ret_is = false; }
				const double_t ini_data[2] = { 0.5, 0.0 };
				[[nodiscard]] inline double_t Expression(variables const & u)noexcept final{
					B[1] = q_N * u.dp_N * Hill(1.0, u.hel, c_H, 2) * u.nec;
					B[2] = q_A_n * u.dp_N * Hill(1.0, u.acu_c, c_A, 5);
					B[3] = -q_epN * u.eps_s * u.nec;
					return Sum_B();
				}
			};
		}

		namespace AcuteChanges
		{
			class Equation : public IEquation<variables const&> {
			public:
				Equation() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false;}
				//const double_t ini_data[2] = { 0.5, 0.2890 };
				const double_t ini_data[2] = { 0.5, 0.01 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = (q_A - p_r) * u.dp_A * u.hel * Hill(1.0, u.acu_c, q_A_h, 1);
					B[2] = q_H * u.dp_N * u.hel * u.acu_c;
					B[3] = -q_A_n * u.dp_N * Hill(1.0, u.acu_c, c_A, 5);
					B[4] = -q_epA * u.eps_s * u.acu_c;;
					return Sum_B();
				}
			};
		}

		namespace IntactCells 
		{
			class Equation : public IEquation<variables const&> {
			public:
				Equation() { B.insert(B.end(), 4, 0.0); B.shrink_to_fit(); ret_is = false; }
				//const double_t ini_data[2] = { 0.5, 0.7110 };
				const double_t ini_data[2] = { 0.5, 0.99 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = (-1.0 + p_r) * u.hel * q_A * u.dp_A * Hill(1.0, u.acu_c, q_A_h, 1);
					B[2] = -q_N * u.dp_N * u.nec * Hill(1.0, u.hel, c_H, 2);
					B[3] = -q_H * u.dp_N * u.hel * u.acu_c;
					return Sum_B();
				}
			};
		}

	} // namespace Neurons

	namespace Cytokines
	{
		// chemotaxis
		//static double_t p_Mach = 4.5, C_Ma = 0.4, p_Lmch = 4.5, C_Lm = 0.18, e_ch = 0.18;
		static double_t C_Ma = 9.9713;
		// ordinal
		//static double_t p_Macy = 2.9481, p_Lmcy = 5.6667, p_Lncy = 8.2513, C_Ln = 3.3918, e_cy = 5.577, l1 = 0.75752, l2 = 5.1694, l3 = 1.3089, _A = 0.14514;
		static double_t p_Macy = 3.7386, p_Lmcy = 4.9559, p_Lncy = 6.2865,
			e_cy = 1.166, l1 = 0.75752, l2 = 5.1694, l3 = 1.3089, _A = 0.014514;																					

		namespace Pro_Inflam
		{
			class Equation : public IEquation<variables const&> {
			public:
				Equation() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false;}
				//const double_t ini_data[2] = { 0.5, 0.0906 };
				const double_t ini_data[2] = { 0.5, 0.0 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = CHill(p_Macy, u.mia, C_Ma, l1) * (u.nec + GetPositive(u.acu_c - _A));
					B[2] = p_Lmcy * u.lm * (u.nec + u.acu_c);
					B[3] = p_Lncy * u.ln * (u.nec + u.acu_c);
					B[4] = - e_cy * u.cy;
					
					return Sum_B();
				}
			};
		}
	} // namespace Cytokines

	namespace Adhesion 
	{
		static double_t 
			o_cy_1 = 2.62e-2, o_cy_2 = 3.14e-2, e_adh = 1.31e-2;

		class Equation : public IEquation<variables const&>{
		public:
			Equation() { B.insert(B.end(), 4, 0.0); B.shrink_to_fit(); ret_is = false; }
			const double_t ini_data[2] = { 0.5, 0.0 };
			//const double_t ini_data[2] = { 0.5, 0.2 };
			[[nodiscard]] double_t Expression(variables const& u)noexcept final{
				B[1] = o_cy_1 * u.cy;
				B[2] = o_cy_2 * u.adh * u.cy;
				B[3] = -e_adh * u.adh;
				return Sum_B();
			}
		};

	} // namespace Adhesion

	namespace LeuMacrophags
	{
		static double_t
			c_Lm = 9.63e-2, c_dLm = 0.08, K_Lm = 43.9, d_Lm = 0.139;

		class Equation : public IEquation<variables const&> {
		public:
			Equation() { B.insert(B.end(), 4, 0.0); B.shrink_to_fit(); ret_is = true; }
			const double_t ini_data[2] = { 0.5, 0.0 };
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1]= c_Lm * u.ret.adh_4;
				B[2] = Hill(c_dLm, u.cy, K_Lm, 4) * u.lm;
				B[3] = -d_Lm * u.lm;
				return Sum_B();
			}
		};

	} // namespace LeuMacrophags

	namespace LeuNeutrophils
	{
		static double_t
			c_Ln = 1.47,  d_Ln = 84.0,
			c_dLn1 = 16.2, K_Ln1 = 3.96,
			c_dLn2 = 4.43, K_Ln2 = 9.899;

		class Equation : public IEquation<variables const&>{
		public:
			Equation() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false; }
			//const double_t ini_data[2] = { 0.5, 0.0861204 };
			const double_t ini_data[2] = { 0.5, 0.0 };
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Ln * u.adh * u.ln;
				B[2] = Hill(c_dLn1, u.cy, K_Ln1, 1);
				B[3] = Hill(c_dLn2, u.cy, K_Ln2, 1) * u.ln;
				B[4] = -d_Ln * u.eps_w * u.ln;
				return Sum_B();
			}
		};


	} // namespace LeuNeutrophils

	namespace Microglia
	{
		static double_t
			c_A = 0.06, c_N = 2.1, c_pro = 2.5,
			T_M1 = 60.0,  
			c_dMi = 0.2, K_Mi = 0.1;
	
		namespace Active
		{
			class Equation : public IEquation<variables const&>{
			public:
				Equation() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false; }
				//const double_t ini_data[2] = { 0.5, 0.01 };
				const double_t ini_data[2] = { 0.5, 0.001 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = u.mii * c_A * u.acu_c;
					B[2] = u.mii * c_N * u.nec;
					B[3] = Hill(c_dMi, u.cy, K_Mi, 1) * u.mii;
					B[4] = -c_pro * u.mia / T_M1;
					return Sum_B();
					//return (c_A * u.acu_c + c_N * u.nec) * u.mii + Hill(c_dMi, u.cy, K_Mi, 1) * u.mii - c_pro * u.mia / T_M1;
				}
			};

		}
		namespace Inactive
		{
			class Equation : public IEquation<variables const&>
			{
			public:
				Equation() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false; }
				//const double_t ini_data[2] = { 0.5, 0.99 };
				const double_t ini_data[2] = { 0.5, 0.999 };
				[[nodiscard]] inline double_t Expression(variables const& u) noexcept final {
					B[1] = - u.mii * c_A * u.acu_c;
					B[2] = - u.mii * c_N * u.nec;
					B[3] = - Hill(c_dMi, u.cy, K_Mi, 1) * u.mii;
					B[4] = c_pro* u.mia / T_M1;
					return Sum_B();
					//return  -(c_A * u.acu_c + c_N * u.nec) * u.mii - Hill(c_dMi, u.cy, K_Mi, 1) * u.mii + c_pro * u.mia / T_M1;
				}
			};
		}
	} // namespace Microglia

	template<typename argType>	class ISubVal : public IRightPart<argType> {};

	namespace ToxDamage
	{

		static double_t
			p_q1 = 0.1, p_q2 = 0.4, p_q3 = 0.4, p_q4 = 0.1,
			c_q1 = 13.206, c_q2 = 1.32, c_q3 = 0.1, c_q4 = 0.1,

			D_0 = 0;


		namespace Full
		{
		
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = Hill(p_q1, u.cy, c_q1, 1);
					B[2] = Hill(p_q2, u.ln, c_q2, 1);
					B[3] = Hill(p_q3, u.lm, c_q3, 1);
					B[4] = Hill(p_q4, u.mia, c_q4, 1);
					return Sum_B();
				}
			};

		}

		namespace Initial
		{
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { ret_is = false; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}
			};
		}

		namespace Nec_partial
		{
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { ret_is = false; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}
			};

		}

		namespace Apop_partial
		{
	
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { ret_is = false; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F);
				}
			};
		}

	} // namespace ToxDamage

	namespace Phagocytosis
	{
		static double_t 
			e_Ma = 0.05, 
			e_Mi = 1.0e-6,
			e_Ln = 5.68e-3,
			e_Lm = 0.05;

		namespace Strong
		{
			
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { B.insert(B.end(), 5, 0.0); B.shrink_to_fit(); ret_is = false; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = e_Ma * u.mia;
					B[2] = e_Mi * u.mii;
					B[3] = e_Lm * u.lm;
					B[4] = e_Ln * u.ln;
					return Sum_B();
				}
			};
		}
		namespace Weak
		{
			class SubVal : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					return e_Mi * u.mii + e_Ma * u.mia + e_Lm * u.lm;
				}
			};
		}

	} // namespace Phagocytosis

} // namespace StraightTask