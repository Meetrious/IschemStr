#pragma once
#include <ISCHEM_variables.h>

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
		if (val > 0) return val; else return 0;
	}

	[[nodiscard]] inline double_t indicator(double_t val, double_t border){
		if (val > border) return 1.0; else return 0;
	}

	template<typename argType>
	class IRightPart
	{
	protected:
		double_t Sum_B() {
			B[0] = 0.0;
			for (uint16_t i = 1; i < B.size(); i++) { B[0] += B[i]; }
			return B[0];
		}
	
	public:
		vector<double_t> B;
		virtual double_t Expression(argType u) = 0;
		virtual double_t Solution(double_t t) { return 0.0; }
	};

	template<typename argType>	class IEquation: public IRightPart<argType>{};


	namespace Neurons
	{
		static double_t
			p_N = 0.9, Rep = 1.0e-4, k_N = 0.56366, k_A = 3.0895, p_R = 0.267460,
			q_N = 1.5, q_A_n = 0.406, c_H = 0.006, c_A = 1.2e-2, q_epN = 0.15, q_epA = 0.001,
			q_A = 1.0, q_A_h = 0.006, p_r = 0.1, q_H = 8.13;


		namespace NecroticCells
		{
			class OrigModel: public IEquation<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return p_N * u.d_ini * u.hel - u.eps_s * u.nec;
				}
			};

			class SecondModel : public IEquation<variables const&>{
			public:
				SecondModel() { B.insert(B.end(), 3, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const & u)noexcept final{
					B[1] = k_N * u.dp_N * u.hel;
					B[2] = -u.eps_s * u.nec;
					return Sum_B();
				}
			};

			class ThirdModel : public IEquation<variables const&>{
			public:
				ThirdModel() { B.insert(B.end(), 4, 0.0); }
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
			class OrigModel : public IEquation<variables const&> {
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return (1-p_N) * u.d_ini * u.hel - u.eps_s * u.ap_e;
				}
			};

			class SecondModel : public IEquation<variables const&>  {
			public:
				SecondModel() { B.insert(B.end(), 4, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = k_A * u.dp_A * u.hel;
					B[2] = -u.acu_c * p_R * k_A * u.dp_A;
					B[3] = -u.acu_c * u.eps_s;
					return Sum_B();
				}
			};
			
			class ThirdModel : public IEquation<variables const&> {
			public:
				ThirdModel() { B.insert(B.end(), 5, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = (q_A - p_r) * u.dp_A * u.hel * Hill(1.0, u.acu_c, q_A_h, 1);
					B[2] = q_H * u.dp_N * u.hel * u.acu_c;
					B[3] = -q_A_n * u.dp_N * Hill(1.0, u.acu_c, c_A, 5);
					B[4] = -q_epA * u.eps_s * u.acu_c;;
					return Sum_B();
				}
			};
		}

		namespace Apoptosis
		{
			namespace Started
			{
				class OrigModel : public IEquation<variables const&> {
				public:
					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						return (1.0 - p_N) * u.d_ini * u.hel - (1.0 - p_N) * u.ret.d_12 * u.ret.hel_12;
					}
				};
			}

			namespace Ended
			{
				class OrigModel : public IEquation<variables const&> {
				public:
					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						return (1.0 - p_N) * u.ret.d_12 * u.ret.hel_12 - u.eps_s * u.ap_e;
					}
				};
			}
		}

		namespace IntactCells 
		{
			class OrigModel : public IEquation<variables const&> {
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return -u.d_ini * u.hel;
				}
			};

			class SecondModel : public IEquation<variables const&> {
			public:
				SecondModel() { B.insert(B.end(), 4, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = -u.hel * k_N * u.dp_N;
					B[2] = -u.hel * k_A * u.dp_A;
					B[3] = p_R * k_A * u.dp_A * u.acu_c;
					return Sum_B();
				}
			};
			class ThirdModel : public IEquation<variables const&> {
			public:
				ThirdModel() { B.insert(B.end(), 4, 0.0); }
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
		static double_t p_Mach = 4.5, C_Ma = 9.9713, p_Lmch = 4.5, C_Lm = 8.1823, e_ch = 0.18;
		// ordinal
		//static double_t p_Macy = 2.9481, p_Lmcy = 5.6667, p_Lncy = 8.2513, C_Ln = 3.3918, e_cy = 5.577, l1 = 0.75752, l2 = 5.1694, l3 = 1.3089, _A = 0.14514;
		static double_t p_Macy = 3.7386, p_Lmcy = 4.9559, p_Lncy = 6.2865, C_Ln = 5.0366, 
			e_cy = 1.166, l1 = 0.75752, l2 = 5.1694, l3 = 1.3089, _A = 0.14514;
		// psy
		static double_t p_xcy0 = 0.5, cy_max = 1.0, T_ = 1.0, t_0 = 8.0;
// 																						

		class Psy : public IEquation<variables const &>{
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				return p_xcy0 * ((cy_max - u.cy) / (T_ + (u.tj / t_0)));
			}
		};

		namespace Pro_Inflam
		{
			class OrigModel : public IEquation<variables const&> {
			public:
				OrigModel() { B.insert(B.end(), 6, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = Hill(u.psy, u.mia, C_Ma, 1) * u.nec;
					B[2] = Hill(u.psy, u.mia, C_Ma, 1) * u.ap_e;
					B[3] = Hill(u.lm, u.lm, C_Lm, 1) * u.nec;
					B[4] = Hill(u.lm, u.lm, C_Lm, 1) * u.ap_e;
					B[5] = -e_cy * u.cy;
					return Sum_B();
				}
			};

			class SecondModel : public IEquation<variables const&> {
			public:
				SecondModel() { B.insert(B.end(), 5, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = Hill(u.psy, u.mia, C_Ma, 1);
					B[2] = Hill(u.psy, u.lm, C_Lm, 1);
					B[3] = Hill(2.0, u.ln, C_Ln, 1);
					B[4] = -e_cy * u.cy;
					return Sum_B();
				}
			};
			class ThirdModel : public IEquation<variables const&> {
			public:
				ThirdModel() { B.insert(B.end(), 5, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					
					B[1] = CHill(p_Macy, u.mia, C_Ma, l1) * (u.nec + GetPositive(u.acu_c - _A));
					B[2] = p_Lmcy * u.lm * (u.nec + u.acu_c);
					B[3] = p_Lncy * u.ln * (u.nec + u.acu_c);
					B[4] = - e_cy * u.cy;
					
					return Sum_B();
				}
			};
		}
		namespace Chemotaxis
		{
			class OrigModel : public IEquation<variables const&> {
			public:
				OrigModel() { B.insert(B.end(), 6, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = Hill(p_Mach, u.mia, C_Ma, 1) * u.nec;
					B[2] = Hill(p_Mach, u.mia, C_Ma, 1) * u.ap_e;
					B[3] = Hill(p_Lmch, u.lm, C_Lm, 1) * u.nec;
					B[4] = Hill(p_Lmch, u.lm, C_Lm, 1) * u.ap_e;
					B[5] = -e_ch * u.ch;
					return Sum_B();
					//return (Hill(p_Mach, u.mia, C_Ma, 1) + Hill(p_Lmch, u.lm, C_Lm, 1)) * (u.nec + u.ap_e) - e_ch * u.ch;
					
				}
			};
		}
	} // namespace Cytokines

	namespace Adhesion 
	{
		static double_t o_cy_1 = 2.62e-2,
					o_cy_2 = 3.14e-2,
					e_adh = 1.31e-2;

		class OrigModel : public IEquation<variables const&>{
		public:
			OrigModel() { B.insert(B.end(), 4, 0.0); }
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
			c_Lm = 9.63e-2, p_dLm = 3.6, T_Lm = 1.0,
			c_dLm = 0.08, K_Lm = 43.9, d_Lm = 0.139;

		class OrigModel : public IEquation<variables const&>{
		public:
			OrigModel() { B.insert(B.end(), 3, 0.0); }
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Lm * u.ret.adh_24;
				B[2] = -p_dLm * (u.lm / T_Lm);
				return Sum_B();
			}
		};
		class SecondModel : public IEquation<variables const&>{
		public:
			SecondModel() { B.insert(B.end(), 4, 0.0); }
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Lm * u.ret.adh_4;
				B[2] = Hill(c_dLm, u.cy, K_Lm, 5);
				B[3] = -p_dLm * (u.lm / T_Lm);
				return Sum_B();
			}
		};
		class ThirdModel : public IEquation<variables const&> {
		public:
			ThirdModel() { B.insert(B.end(), 4, 0.0); }
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
			c_Ln = 1.47, p_dLn = 3.2, T_Ln = 1.0,
			c_dLn = 1.0, K_Ln = 1.0, d_Ln = 84,
			c_dLn1 = 16.2, K_Ln1 = 3.96,
			c_dLn2 = 4.43, K_Ln2 = 9.899;

		class OrigModel : public IEquation<variables const&> {
		public:
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
				return c_Ln * u.ret.adh_12 - p_dLn * u.ln/T_Ln;
			}
		};
		

		class ThirdModel : public IEquation<variables const&>{
		public:
			ThirdModel() { B.insert(B.end(), 5, 0.0); }
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
			p_1 = 1.0, c_A = 0.06, c_N = 2.1, c_pro = 2.5,
			T_M1 = 60.0, T_M2 = 60.0, c_Mi1 = 0.3, c_Mi2 = 0.3,
			c_dMi = 0.2, c_Mi = 1.0, K_Mi = 0.1;
	
		namespace Active
		{

			class OrigModel : public IEquation<variables const&> {
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return p_1*(c_A * u.ap_e + c_N * u.nec) - c_pro * u.mia/T_M1;
				}
			};

			class SecondModel : public IEquation<variables const&>	{
			public:
				SecondModel() { B.insert(B.end(), 5, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = p_1 * u.mii * c_A * u.acu_c;
					B[2] = p_1 * u.mii * c_N * u.nec;
					B[3] = Hill(c_Mi, u.cy, K_Mi, 1) * u.mia;
					B[4] = - c_pro * (u.mia / T_M1);
					return Sum_B();
					//return p_1 * (c_A * u.acu_c + c_N * u.nec) * u.mii + Hill(c_Mi, u.cy, K_Mi, 1) * u.mia - c_pro * (u.mia / T_M1);
				}
			};
			class ThirdModel : public IEquation<variables const&>{
			public:
				ThirdModel() { B.insert(B.end(), 5, 0.0); }
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
			class OrigModel : public IEquation<variables const&> {
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return - (c_A * u.ap_e + c_N * u.nec) + c_pro * u.mia / T_M1 
						- c_dMi * u.mii + (c_Mi1 * u.mii - c_Mi2 * u.mii * u.mii) * indicator(u.tj, T_M2);
				}
			};

			class SecondModel : public IEquation<variables const&>
			{
				SecondModel() { B.insert(B.end(), 5, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u) noexcept final {
					B[1] = c_pro * (u.mia / T_M1);
					B[2] = -u.mii * c_A * u.acu_c;
					B[3] = -u.mii * c_N * u.nec;
					B[4] = -Hill(c_Mi, u.cy, K_Mi, 1) * u.mii;
					return Sum_B();
					//return c_pro * (u.mia / T_M1) - (c_A * u.acu_c + c_N * u.nec) * u.mii - Hill(c_Mi, u.cy, K_Mi, 1) * u.mii;
				}
			};
			class ThirdModel : public IEquation<variables const&>
			{
			public:
				ThirdModel() { B.insert(B.end(), 5, 0.0); }
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
			p_ncy = 17.83860, p_Ln = 27.078500, C_DLn = 40.5670, P_nn = 4.820900,
			p_Lm = 34.70, C_DLm = 8.6732, C_D = 1.0, C_Dcy = 9.7557,

			p_q1 = 0.1, p_q2 = 0.4, p_q3 = 0.4, p_q4 = 0.1,
			c_q1 = 13.206, c_q2 = 1.32, c_q3 = 0.1, c_q4 = 0.1,

			D_0 = 0, p_D = 0.28817;

		[[nodiscard]] inline double_t State_dpN(const variables& u)noexcept
		{
			if (u.d_F <= D_0) return p_D * u.d_F;
			else return u.d_F - (D_0 - p_D * u.d_F);
		}
		[[nodiscard]] inline double_t State_dpA(const variables& u)noexcept
		{
			if (u.d_F <= D_0) return u.d_F - p_D * u.d_F;
			else return D_0 - p_D * u.d_F;
		}
		[[nodiscard]] inline double_t St_dpA(double_t d_F, double_t D_0)noexcept
		{
			if (d_F <= D_0) return d_F;
			else return 0;
		}

		namespace Full
		{
			class SecondModel : public ISubVal<variables const&> {
			public:
				SecondModel() { B.insert(B.end(), 7, 0.0); }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = Hill(p_ncy, u.cy, C_Dcy, 1);
					B[2] = u.nec * Hill(p_Lm, u.lm, C_DLm, 1) / C_D;
					B[3] = u.nec * Hill(p_Ln, u.ln, C_DLn, 1) / C_D;
					B[4] = u.acu_c * Hill(p_Lm, u.lm, C_DLm, 1) / C_D;
					B[5] = u.acu_c * Hill(p_Ln, u.ln, C_DLn, 1) / C_D;
					B[6] = P_nn * u.nec / C_D;
					return Sum_B();
					//return Hill(p_ncy, u.cy, C_Dcy, 1) + (Hill(p_Lm, u.lm, C_DLm, 1) + Hill(p_Ln, u.ln, C_DLn, 1)) * (u.nec + u.acu_c) / C_D 	+ P_nn * u.nec / C_D;
				}
			};
			class ThirdModel : public ISubVal<variables const&>{
			public:
				ThirdModel() { B.insert(B.end(), 5, 0.0); }
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
			class OrigModel : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}
			};
		}

		namespace Nec_partial
		{
			class SecondModel : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(State_dpN(u));
				}
			};
			class ThirdModel : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}
			};

		}

		namespace Apop_partial
		{
			class SecondModel : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(State_dpA(u));
				}
			};
			class ThirdModel : public ISubVal<variables const&>{
			public:
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
			
			class OrigModel : public ISubVal<variables const&>{
			public:
				OrigModel() { B.insert(B.end(), 5, 0.0); }
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
			class OrigModel : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					return e_Mi * u.mii + e_Ma * u.mia + e_Lm * u.lm;
				}
			};
		}

	} // namespace Phagocytosis

	namespace Test
	{
		namespace OneDim 
		{
			static double_t k = 2;

			size_t factorial(size_t n) {
				if (n > 1) return n * factorial(n - 1);	else return 1;
			}

			class Exp_ret :public IEquation<variables const&>{
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					return u.ret.x_1;
				}
			public:
				double_t Solution(double_t t) noexcept final {
					double_t value = 1;
					for (size_t n = 1; n < floor(t); n++) {
						value += std::pow((t - (n - 1)), n) / factorial(n);
					}
					return value;
				}
				
			};
			class Exp :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return u.x;
				}
			public:
				double_t Solution(double_t t) noexcept final { return exp(t); }
			};	

			class Sin :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return -k * u.ret.x_pi2 ;
				}
			public:
				double_t Solution(double_t t) noexcept final { return k * sin(t); }
			};


		}
		namespace ThreeDim
		{

			class Ox_ret :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return (2.0 / pi) * (u.x + u.ret.y_pi2) - u.ret.x_pi2 - (0.5 * pi) * (u.y / u.z);
				}
			public:
				double_t Solution(double_t t)noexcept final { return t * cos(t); }
			};
			class Oy_ret :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return (2.0 / pi) * (u.y - u.ret.x_pi2) - u.ret.y_pi2 + (0.5 * pi) * (u.x / u.z);
				}
			public:
				double_t Solution(double_t t)noexcept final { return t * sin(t); }
			};
			class Oz_ret :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return std::pow(u.ret.x_pi2 * u.ret.x_pi2 + u.ret.y_pi2 * u.ret.y_pi2, 0.5) / u.ret.z_pi2 ;
				}
			public:
				double_t Solution(double_t t)noexcept final { return t; }
			};

			class Ox :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return u.x/u.z - u.y;
				}
			public:
				double_t Solution(double_t t)noexcept final{ return t * cos(t);	}
			};
			class Oy :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return u.y/u.z + u.x;
				}
			public:
				double_t Solution(double_t t)noexcept final { return t * sin(t); }
			};
			class Oz :public IEquation<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return (u.x * u.x + u.y * u.y) / (u.z * u.z);
				}
			public:
				double_t Solution(double_t t)noexcept final { return t; }
			};
		}
	}

} // namespace StraightTask