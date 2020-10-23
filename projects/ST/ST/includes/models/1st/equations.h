/*This header constains the right parts of equations of ODE_system of the 1st-model //*/

#pragma once
#include <models/1st/variables.h>
#include <base/Eq_base.h>
#include <cmath>

namespace StraightTask
{
	[[nodiscard]] inline double_t Hill(double_t C, double_t ed, double_t K, int32_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	[[nodiscard]] inline double_t CHill(double_t C, double_t ed, double_t K, double_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	[[nodiscard]] inline double_t GetPositive(double_t val) {
		if (val > 0.0) return val; else return 0.0;
	}

	[[nodiscard]] inline double_t indicator(double_t val, double_t border) {
		if (val > border) return 1.0; else return 0.0;
	}


	template<typename argType> class IEqMember: public IRightPart<argType>{};

	namespace Neurons
	{
		static double_t
			p_N = 0.5;

		namespace NecroticCells
		{
			class RightPart: public IEqMember<variables const&>{
			public:
				RightPart() { comp_amount = 2; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.4 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = p_N * u.d_ini * u.hel;
					B[2] = -u.eps_s * u.nec;
					return Sum_B();
				}
				
			};
		}

		namespace AcuteChanges
		{
			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 2; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.0 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = (1.0 - p_N) * u.d_ini * u.hel;
					B[2] = -u.eps_s * u.ap_e;
					return Sum_B();
				}
			};
		}

		namespace Apoptosis
		{
			namespace Started
			{
				class RightPart : public IEqMember<variables const&> {
				public:
					RightPart() { comp_amount = 2; ret_is = true; }
					const double_t ini_data[2] = { 0.0, 0.0 };
					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						B[1] = (1.0 - p_N) * u.d_ini * u.hel;
						B[2] = -(1.0 - p_N) * u.ret.d_12 * u.ret.hel_12;
						return Sum_B();
					}
				};
			}

			namespace Ended
			{
				class RightPart : public IEqMember<variables const&> {
				public:
					RightPart() { comp_amount = 2;  ret_is = true; }
					const double_t ini_data[2] = { 0.0, 0.0 };
					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						B[1] = (1.0 - p_N) * u.ret.d_12 * u.ret.hel_12;
						B[2] = -u.eps_s * u.ap_e;
						return Sum_B();
					}
				};
			}
		}

		namespace IntactCells 
		{
			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 1; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.6 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = -u.d_ini * u.hel;
					return Sum_B();
				}
			};
		}

	} // namespace Neurons

	namespace Cytokines
	{
		static double_t p_Mach = 4.5, C_Ma = 0.9, p_Lmch = 4.5, C_Lm = 0.1, e_cy = 0.1, e_ch = 0.18;

		// psy
		static double_t p_xcy0 = 10.0, cy_max = 1.0, T_ = 1.0, t0 = 8.0;
// 																						

		class Psy : public IEqMember<variables const &>{
		public:
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				return p_xcy0 * ((cy_max - u.cy) / (T_ + (u.tj / t0)));
			}
		};

		namespace Pro_Inflam
		{
			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 5; ret_is = false;}
				const double_t ini_data[2] = { 0.0, 0.0 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = Hill(u.psy, u.mia, C_Ma, 1) * u.nec;
					B[2] = Hill(u.psy, u.mia, C_Ma, 1) * u.ap_e;
					B[3] = Hill(u.psy, u.lm, C_Lm, 1) * u.nec;
					B[4] = Hill(u.psy, u.lm, C_Lm, 1) * u.ap_e;
					B[5] = -e_cy * u.cy;
					return Sum_B();
				}
			};
		}
		namespace Chemotaxis
		{
			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 5; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.0 };
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
		static double_t o_cy_1 = 5.0,
					o_cy_2 = 5.0,
					e_adh = 0.1;

		class RightPart : public IEqMember<variables const&>{
		public:
			RightPart() { comp_amount = 3;  ret_is = false; }
			const double_t ini_data[2] = { 0.0, 0.0 };
			[[nodiscard]] double_t Expression(variables const& u)noexcept final{
				B[1] = o_cy_1 * u.cy;
				B[2] = - o_cy_2 * u.adh * u.cy;
				B[3] = - e_adh * u.adh;
				return Sum_B();
			}
		};

	} // namespace Adhesion

	namespace LeuMacrophags
	{
		static double_t
			c_Lm = 2.4, p_dLm = 3.6, T_Lm = 1.0;
			
		class RightPart : public IEqMember<variables const&>{
		public:
			RightPart() { comp_amount = 2;  ret_is = true; }
			const double_t ini_data[2] = { 0.0, 0.0 };
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Lm * u.ret.adh_24;
				B[2] = - p_dLm * (u.lm / T_Lm);
				return Sum_B();
			}
		};
	} // namespace LeuMacrophags

	namespace LeuNeutrophils
	{
		static double_t
			c_Ln = 2.8, p_dLn = 3.2, T_Ln = 1.0;

		class RightPart : public IEqMember<variables const&> {
		public:
			RightPart() { comp_amount = 2; ret_is = true; }
			const double_t ini_data[2] = { 0.0, 0.0 };
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
				B[1] = c_Ln * u.ret.adh_12;
				B[2] = -p_dLn * (u.ln / T_Ln);
				return Sum_B();
			}
		};
	} // namespace LeuNeutrophils

	namespace Microglia
	{
		static double_t
			p_1 = 0.55, c_A = 2.1, c_N = 2.1, c_pro = 1.0,
			T_M1 = 60.0, T_M2 = 18.0, c_Mi1 = 0.3, c_Mi2 = 0.3,	c_dMi = 0.2;
	
		namespace Active
		{

			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 3; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.0 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = p_1 * c_A * u.ap_e * u.mii;
					B[2] = p_1 * c_N * u.nec * u.mii;
					B[3] = -c_pro * (u.mia / T_M1);
					return Sum_B();
				}
			};
		}
		namespace Inactive
		{
			class RightPart : public IEqMember<variables const&> {
			public:
				RightPart() { comp_amount = 5; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 1.0 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = -c_A * u.ap_e * u.mii;
					B[2] = -c_N * u.nec * u.mii;
					B[3] = c_pro * (u.mia / T_M1);
					B[4] = -c_dMi * u.mii;
					B[5] = u.mii * (c_Mi1 - c_Mi2 * u.mii) * indicator(u.tj, T_M2);

					return Sum_B();
				}
			};
		}
	} // namespace Microglia

	template<typename argType>	class ISubVal : public IRightPart<argType> {};

	namespace ToxDamage
	{

		static double_t
			p_ncy = 0.1, p_Ln = 0.4, C_DLn = 0.6, P_nn = 0.05, C_D = 1.0,
			D_0 = 0.02;

		namespace Full
		{
			class SubVal : public ISubVal<variables const&> {
			public:
				SubVal() { comp_amount = 4; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = p_ncy * u.cy;
					B[2] = (u.nec / C_D) * Hill(p_Ln, u.ln, C_DLn, 1);
					B[3] = (u.ap_e / C_D) * Hill(p_Ln, u.ln, C_DLn, 1);
					B[4] = P_nn * (u.nec / C_D);
					return Sum_B();
				}
			};
		}

		namespace Initial
		{
			class SubVal : public ISubVal<variables const&>{
			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}
			};
		}

	} // namespace ToxDamage

	namespace Phagocytosis
	{
		static double_t 
			e_Ma = 0.1, 
			e_Mi = 0.025,
			e_Ln = 0.025,
			e_Lm = 0.1;

		namespace Strong
		{
			class SubVal : public ISubVal<variables const&>{
			public:
				SubVal() { comp_amount = 4; }
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = e_Ma * u.mia;
					B[2] = e_Mi * u.mii;
					B[3] = e_Lm * u.lm;
					B[4] = e_Ln * u.ln;
					return Sum_B();
				}
			};
		}
	} // namespace Phagocytosis

} // namespace StraightTask