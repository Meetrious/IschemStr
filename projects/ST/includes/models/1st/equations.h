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


	template<typename argType>
	class IEqMember: public IRightPart<argType>{};

	namespace Neurons
	{
		class constants {
		public:
			constants() :
				p_N("p_N", 0.5)
			{}
			IConstant p_N;
		};

		namespace NecroticCells
		{
			class RightPart: public IEqMember<variables const&>{
				double_t p_N;
			public:
				RightPart() { comp_amount = 2; ret_is = false; }

				const double_t ini_data[2] = { 0.0, 0.4 };

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = p_N * u.d_ini * u.hel;
					B[2] = -u.eps_s * u.nec;
					return Sum_B();
				}
				void SynchronizeCoefs(constants& c)noexcept {
					c.p_N.sync(p_N);
				}
				
			};
		}

		namespace AcuteChanges
		{
			class RightPart : public IEqMember<variables const&> {
				double_t p_N;

			public:
				RightPart() { comp_amount = 2; ret_is = false; }

				const double_t ini_data[2] = { 0.0, 0.0 };

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = (1.0 - p_N) * u.d_ini * u.hel;
					B[2] = -u.eps_s * u.ap_e;
					return Sum_B();
				}

				void SynchronizeCoefs(constants& c)noexcept {
					c.p_N.sync(p_N);
				}
			};
		}

		namespace Apoptosis
		{
			namespace Started
			{
				class RightPart : public IEqMember<variables const&> {
					double_t p_N;

				public:
					RightPart() { comp_amount = 2; ret_is = true; }

					const double_t ini_data[2] = { 0.0, 0.0 };

					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						B[1] = (1.0 - p_N) * u.d_ini * u.hel;
						B[2] = -(1.0 - p_N) * u.ret.d_12 * u.ret.hel_12;
						return Sum_B();
					}
					void SynchronizeCoefs(constants& c)noexcept {
						c.p_N.sync(p_N);
					}
				};
			}

			namespace Ended
			{
				class RightPart : public IEqMember<variables const&> {
					double_t p_N;

				public:
					RightPart() { comp_amount = 2;  ret_is = true; }

					const double_t ini_data[2] = { 0.0, 0.0 };

					[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
						B[1] = (1.0 - p_N) * u.ret.d_12 * u.ret.hel_12;
						B[2] = -u.eps_s * u.ap_e;
						return Sum_B();
					}
					void SynchronizeCoefs(constants& c)noexcept {
						c.p_N.sync(p_N);
					}
				};
			}
		}

		namespace IntactCells 
		{
			class RightPart : public IEqMember<variables const&> {
				double_t p_N;
			public:
				RightPart() { comp_amount = 1; ret_is = false; }
				const double_t ini_data[2] = { 0.0, 0.6 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = -u.d_ini * u.hel;
					return Sum_B();
				}
				void SynchronizeCoefs(constants& c)noexcept {
					c.p_N.sync(p_N);
				}
			};
		}

	} // namespace Neurons

	namespace Cytokines
	{
		class constants {
		public:
			constants() :
				p_Mach("p_{Ma_{ch}}", 4.5),
				C_Ma("C_{Ma}", 0.9),
				p_Lmch("p_{Lm_{cy}}", 4.5),
				C_Lm("C_{Lm}", 0.1),
				e_cy("e_{cy}", 0.1),
				e_ch("e_{ch}", 0.18),
				p_xcy0("p_{xcy_0}", 10.0),
				cy_max("cy_{max}", 1.0),
				T_("T_trh", 1.0),
				t0("tau_0", 8.0)
			{}

			IConstant p_Mach, C_Ma, p_Lmch, C_Lm, e_cy, e_ch,
				p_xcy0, cy_max, T_, t0;
		};

		
// 																						

		class Psy : public IEqMember<variables const &>{
			double_t p_xcy0, cy_max, T_, t0;

		public:
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				return p_xcy0 * ((cy_max - u.cy) / (T_ + (u.tj / t0)));
			}
			void SynchronizeCoefs(constants& c)noexcept {
				c.p_xcy0.sync(p_xcy0);
				c.cy_max.sync(cy_max);
				c.T_.sync(T_);
				c.t0.sync(t0);
			}
		};

		namespace Pro_Inflam
		{
			class RightPart : public IEqMember<variables const&> {
				double_t C_Ma, C_Lm, e_cy;\

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

				void SynchronizeCoefs(constants& c)noexcept {
					c.C_Ma.sync(C_Ma);
					c.C_Lm.sync(C_Lm);
					c.e_cy.sync(e_cy);
				}
			};
		}
		namespace Chemotaxis
		{
			class RightPart : public IEqMember<variables const&> {
				double_t p_Mach, p_Lmch, C_Ma, C_Lm, e_ch;

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

				void SynchronizeCoefs(constants& c)noexcept {
					c.p_Mach.sync(p_Mach);
					c.p_Lmch.sync(p_Lmch);
					c.C_Ma.sync(C_Ma);
					c.C_Lm.sync(C_Lm);
					c.e_ch.sync(e_ch);
				}
			};
		}
	} // namespace Cytokines

	namespace Adhesion 
	{
		class constants {
		public:
			constants() :
				o_cy_1("o_{cy_1}", 5.0),
				o_cy_2("o_{cy_2}", 5.0),
				e_adh("e_{adh}", 0.1)
			{}

			IConstant o_cy_1, o_cy_2, e_adh;
		};

		class RightPart : public IEqMember<variables const&>{
			double_t o_cy_1, o_cy_2, e_adh;
		public:
			RightPart() { comp_amount = 3;  ret_is = false; }

			const double_t ini_data[2] = { 0.0, 0.0 };

			[[nodiscard]] double_t Expression(variables const& u)noexcept final{
				B[1] = o_cy_1 * u.cy;
				B[2] = - o_cy_2 * u.adh * u.cy;
				B[3] = - e_adh * u.adh;
				return Sum_B();
			}
			void SynchronizeCoefs(constants& c)noexcept {
				c.o_cy_1.sync(o_cy_1);
				c.o_cy_2.sync(o_cy_2);
				c.e_adh.sync(e_adh);
			}
		};

	} // namespace Adhesion

	namespace LeuMacrophags
	{
		class constants {
		public:
			constants() :
				c_Lm("c_{Lm}", 2.4),
				p_dLm("p_{d_{Lm}}", 3.6),
				T_Lm("T_{Lm}", 1.0)
			{}
			IConstant c_Lm, p_dLm, T_Lm; //*/
		};
			
		class RightPart : public IEqMember<variables const&>{
			double_t c_Lm, p_dLm, T_Lm;

		public:
			RightPart() { comp_amount = 2;  ret_is = true; }

			const double_t ini_data[2] = { 0.0, 0.0 };

			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Lm * u.ret.adh_24;
				B[2] = - p_dLm * (u.lm / T_Lm);
				return Sum_B();
			}

			void SynchronizeCoefs(constants& c) noexcept {
				c.c_Lm.sync(c_Lm);
				c.p_dLm.sync(p_dLm);
				c.T_Lm.sync(T_Lm);
			}
		};
	} // namespace LeuMacrophags

	namespace LeuNeutrophils
	{
		class constants {
		public:
			constants() :
				c_Ln("c_{Ln}", 2.8),
				p_dLn("p_{d_{Ln}}}", 3.2),
				T_Ln("T_{Ln}}", 1.0)
			{}
			IConstant c_Ln, p_dLn, T_Ln;
		};

		class RightPart : public IEqMember<variables const&> {
			double_t c_Ln, p_dLn, T_Ln;

		public:
			RightPart() { comp_amount = 2; ret_is = true; }

			const double_t ini_data[2] = { 0.0, 0.0 };

			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
				B[1] = c_Ln * u.ret.adh_12;
				B[2] = -p_dLn * (u.ln / T_Ln);
				return Sum_B();
			}

			void SynchronizeCoefs(constants& c) noexcept {
				c.c_Ln.sync(c_Ln);
				c.p_dLn.sync(p_dLn);
				c.T_Ln.sync(T_Ln);
			}
		};
	} // namespace LeuNeutrophils

	namespace Microglia
	{
		class constants {
		public:
			constants() :
				p_1("p_1", 0.55),
				c_m_A("c_{m_A}", 2.1),
				c_m_N("c_{m_N}", 2.1),
				c_pro("c_{pro}", 1.0),
				T_M1("T_{M1}", 60.0),
				T_M2("T_{M2}", 18.0),
				c_dMi("c_{d_{Mi}}", 0.2),
				c_Mi1("c_{Mi_1}", 0.3),
				c_Mi2("c_{Mi_2}", 0.3)
			{}

			IConstant p_1, c_m_A, c_m_N, c_pro,
				T_M1, T_M2, c_Mi1, c_Mi2, c_dMi;
		};
	
		namespace Active
		{

			class RightPart : public IEqMember<variables const&> {
				double_t p_1, c_m_A, c_m_N, c_pro, T_M1;

			public:
				RightPart() { comp_amount = 3; ret_is = false; }

				const double_t ini_data[2] = { 0.0, 0.0 };

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = p_1 * c_m_A * u.ap_e * u.mii;
					B[2] = p_1 * c_m_N * u.nec * u.mii;
					B[3] = -c_pro * (u.mia / T_M1);
					return Sum_B();
				}
				void SynchronizeCoefs(constants& c) noexcept {
					c.p_1.sync(p_1);
					c.c_m_A.sync(c_m_A);
					c.c_m_N.sync(c_m_N);
					c.c_pro.sync(c_pro);
					c.T_M1.sync(T_M1);
				}
			};
		}
		namespace Inactive
		{
			class RightPart : public IEqMember<variables const&> {
				double_t c_m_A, c_m_N, c_dMi, c_Mi1, c_Mi2, c_pro, T_M1, T_M2;

			public:
				RightPart() { comp_amount = 5; ret_is = false; }

				const double_t ini_data[2] = { 0.0, 1.0 };

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = -c_m_A * u.ap_e * u.mii;
					B[2] = -c_m_N * u.nec * u.mii;
					B[3] = c_pro * (u.mia / T_M1);
					B[4] = -c_dMi * u.mii;
					B[5] = u.mii * (c_Mi1 - c_Mi2 * u.mii) * indicator(u.tj, T_M2);
					return Sum_B();
				}

				void SynchronizeCoefs(constants& c) noexcept {
					c.c_m_A.sync(c_m_A);
					c.c_m_N.sync(c_m_N);
					c.c_dMi.sync(c_dMi);
					c.c_Mi1.sync(c_Mi1);
					c.c_Mi2.sync(c_Mi2);
					c.c_pro.sync(c_pro);
					c.T_M1.sync(T_M1);
					c.T_M2.sync(T_M2);
				}
			};
		}
	} // namespace Microglia

	template<typename argType>	class ISubVal : public IRightPart<argType> {};

	namespace ToxDamage
	{
		class constants {
		public:
			constants() :
				p_ncy("p_{n_{cy}}", 0.1),
				p_Ln("p_{L_n}", 0.4),
				C_DLn("C_{D_{Ln}}", 0.6),
				P_nn("P_{nn}", 0.05),

				C_D("C_D", 1.0),
				D_0("D_0", 0.02)
			{}

			IConstant p_ncy, p_Ln, C_DLn, P_nn,
				C_D, D_0;

		};

		namespace Full
		{
			class SubVal : public ISubVal<variables const&> {
				double_t p_ncy, p_Ln, C_DLn, P_nn, C_D, D_0;

			public:
				SubVal() { comp_amount = 4; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {

					B[1] = p_ncy * u.cy;
					B[2] = (u.nec / C_D) * Hill(p_Ln, u.ln, C_DLn, 1);
					B[3] = (u.ap_e / C_D) * Hill(p_Ln, u.ln, C_DLn, 1);
					B[4] = P_nn * (u.nec / C_D);
					return Sum_B();
				}

				void SynchronizeCoefs(constants& c) noexcept {
					c.p_ncy.sync(p_ncy);
					c.p_Ln.sync(p_Ln);
					c.C_DLn.sync(C_DLn);
					c.P_nn.sync(P_nn);

					c.C_D.sync(C_D);
					c.D_0.sync(D_0);
				}
			};
		}

		namespace Initial
		{
			class SubVal : public ISubVal<variables const&>{
				double_t D_0;

			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {

					return GetPositive(u.d_F - D_0);
				}

				void SynchronizeCoefs(constants& c) noexcept {
					c.D_0.sync(D_0);
				}
			};
		}

	} // namespace ToxDamage

	namespace Phagocytosis
	{
		class constants {
		public:
			constants() :
				e_Ma("e_{Ma}", 0.1),
				e_Mi("e_{Mi}", 0.025),
				e_Ln("e_{Ln}", 0.025),
				e_Lm("e_{Lm}", 0.1)
			{}

			IConstant e_Ma, e_Mi, e_Lm, e_Ln;
		};
		
		namespace Strong
		{
			class SubVal : public ISubVal<variables const&>{
				double_t e_Ma, e_Mi, e_Ln, e_Lm;

			public:
				SubVal() { comp_amount = 4; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = e_Ma * u.mia;
					B[2] = e_Mi * u.mii;
					B[3] = e_Lm * u.lm;
					B[4] = e_Ln * u.ln;
					return Sum_B();
				}

				void SynchronizeCoefs(constants& c) noexcept {
					c.e_Ma.sync(e_Ma);
					c.e_Mi.sync(e_Mi);
					c.e_Lm.sync(e_Lm);
					c.e_Ln.sync(e_Ln);
				}
			};
		}
	} // namespace Phagocytosis

} // namespace StraightTask