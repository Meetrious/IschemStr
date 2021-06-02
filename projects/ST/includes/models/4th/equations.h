/* This header constains the right parts for equations from 4th-model ODE_system //*/

#pragma once
#include <models/4th/variables.h>
#include <base/Eq_base.h>

#include <cmath>


namespace StraightTask
{

	[[nodiscard]] inline double_t Hill (double_t C, double_t ed, double_t K, int32_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	[[nodiscard]] inline double_t CHill (double_t C, double_t ed, double_t K, double_t n) {
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}


	[[nodiscard]] double_t GetPositive (double_t val) {
		if (val > 0.0) return val; else return 0;
	}

	[[nodiscard]] double_t indicator(double_t val, double_t border){
		if (val > border) return 1.0; else return 0;
	}


	template<typename argType>
	class IEqMember: public IRightPart<argType>{};


	namespace Neurons
	{
		class constants {
		public:
			constants() :
				q_N("q_N", 2.0473),
				q_A_n("q_{A_n}", 2.1281),
				c_H("c_H", 2.5482),
				c_A("c_A", 1.4925),
				q_epN("q_{ep_N}", 0.072613),
				q_epA("q_{ep_A}", 0.37173),
				q_A("q_A", 1.50756),
				c_A_h("c_{A_h}", 1.34878),
				p_r("p_r", 0.54062),
				q_H("q_H", 14.081),
				l1("l_1", 2.0),
				l2("l_2", 2.5),
				l3("l_3", 1.0)
			{}
			IConstant q_N, q_A_n, c_H, c_A,
				q_epN, q_epA,
				q_A, c_A_h, p_r, q_H, l1, l2, l3;
		};


		namespace NecroticCells
		{
			class Equation : public IEqMember<variables const&>{
				double_t q_N, q_A_n, c_H, c_A, q_epN, l1, l2;
			public:
				Equation()				
				{ comp_amount = 3; ret_is = false; }

				const double_t ini_data[2] = { 0.5 , 0.0 };
				[[nodiscard]] inline double_t Expression(variables const & u)noexcept final{
					B[1] = q_N * u.dp_N * CHill(1.0, u.hel, c_H, l1) * u.nec;
					B[2] = q_A_n * u.dp_N * CHill(1.0, u.acu_c, c_A, l2);
					B[3] = -q_epN * u.eps_s * u.nec;
					return Sum_B();
				}
				
				void SynchronizeCoefs(constants & c)noexcept{
					c.q_N.sync(q_N);
					c.c_H.sync(c_H);
					c.q_A_n.sync(q_A_n);
					c.c_A.sync(c_A);
					c.l1.sync(l1);
					c.l2.sync(l2);
					c.q_epN.sync(q_epN);
				}
			};
		}

		namespace AcuteChanges
		{
			class Equation : public IEqMember<variables const&> {
				double_t p_r, q_A, c_A_h,
					q_H, q_A_n, c_A, l2, l3,
					q_epA; //*/
					 
			public:
				Equation()
				{ comp_amount = 4; ret_is = false; }

				const double_t ini_data[2] = { 0.5, 0.2890 }; //experiment

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = (1.0 - p_r) * q_A * u.dp_A * u.hel * CHill(1.0 , u.acu_c , c_A_h , l3);
					B[2] = q_H * u.dp_N * u.hel * u.acu_c;
					B[3] = -q_A_n * u.dp_N * CHill(1.0, u.acu_c, c_A, l2);
					B[4] = -q_epA * u.eps_s * u.acu_c;
					return Sum_B();
				}

				void SynchronizeCoefs(constants & c)noexcept {
					c.p_r.sync(p_r);
					c.q_A.sync(q_A);
					c.c_A_h.sync(c_A_h);
					c.q_H.sync(q_H);
					c.q_A_n.sync(q_A_n);
					c.c_A.sync(c_A);
					c.l2.sync(l2);
					c.l3.sync(l3);
					c.q_epA.sync(q_epA);
				}
			};
		}

		namespace IntactCells 
		{
			class Equation : public IEqMember<variables const&> {

				double_t p_r, q_A, c_A_h,
					q_N, c_H, q_H, l1, l3;

			public:
				Equation()
				{ comp_amount = 3; ret_is = false; }

				const double_t ini_data[2] = { 0.5, 0.7110 }; // experiment
				
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = (-1.0 + p_r)  * q_A * u.dp_A * u.hel * CHill(1.0, u.acu_c, c_A_h, l3);
					B[2] = -q_N * u.dp_N * u.nec * CHill(1.0, u.hel, c_H, l1);
					B[3] = -q_H * u.dp_N * u.hel * u.acu_c;
					return Sum_B();
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.p_r.sync(p_r);
					c.q_A.sync(q_A); 
					c.c_A_h.sync(c_A_h);
					c.q_N.sync(q_N);
					c.c_H.sync(c_H);
					c.q_H.sync(q_H);
					c.l1.sync(l1);
					c.l3.sync(l3);
				}
			};
		}

	} // namespace Neurons

	namespace Cytokines
	{

		class constants {
		public:
			constants() :
				p_Macy("p_{Ma_{cy}}", 0.47279),
				C_Ma("C_{Ma}", 4.1228),
				//l1("l_1", 1.178581),
				l1("l_1", 1.0),
				_A("A_{trh}", 0.18645),
				p_Lmcy("p_{Lm_{cy}}", 6.096),
				C_Lm("C_{Lm}", 2.9515),
				l2("l_2", 3.5),
				p_Lncy("p_{Ln_{cy}}", 7.9203),
				C_Ln("C_{Ln}", 0.85552),
				//l3("l_3", 0.882709),
				l3("l_3", 1.0),
				e_cy("e_{cy}", 3.2248)// */
			{}

			IConstant p_Macy, C_Ma, l1, _A,
				p_Lmcy, C_Lm, l2, 
				p_Lncy, C_Ln, l3,
				e_cy; // */
		};

		namespace Pro_Inflam
		{
			class Equation : public IEqMember<variables const&> {
				double_t p_Macy, C_Ma, l1, _A,
					p_Lmcy, C_Lm, l2,
					p_Lncy, C_Ln, l3,
					e_cy; // */
			public:
				Equation()
				{ comp_amount = 4; ret_is = false;}
				const double_t ini_data[2] = { 0.5, 0.0906 };
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					
					B[1] = CHill(p_Macy, u.mia, C_Ma, l1) * (u.nec + GetPositive(u.acu_c - _A));
					B[2] = CHill(p_Lmcy, u.lm, C_Lm, l2) * (u.nec + u.acu_c);
					B[3] = CHill(p_Lncy, u.ln, C_Ln, l3) * (u.nec + u.acu_c);
					B[4] = - e_cy * u.cy;
					

					return Sum_B();
				}
				void SynchronizeCoefs(constants & c)noexcept {
					c.p_Macy.sync(p_Macy);
					c.C_Ma.sync(C_Ma);
					c.l1.sync(l1);
					c._A.sync(_A);
					c.p_Lmcy.sync(p_Lmcy);
					c.C_Lm.sync(C_Lm);
					c.l2.sync(l2);
					c.p_Lncy.sync(p_Lncy);
					c.C_Ln.sync(C_Ln);
					c.l3.sync(l3);
					c.e_cy.sync(e_cy);
				}
			};
		}
	} // namespace Cytokines

	namespace Adhesion 
	{
		class constants {
		public:
			constants() :
				o_cy_1("o_{cy_1}", 2.62e-2),
				o_cy_2("o_{cy_2}", 3.14e-2),
				e_adh("e_{adh}", 1.31e-2)
			{}

			IConstant o_cy_1, o_cy_2, e_adh;
		};

		class Equation : public IEqMember<variables const&>{
			double_t o_cy_1, o_cy_2, e_adh;

		public:
			Equation():
				o_cy_1(2.62e-2),
				o_cy_2(3.14e-2),
				e_adh(1.31e-2)

			{ comp_amount = 3; ret_is = false; }

			const double_t ini_data[2] = { 0.5, 0.2 };

			[[nodiscard]] double_t Expression(variables const& u)noexcept final{
				B[1] = o_cy_1 * u.cy;
				B[2] = -o_cy_2 * u.adh * u.cy;
				B[3] = -e_adh * u.adh;
				return Sum_B();
			}

			void SynchronizeCoefs(constants & c)noexcept {
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
				c_Lm("c_{Lm}", 0.0920115),
				c_dLm("c_{d_{Lm}}", 0.0145086),
				K_Lm("K_{Lm}", 10.8869),
				d_Lm("d_{Lm}", 0.079728),
				l1("l_1", 1.0)
			{}
			IConstant c_Lm, c_dLm, K_Lm, d_Lm, l1; //*/
		};

		class Equation : public IEqMember<variables const&> {
			double_t c_Lm, c_dLm, K_Lm, d_Lm, l1; //*/
		public:
			Equation():
				c_Lm(9.63e-2), c_dLm(0.08), K_Lm(43.9), d_Lm(0.139), l1(4.0)
			{ comp_amount = 3; ret_is = true; }

			const double_t ini_data[2] = { 0.5, 0.0 };
			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1]= c_Lm * u.ret.adh_4;
				// B[2] = CHill(c_dLm, u.cy, K_Lm, l1) * u.lm;
				B[2] = CHill(c_dLm, u.cy, K_Lm, l1);
				B[3] = -d_Lm * u.lm;
				return Sum_B();
			}
			void SynchronizeCoefs(constants & c) noexcept {
				c.c_Lm.sync(c_Lm); 
				c.c_dLm.sync(c_dLm);
				c.K_Lm.sync(K_Lm);
				c.d_Lm.sync(d_Lm);
				c.l1.sync(l1);
			}
		};

	} // namespace LeuMacrophags

	namespace LeuNeutrophils
	{
		class constants {
		public:
			constants():
				//c_Ln("c_{Ln}", 0.0294695),
				c_Ln("c_{Ln}", 0.1095),
				d_Ln("d_{Ln}", 14.263),
				c_dLn1("c_{d_{Ln_1}}", 3.9366),
				K_Ln1("K_{Ln_1}", 6.2388),
				l1("l_1", 1.5)
			{}
			IConstant c_Ln, d_Ln, c_dLn1, K_Ln1, l1;
		};
		
		class Equation : public IEqMember<variables const&>{	
			double_t c_Ln, d_Ln, c_dLn1, K_Ln1, l1;
		public:
			Equation():
				c_Ln(1.47), d_Ln(2.0),
				c_dLn1(16.2), K_Ln1(3.96)
			{ comp_amount = 3; ret_is = false; }

			const double_t ini_data[2] = { 0.5, 0.0861204 };
			//const double_t ini_data[2] = { 0.5, 0.0 };

			[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
				B[1] = c_Ln * u.adh * u.cy;
				B[2] = CHill(c_dLn1, u.cy, K_Ln1, l1);			
				B[3] = - d_Ln * u.eps_w * u.ln;
				return Sum_B();
			}
			void SynchronizeCoefs(constants & c) noexcept {
				c.c_Ln.sync(c_Ln);
				c.d_Ln.sync(d_Ln);
				c.c_dLn1.sync(c_dLn1);
				c.K_Ln1.sync(K_Ln1);
				c.l1.sync(l1);
			}
		};
	} // namespace LeuNeutrophils

	namespace Microglia
	{
		class constants {
		public:
			constants() : 
				c_m_A("c_{m_A}", 0.06),
				c_m_N("c_{m_N}", 2.1),
				c_pro("c_{pro}", 2.5),
				T_M1("T_{M1}", 60.0),
				c_dMi("c_{d_{Mi}}", 0.2),
				K_Mi("K_{Mi}", 0.1)
				{}

			IConstant c_m_A, c_m_N, c_pro,
				T_M1, c_dMi, K_Mi ;
		};

		namespace Active
		{
			class Equation : public IEqMember<variables const&>{
				double_t c_m_A, c_m_N, c_dMi, K_Mi, c_pro, T_M1;

			public:
				Equation()
				{ comp_amount = 4; ret_is = false; }
				const double_t ini_data[2] = { 0.5, 0.01 };

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = u.mii * c_m_A * u.acu_c;
					B[2] = u.mii * c_m_N * u.nec;
					B[3] = Hill(c_dMi, u.cy, K_Mi, 1) * u.mii;
					B[4] = -c_pro * u.mia / T_M1;
					return Sum_B();
				}
				void SynchronizeCoefs(constants & c) noexcept {
					c.c_m_A.sync(c_m_A); 
					c.c_m_N.sync(c_m_N);
					c.c_dMi.sync(c_dMi);
					c.K_Mi.sync(K_Mi);
					c.c_pro.sync(c_pro); 
					c.T_M1.sync(T_M1);
				}
			};

		}
		namespace Inactive
		{
			class Equation : public IEqMember<variables const&> {
				double_t c_m_A, c_m_N, c_dMi, K_Mi, c_pro, T_M1;

			public:
				Equation():
					c_m_A(0.06), c_m_N(2.1), c_dMi(0.2), K_Mi(0.1),
					c_pro(2.5), T_M1(60.0)
				{ comp_amount = 4; ret_is = false; }

				const double_t ini_data[2] = { 0.5, 0.99 };

				[[nodiscard]] inline double_t Expression(variables const& u) noexcept final {
					B[1] = - u.mii * c_m_A * u.acu_c;
					B[2] = - u.mii * c_m_N * u.nec;
					B[3] = - Hill(c_dMi, u.cy, K_Mi, 1) * u.mii;
					B[4] = c_pro * u.mia / T_M1;
					return Sum_B();
				}
				void SynchronizeCoefs(constants & c) noexcept {
					c.c_m_A.sync(c_m_A); 
					c.c_m_N.sync(c_m_N);
					c.c_dMi.sync(c_dMi);
					c.K_Mi.sync(K_Mi);
					c.c_pro.sync(c_pro);
					c.T_M1.sync(T_M1);
				}
			};
		}
	} // namespace Microglia

	template<typename argType>	class ISubVal : public IRightPart<argType> {};

	namespace ToxDamage
	{
		class constants{
		public:
			constants() : 
				p_q1("p_{q_1}", 0.1),
				p_q2("p_{q_2}", 0.4),
				p_q3("p_{q_3}", 0.6),
				p_q4("p_{q_4}", 0.1),

				c_q1 ("c_{q_1}", 13.206),
				c_q2("c_{q_2}", 1.32),
				c_q3("c_{q_3}", 0.6),
				c_q4("c_{q_4}", 0.1),

				D_0("D_0", 0.033 * 1.01)
			{}

		IConstant p_q1, p_q2 , p_q3 , p_q4 ,
		c_q1 , c_q2 , c_q3 , c_q4 ,	D_0; 

		};

		namespace Full
		{
		
			class SubVal : public ISubVal<variables const&>{
				double_t p_q1, p_q2, p_q3, p_q4,
					c_q1, c_q2, c_q3, c_q4;

			public:
				SubVal()
				{ comp_amount = 4; ret_is = false; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = Hill(p_q1, u.cy, c_q1, 1);
					B[2] = Hill(p_q2, u.ln, c_q2, 1);
					B[3] = Hill(p_q3, u.lm, c_q3, 1);
					B[4] = Hill(p_q4, u.mia, c_q4, 1);
					return Sum_B();
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.p_q1.sync(p_q1);
					c.p_q2.sync(p_q2);
					c.p_q3.sync(p_q3);
					c.p_q4.sync(p_q4);

					c.c_q1.sync(c_q1);
					c.c_q2.sync(c_q2);
					c.c_q3.sync(c_q3);
					c.c_q4.sync(c_q4);
				}
			};

		}

		namespace Initial
		{
			class SubVal : public ISubVal<variables const&>{
				double_t D_0;
			public:
				SubVal():D_0(0.033 * 1.01)
				{ ret_is = false; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return GetPositive(u.d_F - D_0);
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.D_0.sync(D_0);
				}
			};
		}

		namespace Nec_partial
		{
			class SubVal : public ISubVal<variables const&>{
				double_t D_0;

			public:
				SubVal():D_0(0.033 * 1.01)				
				{ comp_amount = 1; ret_is = false; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = GetPositive(u.d_F - D_0);
					return B[1];
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.D_0.sync(D_0);
				}
			};

		}

		namespace Apop_partial
		{
	
			class SubVal : public ISubVal<variables const&>{
				double_t D_0;
			public:

				SubVal(): D_0(0.033 * 1.01)
				{ comp_amount = 1; ret_is = false; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					B[1] = u.d_F * indicator(D_0, u.d_F);
					return B[1];
				}
				void SynchronizeCoefs(constants & c) noexcept {
					c.D_0.sync(D_0);
				}
			};
		}

	} // namespace ToxDamage

	namespace Phagocytosis
	{
		class constants {
		public:
			constants():
				e_Ma("e_{Ma}", 0.05),
				e_Mi ("e_{Mi}", 1e-6),
				e_Lm("e_{Lm}", 0.0332),
				e_Ln ("e_{Ln}", 5.68e-3)
			{}

			IConstant e_Ma, e_Mi, e_Lm, e_Ln;
		};
		


		namespace Strong
		{
			
			class SubVal : public ISubVal<variables const&>{
				double_t e_Ma, e_Mi, e_Ln, e_Lm;

			public:
				SubVal() { comp_amount = 4; ret_is = false; }

				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					B[1] = e_Ma * u.mia;
					B[2] = e_Mi * u.mii;
					B[3] = e_Lm * u.lm;
					B[4] = e_Ln * u.ln;
					return Sum_B();
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.e_Ma.sync(e_Ma);
					c.e_Mi.sync(e_Mi); 
					c.e_Lm.sync(e_Lm);
					c.e_Ln.sync(e_Ln);
				}
			};
		}
		namespace Weak
		{
			class SubVal : public ISubVal<variables const&>{
				double_t e_Ma, e_Mi, e_Lm;

			public:
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					return e_Mi * u.mii + e_Ma * u.mia + e_Lm * u.lm;
				}

				void SynchronizeCoefs(constants & c) noexcept {
					c.e_Ma.sync(e_Ma);
					c.e_Mi.sync(e_Mi);
					c.e_Lm.sync(e_Lm);
				}

			};
		}

	} // namespace Phagocytosis

} // namespace StraightTask