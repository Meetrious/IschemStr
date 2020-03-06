#pragma once
#include <functional>	// полиморфные обёртки функций

#include "ISCHEM_variables.h"



namespace StraightTask
{
	struct Neurons
	{

		Neurons()
		{}

		~Neurons()
		{
			if (OutNec.is_open()) OutNec.close();
			if (OutAcu.is_open()) OutAcu.close();
			if (OutAp_s.is_open()) OutAp_s.close();
			if (OutAp_e.is_open()) OutAp_e.close();
			if (OutHel.is_open()) OutHel.close();
		}



		// уравнения динамики некротической плотности
		static const std::array <std::function<double_t(variables const &)>, 6>  necrotic_cells;
		static const std::array <std::function<double_t(variables const &)>, 6>  acute_changes;
		static const std::array <std::function<double_t(variables const &)>, 3>  apoptose_started;
		static const std::array <std::function<double_t(variables const &)>, 3>  apoptose_ended;
		static const std::array <std::function<double_t(variables const &)>, 6>  intact_cells;

		std::ofstream OutNec;
		std::ofstream OutAcu;
		std::ofstream OutAp_s;
		std::ofstream OutAp_e;
		std::ofstream OutHel;

		
		static double_t p_N, Rep, k_N, k_A, p_R;
		static double_t q_N, q_A_n, c_H, c_A, q_epN, q_epA;
		static double_t q_A, q_A_h, p_R3, q_H;
		static double_t* D_0;

	};

	struct Microglia
	{
		Microglia() {}
		~Microglia()
		{
			if (OutMia.is_open()) OutMia.close();
			if (OutMii.is_open()) OutMii.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 4> mi_active;
		static const std::array <std::function<double_t(variables const &)>, 5> mi_inactive;

		std::ofstream OutMia;
		std::ofstream OutMii;


		static double_t p_1, c_A, c_N, c_pro;
		static double_t T_M1, T_M2, c_Mi1, c_Mi2;
		static double_t c_dMi, c_Mi, K_Mi;

	};

	struct Cytokines
	{
		Cytokines() {}

		void CollectData(std::string const & data_dir)
		{
			std::ifstream in(data_dir);
			if (!in)
			{
				// ПИШИ ОБРАБОТКУ!
			}
			double_t tmp;
			while (!in.eof())
			{
				in >> tmp; this->data[0].emplace_back(tmp);
				in >> tmp; this->data[1].emplace_back(tmp);
			}
			in.close();
		}

		~Cytokines()
		{
			if (OutCY.is_open()) OutCY.close();
			if (OutPsy.is_open()) OutPsy.close();
			if (OutCH.is_open()) OutCH.close();
		}


		static const double_t GetPsy(variables const & u)noexcept;

		static const std::array <std::function<double_t(variables const &)>, 2> chemokines;
		static const std::array <std::function<double_t(variables const &)>, 3> cytokines;

		std::ofstream OutCY;
		std::ofstream OutPsy;
		std::ofstream OutCH;


		// хемотаксические
		static double_t p_Mach, C_Ma, p_Lmch, C_Lm, e_ch;
		// нормальные
		static double_t e_cy, C_Ln;
		// psy
		static double_t p_xcy0, cy_max, T_, t_0;

		std::array<std::vector<double_t>, 2> data;

	};

	struct Adhension
	{

		Adhension() {}
		~Adhension()
		{
			if (OutAdh.is_open()) OutAdh.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 2> adhension;

		std::ofstream OutAdh;


		static double_t o_cy_1, o_cy_2, e_adh;

	};

	struct LeuMacrophags
	{
		LeuMacrophags() {}
		~LeuMacrophags()
		{
			if (OutLm.is_open()) OutLm.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 5> macrophags;

		std::ofstream OutLm;


		static double_t c_Lm, p_dLm, T_Lm;
		static double_t c_dLm, K_Lm, d_Lm;
	};

	struct LeuNeutrophils
	{
		LeuNeutrophils() {}
		~LeuNeutrophils()
		{
			if (OutLn.is_open()) OutLn.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 5> neutrophils;

		std::ofstream OutLn;


		static double_t c_Ln, p_dLn, T_Ln;
		static double_t c_dLn, K_Ln, d_Ln;
		static double_t c_dLn1, K_Ln1;
		static double_t c_dLn2, K_Ln2;
	};

	struct ToxDamage
	{
		ToxDamage() {}
		~ToxDamage()
		{
			if (OutD_full.is_open()) OutD_full.close();
			if (OutD_ini.is_open()) OutD_ini.close();
			if (OutDPN.is_open()) OutDPN.close();
			if (OutDPA.is_open()) OutDPA.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 4> full;
		static const std::array <std::function<double_t(variables const &)>, 1> initial;
		static const std::array <std::function<double_t(variables const &)>, 2> nec_partial;
		static const std::array <std::function<double_t(variables const &)>, 2> A_partial;

		std::ofstream OutD_full;
		std::ofstream OutD_ini;
		std::ofstream OutDPN;
		std::ofstream OutDPA;

	//private:

		[[nodiscard]] static inline double_t State_dpN(const variables & u)
		{
			if (u.d_F <= D_0) return p_D * u.d_F;
			else return u.d_F - (D_0 - p_D * u.d_F);
		}
		[[nodiscard]] static inline double_t State_dpA(const variables & u)
		{
			if (u.d_F <= D_0) return u.d_F - p_D * u.d_F;
			else return D_0 - p_D * u.d_F;
		}
		[[nodiscard]] static inline double_t GetPositive(double_t const & val)
		{
			if (val > 0) return val;
			else return 0;
		}
		[[nodiscard]] static inline double_t St_dpA(double_t d_F, double_t D_0)
		{
			if (d_F <= D_0) return d_F;
			else return 0;
		}

		static double_t p_ncy, p_Ln, C_DLn, P_nn;
		static double_t p_Lm, C_DLm, C_D, C_Dcy;

		static double_t p_q1, p_q2, p_q3, p_q4;
		static double_t c_q1, c_q2, c_q3, c_q4;


		static double_t D_0, p_D;
	};

	struct Phagocytosis
	{
		Phagocytosis() {}
		~Phagocytosis()
		{
			if (OutEps.is_open()) OutEps.close();
		}

		static const std::array <std::function<double_t(variables const &)>, 2> phagocytosis;

		std::ofstream OutEps;
		std::ofstream OutePS;
	
	//private:

		static double_t e_Ma, e_Ln, e_Lm, e_Mi;
	};

#include "DefaultCoefValues.h"

	[[nodiscard]] inline double_t Hill(double_t C, double_t ed, double_t K, int32_t n)
	{
		return C * std::pow(ed, n) / (K + std::pow(ed, n));
	}

	template <std::size_t SIZE> using RPartsDef =
		std::array <std::function<double_t(const variables &u)>, SIZE>;

#define FUN(EQUAT) [](const variables & u) noexcept -> double_t {return EQUAT;}

	inline const RPartsDef<6> Neurons::necrotic_cells = {
		/*0*/	FUN(Splines::GetValue(u.tim, u.SplineData[1][0], u.SplineData[1][1], u.SplineData[1][2], u.SplineData[1][3], u.SplineData[1][4])),

		/*1*/	FUN(p_N * u.d_ini * u.hel - u.eps * u.nec),
		/*2*/	FUN(u.d_ini * u.hel - u.eps * u.nec),
		/*3*/	FUN(k_N * u.dp_N * u.hel - u.eps * u.nec),
		/*4*/	FUN(q_N * u.dp_N * Hill(1, u.hel, c_H, 2) * u.nec + q_A_n * u.dp_N * Hill(1, u.acu_c, c_A, 5) - q_epN * u.eps * u.nec ), // модель 3
		/*5*/	FUN(0)
	};
	inline const RPartsDef<6> Neurons::acute_changes{
/*0*/	FUN(Splines::GetValue(u.tim, u.SplineData[2][0], u.SplineData[2][1], u.SplineData[2][2], u.SplineData[2][3], u.SplineData[2][4])),

/*1*/	FUN((1.0 - p_N) * u.d_ini * u.hel - u.eps * u.acu_c),
/*2*/	FUN(u.d_ini * u.hel - u.acu_c * (Rep + u.eps)),
/*3*/	FUN(k_A * u.dp_A * u.hel - u.acu_c * (p_R * k_A * u.dp_A + u.eps)),
/*4*/	FUN(q_A * u.dp_N * u.hel * ( 1-p_R3 ) * Hill(1, u.acu_c, q_A_h, 1) + q_H * u.dp_N * u.hel * u.acu_c - q_A_n * u.dp_N * Hill(1, u.acu_c, c_A, 5) - q_epA * u.eps * u.acu_c), // модель 3
/*5*/	FUN(0)
	};
	inline const RPartsDef<3> Neurons::apoptose_started = {
/*0*/	FUN(0),

/*1*/	FUN((1.0 - p_N) * u.d_ini * u.hel - (1.0 - p_N) * u.d_12 * u.hel_12),
/*2*/	FUN(u.d_ini * u.hel - u.d_12 * u.a_12 - Rep * u.ap_s)
	};
	inline const RPartsDef<3> Neurons::apoptose_ended = {
/*0*/	FUN(0),

/*1*/	FUN((1.0 - p_N) * u.d_12 * u.hel_12 - u.eps * u.ap_e),
/*2*/	FUN(u.d_12 * u.a_12 - u.eps * u.ap_s)
	};
	inline const RPartsDef<6> Neurons::intact_cells = {
		/*0*/	FUN(Splines::GetValue(u.tim, u.SplineData[5][0], u.SplineData[5][1], u.SplineData[5][2], u.SplineData[5][3], u.SplineData[5][4])),

		/*1*/	FUN(-u.hel * u.d_ini),
		/*2*/	FUN(-u.hel * (*D_0 + u.d_ini) + Rep * u.ap_s),
		/*3*/	FUN(-u.hel * (k_N * u.dp_N + k_A * u.dp_A) + p_R * k_A * u.dp_A * u.acu_c),
		/*4*/	FUN( u.hel * q_A * u.dp_A * (-1 + p_R3 ) * Hill(1, u.acu_c, q_A_h, 1) - q_N * u.dp_N * Hill(1, u.hel, c_H, 2) * u.nec - q_H * u.dp_N * u.hel * u.acu_c), //	модель 3
		/*5*/	FUN(0)
	};

	inline const RPartsDef<3> Cytokines::cytokines = {
/*0*/	FUN(Splines::GetValue(u.tim, u.SplineData[6][0], u.SplineData[6][1], u.SplineData[6][2], u.SplineData[6][3], u.SplineData[6][4])),

/*1*/	FUN((Hill(u.psy, u.mia, C_Ma, 1) + Hill(u.psy, u.lm, C_Lm, 1)) * (u.nec + u.ap_e) - e_cy * u.cy),
/*2*/	FUN((Hill(u.psy, u.mia, C_Ma, 1) + Hill(u.psy, u.lm, C_Lm, 1) + Hill(2, u.ln, C_Ln, 1)) - e_cy * u.cy)
	};
	inline const RPartsDef<2> Cytokines::chemokines = {
/*0*/	FUN(0),

/*1*/	FUN((Hill(p_Mach, u.mia, C_Ma, 1) + Hill(p_Lmch, u.lm, C_Lm, 1)) * (u.nec + u.ap_e) - e_ch * u.ch),
	};
	inline const RPartsDef<2> Adhension::adhension = {
/*0*/	FUN(0),

/*1*/	FUN(o_cy_1 * u.cy - o_cy_2 * u.adh * u.cy - e_adh * u.adh)
	};

	inline const double_t Cytokines::GetPsy(variables const & u) noexcept
	{
		return p_xcy0 * ((cy_max - u.cy) / (T_ + (u.tim / t_0)));
	}

	inline const RPartsDef<4> Microglia::mi_active = {
		/*0*/	FUN(0),

		/*1*/	FUN( (c_A * u.ap_e + c_N * u.nec) * u.mii - c_pro * (u.mia / T_M1)),
		/*2*/	FUN( p_1 * (c_A * u.acu_c + c_N * u.nec) * u.mii + Hill(c_Mi, u.cy, K_Mi, 1) * u.mia - c_pro * (u.mia / T_M1)),
		/*3*/	FUN( (c_A * u.acu_c + c_N * u.nec) * u.mii + Hill(c_dMi, u.cy, K_Mi, 1) * u.mii - c_pro * u.mia / T_M1) // модель 3
	};
	inline const RPartsDef<5> Microglia::mi_inactive = {
		/*0*/	FUN(0),

		/*1*/	FUN(-(c_A * u.ap_e + c_N * u.nec) * u.mii + c_pro * (u.mia / T_M1) - c_dMi * u.mii),
		/*2*/	FUN(-p_1 * (c_A * u.ap_e + c_N * u.nec) * u.mii + c_pro * (u.mia / T_M1) - c_dMi * u.mii),
		/*3*/	FUN(c_pro * (u.mia / T_M1) - (c_A * u.acu_c + c_N * u.nec) * u.mii - Hill(c_Mi, u.cy, K_Mi, 1) * u.mii),
		/*4*/	FUN( -(c_A * u.acu_c + c_N * u.nec) * u.mii + c_pro * u.mia / T_M1 - Hill(c_dMi, u.cy, K_Mi, 1) * u.mii) // модель 3
	};
	inline const RPartsDef<5> LeuMacrophags::macrophags = {
		/*0*/	FUN(Splines::GetValue(u.tim, u.SplineData[11][0], u.SplineData[11][1], u.SplineData[11][2], u.SplineData[11][3], u.SplineData[11][4])),

		/*1*/	FUN( c_Lm * u.adh_24 - p_dLm * (u.lm / T_Lm)),
		/*2*/	FUN( c_Lm * u.adh_4 - p_dLm * (u.lm / T_Lm)),
		/*3*/	FUN( c_Lm * u.adh_4 + Hill(c_dLm, u.cy, K_Lm, 5) - p_dLm * (u.lm / T_Lm)), /* + Hill(C[cLmN], Nec, C[KN1], 1))*/
		/*4*/	FUN( c_Lm * u.adh_4 + (Hill(c_dLm, u.cy, K_Lm, 4) - d_Lm) * u.lm) // модель 3
	};
	inline const RPartsDef<5> LeuNeutrophils::neutrophils = {
		/*0*/	FUN( Splines::GetValue(u.tim, u.SplineData[12][0], u.SplineData[12][1], u.SplineData[12][2], u.SplineData[12][3], u.SplineData[12][4])),

		/*1*/	FUN( c_Ln * u.adh_12 - p_dLn * (u.ln / T_Ln)),
		/*2*/	FUN( c_Ln * u.adh - p_dLn * (u.ln / T_Ln)),
		/*3*/	FUN( (1. / 4.) * c_Ln * u.adh + Hill(c_dLn, u.cy, K_Ln, 1) - u.eps * u.ln * 320.),
		/*4*/	FUN( c_Ln * u.adh * u.ln + Hill( c_dLn1, u.cy, K_Ln1, 1 ) + Hill( c_dLn2, u.cy, K_Ln2, 1 ) * u.ln - d_Ln * u.Eps * u.ln  ) // модель 3
	};

	inline const RPartsDef<4> ToxDamage::full = {
		/*0*/	FUN(p_ncy * u.cy + Hill(p_Ln, u.ln, C_DLn, 1) * (u.nec + u.ap_e) + P_nn * u.nec),

		/*1*/	FUN(p_ncy * u.cy + Hill(p_Ln, u.ln, C_DLn, 1) * (u.nec + u.acu_c) + P_nn * u.nec),
		/*2*/	FUN(Hill(p_ncy, u.cy, C_Dcy, 1) + (Hill(p_Lm, u.lm, C_DLm, 1) + Hill(p_Ln, u.ln, C_DLn, 1))  * (u.nec + u.acu_c) / C_D + P_nn * u.nec / C_D),
		/*3*/	FUN(Hill(p_q1,u.cy,c_q1,1) + Hill(p_q2,u.ln,c_q2,1) + Hill(p_q3,u.lm ,c_q3,1) + Hill(p_q4,u.mia,c_q4,1) ) // модель 3
	};
	inline const RPartsDef<1> ToxDamage::initial = {
		FUN(GetPositive(u.d_F - D_0))
	};
	inline const RPartsDef<2> ToxDamage::nec_partial = {
		FUN(GetPositive(State_dpN(u))),
		FUN(GetPositive(u.d_F - D_0))
	};
	inline const RPartsDef<2> ToxDamage::A_partial = {
		FUN(GetPositive(State_dpA(u))),
		FUN(St_dpA(u.d_F, D_0))
	};

	inline const RPartsDef<2> Phagocytosis::phagocytosis = {
		FUN(e_Mi * u.mii + e_Ma * u.mia + e_Lm * u.lm + e_Ln * u.ln ),
		FUN(e_Mi * u.mii + e_Ma * u.mia + e_Lm * u.lm)
	};

#undef FUN
}