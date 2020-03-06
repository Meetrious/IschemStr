#pragma once
#include "Equations.h"

namespace StraightTask
{

	struct Methods :
		Neurons,
		Cytokines, Adhension,
		LeuMacrophags, LeuNeutrophils, Microglia,
		ToxDamage, Phagocytosis
	{

		Methods() {}

		void AllocateOutputStream()
		{
#define OpenCheck(OUTSTREAM, ENUM) \
		OUTSTREAM.open(working_directory + ST_sol + Parameters::sol_names[ENUM], 2);\
		if(!OUTSTREAM) \
		{ std::cout << "\n" << #ENUM << " -- stream was not allocated for some reason. Care to attend!"; system("pause"); }
			OpenCheck(OutNec, NEC)
			//OpenCheck(OutAp_s,AP_S)
			//OpenCheck(OutAp_e, AP_E)
			OpenCheck(OutAcu, ACU)
			OpenCheck(OutHel, HEL)

			OpenCheck(OutCY, CY)
			//OpenCheck(OutCH,CH)
			OpenCheck(OutAdh, ADH)

			OpenCheck(OutMia, MIA)
			OpenCheck(OutMii, MII)
			OpenCheck(OutLm, LM)
			OpenCheck(OutLn, LN)

			OpenCheck(OutD_full, D_F)
//			OpenCheck(OutD_ini, D_INI)
			OpenCheck(OutDPN, DP_N)
			OpenCheck(OutDPA, DP_A)
			OpenCheck(OutePS, ePS)
			OpenCheck(OutEps, EPS)
//			OpenCheck(OutPsy, PSY)
#undef OpenCheck
		}

		~Methods() = default;

		static double_t H;

		// метод для изменения значений коэффициентов в уравнениях
		static void setCoefs(double_t* CFC[], int16_t CFC_cap, std::string way)
		{
			std::vector<double> inCoefsVals; double_t tmp;
			std::ifstream in(way);
			while (!in.eof())
			{
				in >> tmp;
				inCoefsVals.emplace_back(tmp);
			}
			if (inCoefsVals.size() != CFC_cap)
			{	
				std::cout << "\n stop this MADNES RIGHT HERE, criminal scum!, or I'm throwing the exception!";
				for (size_t i = 0; i < CFC_cap; i++)
					*CFC[i] = 0;
				return;
			}
			for (size_t i = 0; i < CFC_cap; i++)
				*CFC[i] = inCoefsVals[i];
			return;
		}

		static void direct(double_t & U, const variables & X, std::function<double_t(variables const &)> const &Eq)
		{
			U = Eq(X);
		}

		// кол-во уравнений в системе = N_eq - 1
		constexpr static size_t N_eq = 13;

		// массив вызовов сплайнов 
		static std::array<std::function<void(variables &)>, Methods::N_eq> Splines;

		// массив вызовов вспомагательных величин
		static std::array<std::function<void(variables &)>, 5> SubValues;

		void OutputCurSol(std::ofstream & out, size_t Nj, double_t tim, double_t sol_val)
		{
			out << Nj << "\t\t\t" << std::setprecision(5) << tim
				<< "\t\t\t" << std::setprecision(15) << sol_val << "\n";
		}

#define OUTP(STREAM, VALUE) [this](variables const & X, size_t Nj) -> void {OutputCurSol(STREAM, Nj, X.tim, VALUE); }
		std::array<std::function<void(variables const &, size_t)>, Methods::N_eq + 5 > Output = {
			[](variables const & X, size_t Nj) -> void {},
			OUTP(OutNec, X.nec),
			OUTP(OutAcu, X.acu_c),
			[](variables const & X, size_t Nj) -> void {}, //OUTP(OutAp_s, X.ap_s),
			[](variables const & X, size_t Nj) -> void {}, //OUTP(OutAp_e, X.ap_e),
			OUTP(OutHel, X.hel),

			OUTP(OutCY, X.cy),
			[](variables const & X, size_t Nj) -> void {}, //OUTP(OutCH, X.ch),
			OUTP(OutAdh, X.adh),

			OUTP(OutMia, X.mia),
			OUTP(OutMii, X.mii),
			OUTP(OutLm, X.lm),
			OUTP(OutLn, X.ln),

			OUTP(OutD_full, X.d_F),
			// [](variables const & X, size_t Nj) -> void {}, //OUTP(OutD_ini, X.d_ini),
			OUTP(OutDPN, X.dp_N),
			OUTP(OutDPA, X.dp_A),
			OUTP(OutEps, X.Eps),
			OUTP(OutePS,X.eps),
			//[](variables const & X, size_t Nj) -> void {}, //OUTP(OutPsy,X.psy),
		};
#undef OUTP


		struct PredCor
		{
			PredCor(Parameters const & CurSTPar)
			{
				// определяем исходное состояние запаздываний
				rets::SetInitialState(CurSTPar.N);

				X_init.MakeVarsInitial(CurSTPar); // инициализируем начальные условия в X_init
				X_init.AllocCurRets(0, CurSTPar.N); // инициализируем значения запаздывающих аргументов

				// инициализируем значения вспомогательных величин от выше определённого
				for (size_t i = 0; i < 5; i++) Methods::SubValues[i](X_init);

				X_prev = X_init;
				X_pred = X_prev;
				X_cor = X_prev;
			}

			void ResetToInitial(Parameters const & CurSTPar)
			{
				// определяем исходное состояние запаздываний

				//rets::SetInitialState(CurSTPar.N);

				X_init.MakeVarsInitial(CurSTPar); // инициализируем начальные условия в X_init
				
				// X_init.AllocCurRets(0, CurSTPar.N); // инициализируем значения запаздывающих аргументов

				// инициализируем значения вспомогательных величин от выше определённого
				for (size_t i = 0; i < 6; i++) Methods::SubValues[i](X_init);

				X_prev = X_init;
				X_pred = X_prev;
				X_cor = X_prev;
			}

			~PredCor() = default;

			variables X_init, X_prev, X_pred, X_cor;


			// массив вызовов предиктора
			static std::array<std::function<void(variables &, const variables &)>, Methods::N_eq> Predictions;

			// массив вызовов корректора
			static std::array<std::function<void(variables &, const variables &, const variables &)>, Methods::N_eq> Corrections;

		private:

			static void predictor
			(double_t &pred_U,
				const double_t prev_U,
				const variables &prev_X,
				std::function<double_t(variables const &)> const &Eq
			)
			{
				pred_U = prev_U + H * Eq(prev_X);
			}

			static void corrector
			(double_t &cor_U,
				const double_t prev_U,
				const variables &prev_X,
				const variables &pred_X,
				std::function<double_t(variables const &)> const &Eq
			)
			{
				cor_U = prev_U + H * (Eq(prev_X) + Eq(pred_X)) * 0.5;
			}

		};
	};

	double_t Methods::H = 0.016;

	std::array<std::function<void(variables &)>, Methods::N_eq> Methods::Splines = {
		[](variables & X) -> void {},

		[](variables & X) -> void { direct(X.nec, X,Neurons::necrotic_cells[0]); },
		[](variables & X) -> void { direct(X.acu_c, X,Neurons::acute_changes[0]); },
		[](variables & X) -> void { X.ap_s = 0; },
		[](variables & X) -> void { X.ap_e = 0; },
		[](variables & X) -> void { direct(X.hel, X, Neurons::intact_cells[0]); },

		
		[](variables & X) -> void { direct(X.cy, X, Cytokines::cytokines[0]); },
		[](variables & X) -> void { X.ch = 0; },
		[](variables & X) -> void { X.adh = 0; },

		[](variables & X) -> void { X.mia = 0; },
		[](variables & X) -> void { X.mii = 0; },
		[](variables & X) -> void { direct(X.lm, X, LeuMacrophags::macrophags[0]); },
		[](variables & X) -> void { direct(X.ln, X, LeuNeutrophils::neutrophils[0]); }
	};


#define PRED(VAL) [](variables &pred_X, const variables &prev_X) -> void {VAL;}
	std::array <std::function<void(variables &, const variables &)>, Methods::N_eq> Methods::PredCor::Predictions = {
		PRED(0),

		PRED(predictor(pred_X.nec, prev_X.nec, prev_X, Neurons::necrotic_cells[4])),
		PRED(predictor(pred_X.acu_c, prev_X.acu_c, prev_X, Neurons::acute_changes[4])),
		PRED(0), //PRED(predictor(pred_X.ap_s, prev_X.ap_s, prev_X, Neurons::apoptose_started[1])),
		PRED(0), //PRED(predictor(pred_X.ap_e, prev_X.ap_e, prev_X, Neurons::apoptose_ended[1])),
		PRED(predictor(pred_X.hel, prev_X.hel, prev_X, Neurons::intact_cells[4])),

		PRED(0), //PRED(predictor(pred_X.cy, prev_X.cy, prev_X, Cytokines::cytokines[2])),
		PRED(0), //PRED(predictor(pred_X.ch, prev_X.ch, prev_X, Cytokines::chemokines[1])),
		PRED(predictor(pred_X.adh, prev_X.adh, prev_X, Adhension::adhension[1])),

		PRED(predictor(pred_X.mia, prev_X.mia, prev_X, Microglia::mi_active[3])),
		PRED(predictor(pred_X.mii, prev_X.mii, prev_X, Microglia::mi_inactive[4])),
		PRED(predictor(pred_X.lm, prev_X.lm, prev_X, LeuMacrophags::macrophags[4])),
		PRED(predictor(pred_X.ln, prev_X.ln, prev_X, LeuNeutrophils::neutrophils[4]))
	};
#undef PRED

#define COR(VAL) [](variables &cor_X, const  variables &prev_X, const variables &pred_X ) -> void {VAL;}
	std::array <std::function<void(variables &, const variables &, const variables &)>, Methods::N_eq> Methods::PredCor::Corrections = {
		COR(0),

		COR(corrector(cor_X.nec, prev_X.nec, prev_X, pred_X, Neurons::necrotic_cells[4])),
		COR(corrector(cor_X.acu_c, prev_X.acu_c, prev_X, pred_X, Neurons::acute_changes[4])),
		COR(0), //COR(corrector(cor_X.ap_s, prev_X.ap_s, prev_X, pred_X, Neurons::apoptose_started[1])),
		COR(0), //COR(corrector(cor_X.ap_e, prev_X.ap_e, prev_X, pred_X, Neurons::apoptose_ended[1])),
		COR(corrector(cor_X.hel, prev_X.hel, prev_X, pred_X, Neurons::intact_cells[4])),

		COR(0), //COR(corrector(cor_X.cy, prev_X.cy, prev_X, pred_X, Cytokines::cytokines[2])),
		COR(0),	//COR(corrector(cor_X.ch, prev_X.ch, prev_X, pred_X, Cytokines::chemokines[1])),
		COR(corrector(cor_X.adh, prev_X.adh, prev_X, pred_X, Adhension::adhension[1])),

		COR(corrector(cor_X.mia, prev_X.mia, prev_X, pred_X, Microglia::mi_active[3])),
		COR(corrector(cor_X.mii, prev_X.mii, prev_X, pred_X, Microglia::mi_inactive[4])),
		COR(corrector(cor_X.lm, prev_X.lm, prev_X, pred_X, LeuMacrophags::macrophags[4])),
		COR(corrector(cor_X.ln, prev_X.ln, prev_X, pred_X, LeuNeutrophils::neutrophils[4])),
	};
#undef COR

	std::array<std::function<void(variables &)>, 5> Methods::SubValues = {
		[](variables & X) -> void { direct(X.Eps, X, Phagocytosis::phagocytosis[1]); },
		[](variables & X) -> void { direct(X.eps, X, Phagocytosis::phagocytosis[0]); },
		[](variables & X) -> void { direct(X.d_F, X, ToxDamage::full[3]); },
		[](variables & X) -> void { direct(X.dp_A, X, ToxDamage::A_partial[1]); },
		[](variables & X) -> void { direct(X.dp_N, X, ToxDamage::nec_partial[1]); },
	};



	void launchPredCor(Methods::PredCor & TL)
	{
		// предикция/сплайны
		for (size_t i = 1; i < Methods::N_eq; i++)
		{
			if (TL.X_pred.SIS[i]) Methods::Splines[i](TL.X_pred);
			else Methods::PredCor::Predictions[i](TL.X_pred, TL.X_prev);
		}

		// вспомогательные величины от вычисленных предварительных значений
		for (size_t i = 0; i < 5; i++) Methods::SubValues[i](TL.X_pred);

		// коррекция/сплайны
		for (size_t i = 1; i < Methods::N_eq; i++)
		{
			if (TL.X_pred.SIS[i]) TL.X_cor.SetMatchedVar(i, (*TL.X_pred.GetMatchedVar(i))); // Methods::Splines[i](cor_X);
			else Methods::PredCor::Corrections[i](TL.X_cor, TL.X_prev, TL.X_pred);
		}

		// вспомогательные величины от коррекции
		for (size_t i = 0; i < 5; i++) Methods::SubValues[i](TL.X_cor);
	}
}