#pragma once
#include <functional>	// поддержка полиморфных обёрток функций
#include <Settings.h>

namespace StraightTask
{
	class PhaseTrajOutput
	{
	public:
		std::ofstream str;
		PhaseTrajOutput(const char* name) { AllocateOutputStream(name); }

		void AllocateOutputStream(const char* name)
		{
			str.open(output_dir + "ST/SOL/PTraj/" + name + ".txt");
			if (!str)
			{
				std::cout << "\n For some reason budgets output stream for <" << name << "> wasn't allocated \n";
				getchar();
				return;
			}
		}
		void OutputPhaseTraj(uint32_t Nj, double_t Tj, vector<double_t> Vals) {
			str << Nj << "\t\t\t" << std::setprecision(5) << Tj << "\t\t\t" << std::setprecision(15);
			for (auto const& cur : Vals)
				str << cur << "\t\t\t";
			str << std::endl;
		}

		~PhaseTrajOutput() { if (str.is_open()) str.close(); }
	};

	class IAggregate
	{
	public:

		// list of equations in the processed system
		Neurons::NecroticCells::M_ODE NEC;
		Neurons::AcuteChanges::M_ODE AC;
		Neurons::Apoptosis::Started::M_ODE AS;
		Neurons::Apoptosis::Started::M_ODE AE;
		Neurons::IntactCells::M_ODE HEL;

		Cytokines::Pro_Inflam::M_ODE CY;
		Adhesion::M_ODE ADH;
		LeuMacrophags::M_ODE LM;
		LeuNeutrophils::M_ODE LN;

		Microglia::Active::M_ODE MIA;
		Microglia::Inactive::M_ODE MII;

		Test::OneDim::M_ODE EXP;

		Test::ThreeDim::X_m_ODE TCOS;
		Test::ThreeDim::Y_m_ODE TSIN;
		Test::ThreeDim::Z_m_ODE T;


		// subbordinate values relevant for a system 
		ToxDamage::Full::M_Sub DF;
		ToxDamage::Nec_partial::M_Sub DPN;
		ToxDamage::Apop_partial::M_Sub DPA;

		Phagocytosis::Strong::M_Sub EPS_S;
		Phagocytosis::Weak::M_Sub EPS_W;


		std::array<std::function<void()>, 14> DataCollector = {
			[&]()-> void {NEC.CollectSolData(); },
			[&]()-> void {AC.CollectSolData(); },
			[&]()-> void {HEL.CollectSolData(); },

			[&]()-> void {CY.CollectSolData(); },
			[&]()-> void {ADH.CollectSolData(); },
			[&]()-> void {LM.CollectSolData(); },
			[&]()-> void {LN.CollectSolData(); },
			[&]()-> void {MIA.CollectSolData(); },
			[&]()-> void {MII.CollectSolData(); },

			[&]()-> void {DF.CollectSolData(); },
			[&]()-> void {DPN.CollectSolData(); },
			[&]()-> void {DPA.CollectSolData(); },

			[&]()-> void {EPS_S.CollectSolData(); },
			[&]()-> void {EPS_W.CollectSolData(); }
		};

		std::vector<std::function<bool(variables&)>> SplineQueue = {
			[&](variables& X) -> bool { return NEC.QueueRule(X.nec, X.tj, X.spl_gap_counter); },
			[&](variables& X) -> bool { return AC.QueueRule(X.acu_c, X.tj, X.spl_gap_counter); },
			[&](variables& X) -> bool { return HEL.QueueRule(X.hel, X.tj, X.spl_gap_counter); },
			[&](variables& X) -> bool { return CY.QueueRule(X.cy, X.tj, X.spl_gap_counter); },
			[&](variables& X) -> bool { return LM.QueueRule(X.lm, X.tj, X.spl_gap_counter); },
			[&](variables& X) -> bool { return LN.QueueRule(X.ln, X.tj, X.spl_gap_counter); }
		};

		std::array<std::function<void(variables&)>, 9> IniDataInitialiser = {
			/*[&](variables& X) -> void {X.x = TCOS.GetInitialData(); },
			[&](variables& X) -> void {X.y = TSIN.GetInitialData(); },
			[&](variables& X) -> void {X.z = T.GetInitialData(); }// */

			//[&](variables& X) -> void {X.x = EXP.GetInitialData(); }

			[&](variables & X) -> void {X.nec = NEC.GetInitialData(); },
			[&](variables & X) -> void {X.acu_c = AC.GetInitialData(); },
			[&](variables & X) -> void {X.hel = HEL.GetInitialData(); },
			[&](variables & X) -> void {X.cy = CY.GetInitialData(); },
			[&](variables & X) -> void {X.adh = ADH.GetInitialData(); },
			[&](variables & X) -> void {X.lm = LM.GetInitialData(); },
			[&](variables & X) -> void {X.ln = LN.GetInitialData(); },
			[&](variables & X) -> void {X.mia = MIA.GetInitialData(); },
			[&](variables & X) -> void {X.mii = MII.GetInitialData(); }// 9шт */
		};

		std::array<std::function<void(variables&)>, 3> IniRetInitialiser = {
			/*[&](variables& X) -> void {X.ret.x_pi2 = TCOS.GetRetValue(pi * 0.5, 0); },
			[&](variables& X) -> void {X.ret.y_pi2 = TSIN.GetRetValue(pi * 0.5, 0); },
			[&](variables& X) -> void {X.ret.z_pi2 = T.GetRetValue(pi * 0.5, 0); } // */

			//[&](variables& X) -> void {X.ret.x_1 = EXP.GetRetValue(1.0, 0); }

			//[&](variables& X) -> void {X.ret.x_pi2 = EXP.GetRetValue(pi / 2, 0); }

			[&](variables& X) -> void {X.ret.hel_12 = HEL.GetRetValue(12.0, 0); },
			[&](variables& X) -> void {X.ret.adh_12 = ADH.GetRetValue(12.0, 0); },
			[&](variables& X) -> void {X.ret.adh_4 = ADH.GetRetValue(4.0, 0); } // */

		};

		std::array<std::function<void(uint16_t, uint32_t, uint32_t, variables&)>, 14> SolDataInitialiser = {
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.nec = NEC.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.acu_c = AC.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.hel = HEL.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.cy = CY.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.adh = ADH.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.lm = LM.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.ln = LN.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.mia = MIA.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.mii = MII.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.d_F = DF.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.dp_N = DPN.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.dp_A = DPA.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.eps_s = EPS_S.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.eps_w = EPS_W.GetSolData(day, Nj, N); }
		};

		std::array<std::function<void(variables&)>, 5> ExpressSubValues = {
			[&](variables& X) -> void { X.d_F = DF.RP.Expression(X); },
			[&](variables& X) -> void { X.dp_N = DPN.RP.Expression(X); },
			[&](variables& X) -> void { X.dp_A = DPA.RP.Expression(X); },
			[&](variables& X) -> void { X.eps_s = EPS_S.RP.Expression(X); },
			[&](variables& X) -> void { X.eps_w = EPS_W.RP.Expression(X); }
		};

		std::array<std::function<void()>, 14> OutStreamAllocator = {
			/*[&]() -> void {TCOS.AllocateOutputStreams(); },
			[&]() -> void {TSIN.AllocateOutputStreams(); },
			[&]() -> void {T.AllocateOutputStreams(); } // */


			//[&]() -> void {EXP.AllocateOutputStreams(); }

			[&]() -> void {NEC.AllocateOutputStreams(); },
			[&]() -> void {AC.AllocateOutputStreams(); },
			[&]() -> void {HEL.AllocateOutputStreams(); },

			[&]() -> void {CY.AllocateOutputStreams(); },
			[&]() -> void {ADH.AllocateOutputStreams(); },
			[&]() -> void {LM.AllocateOutputStreams(); },
			[&]() -> void {LN.AllocateOutputStreams(); },
			[&]() -> void {MIA.AllocateOutputStreams(); },
			[&]() -> void {MII.AllocateOutputStreams(); },

			[&]() -> void {DF.AllocateOutputStreams(); },
			[&]() -> void {DPA.AllocateOutputStreams(); },
			[&]() -> void {DPN.AllocateOutputStreams(); },
			[&]() -> void {EPS_S.AllocateOutputStreams(); },
			[&]() -> void {EPS_W.AllocateOutputStreams(); }// 14шт */
		};

		std::array<std::function<void(uint32_t, double_t, variables&)>, 14> SolutionOutputter = {
			/*[&](uint32_t Nj, double_t Tj, variables& X) -> void {TCOS.OutputSol(Nj, Tj, X.x); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {TSIN.OutputSol(Nj, Tj, X.x); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {T.OutputSol(Nj, Tj, X.x); }// */

			//[&](uint32_t Nj, double_t Tj, variables& X) -> void {EXP.OutputSol(Nj, Tj, X.x); }

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {NEC.OutputSol(Nj, Tj, X.nec); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {AC.OutputSol(Nj, Tj, X.acu_c); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {HEL.OutputSol(Nj, Tj, X.hel); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {CY.OutputSol(Nj, Tj, X.cy); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {ADH.OutputSol(Nj, Tj, X.adh); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {LM.OutputSol(Nj, Tj, X.lm); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {LN.OutputSol(Nj, Tj, X.ln); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {MIA.OutputSol(Nj, Tj, X.mia); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {MII.OutputSol(Nj, Tj, X.mii); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {DF.OutputSol(Nj, Tj, X.d_F); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {DPA.OutputSol(Nj, Tj, X.dp_N); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {DPN.OutputSol(Nj, Tj, X.dp_A); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {EPS_S.OutputSol(Nj, Tj, X.eps_s); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {EPS_W.OutputSol(Nj, Tj, X.eps_w); } // 14шт */
		};

		std::array<std::function<void(uint32_t, double_t)>, 9> BudgetOutputter = {
			[&](uint32_t Nj, double_t Tj) -> void {NEC.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {AC.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {HEL.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {CY.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {ADH.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {LM.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {LN.OutputBuds(Nj, Tj); },

			//[&](uint32_t Nj, double_t Tj) -> void {miact.OutputBuds(Nj, Tj); },
			//[&](uint32_t Nj, double_t Tj) -> void {minact.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {DF.OutputBuds(Nj, Tj); },
			//[&](uint32_t Nj, double_t Tj) -> void {DPA.OutputBuds(Nj, Tj); },
			//[&](uint32_t Nj, double_t Tj) -> void {DPN.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {EPS_S.OutputBuds(Nj, Tj); }
			//[&](uint32_t Nj, double_t Tj) -> void {eps_w.OutputBuds(Nj, Tj); }
		};

		std::array<std::function<void(uint32_t, variables&)>, 1> RetUploader = {
			[&](uint32_t Nj, variables& X)->void { X.ret.x_1 = EXP.GetRetValue(1.0, Nj); }
			/*[&](uint32_t Nj, variables & X)->void { X.ret.hel_12 = HEL.GetRetValue(12.0, Nj); },
			[&](uint32_t Nj, variables & X)->void { X.ret.adh_12 = ADH.GetRetValue(12.0, Nj); },
			[&](uint32_t Nj, variables & X)->void { X.ret.adh_4 = ADH.GetRetValue(4.0, Nj); }// */
		};

		std::array<std::function<void(uint32_t, variables&)>, 3> RetDataUpdater = {
			[&](uint32_t Nj, variables& X) -> void { HEL.ShiftRets(Nj, X.hel); },
			[&](uint32_t Nj, variables& X) -> void { DF.ShiftRets(Nj, X.d_F); },
			[&](uint32_t Nj, variables& X) -> void { ADH.ShiftRets(Nj, X.adh); }
		};

		std::array<std::function<void(uint32_t, double_t, double_t)>, 3> RetInitialiser = {
			//[&](uint32_t N, double_t t0, double_t gapWidth) -> void {EXP.SetInitialRet(N, t0, gapWidth); }

			/*[&](uint32_t N, double_t t0, double_t gapWidth) -> void {TCOS.SetInitialRet(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {TSIN.SetInitialRet(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {T.SetInitialRet(N, t0, gapWidth); } // */

			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {ADH.SetInitialRet(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {HEL.SetInitialRet(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {DF.SetInitialRet(N, t0, gapWidth); }// */
		};

	};

	template<typename Method>
	class ISolver : public IAggregate {
	public:
		Method ST; // ST is for Straight Task

		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };

		void RetUpload(uint32_t Nj) {

			//ST.X_pred.ret.x_pi2 = EXP.GetRetValue(1.0, Nj);

			//ST.X_pred.ret.x_pi2 = EXP.GetRetValue(pi / 2.0, Nj);

			/*ST.X_pred.ret.x_pi2 = TCOS.GetRetValue(pi * 0.5, Nj);
			ST.X_pred.ret.y_pi2 = TSIN.GetRetValue(pi * 0.5, Nj);
			ST.X_pred.ret.z_pi2 = T.GetRetValue(pi * 0.5, Nj); //*/


			ST.X_pred.ret.hel_12 = HEL.GetRetValue(12.0, Nj);
			ST.X_pred.ret.adh_12 = ADH.GetRetValue(12.0, Nj);
			ST.X_pred.ret.adh_4 = ADH.GetRetValue(4.0, Nj); // */
		}

		virtual void NodeShift() { ST.X_prev = *ST.X_sol; }

		void RetDataUpdate(uint32_t Nj) {
			//EXP.ShiftRets(Nj, ST.X_prev.x);

			/*TCOS.ShiftRets(Nj, ST.X_prev.x);
			TSIN.ShiftRets(Nj, ST.X_prev.y);
			T.ShiftRets(Nj, ST.X_prev.z);//*/

			HEL.ShiftRets(Nj, ST.X_prev.hel);
			DF.ShiftRets(Nj, ST.X_prev.d_F);
			ADH.ShiftRets(Nj, ST.X_prev.adh);// */
		}

		virtual void ApplyMethod() = 0;

		void ErrSaving(double_t Tj) {
			EXP.SaveError(Tj, IErr::GetDif(ST.X_sol->x, EXP.RP.Solution(Tj), 1));

			/*TCOS.SaveError(Tj, IErr::GetDif(ST.X_prev.x, TCOS.RP.Solution(Tj), 1));
			TSIN.SaveError(Tj, IErr::GetDif(ST.X_prev.y, TSIN.RP.Solution(Tj), 1));
			T.SaveError(Tj, IErr::GetDif(ST.X_prev.z, T.RP.Solution(Tj), 1)); // */
		}

		void ErrOutput(double_t Tj) {
			EXP.OutputError(Tj, IErr::GetDif(ST.X_sol->x, EXP.RP.Solution(Tj), 1));

			/*TCOS.OutputError(Tj, IErr::GetDif(ST.X_prev.x, TCOS.RP.Solution(Tj), 1));
			TSIN.OutputError(Tj, IErr::GetDif(ST.X_prev.y, TSIN.RP.Solution(Tj), 1));
			T.OutputError(Tj, IErr::GetDif(ST.X_prev.z, T.RP.Solution(Tj), 1)); // */

		}

	private:
		virtual void ApplyPrepMethod() { return; }

	};

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final
		{
			//ST.X_pred.x = ST.Predictor(ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_prev.z, T.RP);// */

			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}

	};

	class PredCor : public ISolver<Methods::ModEuler>
	{
	public:
		PredCor() { ST.X_sol = &ST.X_cor; }
		void ApplyMethod() final
		{
			ApplyEuler();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			// ST.X_cor.x = ST.Corrector(ST.X_prev.x, EXP.RP);

			/*ST.X_cor.x = ST.Corrector(ST.X_prev.x, TCOS.RP);
			ST.X_cor.y = ST.Corrector(ST.X_prev.y, TSIN.RP);
			ST.X_cor.z = ST.Corrector(ST.X_prev.z, T.RP);//*/

			ST.X_cor.nec = ST.Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.Corrector(ST.X_prev.acu_c, AC.RP);
			ST.X_cor.hel = ST.Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_cor); }// */
		}

	private:
		void ApplyEuler() {
			//ST.X_pred.x = ST.Predictor(ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_prev.z, T.RP); //*/
			
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}
	};

	class RunKut : public ISolver<Methods::RunKut4> {
	public:
		RunKut() { *ST.X_sol = ST.X_pred; }
		void ApplyMethod() final
		{
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_sub.y, ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_sub.z, ST.X_prev.z, T.RP);// */


			/*ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); } // */
		}
	};

	class Gear : public ISolver<Methods::RKGear> {

	public:
		Gear() { *ST.X_sol = ST.X_cor; }
		void ApplyPrepStep(uint32_t& Nj, double_t& Tj) final {
			for (; Nj < ST.X.size(); Nj++)
			{
				// сдвиг на следующий шаг по времени
				Tj = ST.X_pred.tj += ST.H;

				//назначаем соответствующие запаздывания
				RetUpload(Nj);

				ApplyPrepMethod();

				ST.X[Nj] = ST.X_pred;

				//вывод решения
				for (auto const& cur : SolutionOutputter) { cur(Nj, Tj, ST.X_pred); }

				// сдвигаем запаздывания (если они есть, лол)
				RetDataUpdate(Nj);

				// обновляем предшествующий временной ряд для перехода на следующий шаг
				ST.X_prev = ST.X_pred;
			}
		}

		void ApplyMethod() final
		{
			ApplyPred();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.x = ST.Corrector(ST.X[1].x, ST.X[2].x, ST.X[3].x, ST.X_prev.x, EXP.RP);

			/*ST.X_cor.x = ST.Corrector(ST.X[1].x, ST.X[2].x, ST.X[3].x, ST.X_prev.x, TCOS.RP);
			ST.X_cor.y = ST.Corrector(ST.X[1].y, ST.X[2].y, ST.X[3].y, ST.X_prev.y, TSIN.RP);
			ST.X_cor.z = ST.Corrector(ST.X[1].z, ST.X[2].z, ST.X[3].z, ST.X_prev.z, T.RP);// */


			/*ST.X_cor.nec = ST.Corrector(ST.X_1.nec, ST.X_2.nec, ST.X_3.nec, ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.Corrector(ST.X_1.acu_c, ST.X_2.acu_c, ST.X_3.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_cor.hel = ST.Corrector(ST.X_1.hel, ST.X_2.hel, ST.X_3.hel, ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X_1.cy, ST.X_2.cy, ST.X_3.cy, ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X_1.adh, ST.X_2.adh, ST.X_3.adh, ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X_1.lm, ST.X_2.lm, ST.X_3.lm, ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X_1.ln, ST.X_2.ln, ST.X_3.ln, ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X_1.mia, ST.X_2.mia, ST.X_3.mia, ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X_1.mii, ST.X_2.mii, ST.X_3.mii, ST.X_prev.mii, MII.RP);

			ST.X_cor.ret = ST.X_pred.ret; // in case of emergency to update ret-values on X_cor, but what no earth for

			for (auto const& cur : ExpressSubValues) { cur(ST.X_cor); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPred()
		{
			ST.X_pred.x = ST.GPred(ST.X[1].x, ST.X[2].x, ST.X[3].x, ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.EulPred(ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.EulPred(ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.EulPred(ST.X_prev.z, T.RP);// */



			/*ST.X_pred.nec = ST.EulPred(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.EulPred(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.EulPred(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.EulPred(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.EulPred(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.EulPred(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.EulPred(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.EulPred(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.EulPred(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}

		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_sub.y, ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_sub.z, ST.X_prev.z, T.RP);// */

			/*ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}
	};

	class Adams : public ISolver<Methods::Adams> {
	public:
		Adams() { ST.X_sol = &ST.X_pred; }

		void ApplyPrepStep(uint32_t& Nj, double_t& Tj) final {
			for (; Nj < ST.X.size(); Nj++)
			{
				// сдвиг на следующий шаг по времени
				Tj = ST.X_pred.tj += ST.H;

				//назначаем соответствующие запаздывания
				RetUpload(Nj);

				ApplyPrepMethod();

				ST.X[Nj] = ST.X_pred;
				//вывод, когда предиктор или констатация
				for (auto const& cur : SolutionOutputter) { cur(Nj, Tj, ST.X_pred); }

				// сдвигаем запаздывания (если они есть, лол)
				RetDataUpdate(Nj);

				// обновляем предшествующий временной ряд для перехода на следующий шаг
				ST.X_prev = ST.X_pred;
			}
		}

		void ApplyMethod()
		{
			ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.A_Predictor(ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.A_Predictor(ST.X_prev.z, T.RP); // */

			/*ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.A_Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_sub.y, ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_sub.z, ST.X_prev.z, T.RP);// */

			/*ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}

	};

	class ABM : public ISolver<Methods::ABM> {
	public:
		ABM() { ST.X_sol = &ST.X_cor; }

		void ApplyPrepStep(uint32_t& Nj, double_t& Tj) final {
			for (; Nj < ST.X.size(); Nj++)
			{
				// сдвиг на следующий шаг по времени
				Tj = ST.X_pred.tj += ST.H;

				//назначаем соответствующие запаздывания
				RetUpload(Nj);

				ApplyPrepMethod();

				ST.X[Nj] = ST.X_pred;
				//вывод, когда предиктор или констатация
				for (auto const& cur : SolutionOutputter) { cur(Nj, Tj, ST.X_pred); }

				// сдвигаем запаздывания (если они есть, лол)
				RetDataUpdate(Nj);

				// обновляем предшествующий временной ряд для перехода на следующий шаг
				ST.X_prev = ST.X_pred;
			}
		}

		void ApplyMethod() final
		{
			ApplyAdams();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.x = ST.A_Corrector(ST.X_prev.x, EXP.RP);

			/*ST.X_cor.x = ST.A_Corrector(ST.X_prev.x, TCOS.RP);
			ST.X_cor.y = ST.A_Corrector(ST.X_prev.y, TSIN.RP);
			ST.X_cor.z = ST.A_Corrector(ST.X_prev.z, T.RP); // */

			/*ST.X_cor.nec = ST.A_Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.A_Corrector(ST.X_prev.acu_c, AC.RP);
			ST.X_cor.hel = ST.A_Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.A_Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.A_Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.A_Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.A_Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.A_Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.A_Corrector(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_cor); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyAdams()
		{
			ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.A_Predictor(ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.A_Predictor(ST.X_prev.z, T.RP); // */

			/*ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.A_Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}

		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);

			/*ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, TCOS.RP);
			ST.X_pred.y = ST.Predictor(ST.X_sub.y, ST.X_prev.y, TSIN.RP);
			ST.X_pred.z = ST.Predictor(ST.X_sub.z, ST.X_prev.z, T.RP);// */

			/*ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }// */
		}



	};


//______________________________________________________________________________
//==============================================================================


	class INI_data {
	public:
		INI_data() {}

		INI_data(uint32_t N, double_t t_0, double_t tau, uint16_t Ns) :
			m_N(N), m_t_0(t_0), m_tau(tau), m_Ns(Ns) {}

		~INI_data() = default;

		uint32_t m_N = 200; // размер сетки в итерации метода последовательного интегрирования (МПИ)
		double_t m_t_0 = 1; // начальный момент времени
		double_t m_tau = 0; // величина шага МПИ
		uint16_t m_Ns = 2; // число шагов МПИ
	};

	// решает диффур для вывода результата
	void ODE_solver(INI_data& IDT)
	{
		// choosing num-method from StraightTask
		Euler SYS;

		//PhaseTrajOutput PTO("APPROX");
		// здесь определяем параметры численного метода
		SYS.ST.Set(IDT.m_t_0, IDT.m_N, IDT.m_tau, IDT.m_Ns);

		// Инициализируем накопитель запаздывающих аргументов
		for (auto const& cur : SYS.RetInitialiser) { cur(SYS.ST.N, SYS.ST.t_0, SYS.ST.gap_width); }

		// заполняем массивы из уже существующих решений на случай,
		// если надо зафиксировать систему вдоль некоторого решения компоненты системы
		//for (auto const& cur : SYS.DataCollector) { cur(); }

		// и определяем начальные условия из предзаписанного решения
		//SYS.ST.t_0 = SYS.ST.X_init.tim = 0.5;
		//for (auto const& cur : SYS.SolDataInitialiser) { cur(0, 0, SYS.ST.X_init); }


		// и определяем начальные условия
		SYS.ST.X_init.tj = SYS.ST.t_0;
		for (auto const& cur : SYS.IniDataInitialiser) { cur(SYS.ST.X_init); }
		for (auto const& cur : SYS.IniRetInitialiser) { cur(SYS.ST.X_init); }
		//for (auto const& cur : SYS.ExpressSubValues) { cur(SYS.ST.X_init); } // */

		// открываем потоки для вывода численного решения 
		for (auto const& cur : SYS.OutStreamAllocator) { cur(); }

		// выводим данные начальных условий во внешний файл с решением
		for (auto const& cur : SYS.SolutionOutputter) { cur(0, SYS.ST.X_init.tj, SYS.ST.X_init); }


		uint16_t current_gap = 0;
		double_t Tj;

		while (current_gap < SYS.ST.full_amount_of_gaps)
		{

			std::cout << SYS.ST.full_amount_of_gaps - current_gap << " "; // вывод на экран "кол-во промежутков, которые надо просчитать"
			SYS.ST.X_prev = SYS.ST.X_init; // даём начальное условие на предшествующий вектор
			Tj = SYS.ST.X_pred.tj = SYS.ST.X_prev.tj;
			uint32_t Nj = 1;

			if (current_gap == 0) SYS.ApplyPrepStep(Nj, Tj); // 1st approximations for multistep-methods

			// цикл на следующие 24 часа
			for (; Nj <= SYS.ST.N; Nj++)
			{
				// сдвиг на следующий шаг по времени
				Tj = SYS.ST.X_pred.tj += SYS.ST.H;

				//назначаем соответствующие запаздывания
				//for (auto const& cur : SYS.RetUploader) { cur(Nj, SYS.ST.X_pred); }
				SYS.RetUpload(Nj);

				// контролируем промежуток интерполяции
				//ST.X_pred.CheckShiftInterpGap(STpar);

				// предзаписанное решение
				//for (auto const& cur : SYS.SolDataInitialiser) { cur(current_day, Nj, SYS.ST.X_pred); }			

				//SYS.SolDataInitialiser[3](current_day, Nj, SYS.ST.X_pred); // цитокины

				SYS.ApplyMethod();

				//вывод, когда предиктор или констатация
				for (auto const& cur : SYS.SolutionOutputter) { cur(Nj, Tj, *SYS.ST.X_sol); }

				// вывод фазового портрета решения
				//PTO.OutputPhaseTraj(Nj, Tj, { (*SYS.ST.X_sol).x, (*SYS.ST.X_sol).y, (*SYS.ST.X_sol).z });


				//SYS.ErrSaving(Tj); // сохранение локальной ошибки в массив
				//SYS.ErrOutput(Tj); // вывод локальной ошибки во внешний файл

				//for (auto const& cur : SYS.BudgetOutputter) { cur(Nj, Tj); }

				// сдвигаем запаздывания (если они есть, лол)
				SYS.RetDataUpdate(Nj);

				// обновляем предшествующий временной ряд для перехода на следующий шаг
				SYS.NodeShift();


			} // конец цикла рассчётов на текущий день for(Nj: 1->N)

			current_gap++; std::cout << " ; ";

			// Проверка, «а надо ли решать дальше?» 
			if (current_gap == SYS.ST.full_amount_of_gaps)
			{
				std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
				//TRIG = false; // условие на выход из цикла while
				// else не нужен, ведь если !ID, то внешнее условие while его поймает
			}
			else
			{
				SYS.ST.X_init = *SYS.ST.X_sol;
				// SYS.ST.X_init = SYS.ST.X_pred;
				Nj = 1;
			}
		}


		//SYS.TCOS.GetOrdOfMaxError();
		//SYS.TSIN.GetOrdOfMaxError();
		//SYS.T.GetOrdOfMaxError();

		//std::cout << std::endl;
		/*std::cout << "( " << SYS.TCOS.Error[0] <<"\t" <<  SYS.TSIN.Error[0] << "\t" << SYS.T.Error[0] <<  " ) " << "\n"
			<< SYS.TCOS.Error[1] + SYS.TSIN.Error[1] + SYS.T.Error[1] << std::endl;// */
			//std::cout << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][0] << "  " << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][1] << std::endl;
	}

	void SolutionOutput()
	{
		Test::ThreeDim::Ox_ret X;
		Test::ThreeDim::Oy_ret Y;
		Test::ThreeDim::Oz_ret Z;

		PhaseTrajOutput PTO("SOL");

		uint16_t GA = 10;
		uint32_t N = 2000;
		double_t H = pi / N;

		//std::ofstream sol(output_dir + "ST/SOL/PTraj/" + "SOL" + ".txt");
		for (uint32_t i = 0; i < GA * N; i++) {
			PTO.OutputPhaseTraj(i, pi + i * H, { X.Solution(pi + i * H), Y.Solution(pi + i * H), Z.Solution(pi + i * H) });
		}
		//sol.close();
		return;
	}


}