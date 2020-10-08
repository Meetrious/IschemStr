/* This header serves as pre-surface shell for 1st model pipe realisation.
 It contains the class IAggregate that lists members of ODE system
 and the class ISolver publicly inherited in class <Method>, where presented the list of all functions required in algorythm stated in Solver.h // */

// This header ends with #include<Solver.h> as for other header located on top of ST/include directory

#pragma once
#include <functional> 
#include <array>

#include <models/1st/settings.h>

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
		Neurons::Apoptosis::Started::M_ODE AS;
		Neurons::Apoptosis::Ended::M_ODE AE;
		Neurons::IntactCells::M_ODE HEL;

		Cytokines::Pro_Inflam::M_ODE CY;
		Adhesion::M_ODE ADH;
		LeuMacrophags::M_ODE LM;
		LeuNeutrophils::M_ODE LN;

		Microglia::Active::M_ODE MIA;
		Microglia::Inactive::M_ODE MII;


		// subbordinate values relevant for a system 
		Cytokines::M_Sub PSY;
		ToxDamage::Full::M_Sub DF;
		ToxDamage::Initial::M_Sub D_INI;

		Phagocytosis::Strong::M_Sub EPS_S;



		std::array<std::function<void()>, 15> DataCollector = {
			[&]()-> void {NEC.CollectSolData(); },
			[&]()-> void {AE.CollectSolData(); },
			[&]()-> void {HEL.CollectSolData(); },

			[&]()-> void {CY.CollectSolData(); },
			[&]()-> void {ADH.CollectSolData(); },
			[&]()-> void {LM.CollectSolData(); },
			[&]()-> void {LN.CollectSolData(); },
			[&]()-> void {MIA.CollectSolData(); },
			[&]()-> void {MII.CollectSolData(); },

			[&]()-> void {DF.CollectSolData(); },
			[&]()-> void {D_INI.CollectSolData(); },

			[&]()-> void {EPS_S.CollectSolData(); },
			
		};

		std::array<std::function<void(variables&)>, 10> IniDataInitialiser_FO = {
			[&](variables& X) -> void {NEC.GetInitialDataFromOutside(X.tj, X.nec); },
			[&](variables& X) -> void {AS.GetInitialDataFromOutside(X.tj,X.ap_s); },
			[&](variables& X) -> void {AE.GetInitialDataFromOutside(X.tj,X.ap_e); },
			[&](variables& X) -> void {HEL.GetInitialDataFromOutside(X.tj,X.hel); },
			[&](variables& X) -> void {CY.GetInitialDataFromOutside(X.tj,X.cy); },
			[&](variables& X) -> void {ADH.GetInitialDataFromOutside(X.tj,X.adh); },
			[&](variables& X) -> void {LM.GetInitialDataFromOutside(X.tj,X.lm); },
			[&](variables& X) -> void {LN.GetInitialDataFromOutside(X.tj,X.ln); },
			[&](variables& X) -> void {MIA.GetInitialDataFromOutside(X.tj,X.mia); },
			[&](variables& X) -> void {MII.GetInitialDataFromOutside(X.tj,X.mii); }// 9шт */
		};
		std::array<std::function<void(variables&)>, 10> IniDataInitialiser = {
			[&](variables& X) -> void {NEC.GetInitialData(X.tj, X.nec); },
			[&](variables& X) -> void {AS.GetInitialData(X.tj,X.ap_s); },
			[&](variables& X) -> void {AE.GetInitialData(X.tj,X.ap_e); },
			[&](variables& X) -> void {HEL.GetInitialData(X.tj,X.hel); },
			[&](variables& X) -> void {CY.GetInitialData(X.tj,X.cy); },
			[&](variables& X) -> void {ADH.GetInitialData(X.tj,X.adh); },
			[&](variables& X) -> void {LM.GetInitialData(X.tj,X.lm); },
			[&](variables& X) -> void {LN.GetInitialData(X.tj,X.ln); },
			[&](variables& X) -> void {MIA.GetInitialData(X.tj,X.mia); },
			[&](variables& X) -> void {MII.GetInitialData(X.tj,X.mii); }// 9шт */
		};

		std::array<std::function<void(variables&)>, 4> IniRetInitialiser = {
			[&](variables& X) -> void {X.ret.hel_12 = HEL.GetRetValue(12.0, 0); },
			[&](variables& X) -> void {X.ret.adh_24 = ADH.GetRetValue(24.0, 0); },
			[&](variables& X) -> void {X.ret.d_12 = D_INI.GetRetValue(12.0, 0); },
			[&](variables& X) -> void {X.ret.adh_12 = ADH.GetRetValue(12.0, 0); } // */

		};

		std::array<std::function<void(uint16_t, uint32_t, uint32_t, variables&)>, 13> SolDataGetter = {
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.nec = NEC.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.ap_s = AS.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.ap_e = AE.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.hel = HEL.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.cy = CY.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.adh = ADH.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.lm = LM.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.ln = LN.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.mia = MIA.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.mii = MII.GetSolData(day, Nj, N); },

			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.d_F = DF.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.d_ini = D_INI.GetSolData(day, Nj, N); },
			[&](uint16_t day, uint32_t Nj, uint32_t N, variables& X) -> void {X.eps_s = EPS_S.GetSolData(day, Nj, N); }
			
		};

		std::array<std::function<void(variables&)>, 4> SubValuesExpressor = {
			[&](variables& X) -> void { X.d_F = DF.RP.Expression(X); },
			[&](variables& X) -> void { X.d_ini = D_INI.RP.Expression(X); },
			[&](variables& X) -> void { X.psy = PSY.RP.Expression(X); },
			[&](variables& X) -> void { X.eps_s = EPS_S.RP.Expression(X); },
		};

		std::array<std::function<void()>, 13> OutStreamAllocator = {
			[&]() -> void {NEC.AllocateOutputStreams(); },
			[&]() -> void {AS.AllocateOutputStreams(); },
			[&]() -> void {AE.AllocateOutputStreams(); },
			[&]() -> void {HEL.AllocateOutputStreams(); },

			[&]() -> void {CY.AllocateOutputStreams(); },
			[&]() -> void {ADH.AllocateOutputStreams(); },
			[&]() -> void {LM.AllocateOutputStreams(); },
			[&]() -> void {LN.AllocateOutputStreams(); },
			[&]() -> void {MIA.AllocateOutputStreams(); },
			[&]() -> void {MII.AllocateOutputStreams(); },

			[&]() -> void {DF.AllocateOutputStreams(); },
			[&]() -> void {D_INI.AllocateOutputStreams(); },
			[&]() -> void {EPS_S.AllocateOutputStreams(); } // 13 */
		};

		//why wouldn't you put this output realisation as a method in variables class?
		std::array<std::function<void(uint32_t, double_t, variables&)>, 13> SolutionOutputter = {

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {NEC.OutputSol(Nj, Tj, X.nec); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {AS.OutputSol(Nj, Tj, X.ap_s); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {AE.OutputSol(Nj, Tj, X.ap_e); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {HEL.OutputSol(Nj, Tj, X.hel); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {CY.OutputSol(Nj, Tj, X.cy); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {ADH.OutputSol(Nj, Tj, X.adh); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {LM.OutputSol(Nj, Tj, X.lm); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {LN.OutputSol(Nj, Tj, X.ln); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {MIA.OutputSol(Nj, Tj, X.mia); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {MII.OutputSol(Nj, Tj, X.mii); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {DF.OutputSol(Nj, Tj, X.d_F); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {D_INI.OutputSol(Nj, Tj, X.d_ini); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {EPS_S.OutputSol(Nj, Tj, X.eps_s); },
			
		};

		std::array<std::function<void(uint32_t, double_t)>, 12> BudgetOutputter = {
			[&](uint32_t Nj, double_t Tj) -> void {NEC.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {AS.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {AE.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {HEL.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {CY.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {ADH.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {LM.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {LN.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {MIA.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {MII.OutputBuds(Nj, Tj); },

			[&](uint32_t Nj, double_t Tj) -> void {DF.OutputBuds(Nj, Tj); },
			[&](uint32_t Nj, double_t Tj) -> void {EPS_S.OutputBuds(Nj, Tj); }
		};

		std::array<std::function<void(uint32_t, double_t, double_t)>, 3> RetInitialiser = {
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {ADH.StateRetArray(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {HEL.StateRetArray(N, t0, gapWidth); },
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {D_INI.StateRetArray(N, t0, gapWidth); }// */
		};

	};

	template<typename Method>
	class ISolver : public IAggregate {
	public:
		Method ST; // ST is for Straight Task

		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };

		void InitialiseIniData() {
			vector<double_t> T0s;
			for (auto const& cur : IniDataInitialiser) {
				cur(ST.X_init); T0s.emplace_back(ST.X_init.tj);
			}
			double_t mean = 0;
			for (size_t i = 0; i < T0s.size(); i++) {
				mean += T0s[0];
			}
			mean = mean / T0s.size(); 
			
			ST.X_init.tj = ST.t_0 = mean;
		
			bool trg = false;

			// making sure, that each initial datum is set in the same t0 
			for (size_t i = 0; i < T0s.size(); i++) {
				if (mean == T0s[i])
					continue;
				else trg = true;
			}
			if (trg)
			{
				std::cerr << "\n WARNING: initial data for the ODE system is not consistent:"
					<< "\n initial time moments are not the same: \n {_";
				for (auto const& cur : T0s) std::cerr << cur << '_';
				std::cerr << "} \n t_0 will be defined as mean."
					<< "\n Do you wish to proceed? :\n  0. NO; \n 1. YES;" << std::endl;
				std::cin >> trg;
				if (!trg) throw(" Initial data is inacceptable. ");
				else std::cerr << "Approved \n Proceeding. \n";
			}
			for (auto const& cur : SubValuesExpressor) { cur(ST.X_init); }// */
			return;
		}

		bool is_SYS_deflecting() {
			bool if_ret = NEC.RP.ret_is + AS.RP.ret_is + AE.RP.ret_is + HEL.RP.ret_is +
				CY.RP.ret_is + ADH.RP.ret_is + LM.RP.ret_is + LN.RP.ret_is +
				MIA.RP.ret_is + MII.RP.ret_is;
			return if_ret;
		}

		void PrepairTheTask()
		{
			try {
				InitialiseIniData();

				// collecting presolved solution to freeze the system relatively given behaviour
				// for (auto const& cur : DataCollector) { cur(); } // full
				//DataCollector[10]();

				// setting initial data from presolved_solution_data
				// for (auto const& cur : SolDataGetter) { cur(0, 0, ST.X_init); } // full
				//SolDataGetter[11](0, 0, ST.N, ST.X_init);

				if (is_SYS_deflecting())
				{
					// initialise ret-value storage
					for (auto const& cur : RetInitialiser) { cur(ST.N, ST.t_0, ST.gap_width); }

					// setting first ret-values required on the first step of calculation
					for (auto const& cur : IniRetInitialiser) { cur(ST.X_init); }
				}

			}
			catch (const char* exception) {
				std::cerr << "WARNING:" << exception << "\n Terminating.";
				throw(exception);
			}

		}

		void PrepairTheOutput()
		{
			try{

			// streams for output
			for (auto const& cur : OutStreamAllocator) { cur(); }

			// outputting initial solution data
			for (auto const& cur : SolutionOutputter) { cur(0, ST.X_init.tj, ST.X_init); }
			
			}
			catch (const char* exception){
				std::cerr << "WARNING:" << exception << "\n Terminating.";
				throw(exception);
			}
		}

		void RetUpload(uint32_t Nj) {
			
			ST.X_pred.ret.hel_12 = HEL.GetRetValue(12.0, Nj);
			ST.X_pred.ret.adh_12 = ADH.GetRetValue(12.0, Nj);
			ST.X_pred.ret.adh_24 = ADH.GetRetValue(24.0, Nj);
			ST.X_pred.ret.d_12 = D_INI.GetRetValue(12.0, Nj);
			// */
		}

		virtual void NodeShift() { ST.X_prev = *ST.X_sol; }

		void RetDataUpdate(uint32_t Nj) {
			HEL.ShiftRets(Nj, ST.X_prev.hel);
			D_INI.ShiftRets(Nj, ST.X_prev.d_ini);
			ADH.ShiftRets(Nj, ST.X_prev.adh);// */
		}

		virtual void ApplyMethod() = 0;

	private:
		virtual void ApplyPrepMethod() { return; }

	};

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
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

			ST.X_cor.nec = ST.Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.ap_s = ST.Corrector(ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.Corrector(ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_cor); }// */
		}

	private:
		void ApplyEuler() {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
		}
	};

	class RunKut : public ISolver<Methods::RunKut4> {
	public:
		RunKut() { *ST.X_sol = ST.X_pred; }
		void ApplyMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); } // */
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

			ST.X_cor.nec = ST.Corrector(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_cor.ap_s = ST.Corrector(ST.X[1].ap_s, ST.X[2].ap_s, ST.X[3].ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.Corrector(ST.X[1].ap_e, ST.X[2].ap_e, ST.X[3].ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.Corrector(ST.X[1].hel, ST.X[2].hel, ST.X[3].hel, ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X[1].cy, ST.X[2].cy, ST.X[3].cy, ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X[1].adh, ST.X[2].adh, ST.X[3].adh, ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X[1].lm, ST.X[2].lm, ST.X[3].lm, ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X[1].ln, ST.X[2].ln, ST.X[3].ln, ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X[1].mia, ST.X[2].mia, ST.X[3].mia, ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X[1].mii, ST.X[2].mii, ST.X[3].mii, ST.X_prev.mii, MII.RP);

			ST.X_cor.ret = ST.X_pred.ret; // in case of emergency to update ret-values on X_cor, but what no earth for

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_cor); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPred() {
			ST.X_pred.nec = ST.GPred(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.GPred(ST.X[1].ap_s, ST.X[2].ap_s, ST.X[3].ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.GPred(ST.X[1].ap_e, ST.X[2].ap_e, ST.X[3].ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.GPred(ST.X[1].hel, ST.X[2].hel, ST.X[3].hel, ST.X_prev.hel, HEL.RP);
			ST.X_pred.cy = ST.GPred(ST.X[1].cy, ST.X[2].cy, ST.X[3].cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.GPred(ST.X[1].adh, ST.X[2].adh, ST.X[3].adh, ST.X_prev.adh, ADH.RP);
			ST.X_pred.lm = ST.GPred(ST.X[1].lm, ST.X[2].lm, ST.X[3].lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.GPred(ST.X[1].ln, ST.X[2].ln, ST.X[3].ln, ST.X_prev.ln, LN.RP);
			ST.X_pred.mia = ST.GPred(ST.X[1].mia, ST.X[2].mia, ST.X[3].mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.GPred(ST.X[1].mii, ST.X[2].mii, ST.X[3].mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
		}

		void ApplyPrepMethod() final {
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
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

		void ApplyMethod() {
			ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.A_Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.A_Predictor(ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
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

			ST.X_cor.nec = ST.A_Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.ap_s = ST.A_Corrector(ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.A_Corrector(ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.A_Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.A_Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.A_Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.A_Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.A_Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.A_Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.A_Corrector(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_cor); }// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyAdams() {
			ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.A_Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.A_Predictor(ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
		}

		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			for (auto const& cur : SubValuesExpressor) { cur(ST.X_pred); }// */
		}
	};
}
#include <base/Solver.h>
