#pragma once
#include <functional>
#include <array>
#include <tests/settings.h>

namespace StraightTask
{
	class IAggregate
	{
	public:

		// list of equations in the processed system
		Test::OneDim::M_ODE EXP;


		std::array<std::function<void(variables&)>, 1> IniDataInitialiser = {
			[&](variables& X) -> void {EXP.GetInitialData(X.tj, X.x); }
		};

		std::array<std::function<void(variables&)>, 1> IniRetInitialiser = {
			[&](variables& X) -> void {X.ret.x_1 = EXP.GetRetValue(1.0, 0); }
		};


		std::array<std::function<void()>, 1> OutStreamAllocator = {
			[&]() -> void {EXP.AllocateOutputStreams(); }
		};

		std::array<std::function<void(uint32_t, double_t, variables&)>, 1> SolutionOutputter = {
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {EXP.OutputSol(Nj, Tj, X.x); }
		};


		std::array<std::function<void(uint32_t, variables&)>, 1> RetUploader = {
			[&](uint32_t Nj, variables& X)->void { X.ret.x_1 = EXP.GetRetValue(1.0, Nj); }
		};


		std::array<std::function<void(uint32_t, double_t, double_t)>, 1> RetInitialiser = {
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {EXP.StateRetArray(N, t0, gapWidth); }
		};

	};

	template<typename Method>
	class ISolver : public IAggregate {
	public:
		Method ST; // ST is for Straight Task

		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };

		bool is_SYS_deflecting() { return EXP.RP.ret_is; }

		void InitialiseIniData() throw(const char*) {
			vector<double_t> T0s;
			for (auto const& cur : IniDataInitialiser){
				cur(ST.X_init); T0s.emplace_back(ST.X_init.tj);
			}
			double_t mean = 0;
			for (size_t i = 0; i < T0s.size(); i++){
				mean += T0s[0];
			}
			mean = mean / T0s.size(); ST.X_init.tj = mean;

			bool trg = false;

			// making sure, that each initial data is set in the same t0 
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
			return;
		}

		void PrepairTheTask()
		{
			try {
				SYS.InitialiseIniData();

				// collecting presolved solution to freeze the system relatively given behaviour
				// for (auto const& cur : DataCollector) { cur(); } // full

				// setting initial data from presolved_solution_data
				// for (auto const& cur : SolDataGetter) { cur(0, 0, SYS.ST.X_init); } // full

				if (is_SYS_deflecting())
				{
					// initialise ret-value storage
					for (auto const& cur : RetInitialiser) { cur(ST.N, ST.t_0, ST.gap_width); }

					// setting first ret-values required on the first step of calculation
					for (auto const& cur : IniRetInitialiser) { cur(ST.X_init); }
				}

				// streams for output
				for (auto const& cur : OutStreamAllocator) { cur(); }

				// outputting initial solution data
				for (auto const& cur : SolutionOutputter) { cur(0, ST.X_init.tj, ST.X_init); }
			}
			catch (const char* exception) {
				std::cerr << "WARNING:" << exception << "\n Terminating.";
				throw(exception);
			}

		}


		

		void RetUpload(uint32_t Nj) {
			ST.X_pred.ret.x_1 = EXP.GetRetValue(ST.gap_width, Nj);
		}

		virtual void NodeShift() { ST.X_prev = *ST.X_sol; }

		void RetDataUpdate(uint32_t Nj) {
			EXP.ShiftRets(Nj, ST.X_prev.x);
		}

		virtual void ApplyMethod() = 0;

		void ErrSaving(double_t Tj) {
			EXP.SaveError(Tj, IErr::GetDif(ST.X_sol->x, EXP.RP.Solution(Tj), 1));
		}

		void ErrOutput(double_t Tj) {
			EXP.OutputError(Tj, IErr::GetDif(ST.X_sol->x, EXP.RP.Solution(Tj), 1));
		}

	private:
		virtual void ApplyPrepMethod() { return; }
	};

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final {
			ST.X_pred.x = ST.Predictor(ST.X_prev.x, EXP.RP);
		}

	};

	class PredCor : public ISolver<Methods::ModEuler>
	{
	public:
		PredCor() { ST.X_sol = &ST.X_cor; }
		void ApplyMethod() final {
			ApplyEuler();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.x = ST.Corrector(ST.X_prev.x, EXP.RP);

		}

	private:
		void ApplyEuler() {
			ST.X_pred.x = ST.Predictor(ST.X_prev.x, EXP.RP);
		}
	};

	class RunKut : public ISolver<Methods::RunKut4> {
	public:
		RunKut() { *ST.X_sol = ST.X_pred; }
		void ApplyMethod() final
		{
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);
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
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPred(){
			ST.X_pred.x = ST.GPred(ST.X[1].x, ST.X[2].x, ST.X[3].x, ST.X_prev.x, EXP.RP);
		}

		void ApplyPrepMethod() final {
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);
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
			ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, EXP.RP);
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPrepMethod() final{
			ST.X_sub = ST.X_prev;

			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);
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
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyAdams() {
			ST.X_pred.x = ST.A_Predictor(ST.X_prev.x, EXP.RP);
		}

		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;
			ST.X_pred.x = ST.Predictor(ST.X_sub.x, ST.X_prev.x, EXP.RP);
		}

	};


//______________________________________________________________________________
//==============================================================================

}

#include <base/Solver.h>