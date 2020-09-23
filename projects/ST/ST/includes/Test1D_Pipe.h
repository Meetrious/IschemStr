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
		Test::OneDim::M_ODE EXP;


		std::array<std::function<void(variables&)>, 1> IniDataInitialiser = {
			[&](variables& X) -> void {EXP.GetInitialData(X.tj, X.x); }
		};

		std::array<std::function<void(variables&)>, 1> IniRetInitialiser = {
			[&](variables& X) -> void {X.ret.x_pi2 = EXP.GetRetValue(pi / 2, 0); }
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
			[&](uint32_t N, double_t t0, double_t gapWidth) -> void {EXP.SetInitialRet(N, t0, gapWidth); }
		};

	};

	template<typename Method>
	class ISolver : public IAggregate {
	public:
		Method ST; // ST is for Straight Task

		virtual void ApplyPrepStep(uint32_t& Nj, double_t& Tj) { return; };

		void IniDataInitialise(){
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
			for (size_t i = 0; i < T0s.size(); i++) {
				if (mean == T0s[i])
					continue;
				else trg = true;
			}
			if (trg) std::cout << "\n WARNING: initial data for the equation system is not consistent:"
				<< "\n initial time moments are not the same. "
				<< "\n t_0 will be defined as mean."
				<< "\n Do you wish to proceed? :\n  0. NO; \n 1. YES;";
			std::cin >> trg;
			// clarify the exception!!!!!!
		}

		void RetUpload(uint32_t Nj) {
			ST.X_pred.ret.x_pi2 = EXP.GetRetValue(pi / 2.0, Nj);
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

	// решает диффур для вывода результата
	void ODE_solver()
	{
		// choosing num-method from StraightTask
		Euler SYS;

		// здесь определяем параметры численного метода
		SYS.ST.Set(200, pi/2.0, 5);

		// initialise ini-data
		SYS.IniDataInitialise();


		// Инициализируем накопитель запаздывающих аргументов
		for (auto const& cur : SYS.RetInitialiser) { cur(SYS.ST.N, SYS.ST.t_0, SYS.ST.gap_width); }

		for (auto const& cur : SYS.IniRetInitialiser) { cur(SYS.ST.X_init); }
		
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
				SYS.RetUpload(Nj);

				SYS.ApplyMethod();

				//вывод, когда предиктор или констатация
				for (auto const& cur : SYS.SolutionOutputter) { cur(Nj, Tj, *SYS.ST.X_sol); }

				//SYS.ErrSaving(Tj); // сохранение локальной ошибки в массив
				//SYS.ErrOutput(Tj); // вывод локальной ошибки во внешний файл

				// updating in array from X_prev 
				SYS.RetDataUpdate(Nj);

				// updating X_prev for the next step
				SYS.NodeShift();


			} // конец цикла рассчётов на текущий день for(Nj: 1->N)

			current_gap++; std::cout << " ; ";

			// Проверка, «а надо ли решать дальше?» 
			if (current_gap == SYS.ST.full_amount_of_gaps)
				std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			else{
				SYS.ST.X_init = *SYS.ST.X_sol;
				Nj = 1;
			}
		}

		//std::cout << std::endl;
		/*std::cout << "( " << SYS.TCOS.Error[0] <<"\t" <<  SYS.TSIN.Error[0] << "\t" << SYS.T.Error[0] <<  " ) " << "\n"
			<< SYS.TCOS.Error[1] + SYS.TSIN.Error[1] + SYS.T.Error[1] << std::endl;// */
			//std::cout << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][0] << "  " << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][1] << std::endl;
	}
}