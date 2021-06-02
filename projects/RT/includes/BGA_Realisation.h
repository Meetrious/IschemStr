#pragma once
#include <ReverseTaskBase.h>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

namespace StraightTask
{


	class Aggregator
	{
	public:

		Aggregator(uint32_t grid_spliting_amount, uint16_t full_amount_of_days) :
			ST{ grid_spliting_amount , full_amount_of_days }
		{}

		Methods::PredCor ST;
		// перечисляем уравнения, которые будем решать
		SODE_Unit::NecUnit necr;
		SODE_Unit::AcuUnit acu;
		SODE_Unit::HelUnit healt;

		SODE_Unit::CytoUnit cyto;
		SODE_Unit::AdhUnit adhes;
		SODE_Unit::LmUnit lmacr;
		SODE_Unit::LnUnit lneut;
		SODE_Unit::MiaUnit miact;
		SODE_Unit::MiiUnit minact;


		// перечисляем вспомогательные величины
		SODE_Unit::TDF_Unit df;
		SODE_Unit::TDpN_Unit dpN;
		SODE_Unit::TDpA_Unit dpA;

		SODE_Unit::PhS_Unit eps_s;
		SODE_Unit::PhW_Unit eps_w;

		std::array<std::function<void()>, 14> DataCollector = {
			[&]()-> void {necr.CollectSolData(); },
			[&]()-> void {acu.CollectSolData(); },
			[&]()-> void {healt.CollectSolData(); },

			[&]()-> void {cyto.CollectSolData(); },
			[&]()-> void {adhes.CollectSolData(); },
			[&]()-> void {lmacr.CollectSolData(); },
			[&]()-> void {lneut.CollectSolData(); },
			[&]()-> void {miact.CollectSolData(); },
			[&]()-> void {minact.CollectSolData(); },

			[&]()-> void {df.CollectSolData(); },
			[&]()-> void {dpN.CollectSolData(); },
			[&]()-> void {dpA.CollectSolData(); },

			[&]()-> void {eps_s.CollectSolData(); },
			[&]()-> void {eps_w.CollectSolData(); }
		};

		std::array<std::function<void(variables&)>, 9> IniDataInitialiser = {
			[&](variables& X) -> void {X.nec = necr.GetInitialData(); },
			[&](variables& X) -> void {X.acu_c = acu.GetInitialData(); },
			[&](variables& X) -> void {X.hel = healt.GetInitialData(); },

			[&](variables& X) -> void {X.cy = cyto.GetInitialData(); },
			[&](variables& X) -> void {X.adh = adhes.GetInitialData(); },
			[&](variables& X) -> void {X.lm = lmacr.GetInitialData(); },
			[&](variables& X) -> void {X.ln = lneut.GetInitialData(); },
			[&](variables& X) -> void {X.mia = miact.GetInitialData(); },
			[&](variables& X) -> void {X.mii = minact.GetInitialData(); },
		};

		std::array<std::function<void(uint16_t, uint32_t, variables&)>, 14>	SolDataInitialiser = {
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.nec = necr.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.acu_c = acu.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.hel = healt.GetSolData(day, Nj, ST.N); },

			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.cy = cyto.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.adh = adhes.GetSolData(day, Nj, ST.N); },

			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.lm = lmacr.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.ln = lneut.GetSolData(day, Nj, ST.N); },

			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.mia = miact.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.mii = minact.GetSolData(day, Nj, ST.N); },

			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.d_F = df.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.dp_N = dpN.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.dp_A = dpA.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.eps_s = eps_s.GetSolData(day, Nj, ST.N); },
			[&](uint16_t day, uint32_t Nj, variables& X) -> void {X.eps_w = eps_w.GetSolData(day, Nj, ST.N); }
		};

		std::array<std::function<void(variables&)>, 5> ExpressSubValues = {
			[&](variables& X) -> void { X.d_F = df.RP.Expression(X); },
			[&](variables& X) -> void { X.dp_N = dpN.RP.Expression(X); },
			[&](variables& X) -> void { X.dp_A = dpA.RP.Expression(X); },
			[&](variables& X) -> void { X.eps_s = eps_s.RP.Expression(X); },
			[&](variables& X) -> void { X.eps_w = eps_w.RP.Expression(X); }
		};

		std::array<std::function<void()>, 14> OutStreamAllocator = {
			[&]() -> void {necr.AllocateOutputStreams(); },
			[&]() -> void {acu.AllocateOutputStreams(); },
			[&]() -> void {healt.AllocateOutputStreams(); },

			[&]() -> void {cyto.AllocateOutputStreams(); },
			[&]() -> void {adhes.AllocateOutputStreams(); },
			[&]() -> void {lmacr.AllocateOutputStreams(); },
			[&]() -> void {lneut.AllocateOutputStreams(); },
			[&]() -> void {miact.AllocateOutputStreams(); },
			[&]() -> void {minact.AllocateOutputStreams(); },

			[&]() -> void {df.AllocateOutputStreams(); },
			[&]() -> void {dpA.AllocateOutputStreams(); },
			[&]() -> void {dpN.AllocateOutputStreams(); },
			[&]() -> void {eps_s.AllocateOutputStreams(); },
			[&]() -> void {eps_w.AllocateOutputStreams(); }
		};

		std::array<std::function<void(uint32_t, double_t, variables&)>, 14> SolutionOutputter = {
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {necr.OutputSol(Nj, Tj, X.nec); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {acu.OutputSol(Nj, Tj, X.acu_c); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {healt.OutputSol(Nj, Tj, X.hel); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {cyto.OutputSol(Nj, Tj, X.cy); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {adhes.OutputSol(Nj, Tj, X.adh); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {lmacr.OutputSol(Nj, Tj, X.lm); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {lneut.OutputSol(Nj, Tj, X.ln); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {miact.OutputSol(Nj, Tj, X.mia); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {minact.OutputSol(Nj, Tj, X.mii); },

			[&](uint32_t Nj, double_t Tj, variables& X) -> void {df.OutputSol(Nj, Tj, X.d_F); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {dpA.OutputSol(Nj, Tj, X.dp_N); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {dpN.OutputSol(Nj, Tj, X.dp_A); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {eps_s.OutputSol(Nj, Tj, X.eps_s); },
			[&](uint32_t Nj, double_t Tj, variables& X) -> void {eps_w.OutputSol(Nj, Tj, X.eps_w); }
		};

		std::array<std::function<void(uint32_t, double_t)>, 9> BudgetOutputter = {
		[&](uint32_t Nj, double_t Tj) -> void {necr.OutputBuds(Nj, Tj); },
		[&](uint32_t Nj, double_t Tj) -> void {acu.OutputBuds(Nj, Tj); },
		[&](uint32_t Nj, double_t Tj) -> void {healt.OutputBuds(Nj, Tj); },

		[&](uint32_t Nj, double_t Tj) -> void {cyto.OutputBuds(Nj, Tj); },
		[&](uint32_t Nj, double_t Tj) -> void {adhes.OutputBuds(Nj, Tj); },

		[&](uint32_t Nj, double_t Tj) -> void {lmacr.OutputBuds(Nj, Tj); },
		[&](uint32_t Nj, double_t Tj) -> void {lneut.OutputBuds(Nj, Tj); },

		//[&](uint32_t Nj, double_t Tj) -> void {miact.OutputBuds(Nj, Tj); },
		//[&](uint32_t Nj, double_t Tj) -> void {minact.OutputBuds(Nj, Tj); },

		[&](uint32_t Nj, double_t Tj) -> void {df.OutputBuds(Nj, Tj); },
		//[&](uint32_t Nj, double_t Tj) -> void {dpA.OutputBuds(Nj, Tj); },
		//[&](uint32_t Nj, double_t Tj) -> void {dpN.OutputBuds(Nj, Tj); },
		[&](uint32_t Nj, double_t Tj) -> void {eps_s.OutputBuds(Nj, Tj); }
		//[&](uint32_t Nj, double_t Tj) -> void {eps_w.OutputBuds(Nj, Tj); }
		};

		std::array<std::function<void(uint32_t, variables&)>, 3> RetUploader = {
			[&](uint32_t Nj, variables& X)->void { X.ret.hel_12 = healt.GetRetValue(12.0, Nj); },
			[&](uint32_t Nj, variables& X)->void { X.ret.adh_12 = adhes.GetRetValue(12.0, Nj); },
			[&](uint32_t Nj, variables& X)->void { X.ret.adh_4 = adhes.GetRetValue(4.0, Nj); }
		};

		std::array<std::function<void(uint32_t, variables&)>, 3> RetDataUpdater = {
			[&](uint32_t Nj, variables& X) -> void { healt.ShiftRets(Nj, X.hel); },
			[&](uint32_t Nj, variables& X) -> void { df.ShiftRets(Nj, X.d_F); },
			[&](uint32_t Nj, variables& X) -> void { adhes.ShiftRets(Nj, X.adh); }
		};

		std::array<std::function<void()>, 3> RetInitialiser = {
			[&]() -> void {adhes.SetInitialRet(); },
			[&]() -> void {healt.SetInitialRet(); },
			[&]() -> void {df.SetInitialRet(); }
		};

		void ApplyEulersMethod()
		{
			ST.predictor(ST.X_pred.nec, ST.X_prev.nec, necr.RP.GetExpression(ST.X_prev));
			ST.predictor(ST.X_pred.acu_c, ST.X_prev.acu_c, acu.RP.GetExpression(ST.X_prev));
			ST.predictor(ST.X_pred.hel, ST.X_prev.hel, healt.RP.GetExpression(ST.X_prev));

			ST.predictor(ST.X_pred.cy, ST.X_prev.cy, cyto.RP.GetExpression(ST.X_prev));
			ST.predictor(ST.X_pred.adh, ST.X_prev.adh, adhes.RP.GetExpression(ST.X_prev));

			ST.predictor(ST.X_pred.lm, ST.X_prev.lm, lmacr.RP.GetExpression(ST.X_prev));
			ST.predictor(ST.X_pred.ln, ST.X_prev.ln, lneut.RP.GetExpression(ST.X_prev));

			ST.predictor(ST.X_pred.mia, ST.X_prev.mia, miact.RP.GetExpression(ST.X_prev));
			ST.predictor(ST.X_pred.mii, ST.X_prev.mii, minact.RP.GetExpression(ST.X_prev));

			for (auto const& cur : ExpressSubValues) { cur(ST.X_pred); }
		}

		void ApplyPredCorMethod()
		{
			ApplyEulersMethod();

			ST.corrector(ST.X_cor.nec, ST.X_prev.nec, necr.RP.GetExpression(ST.X_prev), necr.RP.GetExpression(ST.X_pred));
			ST.corrector(ST.X_cor.acu_c, ST.X_prev.acu_c, acu.RP.GetExpression(ST.X_prev), acu.RP.GetExpression(ST.X_pred));
			ST.corrector(ST.X_cor.hel, ST.X_prev.hel, healt.RP.GetExpression(ST.X_prev), healt.RP.GetExpression(ST.X_pred));

			ST.corrector(ST.X_cor.cy, ST.X_prev.cy, cyto.RP.GetExpression(ST.X_prev), cyto.RP.GetExpression(ST.X_pred));
			ST.corrector(ST.X_cor.adh, ST.X_prev.adh, adhes.RP.GetExpression(ST.X_prev), adhes.RP.GetExpression(ST.X_pred));

			ST.corrector(ST.X_cor.lm, ST.X_prev.lm, lmacr.RP.GetExpression(ST.X_prev), lmacr.RP.GetExpression(ST.X_pred));
			ST.corrector(ST.X_cor.ln, ST.X_prev.ln, lneut.RP.GetExpression(ST.X_prev), lneut.RP.GetExpression(ST.X_pred));

			ST.corrector(ST.X_cor.mia, ST.X_prev.mia, miact.RP.GetExpression(ST.X_prev), miact.RP.GetExpression(ST.X_pred));
			ST.corrector(ST.X_cor.mii, ST.X_prev.mii, minact.RP.GetExpression(ST.X_prev), minact.RP.GetExpression(ST.X_pred));

			for (auto const& cur : ExpressSubValues) { cur(ST.X_cor); }
		}
	};

	void SolveCurrentIndiv
	(
		uint32_t grid_spliting_amount,
		uint16_t full_amount_of_days,
		ReverseTask::BGA::Generation & indiv,
		ReverseTask::BGA::Parameters & RT
	)
	{
		// здесь, и только здесь мы выбираем метод решения прямой задачи	
		Aggregator SYS(grid_spliting_amount, full_amount_of_days);

		// Инициализируем накопитель запаздывающих аргументов
		for (auto const& cur : SYS.RetInitialiser) { cur(); }

		for (auto const& cur : SYS.DataCollector) { cur(); }

		// вносим текущие значения коэффициентов на конвеер
		for (uint16_t i = 0; i < RT.CoefsForChange.size(); i++)
			(*RT.CoefsForChange[i].second) = indiv.GetCoefVal(i);


		// и определяем начальные условия
		SYS.ST.X_init.tj = 0.5;
		for (auto const& cur : SYS.IniDataInitialiser) { cur(SYS.ST.X_init); }
		for (auto const& cur : SYS.ExpressSubValues) { cur(SYS.ST.X_init); }

		// и определяем начальные условия из предзаписанного решения
		//SYS.ST.t_0 = SYS.ST.X_init.tim = 0.5;
		//for (auto const& cur : SYS.SolDataInitialiser) { cur(0, 0, SYS.ST.X_init); }
		//for (auto const& cur : SYS.ExpressSubValues) { cur(SYS.ST.X_init); }

		// открываем потоки для вывода численного решения 
		for (auto const& cur : SYS.OutStreamAllocator) { cur(); }

		// выводим данные начальных условий во внешний файл с решением
		for (auto const& cur : SYS.SolutionOutputter) { cur(0, SYS.ST.X_init.tj, SYS.ST.X_init); }


		//bool TRIG = true;
		uint16_t current_day = 0;


		while (current_day < SYS.ST.full_amount_of_days)
		{
			std::cout << SYS.ST.full_amount_of_days - current_day << " "; // вывод на экран "кол-во дней, которые надо просчитать"
			SYS.ST.X_prev = SYS.ST.X_init; // даём начальное условие на предшествующий вектор
			double_t Tj = SYS.ST.X_pred.tj = SYS.ST.X_prev.tj;



			// цикл на следующие 24 часа
			for (uint32_t Nj = 1; Nj <= SYS.ST.N; Nj++)
			{
				// сдвиг на следующий шаг по времени
				Tj = SYS.ST.X_pred.tj += SYS.ST.H;

				//назначаем соответствующие запаздывания
				for (auto const& cur : SYS.RetUploader) { cur(Nj, SYS.ST.X_pred); }

				// контролируем промежуток интерполяции
				//ST.X_pred.CheckShiftInterpGap(STpar);

				// когда коррекция
				//ST.X_cor = ST.X_pred;


				// использовать предзаписанное решение
				//for (auto const& cur : SYS.SolDataInitialiser) { cur(current_day, Nj, SYS.ST.X_pred); }
				//SYS.ST.predictor(SYS.ST.X_pred.cy, SYS.ST.X_prev.cy, SYS.cyto.RP.GetExpression(SYS.ST.X_prev));

				// когда предиктор
				SYS.ApplyEulersMethod();

				//вывод, когда предиктор или предзаписанное
				for (auto const& cur : SYS.SolutionOutputter) { cur(Nj, Tj, SYS.ST.X_pred); }

				for (auto const& cur : SYS.BudgetOutputter) { cur(Nj, Tj); }
				

				// сдвигаем запаздывания (если они есть, лол)
				for (auto const& cur : SYS.RetDataUpdater) { cur(Nj, SYS.ST.X_prev); }

				// обновляем предшествующий временной ряд для перехода на следующий шаг

				SYS.ST.X_prev = SYS.ST.X_pred; //когда предиктор
				//ST.X_prev = ST.X_cor; // когда предиктор-корректор

			} // конец цикла рассчётов на текущий день for(Nj: 1->N)

			current_day++; std::cout << " ; ";

			// Проверка, «а надо ли решать дальше?» 
			if (current_day == SYS.ST.full_amount_of_days)
			{
				std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			}
			else SYS.ST.X_init = SYS.ST.X_pred; // когда предиктор
			//else ST.X_init = ST.X_cor; // когда предиктор-корректор
		}
	}
}


namespace ReverseTask
{

	void StartBGA()
	{
		// Определяем параметры метода обратной задачи
		BGA::Parameters RT;


		// Собираем данные из предзаписанных решений, которыми будем фиксировать поведение системы
		StraightTask::Aggregator SYS(3000, 4);
		for (auto const& cur : SYS.DataCollector) { cur(); }
		
		// предварительно: обязательно надо вызвать;
		//StraightTask::Splines::ConfigureSplines();

#define CONF_COEFS(STRUCT, NAME, ID, L, R) \
			RT.CoefsForChange.emplace_back( std::make_pair(NAME, &StraightTask::STRUCT::ID));\
			RT.bonds[0].emplace_back(L); RT.bonds[1].emplace_back(R)

		/*1*/		CONF_COEFS(Cytokines, "p_{Ma_{cy}}", p_Macy, 0.03, 15.0);
		/*2*/		CONF_COEFS(Cytokines, "C_{Ma}", C_Ma, 0.03, 0.9);

		/*3*/		CONF_COEFS(Cytokines, "p_{Lm_{cy}}", p_Lmcy, 1.0, 5.0);
		/*4*/		CONF_COEFS(Cytokines, "C_{Lm}", C_Lm, 0.03, 0.9);

		/*5*/		CONF_COEFS(Cytokines, "p_{Ln_{cy}}", p_Lncy, 0.1, 3.5);
		/*6*/		CONF_COEFS(Cytokines, "C_{Ln}", C_Ln, 0.03, 0.9);

		/*7*/		CONF_COEFS(Cytokines, "e_{cy}", e_cy, 0.03, 3.0);

		/*8*/		CONF_COEFS(Cytokines, "l_1", l1, 0.5, 5.0);
		/*9*/		CONF_COEFS(Cytokines, "l_2", l2, 0.5, 5.0);
		/*10*/		CONF_COEFS(Cytokines, "l_3", l3, 0.5, 5.0);

		/*11*/		CONF_COEFS(Cytokines, "A_{btm}", _A, 0.05, 1.0);
#undef CONF_COEFS

	


		
		RT.ios.Cout.open(output_dir + "RT/current/leaders_C.txt", std::ios_base::out);
		if (!RT.ios.Cout)
		{
			std::cout << "\n !!!stream for coeficients was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		RT.ios.Cout << '#'
			<< "mu = " << RT.Mutation << "\t"
			<< "d = " << RT.Recombination << "\t"
			<< "N_gen" << RT.iterAmount << "\t"
			<< "population : " << RT.p_0 << "->" << RT.p << std::endl;

		RT.ios.Cout << "#";
		for (size_t i = 0; i < RT.CoefsForChange.size(); i++)
			RT.ios.Cout << RT.CoefsForChange[i].first << "\t";
		RT.ios.Cout << std::endl;

		RT.ios.Cout.close();



		RT.ios.Fout.open(output_dir + "RT/current/leaders_F.txt", std::ios_base::out);
		if (!RT.ios.Fout)
		{
			std::cout << "\n !!!stream for F-values was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		RT.ios.Fout.close();




//		 this was included
		vector<BGA::Generation> individuals(RT.p_0 + 1);
		individuals.shrink_to_fit();// */

//===============================================================================================================================================================================
		//// определяем первого индивида дефолтными значениями коэффициентов
		//for (size_t i = 0; i < RT.CoefsForChange.size(); i++)
		//	individuals[1].AddCoefVal(*RT.CoefsForChange[i].second);

		srand((unsigned int)(time(0))); // запускаем генератор случайных чисел
		
		// определяем первичные случайные значения коэффициентов
		for (size_t q = 1; q <= RT.p_0; q++)
			for (size_t i = 0; i < RT.CoefsForChange.size(); i++)
				individuals[q].AddCoefVal(RT.bonds[0][i], RT.bonds[1][i]);



		// да запустится алгоритм!
		for (uint16_t k = 1; k <= RT.iterAmount; k++)
		{
			std::cout << "\n" << k << '/' << RT.iterAmount << "th iteration of BGA: " << std::endl;

			// 1. Строим поколение алгоритма BGA на решениях задачи
			// подготовим выделенных участников к анализу
			for (size_t q = 1; q < individuals.size(); q++)
			{
				// фиксируем индекс генерации и тем самым проверяем, посчитано ли решение для этого участника
				if (!individuals[q].SetGenu(k)) continue;

				// вносим текущие значения коэффициентов на конвеер
				for (uint16_t i = 0; i < RT.CoefsForChange.size(); i++)
					(*RT.CoefsForChange[i].second) = individuals[q].GetCoefVal(i);

				// вычисляем значение функционала->min для текущего inst[q]
				{
					// определяем начальные условия по предзаписанным решениям
					SYS.ST.X_init.tj = SYS.ST.t_0;
					// для следующей строки необходимо, чтобы выше данные о предзаписанных решениях уже были собраны!
					for (auto const& cur : SYS.SolDataInitialiser) { cur(0, 0, SYS.ST.X_init); }

					// контейнер c величинами, которые составляют значение функционала->min
					for (uint16_t i = 0; i < RT.ControlValues.size(); i++)
						individuals[q].budgets.emplace_back(0);
					

					uint16_t M = 0;
					uint16_t current_day = 0;
					

					while (current_day < SYS.ST.full_amount_of_days)
					{
						// вывод на экран "кол-во дней, которые надо просчитать"
						//std::cout << ST.full_amount_of_days - current_day << " "; 

						SYS.ST.X_prev = SYS.ST.X_init; // даём начальное условие на предшествующий вектор
						SYS.ST.X_pred.tj = SYS.ST.X_prev.tj;
						
						// цикл на следующие 24 часа
						for (uint32_t Nj = 1; Nj <= SYS.ST.N; Nj++)
						{
							// сдвиг на следующий шаг по времени
							SYS.ST.X_pred.tj += SYS.ST.H;

							
							for (auto const& cur : SYS.SolDataInitialiser) { cur(current_day, Nj, SYS.ST.X_pred); }
							SYS.ST.predictor(SYS.ST.X_pred.cy, SYS.ST.X_prev.cy, SYS.cyto.RP.GetExpression(SYS.ST.X_prev));

							// здесь считаем вклад в невязку по каждой компоненте
							if (individuals[q].budgets[0] > 100)
									individuals[q].budgets[0] = 100;
							else
							{
								double_t tmp;
								tmp = SYS.ST.X_pred.cy - SYS.cyto.GetSolData(current_day, Nj, SYS.ST.N);
								individuals[q].budgets[0] += tmp * tmp;
							}
							M++;
							

							// обновляем предшествующий временной ряд для перехода на следующий шаг
							SYS.ST.X_prev = SYS.ST.X_pred;

						} // конец цикла рассчётов на текущий день for(Nj: 1->N)

						current_day++; // std::cout << " ; ";
						SYS.ST.X_init = SYS.ST.X_pred;
					}

					for (uint16_t i = 0; i < RT.ControlValues.size(); i++)
						individuals[q].budgets[i] /= (double_t)(M)* RT.ControlValues.size();
					// умножение на размер массива приближаемых величин нужен для осуществления
					// поддержания однородности значений F в ситуациях, когда число приближаемых значений меняется
				}
				individuals[q].SumUpBuds();
			}

			// 2. Сортируем по возрастанию значения F_value построенный массив до 30% доли его носителя
			BGA::Generation::Sort1stFrac(individuals, RT.Sp);

			// 2.5 Проверяем, сошёлся ли функционал невязки
			//if (abs(indiv[1].Get_F() - indiv[CurRTP.SortFraction].Get_F()) < CurRTP.eps) break;

			// 3. отпускаем часть популяции из рассмотрения
			while (individuals.size() >= RT.Sp + 1) individuals.pop_back();


			// 4. Рекомбинация до исходного числа p участников и их мутация
			for (size_t q = RT.Sp; q <= RT.p; q++)
			{
				individuals.emplace_back();

				int parents[] = { (int)(random(1,RT.Sp)), (int)(random(1,RT.Sp)) };

				for (uint16_t j = 0; j < RT.CoefsForChange.size(); j++)
					individuals[q].Recombine(
								individuals[parents[0]].GetCoefVal(j),
								individuals[parents[1]].GetCoefVal(j),
								RT.Recombination);

				for (uint16_t j = 0; j < RT.CoefsForChange.size(); j++)
					individuals[q].Mutate(j,
								RT.bonds[0][j],
								RT.bonds[1][j],
								RT.Mutation);
			}

			RT.ios.WriteResult(individuals[1]);

		}

		RT.ios.WriteBest(individuals[1]);
		RT.ios.WriteStatData(individuals[1]);

		// отпускаем популяцию, порождённую последней итерацией алгоритма BGA
		for (size_t i = individuals.size() - 1; i >= RT.Sp + 1; i--) individuals.pop_back();
		std::cout << "\n Reversed task is solved for now \n";

		StraightTask::SolveCurrentIndiv(3000, 4, individuals[1], RT);
		getchar();
		return;
	}
}