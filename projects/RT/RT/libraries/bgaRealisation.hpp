#pragma once
#include "ReverseTaskParameters.hpp"

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

void CoutProgress(size_t Nj, size_t N)
{
	int i = 2;
	if ((100 * Nj) / N >= 10) i = 3;
	if ((100 * Nj) / N == 100) i = 4;
	switch (i)
	{
	case 2: std::cout << "\b\b" << (100 * Nj) / N << "%"; break;
	case 3: std::cout << "\b\b\b" << (100 * Nj) / N << "%"; break;
	case 4: std::cout << "\b\b\b\b" << (100 * Nj) / N << "%"; break;
	default: std::cout << "\r" << (100 * Nj) / N << "%"; break;
	}
}

void GetSolution(StraightTask::Parameters const & CurSTP, uint16_t faod)
{
	using namespace StraightTask;

	// тут инициализируются все начальные условия, все
	Methods::PredCor TLayers(CurSTP);
	Methods STask; STask.AllocateOutputStream();


	uint16_t days_remain = faod;
	Methods::H = CurSTP.H;
	bool TRIG = true;
	while (TRIG)
	{
		//std::cout << days_remain << " "; // вывод на экран "кол-во дней, которые надо просчитать"
		TLayers.X_prev = TLayers.X_init; // даём начальное условие на предшествующий вектор

		// цикл на следующие 24 часа
		for (size_t Nj = 1; Nj <= CurSTP.N; Nj++)
		{
			// сдвиг на следующий шаг по времени
			TLayers.X_pred.tim += CurSTP.H;

			//назначаем соответствующие запаздывания
			TLayers.X_pred.AllocCurRets(Nj, CurSTP.N);

			// контролируем промежуток интерполяции
			TLayers.X_pred.CheckShiftInterpGap(CurSTP);
			TLayers.X_cor = TLayers.X_pred;

			launchPredCor(TLayers);

			// вывод во внешние *.txt
			for (size_t i = 1; i < Methods::N_eq + 6; i++) STask.Output[i](TLayers.X_cor, Nj);

			// сдвигаем запаздывания (если они есть, лол)
			rets::ShiftRets(TLayers.X_prev, Nj);

			// обновляем предшествующий временной ряд для перехода на следующий шаг
			TLayers.X_prev = TLayers.X_cor;

		} // конец цикла рассчётов на текущий день for(Nj: 1->N)

		days_remain--; //std::cout << " ; ";

		// Проверка, «а надо ли решать дальше?» 
		if (days_remain == 0)
		{
			std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			TRIG = false; // условие на выход из цикла while
			// else не нужен, ведь если !ID, то внешнее условие while его поймает
		}
		else TLayers.X_init = TLayers.X_cor;
	}
}


namespace ReverseTask
{

	void CountFBuds
	// этот метод призван вычислить значение функционала невязки посредсвом получения численного решения прямой задачи
	(
		std::vector<double_t> & buds,
		std::vector<NumOfDerivative> CVals,
		StraightTask::Parameters & CurSTP,
		StraightTask::Methods::PredCor & TL
	)
	{
		using namespace StraightTask;

		// тут инициализируются все начальные условия, все
		TL.ResetToInitial(CurSTP);

		// контейнер c величинами, которые составляют значение функционала->min
		for (uint16_t i = 0; i < CVals.size(); i++)
			buds.emplace_back(0);
		std::vector<bool> ACVals(CVals.size(), true);
		uint16_t M = 0;

		uint16_t days_remain = CurSTP.FAOD;
		bool TRIG = true;
		while (TRIG)
		{
			//std::cout << days_remain << " "; // вывод на экран "кол-во дней, которые надо просчитать"
			TL.X_prev = TL.X_init; // даём начальное условие на предшествующий вектор

			// цикл на следующие 24 часа
			for (size_t Nj = 1; Nj <= CurSTP.N; Nj++)
			{
				// сдвиг на следующий шаг по времени
				TL.X_pred.tim += CurSTP.H;

				//назначаем соответствующие запаздывания
				//TLayers.X_pred.AllocCurRets(Nj, CurSTP.N);

				// контролируем промежуток интерполяции
				TL.X_pred.CheckShiftInterpGap(CurSTP);
				TL.X_cor = TL.X_pred;

				launchPredCor(TL);

				// здесь считаем вклад в невязку по каждой компоненте

				for (size_t i = 0; i < CVals.size(); i++)
				{
					if (buds[i] > 100)
						buds[i] = 100;
					else
						BGA_Parameters::account_investments2(CVals[i], buds[i], TL.X_cor); // вычисляем текущий вклад в функционал->min 
				}
				M++;

				// сдвигаем запаздывания (если они есть, лол)
				//rets::ShiftRets(TLayers.X_prev, Nj);

				// обновляем предшествующий временной ряд для перехода на следующий шаг
				TL.X_prev = TL.X_cor;

			} // конец цикла рассчётов на текущий день for(Nj: 1->N)

			days_remain--; //std::cout << " ; ";

			// Проверка, «а надо ли решать дальше?» 
			if (days_remain == 0)
			{
				//std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
				TRIG = false; // условие на выход из цикла while
				// else не нужен, ведь если !ID, то внешнее условие while его поймает
			}
			else TL.X_init = TL.X_cor;
		}

		for (uint16_t i = 0; i < CVals.size(); i++)
			buds[i] /= (double_t)(M)* CVals.size();
		// умножение на размер массива приближаемых величин нужен для осуществления
		// поддержания однородности значений F в ситуациях, когда число приближаемых значений меняется
	}


	void StartBGA
	(
		std::vector<Generation> & indiv,
		BGA_Parameters & CurRTP
	)
	{
		std::ofstream outGen;
		std::ofstream outBest;
		std::ofstream outLeaders;
		bool isBestFound = false;

		outGen.open(working_directory + bga_data_dir + "general info.txt", std::ios_base::out);
		if (!outGen)
		{
			std::cout << "\n !!!stream for outGen was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		//написать проверку на открытие потока
		CurRTP.outInfo(outGen);
		CurRTP.coutInfo();

		outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out);
		if (!outLeaders)
		{
			std::cout << "\n !!!stream for outLeaders was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		outLeaders << '#'
			<< "mu = " << CurRTP.MutPar << "\t"
			<< "d = " << CurRTP.RecombPar << "\t"
			<< "N_gen" << CurRTP.iterAmount << "\t"
			<< "population : " << CurRTP.p_0 << "->" << CurRTP.p << std::endl;

		outLeaders << "#";
		for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
			outLeaders << CurRTP.CoefsForChange[i].first << "\t\t";
		outLeaders << std::endl;

		outLeaders.close();

		outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out);
		if (!outBest)
		{
			std::cout << "\n !!!stream for outBest was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		outBest.close();

		// определяем первого индивида дефолтными значениями коэффициентов
		for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
			indiv[1].AddCoefVal(*CurRTP.CoefsForChange[i].second);

		srand((unsigned int)(time(0))); // запускаем генератор случайных чисел
		// определяем первичные случайные значения коэффициентов
		for (size_t q = 2; q <= CurRTP.p_0; q++)
			for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				indiv[q].AddCoefVal(CurRTP.bonds[0][i], CurRTP.bonds[1][i]);

		// первичная настройка терминала
		CurRTP.SetInfoPlot(0);
		CurRTP.start_gif_conf(1);

		// ускоряющие процесс подготовки:
		StraightTask::Methods::PredCor TimeLayers(CurRTP); // чтобы вечно их не инициализировать, а каждый раз только обновлять.
		StraightTask::Methods::H = CurRTP.H; // ОБЯЗАТЕЛЬНО, иначе уравнения в Predictions, Corrections массивы будут определены неправильно

		// да запустится алгоритм!
		for (uint16_t k = 1; k <= CurRTP.iterAmount; k++)
		{
			std::cout << "\n" << k << '/' << CurRTP.iterAmount << "th iteration of BGA: " << std::endl;

			// 1. Строим поколение алгоритма BGA на решениях задачи
			// подготовим выделенных участников к анализу
			for (size_t q = 1; q < indiv.size(); q++)
			{
				// фиксируем индекс генерации и тем самым проверяем, посчитано ли решение для этого участника
				if (!indiv[q].SetGenu(k)) continue;

				// вносим текущие значения коэффициентов на конвеер
				for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					(*CurRTP.CoefsForChange[i].second) = indiv[q].GetCoefVal(i);

				// вычисляем значение функционала->min для текущего inst[q]
				CountFBuds(indiv[q].budgets, CurRTP.ControlValues, CurRTP, TimeLayers);

				indiv[q].SumUpBuds();
			}

			// 2. Сортируем по возрастанию значения F_value построенный массив до 30% доли его носителя
			isBestFound = Generation::Sort1stFrac(indiv, CurRTP.SortFraction);

			outGen << "\n end of " << k << '/' << CurRTP.iterAmount << "th iteration. \n";
			outGen << "\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
				<< std::flush;

			// 2.5 Проверяем, сошёлся ли функционал невязки
			if (abs(indiv[1].Get_F() - indiv[CurRTP.SortFraction].Get_F()) < CurRTP.eps) break;

			// 3. отпускаем часть популяции из рассмотрения
			while (indiv.size() >= CurRTP.SortFraction + 1) indiv.pop_back();


			// 4. Рекомбинация до исходного числа p участников и их мутация
			for (int q = CurRTP.SortFraction; q <= CurRTP.p; q++)
			{
				indiv.emplace_back();

				int parents[] = { (int)(random(1,CurRTP.SortFraction)), (int)(random(1,CurRTP.SortFraction)) };

				for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
					indiv[q].Recombination(indiv[parents[0]].GetCoefVal(j), indiv[parents[1]].GetCoefVal(j), CurRTP.RecombPar);

				for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
					indiv[q].Mutation(j, CurRTP.bonds[0][j], CurRTP.bonds[1][j], CurRTP.MutPar);
			}


			// выводим общий результат лидеров
			for (int q = 1; q <= 20 && q < CurRTP.SortFraction; q++)
				indiv[q].bga_info_pin(outGen, CurRTP.CoefsForChange);
			// выводим результат худшего во внешний файл
			indiv[CurRTP.SortFraction].bga_info_pin(outGen, CurRTP.CoefsForChange);

			// вывод значений коэффициентов от наилучшего решения текущей итерации
			outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out | std::ios_base::app);
			for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outLeaders << std::setprecision(5) << indiv[1].GetCoefVal(i) << "\t\t";
			outLeaders << std::endl;
			outLeaders.close();

			//вывод лучшего значения функционала
			outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out | std::ios_base::app);
			outBest << indiv[1].Get_F() << std::endl;
			outBest.close();

			if (isBestFound)
			{
				// вносим текущие значения коэффициентов на конвеер вывода лучшего решения на текущей итерации
				/*for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					(*CurRTP.CoefsForChange[i].second) = indiv[1].GetCoefVal(i);
				GetSolution(CurRTP, 3);*/
				CurRTP.collect_in_gif(k, indiv[1].Get_F());
				isBestFound = false;
			} //*/



			// отрисовываем текущее положение дел
			CurRTP.ProgressPlot(k, indiv[1].Get_F());
		}

		if (outBest.is_open()) outBest.close();
		if (outLeaders.is_open()) outLeaders.close();

		outGen.close();

		// отпускаем популяцию, порождённую последней итерацией алгоритма BGA
		for (size_t i = indiv.size() - 1; i >= CurRTP.SortFraction + 1; i--) indiv.pop_back();

		/*// блок вывода финального решения
		{
			// вносим текущие значения коэффициентов на конвеер вывода финального решения
			for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				(*CurRTP.CoefsForChange[i].second) = indiv[1].GetCoefVal(i);

			// решаем прямую задачу для финального вывода
			GetSolution(CurRTP, 3);
		} // */

		// блок вывода лучшего индивида за процесс оптимизации
		{
			outBest.open(working_directory + bga_data_dir + "bestIndiv.txt", std::ios_base::out);
			if (!outBest)
			{
				std::cout << "\n outBest wasn't initialised for some reason \n ";
				system("pause");
				return;
			}

			for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outBest << indiv[1].GetCoefVal(i) << "\n";
			outBest.close();
		}

		// блок сохранения лучшего индивида в набор из неоднородной статистики
		{
			outLeaders.open(working_directory + bga_data_dir + "accumul/stat.txt", std::ios_base::out | std::ios_base::app);
			if (!outLeaders)
			{
				std::cout << "\n outLeaders wasn't initialised for some reason \n ";
				system("pause");
				return;
			}

			outLeaders << "#";
			for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outLeaders << CurRTP.CoefsForChange[i].first << "\t\t";
			outLeaders << "F_min" << std::endl;
			for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outLeaders << indiv[1].GetCoefVal(i) << "\t\t";
			outLeaders << indiv[1].Get_F() << std::endl;
			outLeaders.close();
		}

		CurRTP.GifLine("unset output");
		getchar();
		return;
	}



	void StartBGA_For_Stat
	(
		int16_t ExerPow,
		std::vector<Generation> & indiv,
		BGA_Parameters & CurRTP
	)
	{
		std::ofstream outBest;
		std::ofstream outLeaders;

		CurRTP.coutInfo();

		/*outLeaders.open(working_directory + bga_data_dir + "accumul/stat11.txt", std::ios_base::out);
		if (!outLeaders)
		{
			std::cout << "\n outLeaders wasn't initialised for some reason \n ";
			system("pause");
			return;
		}
		outLeaders << "#";
		for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
			outLeaders << CurRTP.CoefsForChange[i].first << "\t\t";
		outLeaders << "F_min" << std::endl;
		outLeaders.close();*/



		outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out);
		if (!outBest)
		{
			std::cout << "\n !!!stream for outBest was not allocated for some reason.\n Care to attend!";
			getchar();
			return;
		}
		outBest.close();

		// ускоряющие процесс подготовки:
		StraightTask::Methods::PredCor TimeLayers(CurRTP); // чтобы вечно их не инициализировать, а каждый раз только обновлять.
		StraightTask::Methods::H = CurRTP.H; // ОБЯЗАТЕЛЬНО, иначе уравнения в Predictions, Corrections будут определены неправильно

		
		srand((unsigned int)(time(0))); // запускаем генератор случайных чисел

		// да запустится сбор статистики!
		while (ExerPow > 0)
		{
			std::cout << std::endl << ExerPow << " runs is left " << std::endl;

			outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out);
			if (!outLeaders)
			{
				std::cout << "\n !!!stream for outLeaders was not allocated for some reason.\n Care to attend!";
				getchar();
				return;
			}
			outLeaders << '#'
				<< "mu = " << CurRTP.MutPar << "\t"
				<< "d = " << CurRTP.RecombPar << "\t"
				<< "N_gen = " << CurRTP.iterAmount << "\t"
				<< "population : " << CurRTP.p_0 << "->" << CurRTP.p << std::endl;

			outLeaders << "#";
			for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outLeaders << CurRTP.CoefsForChange[i].first << "\t\t";
			outLeaders << std::endl;
			outLeaders.close();


			outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out);
			if (!outBest)
			{
				std::cout << "\n !!!stream for outBest was not allocated for some reason.\n Care to attend!";
				getchar();
				return;
			}
			outBest.close();


			// определяем первичные псевдослучайные значения коэффициентов
			for (size_t q = 1; q <= CurRTP.p_0; q++)
				for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					indiv[q].AddCoefVal(CurRTP.bonds[0][i], CurRTP.bonds[1][i]);

			//отрисовка главных параметров экспериментов
			CurRTP.SetInfoPlot(1);


			// да запустится алгоритм!
			for (uint16_t k = 1; k <= CurRTP.iterAmount; k++)
			{
				std::cout << "\n" << k << '/' << CurRTP.iterAmount << "th iteration of BGA: " << std::endl;

				// 1. Строим поколение алгоритма BGA на решениях задачи
				// подготовим выделенных участников к анализу
				for (size_t q = 1; q < indiv.size(); q++)
				{
					// фиксируем индекс генерации и тем самым проверяем, посчитано ли решение для этого участника
					if (!indiv[q].SetGenu(k)) continue;

					// вносим текущие значения коэффициентов на конвеер
					for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
						(*CurRTP.CoefsForChange[i].second) = indiv[q].GetCoefVal(i);

					// вычисляем значение функционала->min для текущего inst[q]
					CountFBuds(indiv[q].budgets, CurRTP.ControlValues, CurRTP, TimeLayers);

					indiv[q].SumUpBuds();
				}

				// 2. Сортируем по возрастанию значения F_value построенный массив до 30% доли его носителя
				Generation::Sort1stFrac(indiv, CurRTP.SortFraction);


				// 2.5 Проверяем, сошёлся ли функционал невязки
				if (abs(indiv[1].Get_F() - indiv[CurRTP.SortFraction].Get_F()) < CurRTP.eps)
					break;

				// 3. отпускаем часть популяции из рассмотрения
				while (indiv.size() >= CurRTP.SortFraction + 1) indiv.pop_back();


				// 4. Рекомбинация до исходного числа p участников и их мутация
				for (int q = CurRTP.SortFraction; q <= CurRTP.p; q++)
				{
					indiv.emplace_back();

					int parents[] = {
						(int)(random(1, CurRTP.SortFraction - 1)),
						(int)(random(1, CurRTP.SortFraction - 1)) 
					};

					for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
						indiv[q].Recombination(indiv[parents[0]].GetCoefVal(j), indiv[parents[1]].GetCoefVal(j), CurRTP.RecombPar);

					for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
						indiv[q].Mutation(j, CurRTP.bonds[0][j], CurRTP.bonds[1][j], CurRTP.MutPar);
				}


				// вывод значений коэффициентов от наилучшего решения текущей итерации
				outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out | std::ios_base::app);
				for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					outLeaders << std::setprecision(5) << indiv[1].GetCoefVal(i) << "\t\t";
				outLeaders << std::endl;
				outLeaders.close();

				//вывод лучшего значения функционала
				outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out | std::ios_base::app);
				outBest << indiv[1].Get_F() << std::endl;
				outBest.close();

				//отрисовываем текущее положение дел
				//CurRTP.ProgressPlot(k, indiv[1].Get_F());
			}

			if (outBest.is_open()) outBest.close();
			if (outLeaders.is_open()) outLeaders.close();

			// отпускаем популяцию, порождённую последней итерацией алгоритма BGA
			for (int i = (int)(indiv.size() - 1); i >= CurRTP.SortFraction + 1; i--)
				indiv.pop_back();

			{
				outLeaders.open(working_directory + bga_data_dir + "accumul/stat11.txt", std::ios_base::out | std::ios_base::app);
				if (!outLeaders)
				{
					std::cout << "\n outLeaders wasn't initialised for some reason \n ";
					system("pause");
					return;
				}
				for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					outLeaders << indiv[1].GetCoefVal(i) << "\t\t";
				outLeaders << indiv[1].Get_F() << std::endl;
				outLeaders.close();
			}

			indiv.clear();
			while (indiv.size() <= CurRTP.p_0) indiv.emplace_back();
			ExerPow--;
		}
		return;
	}
}