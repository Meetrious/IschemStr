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

	// ��� ���������������� ��� ��������� �������, ���
	Methods::PredCor TLayers(CurSTP);
	Methods STask; STask.AllocateOutputStream();


	uint16_t days_remain = faod;
	Methods::H = CurSTP.H;
	bool TRIG = true;
	while (TRIG)
	{
		//std::cout << days_remain << " "; // ����� �� ����� "���-�� ����, ������� ���� ����������"
		TLayers.X_prev = TLayers.X_init; // ��� ��������� ������� �� �������������� ������

		// ���� �� ��������� 24 ����
		for (size_t Nj = 1; Nj <= CurSTP.N; Nj++)
		{
			// ����� �� ��������� ��� �� �������
			TLayers.X_pred.tim += CurSTP.H;

			//��������� ��������������� ������������
			TLayers.X_pred.AllocCurRets(Nj, CurSTP.N);

			// ������������ ���������� ������������
			TLayers.X_pred.CheckShiftInterpGap(CurSTP);
			TLayers.X_cor = TLayers.X_pred;

			launchPredCor(TLayers);

			// ����� �� ������� *.txt
			for (size_t i = 1; i < Methods::N_eq + 6; i++) STask.Output[i](TLayers.X_cor, Nj);

			// �������� ������������ (���� ��� ����, ���)
			rets::ShiftRets(TLayers.X_prev, Nj);

			// ��������� �������������� ��������� ��� ��� �������� �� ��������� ���
			TLayers.X_prev = TLayers.X_cor;

		} // ����� ����� ��������� �� ������� ���� for(Nj: 1->N)

		days_remain--; //std::cout << " ; ";

		// ��������, �� ���� �� ������ ������?� 
		if (days_remain == 0)
		{
			std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			TRIG = false; // ������� �� ����� �� ����� while
			// else �� �����, ���� ���� !ID, �� ������� ������� while ��� �������
		}
		else TLayers.X_init = TLayers.X_cor;
	}
}


namespace ReverseTask
{

	void CountFBuds
	// ���� ����� ������� ��������� �������� ����������� ������� ���������� ��������� ���������� ������� ������ ������
	(
		std::vector<double_t> & buds,
		std::vector<NumOfDerivative> CVals,
		StraightTask::Parameters & CurSTP,
		StraightTask::Methods::PredCor & TL
	)
	{
		using namespace StraightTask;

		// ��� ���������������� ��� ��������� �������, ���
		TL.ResetToInitial(CurSTP);

		// ��������� c ����������, ������� ���������� �������� �����������->min
		for (uint16_t i = 0; i < CVals.size(); i++)
			buds.emplace_back(0);
		std::vector<bool> ACVals(CVals.size(), true);
		uint16_t M = 0;

		uint16_t days_remain = CurSTP.FAOD;
		bool TRIG = true;
		while (TRIG)
		{
			//std::cout << days_remain << " "; // ����� �� ����� "���-�� ����, ������� ���� ����������"
			TL.X_prev = TL.X_init; // ��� ��������� ������� �� �������������� ������

			// ���� �� ��������� 24 ����
			for (size_t Nj = 1; Nj <= CurSTP.N; Nj++)
			{
				// ����� �� ��������� ��� �� �������
				TL.X_pred.tim += CurSTP.H;

				//��������� ��������������� ������������
				//TLayers.X_pred.AllocCurRets(Nj, CurSTP.N);

				// ������������ ���������� ������������
				TL.X_pred.CheckShiftInterpGap(CurSTP);
				TL.X_cor = TL.X_pred;

				launchPredCor(TL);

				// ����� ������� ����� � ������� �� ������ ����������

				for (size_t i = 0; i < CVals.size(); i++)
				{
					if (buds[i] > 100)
						buds[i] = 100;
					else
						BGA_Parameters::account_investments2(CVals[i], buds[i], TL.X_cor); // ��������� ������� ����� � ����������->min 
				}
				M++;

				// �������� ������������ (���� ��� ����, ���)
				//rets::ShiftRets(TLayers.X_prev, Nj);

				// ��������� �������������� ��������� ��� ��� �������� �� ��������� ���
				TL.X_prev = TL.X_cor;

			} // ����� ����� ��������� �� ������� ���� for(Nj: 1->N)

			days_remain--; //std::cout << " ; ";

			// ��������, �� ���� �� ������ ������?� 
			if (days_remain == 0)
			{
				//std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
				TRIG = false; // ������� �� ����� �� ����� while
				// else �� �����, ���� ���� !ID, �� ������� ������� while ��� �������
			}
			else TL.X_init = TL.X_cor;
		}

		for (uint16_t i = 0; i < CVals.size(); i++)
			buds[i] /= (double_t)(M)* CVals.size();
		// ��������� �� ������ ������� ������������ ������� ����� ��� �������������
		// ����������� ������������ �������� F � ���������, ����� ����� ������������ �������� ��������
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
		//�������� �������� �� �������� ������
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

		// ���������� ������� �������� ���������� ���������� �������������
		for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
			indiv[1].AddCoefVal(*CurRTP.CoefsForChange[i].second);

		srand((unsigned int)(time(0))); // ��������� ��������� ��������� �����
		// ���������� ��������� ��������� �������� �������������
		for (size_t q = 2; q <= CurRTP.p_0; q++)
			for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				indiv[q].AddCoefVal(CurRTP.bonds[0][i], CurRTP.bonds[1][i]);

		// ��������� ��������� ���������
		CurRTP.SetInfoPlot(0);
		CurRTP.start_gif_conf(1);

		// ���������� ������� ����������:
		StraightTask::Methods::PredCor TimeLayers(CurRTP); // ����� ����� �� �� ����������������, � ������ ��� ������ ���������.
		StraightTask::Methods::H = CurRTP.H; // �����������, ����� ��������� � Predictions, Corrections ������� ����� ���������� �����������

		// �� ���������� ��������!
		for (uint16_t k = 1; k <= CurRTP.iterAmount; k++)
		{
			std::cout << "\n" << k << '/' << CurRTP.iterAmount << "th iteration of BGA: " << std::endl;

			// 1. ������ ��������� ��������� BGA �� �������� ������
			// ���������� ���������� ���������� � �������
			for (size_t q = 1; q < indiv.size(); q++)
			{
				// ��������� ������ ��������� � ��� ����� ���������, ��������� �� ������� ��� ����� ���������
				if (!indiv[q].SetGenu(k)) continue;

				// ������ ������� �������� ������������� �� �������
				for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					(*CurRTP.CoefsForChange[i].second) = indiv[q].GetCoefVal(i);

				// ��������� �������� �����������->min ��� �������� inst[q]
				CountFBuds(indiv[q].budgets, CurRTP.ControlValues, CurRTP, TimeLayers);

				indiv[q].SumUpBuds();
			}

			// 2. ��������� �� ����������� �������� F_value ����������� ������ �� 30% ���� ��� ��������
			isBestFound = Generation::Sort1stFrac(indiv, CurRTP.SortFraction);

			outGen << "\n end of " << k << '/' << CurRTP.iterAmount << "th iteration. \n";
			outGen << "\n - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - \n"
				<< std::flush;

			// 2.5 ���������, ������� �� ���������� �������
			if (abs(indiv[1].Get_F() - indiv[CurRTP.SortFraction].Get_F()) < CurRTP.eps) break;

			// 3. ��������� ����� ��������� �� ������������
			while (indiv.size() >= CurRTP.SortFraction + 1) indiv.pop_back();


			// 4. ������������ �� ��������� ����� p ���������� � �� �������
			for (int q = CurRTP.SortFraction; q <= CurRTP.p; q++)
			{
				indiv.emplace_back();

				int parents[] = { (int)(random(1,CurRTP.SortFraction)), (int)(random(1,CurRTP.SortFraction)) };

				for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
					indiv[q].Recombination(indiv[parents[0]].GetCoefVal(j), indiv[parents[1]].GetCoefVal(j), CurRTP.RecombPar);

				for (uint16_t j = 0; j < CurRTP.CoefsForChange.size(); j++)
					indiv[q].Mutation(j, CurRTP.bonds[0][j], CurRTP.bonds[1][j], CurRTP.MutPar);
			}


			// ������� ����� ��������� �������
			for (int q = 1; q <= 20 && q < CurRTP.SortFraction; q++)
				indiv[q].bga_info_pin(outGen, CurRTP.CoefsForChange);
			// ������� ��������� ������� �� ������� ����
			indiv[CurRTP.SortFraction].bga_info_pin(outGen, CurRTP.CoefsForChange);

			// ����� �������� ������������� �� ���������� ������� ������� ��������
			outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out | std::ios_base::app);
			for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				outLeaders << std::setprecision(5) << indiv[1].GetCoefVal(i) << "\t\t";
			outLeaders << std::endl;
			outLeaders.close();

			//����� ������� �������� �����������
			outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out | std::ios_base::app);
			outBest << indiv[1].Get_F() << std::endl;
			outBest.close();

			if (isBestFound)
			{
				// ������ ������� �������� ������������� �� ������� ������ ������� ������� �� ������� ��������
				/*for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					(*CurRTP.CoefsForChange[i].second) = indiv[1].GetCoefVal(i);
				GetSolution(CurRTP, 3);*/
				CurRTP.collect_in_gif(k, indiv[1].Get_F());
				isBestFound = false;
			} //*/



			// ������������ ������� ��������� ���
			CurRTP.ProgressPlot(k, indiv[1].Get_F());
		}

		if (outBest.is_open()) outBest.close();
		if (outLeaders.is_open()) outLeaders.close();

		outGen.close();

		// ��������� ���������, ���������� ��������� ��������� ��������� BGA
		for (size_t i = indiv.size() - 1; i >= CurRTP.SortFraction + 1; i--) indiv.pop_back();

		/*// ���� ������ ���������� �������
		{
			// ������ ������� �������� ������������� �� ������� ������ ���������� �������
			for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
				(*CurRTP.CoefsForChange[i].second) = indiv[1].GetCoefVal(i);

			// ������ ������ ������ ��� ���������� ������
			GetSolution(CurRTP, 3);
		} // */

		// ���� ������ ������� �������� �� ������� �����������
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

		// ���� ���������� ������� �������� � ����� �� ������������ ����������
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

		// ���������� ������� ����������:
		StraightTask::Methods::PredCor TimeLayers(CurRTP); // ����� ����� �� �� ����������������, � ������ ��� ������ ���������.
		StraightTask::Methods::H = CurRTP.H; // �����������, ����� ��������� � Predictions, Corrections ����� ���������� �����������

		
		srand((unsigned int)(time(0))); // ��������� ��������� ��������� �����

		// �� ���������� ���� ����������!
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


			// ���������� ��������� ��������������� �������� �������������
			for (size_t q = 1; q <= CurRTP.p_0; q++)
				for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					indiv[q].AddCoefVal(CurRTP.bonds[0][i], CurRTP.bonds[1][i]);

			//��������� ������� ���������� �������������
			CurRTP.SetInfoPlot(1);


			// �� ���������� ��������!
			for (uint16_t k = 1; k <= CurRTP.iterAmount; k++)
			{
				std::cout << "\n" << k << '/' << CurRTP.iterAmount << "th iteration of BGA: " << std::endl;

				// 1. ������ ��������� ��������� BGA �� �������� ������
				// ���������� ���������� ���������� � �������
				for (size_t q = 1; q < indiv.size(); q++)
				{
					// ��������� ������ ��������� � ��� ����� ���������, ��������� �� ������� ��� ����� ���������
					if (!indiv[q].SetGenu(k)) continue;

					// ������ ������� �������� ������������� �� �������
					for (uint16_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
						(*CurRTP.CoefsForChange[i].second) = indiv[q].GetCoefVal(i);

					// ��������� �������� �����������->min ��� �������� inst[q]
					CountFBuds(indiv[q].budgets, CurRTP.ControlValues, CurRTP, TimeLayers);

					indiv[q].SumUpBuds();
				}

				// 2. ��������� �� ����������� �������� F_value ����������� ������ �� 30% ���� ��� ��������
				Generation::Sort1stFrac(indiv, CurRTP.SortFraction);


				// 2.5 ���������, ������� �� ���������� �������
				if (abs(indiv[1].Get_F() - indiv[CurRTP.SortFraction].Get_F()) < CurRTP.eps)
					break;

				// 3. ��������� ����� ��������� �� ������������
				while (indiv.size() >= CurRTP.SortFraction + 1) indiv.pop_back();


				// 4. ������������ �� ��������� ����� p ���������� � �� �������
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


				// ����� �������� ������������� �� ���������� ������� ������� ��������
				outLeaders.open(working_directory + bga_data_dir + "leaders_C.txt", std::ios_base::out | std::ios_base::app);
				for (size_t i = 0; i < CurRTP.CoefsForChange.size(); i++)
					outLeaders << std::setprecision(5) << indiv[1].GetCoefVal(i) << "\t\t";
				outLeaders << std::endl;
				outLeaders.close();

				//����� ������� �������� �����������
				outBest.open(working_directory + bga_data_dir + "leaders_F.txt", std::ios_base::out | std::ios_base::app);
				outBest << indiv[1].Get_F() << std::endl;
				outBest.close();

				//������������ ������� ��������� ���
				//CurRTP.ProgressPlot(k, indiv[1].Get_F());
			}

			if (outBest.is_open()) outBest.close();
			if (outLeaders.is_open()) outLeaders.close();

			// ��������� ���������, ���������� ��������� ��������� ��������� BGA
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