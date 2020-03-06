#include "libraries/bgaRealisation.hpp"

// вычисл€ет значение функционала нев€зки на решении со значени€ми коэффициентов по умолчанию
double_t GetCurrentF()
{
	using namespace ReverseTask;
	StraightTask::Splines::ConfigureSplines();
	BGA_Parameters::SetCoefsForChange();
	BGA_Parameters CurRTP;
	Generation individ;

	double_t* CFC[12] = {
		&StraightTask::Neurons::k_N,
		&StraightTask::Neurons::k_A,
		&StraightTask::Neurons::p_R,

		&StraightTask::ToxDamage::p_ncy,
		&StraightTask::ToxDamage::C_Dcy,
		&StraightTask::ToxDamage::C_DLn,
		&StraightTask::ToxDamage::C_DLm,
		&StraightTask::ToxDamage::P_nn,
		&StraightTask::ToxDamage::p_Lm,
		&StraightTask::ToxDamage::p_Ln,
		&StraightTask::ToxDamage::D_0,
		&StraightTask::ToxDamage::p_D
	};
	StraightTask::Methods::setCoefs(CFC, 12, working_directory+"mstat/best.txt");
	
	// ускор€ющие процесс подготовки:
	StraightTask::Methods::PredCor TimeLayers(CurRTP);
	StraightTask::Methods::H = CurRTP.H; // ќЅя«ј“≈Ћ№Ќќ, иначе уравнени€ в Predictions, Corrections массивы будут определены неправильно

	// вычисл€ем значение функционала->min дл€ текущего inst[q]
	CountFBuds(individ.budgets, CurRTP.ControlValues, CurRTP, TimeLayers);
	individ.SumUpBuds();

	return individ.Get_F();
}

void GetSols()
{
	using namespace ReverseTask;
	StraightTask::Splines::ConfigureSplines();
	BGA_Parameters::SetCoefsForChange();

	StraightTask::Parameters CurSTP(2, 1500);


	std::ifstream in(working_directory + bga_data_dir + "/accumul/stat2.txt");
	if (!in)
	{
		std::cout << "\n file was not found or something.(( Probably there is some typo in the direction of a *.txt file \n";
		system("pause");
		return;
	}

	std::vector<Generation>indivs;
	for (size_t i = 0; !in.eof(); i++)
	{
		double_t tmp;
		indivs.emplace_back();
		for (size_t j = 0; j < BGA_Parameters::CoefsForChange.size(); j++)
		{
			in >> tmp;
			indivs[i].AddCoefVal(tmp);
		}
		in >> tmp;
	}
	in.close();


	Gnuplot nec_line, acu_line, hel_line;
	
	nec_line.cd(working_directory + ST_sol + "plot_scripts/for_comparison");
	acu_line.cd(working_directory + ST_sol + "plot_scripts/for_comparison");
	hel_line.cd(working_directory + ST_sol + "plot_scripts/for_comparison");

	for (size_t i = 0; i < indivs.size(); i++)
	{
		nec_line.png(4700, 4000, 180, working_directory + "BGA/report/current/gnu_output/", "NCD_" + std::to_string(i));
		acu_line.png(4700, 4000, 180, working_directory + "BGA/report/current/gnu_output/", "DCD_" + std::to_string(i));
		hel_line.png(4700, 4000, 180, working_directory + "BGA/report/current/gnu_output/", "ICD_" + std::to_string(i));

		for (size_t j = 0; j < BGA_Parameters::CoefsForChange.size(); j++)
			(*BGA_Parameters::CoefsForChange[j].second) = indivs[i].GetCoefVal(j);

		GetSolution(CurSTP, CurSTP.FAOD);

		nec_line("load \'NEC.plt\'");
		acu_line("load \'ACU.plt\'");
		hel_line("load \'HEL.plt\'");
		getchar();
	}
	nec_line("unset output"); 
	acu_line("unset output");
	hel_line("unset output");
}

void GetSol()
{
	using namespace ReverseTask;

	StraightTask::Splines::ConfigureSplines();
	BGA_Parameters::SetCoefsForChange();

	StraightTask::Parameters STParameters(3, 1500);
	
	std::ifstream in(working_directory + bga_data_dir + "bestIndiv.txt");
	if (!in)
	{
		std::cout << "\n file was not found or something.(( Probably there is some typo in the direction of a *.txt file \n";
		system("pause");
		return;
	}

	for (int i = 0; i < BGA_Parameters::CoefsForChange.size(); i++)
		in >> (*BGA_Parameters::CoefsForChange[i].second);

	in.close(); //*/

	GetSolution(STParameters, STParameters.FAOD);

	Gnuplot line;
	line.cd(working_directory + "SOL/plot_scripts");

	line("load \"I_block.plt\"");
	line("load \"II_block.plt\"");
	line("load \"DPSYEPS.plt\"");
	
	getchar();

}

void BGA()
{
	using namespace ReverseTask;

	// предварительно: об€зательно надо вызвать;
	StraightTask::Splines::ConfigureSplines();
	BGA_Parameters::SetCoefsForChange();
	//без них BGA не будет знать, что делать.

	// конструктор класса, который хранит параметры генетического алгоритма и пр€мой задачи
	BGA_Parameters RTParameters;

	// конфигурируем скриптовый файл gnuplot
	BGA_Parameters::ConfigurePlotScript();


	std::vector <Generation> individuals(RTParameters.p_0 + 1);
	individuals.shrink_to_fit();

	//StartBGA_1(individuals, RTParameters, BGA_Parameters::ErrorFunctions[0]);
	StartBGA(individuals, RTParameters);
	/*bool trig;
	std::cout << "\n\n do you wish to see how it looks?\n 0. NO \n 1. YEAH \n ans := " << std::endl;
	std::cin >> trig;

	if (trig) GetSol(RTParameters, 3, individuals.front());*/

}

void Stat_BGA(int16_t Exerct_power)
{
	using namespace ReverseTask;

	// предварительно: об€зательно надо вызвать;
	StraightTask::Splines::ConfigureSplines();
	BGA_Parameters::SetCoefsForChange();
	//без них BGA не будет знать, что делать.

	// конструктор класса, который хранит параметры генетического алгоритма и пр€мой задачи
	BGA_Parameters RTParameters;

	// конфигурируем скриптовый файл gnuplot
	BGA_Parameters::ConfigurePlotScript();


	std::vector <Generation> individuals(RTParameters.p_0 + 1);
	individuals.shrink_to_fit();

	StartBGA_For_Stat(Exerct_power, individuals, RTParameters);
}

auto main() -> int
{
	using namespace ReverseTask;
	//BGA_Parameters::SetCoefsForChange(); BGA_Parameters::ConfigurePlotScript();
	//GetSol();
	
	//BGA();
	Stat_BGA(50);

	//std::cout << "\n default F = " << GetCurrentF(); getchar();
	
	//GetSols();
	//getchar();
	return EXIT_SUCCESS;
}