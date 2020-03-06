#pragma once
#include <Methods.h>
#include <GNUplot.h>

namespace ReverseTask
{
	// random value in [a,b]
	double_t random(double_t a, double_t b)
	{
		double_t tmp = (double_t)(rand() % 10001) / 10000.0; // tmp = число из [ 0 , 1 ]
		return a + (b - a)*tmp; // число из [ a , b ]
	}

	double_t GetMax(std::vector<double_t>const & A)
	{
		 double_t tmp = 0;
		 for (uint16_t i = 0; i < A.size(); i++)
			 if (tmp < A[i]) tmp = A[i];
		 return tmp;
	}


	class Generation
	{
	public:

		Generation() {}

			std::vector<double> budgets;

			void bga_info_pin(std::ofstream & outGen, std::vector < std::pair <const char*, double_t* > > const & CFC)
			{
				// возраст, значение F, значение коэффициентов, успехи на доске лидеров
				outGen << "achievements: \t\t" << " : { ";
				for (uint16_t i = 0; i < dashboard_index.size(); i++)
					outGen << dashboard_index[i] << ", ";
				outGen << " }" << "\t";
				outGen << "genu number = " << genu_index << "\n" << F_value << " = ";
				for (uint16_t i = 0; i < budgets.size(); i++)
					outGen << budgets[i] << " + ";
				outGen << "0 \n";

				for (uint16_t i = 0; i < CFC.size(); i++)
					outGen << "\t" << CFC[i].first << " = " << coefs_values[i];
				outGen << std::endl;
			}

			bool SetGenu(uint16_t i)
			{
				if (i == 0)
				{
					std::cout << "\n given i = " << i << " value is inappropriate for genu_index";
					this->genu_index = i + 100;
					return false;
				}
				else if (this->genu_index >= 1)
				{
					//cout << "\n current instance is already have its genu = " << this->genu_index;
					return false;
				}
				this->genu_index = i;
				return true;
			}

			auto Get_F() -> double_t { return this->F_value; }

			void SumUpBuds()
			{
				this->F_value = 0;
				for (size_t i = 0; i < budgets.size(); i++)
					this->F_value += budgets[i];
			}

			auto GetCoefVal(size_t i) -> double_t { return this->coefs_values[i]; }

			void AddCoefVal(double_t a, double_t b)
			{
				this->coefs_values.emplace_back(random(a, b));
				return;
			}

			void AddCoefVal(double_t val)
			{
				this->coefs_values.emplace_back(val);
				return;
			}

			void Recombination(double_t x, double_t y, double_t d)
			{
				double alp = random(-d, 1. + d);
				this->coefs_values.emplace_back(x + alp * (y - x));
			}

			void Mutation(size_t j, double_t L, double_t R, double_t mu)
			{
				bool oper = (int)(random(1., 100.)) % 2; // рандомизация +-
				double_t gamma = random(0., 1.);

				if (oper) coefs_values[j] += mu * (L - R) * pow(2, -16 * gamma);
				else coefs_values[j] -= mu * (L - R) * pow(2, -16 * gamma);

			}

			void AddDbIndex(size_t i)
			{
				(this->dashboard_index).emplace_back(i);
				return;
			}

			void static SwapInVect(std::vector<Generation> & indiv, size_t loser, size_t winner)
			{
				Generation tmp = indiv[loser];
				indiv[loser] = indiv[winner];
				indiv[winner] = tmp;
				return;
			}

			bool static Sort1stFrac(std::vector<Generation> & CurIndivs, uint16_t fraction)
			{
				size_t qmin;
				double_t tmp = 100;
				bool isBestfound = true;
				// цикл поиска минимумов 
				for (size_t i = 1; i <= fraction; i++)
				{
					for (size_t q = i; q < CurIndivs.size(); q++)
					{
						if (CurIndivs[q].Get_F() < tmp)
						{
							tmp = CurIndivs[q].Get_F();
							qmin = q;
						}
					}

					CurIndivs[qmin].AddDbIndex(i); // говорим участнику что он победил
					tmp = 100;

					// если он стоит на своём месте, то не трогаем его
					if (qmin == i)
					{ 
						if (i == 1) isBestfound = false;
						continue;
					}
					// иначе меняем его позицию в контейнере
					else SwapInVect(CurIndivs, i, qmin);

					if (i <= 3) // на экран выводим только первый и последний 
					{
						std::cout 
							<< "\t" << qmin << " ===> " << i << " \t\t F = " << CurIndivs[i].Get_F() << std::endl;
						continue;
					}
					if (i == fraction) 
					{
						std::cout << "\t ... \t ... \n"
							<< "\t" << qmin << " ===> " << i << " \t\t F = " << CurIndivs[i].Get_F() << std::endl;
					}
				}
				return isBestfound;
			}

		~Generation() = default;

	private:

		uint16_t genu_index = 0;
		std::vector <double_t> coefs_values;
		std::vector <size_t> dashboard_index;
		double_t F_value = 0;

	};

	struct BGA_Parameters : StraightTask::Parameters
	{
		BGA_Parameters() : StraightTask::Parameters(1, 1500) // sce, FAOD, GridNum
		{
			p_0 = 7000;
			p = 500;
			SortFraction = (int)((p * 30) / 100);
			RecombPar = 0.1;
			MutPar = 0.0012;
			iterAmount = 200;
			eps = 1.0e-7;

			ControlValues = { NEC, ACU };
			ControlValues.shrink_to_fit();

		}

			void static SetCoefsForChange();

			void outInfo(std::ofstream & outGen)
			{
				outGen << "\n\n"
					<< "full amount of iterations in BGA = " << iterAmount << std::endl
					<< "initial amount of individuals = " << p_0 << std::endl
					<< "regular amount of individuals = " << p << std::endl
					<< "recombination parameter d = " << RecombPar << std::endl
					<< "mutation parameter mu = " << MutPar << std::endl
					<< "fraction of sorted survivers = " << ((double)(SortFraction) / p) * 100.0 << "% = " << SortFraction << std::endl;

				outGen << "\n\n"
					<< "Coefficents for adaption: \n";
				for (uint16_t i = 0; i < CoefsForChange.size(); i++)
				{
					outGen << "\t" << '[' << bonds[0][i] << ", " << bonds[1][i] << ']' << "\t\t"
						<< CoefsForChange[i].first << " ~ " << i << std::endl;
				}

				outGen << "\n"
					<< "Splines to hang on: \n";
				for (uint16_t i = 0; i < ControlValues.size(); i++)
					outGen << "\t" << ControlValues[i] << " ~ " << DerEnumToString(ControlValues[i]) << std::endl;

				outGen << "\n==============================================================================" << std::endl;
			}

			void coutInfo()
			{
				std::cout << "\n\n"
					<< "full amount of iterations in BGA = " << iterAmount << std::endl
					<< "initial amount of individuals = " << p_0 << std::endl
					<< "regular amount of individuals = " << p << std::endl
					<< "recombination parameter d = " << RecombPar << std::endl
					<< "mutation parameter mu = " << MutPar << std::endl
					<< "fraction of sorted survivers = " << ((double)(SortFraction) / p) * 100.0 << "% = " << SortFraction << std::endl;

				std::cout << "\n\n"
					 << "Coefficents for adaption: \n";
				for (uint16_t i = 0; i < CoefsForChange.size(); i++)
				{
					std::cout << "\t" << '[' << bonds[0][i] << ", " << bonds[1][i] << ']' << "\t\t"
						<< CoefsForChange[i].first << " ~ " << i << std::endl;
				}

				std::cout << "\n"
					<< "Splines to hang on: \n";
				for (uint16_t i = 0; i < ControlValues.size(); i++)
					std::cout << "\t" << ControlValues[i] << " ~ " << DerEnumToString(ControlValues[i]) << std::endl;

				std::cout << "\n==============================================================================" << std::endl;
			}
			
			void static ConfigurePlotScript();

			void SetInfoPlot(int16_t ID)
			{
				//line("mu = " + std::to_string(MutPar));
				//line("d = " + std::to_string(RecombPar));
				//line("p0 = " + std::to_string(p_0));
				//line("p = " + std::to_string(p));
				//line("sfract = " + std::to_string(SortFraction));
				
				Progressline.wxt(ID, 1360,768 , 12);

				//line("set label 1 sprintf(\"{/Symbol m} = %3.5f;\", mu) at 1, 225 font \"arialbd,15\"");
				//line("set label 2 sprintf(\"d = %3.5f;\", d) at 1, 213 font \"arialbd,15\"");
				//line("set label 3 sprintf(\"initial population = %d;\", p0) at 12, 225 font \"arialbd,15\"");
				//line("set label 4 sprintf(\"regular population = %d;\", p) at 12, 213 font \"arialbd,15\"");
				//line("set label 5 sprintf(\"survived sorted fraction = %d;\", sfract) at 33, 213 font \"arialbd,15\"");

				Progressline("set key font \"arial,12\"");
				Progressline.cd(plotscripts_dir);
				Progressline("set yrange [0:" + std::to_string(iterAmount) + ']');
			}

			void ProgressPlot(uint16_t iteration, double_t fmin)
			{
				Progressline("fmin = " + std::to_string(fmin));
				Progressline("load \'leaders_evo.plt\'");
			}

			void start_gif_conf(int16_t ID)
			{
				GifLine.gif (ID, 1360, 768, 12, 30);
				GifLine.cd(plotscripts_dir);
				GifLine("set output \"" + working_directory + "BGA/report/current/gnu_output/animated.gif\"");
				GifLine("set yrange [0:" + std::to_string(iterAmount) + ']');
				GifLine("set key font \"arialbd,12\"");
			}

			void collect_in_gif(uint16_t iteration, double_t fmin)
			{
				GifLine("fmin = " + std::to_string(fmin));
				GifLine("load \'leaders_evo.plt\'");
			}

			static void account_investments
			(
				uint16_t i,
				std::vector<double_t> & buds,
				StraightTask::variables & sol,
				std::vector<NumOfDerivative> & CVals,
				std::vector<bool> & ACVals
			)
			{
				// проверяем, можно ли определить невязку
				for (size_t i = 0; i < CVals.size(); i++)
				{
					if (sol.gap_counter[CVals[i]] > StraightTask::Splines::max_gap_amount[CVals[i]]) ACVals[i] = false;
					if (buds[i] > 200)	{ buds[i] = 200; ACVals[i] = false; }
				}
				
				for (size_t i = 0; i < CVals.size(); i++)
				{
					if (!ACVals[i]) continue;
					buds[i] += count_instan_bud[1](*sol.GetMatchedVar(i), sol.GetSplineValue(i));
				}

			}

			static const std::array<std::function<double_t(double_t, double_t)>, 2> count_instan_bud;
			
			static void account_investments2(
				uint16_t i, double_t & bud,
				StraightTask::variables & cor_X
			)
			{
				double_t tmp = (*cor_X.GetMatchedVar(i) - cor_X.GetSplineValue(i));
				bud += tmp * tmp;
			}


		~BGA_Parameters() = default;

		static std::vector < std::pair <const char*, double_t* > > CoefsForChange;
		static std::array<std::vector<double_t>, 2> bonds;
		std::vector< NumOfDerivative > ControlValues;

		uint16_t p_0;  // исходное число участников одной генерации  | 0 < p_0 < 32767
		uint16_t p; // регулярное число популяции | 0 < p <= p_0 < 32767
		uint16_t SortFraction; // доля отсортированных счастливчиков | 0 < SortFraction < 32767
		double_t RecombPar; // параметр рекомбинации
		double_t MutPar; // мутационный параметр
		double_t eps;
		uint16_t iterAmount; // кол-во итераций алгоритма BGA | 0 < iterAmount < 32767

		Gnuplot Progressline; // строка исполнения команд в gnuplot для описания прогресса
		Gnuplot GifLine; // строка исполнения команд в gnuplot для компоновки gif
	};

	std::vector < std::pair < const char*, double_t* > > BGA_Parameters::CoefsForChange = {};
	std::array < std::vector<double_t>, 2 > BGA_Parameters::bonds = {};

	inline const std::array<std::function<double_t(double_t,double_t)>,2> BGA_Parameters::count_instan_bud = {
		[](double_t a, double_t b) -> double_t { return abs(a - b); },
		[](double_t a, double_t b) -> double_t { double_t tmp = (a - b);  return tmp*tmp; }
	};

	void BGA_Parameters::SetCoefsForChange()
	{
		if (CoefsForChange.size() != 0) return; // нельзя допустить повторной инициализации!!!

#define CONF_COEFS(STRUCT, NAME, ID, L, R) BGA_Parameters::\
			CoefsForChange.emplace_back( std::make_pair(NAME, &StraightTask::STRUCT::ID));\
			bonds[0].emplace_back(L); bonds[1].emplace_back(R);
	
		/*1*/		CONF_COEFS(Neurons, "k_N", k_N, 0.03, 5.0)
		/*2*/		CONF_COEFS(Neurons, "k_A", k_A, 1.0, 6.0)
		/*3*/		CONF_COEFS(Neurons, "p_R", p_R, 0.03, 2.0)

		/*3*/		CONF_COEFS(ToxDamage, "p_{ncy}", p_ncy, 14.0, 22.0)
		/*4*/		CONF_COEFS(ToxDamage, "C_{Dcy}", C_Dcy, 6.0, 17.0)
		/*5*/		CONF_COEFS(ToxDamage, "C_{DLn}", C_DLn, 30.0, 60.0)
		/*6*/		CONF_COEFS(ToxDamage, "C_{DLm}", C_DLm, 4.0, 15.0)
		/*7*/		CONF_COEFS(ToxDamage, "P_{nn}", P_nn, 1.0, 7.5)
		/*8*/		CONF_COEFS(ToxDamage, "p_{Lm}", p_Lm, 15.0, 36.0)
		/*9*/		CONF_COEFS(ToxDamage, "p_{Ln}", p_Ln, 20.0, 40.0)
		/*10*/		CONF_COEFS(ToxDamage, "D_0", D_0, 0.03, 3.0)
		/*11*/		CONF_COEFS(ToxDamage, "p_D", p_D, 0.03, 1.0)
#undef pair_up
		bonds[0].shrink_to_fit(); bonds[1].shrink_to_fit();
		CoefsForChange.shrink_to_fit();
	}

	void BGA_Parameters::ConfigurePlotScript()
	{
		double_t R = GetMax(bonds[1]);
		std::ofstream PlotScr
		(working_directory + "BGA/report/current/plot_scripts/leaders_evo.plt");
		{
			PlotScr
				<< " set multiplot \n" << std::endl

				<< "set size 0.75 ,0.96" << std::endl
				<< "set origin 0, 0" << std::endl
				<< "load \'C_evo.plt\'\n" << std::endl

				<< "set size 0.27, 0.986" << std::endl
				<< "set origin 0.73, 0" << std::endl
				<< "set x2tics(sprintf(\"min(F) = %3.8f\", fmin)fmin)" << std::endl
				<< "set xrange[0: fmin + 0.05]" << std::endl
				<< "load \'F_evo.plt\'\n" << std::endl

				<< "unset multiplot" << std::flush;
		}
		PlotScr.close();

		PlotScr.open
		(working_directory + "BGA/report/current/plot_scripts/C_evo.plt", std::ios_base::out);
		{
			PlotScr
				<< "set border 9" << std::endl
				<< "set xtics nomirror" << std::endl
				<< "set nox2tics" << std::endl

				<< "set xrange[-0.5:" << std::to_string(R) << "]" << std::endl
				<< "set grid xtics ytics \n" << std::endl
				<< "plot \"" << working_directory << bga_data_dir << "leaders_C.txt\" \\" << std::endl
				<< "u 1:0 with points ls 1 title \"" << CoefsForChange[0].first << "\", \"\" u 1:0 with lines ls 1 notitle,";
			for (size_t i = 2; i <= CoefsForChange.size(); i++)
			{
				PlotScr
					<< "\\" << std::endl
					<< "\"\" u " << i << ":0 with points ls " << i << " title \"" << CoefsForChange[i - 1].first << "\", "
					<< "\"\" u " << i << ":0 with lines ls " << i << " notitle,";
			}
			PlotScr
				<< "\n unset grid \n" << std::flush;
		}
		PlotScr.close();

		PlotScr.open
		(working_directory + "BGA/report/current/plot_scripts/F_evo.plt", std::ios_base::out);
		{
			PlotScr
				
				<< "set x2tics scale 0 \n" << std::endl

				<< "set border 3" << std::endl
				<< "set ytics nomirror" << std::endl
				<< "set grid xtics ytics\n" << std::endl
				<< "plot \"" << working_directory << bga_data_dir << "leaders_F.txt\" \\" << std::endl
				<< "u 1:0 with points lt 8 pt 7 ps 1 t \"F\", \"\" u 1:0 with lines lt 1 notitle \n" << std::endl

				<< "unset grid" << std::flush;
		}
		PlotScr.close();

	}
}
