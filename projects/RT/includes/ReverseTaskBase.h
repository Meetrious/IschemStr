#pragma once
#include <vector>
#include <fstream>

#include <Settings.h>
#include <DirPaths.h>
#include <GNUplot.h>
//#include <Methods.h>


namespace ReverseTask
{
	template <typename T>
	using vector = std::vector<T>;

	template <typename T>
	using matrix = vector<vector<T>>;



	// random value in [a,b]
	double_t random(double_t a, double_t b)
	{
		double_t tmp = (double_t)(rand() % 10001) / 10000.0; // tmp = число из [ 0 , 1 ]
		return a + (b - a) * tmp; // число из [ a , b ]
	}


	double_t GetMax(vector<double_t>const& A)
	{
		double_t tmp = 0;
		for (uint16_t i = 0; i < A.size(); i++)
			if (tmp < A[i]) tmp = A[i];
		return tmp;
	}

	namespace BGA
	{

		class Generation
		{
			class IOs friend;
		public:

			Generation() {}

			vector<double_t> budgets;

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

			void Recombine(double_t x, double_t y, double_t d)
			{
				double alp = random(-d, 1. + d);
				this->coefs_values.emplace_back(x + alp * (y - x));
			}

			void Mutate(size_t j, double_t L, double_t R, double_t mu)
			{
				bool oper = (int)(random(1., 100.)) % 2; // рандомизация +-
				double_t gamma = random(0., 1.);

				if (oper) coefs_values[j] += mu * (L - R) * pow(2, -16 * gamma);
				else coefs_values[j] -= mu * (L - R) * pow(2, -16 * gamma);

			}

			void static SwapInVect(std::vector<Generation>& indiv, size_t loser, size_t winner)
			{
				Generation tmp = indiv[loser];
				indiv[loser] = indiv[winner];
				indiv[winner] = tmp;
				return;
			}

			bool static Sort1stFrac(std::vector<Generation>& CurIndivs, uint16_t fraction)
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
			vector<double_t> coefs_values;
			double_t F_value = 0;

		};

		class IOs 
		{
		public:
			std::ofstream Cout, Fout, Statout;
			std::ifstream Cin;

			void WriteBest(Generation& indiv)
			{
				Cout.open(output_dir + "/RT/current/best.txt", std::ios_base::out);
				for (auto const& cur : indiv.coefs_values) Cout << cur << std::endl;
				Cout.close();
			}

			void WriteResult(Generation& indiv)
			{
				Cout.open(output_dir + "/RT/current/leaders_C.txt", std::ios_base::out | std::ios_base::app);
				for (auto const& cur : indiv.coefs_values) Cout << cur << "\t\t\t";
				Cout << std::endl;
				Cout.close();

				Fout.open(output_dir + "/RT/current/leaders_F.txt", std::ios_base::out | std::ios_base::app);
				Fout << indiv.F_value << std::endl;
				Fout.close();
			}

			void WriteStatData(Generation& best_indiv)
			{
				Statout.open(output_dir + "/RT/accumul/stat_N.txt", std::ios_base::out | std::ios_base::app);
				for (auto const& cur : best_indiv.coefs_values)
					Statout << cur << "\t\t\t";
				Statout << best_indiv.F_value << std::endl;
				Statout.close();
			}

			void ReadBest(Generation& indiv)
			{
				indiv.coefs_values.clear();
				Cin.open(output_dir + "RT/current/best.txt");
				if (!Cin)
				{
					std::cout << "\n stop this maddness right now, and give me proper direction for best RT-solution data!\n";
					return;
				}
				while (!Cin.eof())
				{
					double_t tmp; Cin >> tmp;
					indiv.coefs_values.emplace_back(tmp);
				}
				Cin.close();
			}
			
			
		};

		class Parameters {
		public:


			IOs ios;

			// исходное число участников одной генерации  | 0 < p_0 < 32767
			uint16_t p_0 = 5000;

			// регулярное число популяции | 0 < p <= p_0 < 32767
			uint16_t p = 500; uint16_t Sp;
			
			
			double_t eps = 1.0e-7;
			uint16_t iterAmount = 200; // кол-во итераций алгоритма BGA | 0 < iterAmount < 32767
			
			double_t Recombination = 0.05; // параметр рекомбинации
			double_t Mutation = 0.01; // мутационный параметр
			

			std::vector< std::pair <const char*, double_t * > > CoefsToVariate;

			matrix<double_t> appropriate_boundaries;

			vector<vector<double_t>> ControlValues {CY};
	
			Parameters()
			{
				ControlValues.shrink_to_fit();
				bonds.clear(); bonds.emplace_back(); bonds.emplace_back();
				bonds.shrink_to_fit();
			
				Sp = GetFraction(p);
			}
			

			~Parameters() = default;

		private:
			uint16_t GetFraction(uint16_t general) { return (uint16_t)(SortedFraction * general / 100); }
			uint16_t SortedFraction = 30; // процент отсортированных счастливчиков | 0 < SortFraction < 100 < 
		};
	}
}
