#pragma once
#include <string>

#include <array>
#include <utility>
#include <cmath>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */



#include <base/Solver_base.h>

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

namespace ReverseTask
{
	

	namespace BGA
	{

		// returns pseudorandom value in [a,b]
		[[nodiscard]] double_t random(double_t a, double_t b) noexcept {
			double_t tmp = (double_t)(rand() % 10001) / 10000.0; // tmp \in [ 0 , 1 ]
			return a + (b - a) * tmp; // return -> \in [ a , b ]
		}

		class Species{

			/* all fields are independent */

		public:
			

			// attributes of current individual
			vector<double_t> coefs_values;

			// Aberration value
			double_t F_value = 0;

			bool operator< (Species const& other){
				return this->F_value < other.F_value;
			}

			bool operator< (double_t const& other) {
				return this->F_value < other;
			}

			// add random_number_generated coeficient
			void AddRNGCoef(double_t a, double_t b)
			{
				this->coefs_values.emplace_back(random(a, b));
				return;
			}

			void RandomiseCoef(size_t j, double_t const left, double_t const right) {
				coefs_values.at(j) = random(left, right);
			}

		};
	
		class Parameters{

		/* all fields are independent except SetCTVlist member-function.
		* SetCTVlist is to be defined somewhere further after equation.h inclusion. */

		public:

			Parameters()
				:p_0{ 200 }, p{ 100 }, amount_of_iterations{ 200 }, rc{ 0.05 }, mu{ 0.01 }
			{}

			Parameters(uint16_t iter_amount, uint16_t init_pop_size, uint16_t reg_pop_size, uint16_t Sort_Fract)
				: amount_of_iterations{ iter_amount }, p_0{ init_pop_size }, p{ reg_pop_size }, SortedFraction{ Sort_Fract }{
				
				bounds.clear();
				bounds.emplace_back();
				bounds.emplace_back();
				bounds.shrink_to_fit();

				amount_of_attributes = SetCTVlist();

				Sp = SortedFraction * (p / 100);
				p_rec = 10 * (p / 100);

			}

			~Parameters() = default;

		protected:

			size_t amount_of_attributes;

			// list of coefficients varying in the BGA
			std::array<std::pair<const char*, double_t*>, 30> CoefsToVariate;

			/* list of boundaries within which coefficient values are initially set by RNG :
			* dim = (2 x amount_of_attributes) */
			matrix<double_t> bounds;

			// initial population size | 0 < p_0 < 32767
			uint16_t p_0 = 5000;

			// regular population size | 0 < p <= p_0 < 32767
			uint16_t p = 500;

			// fraction of 
			uint16_t Sp;
			uint16_t p_rec;

			// BGA iteration amount | 0 < iterAmount < 32767
			uint16_t amount_of_iterations = 200; 

			// recombination parameter
			double_t rc = 0.05; 

			// mutation parameter
			double_t mu = 0.01;

			
			double_t eps = 1.0e-7;

		private:

			// defined somewhere further after equation.h inclusion 
			size_t SetCTVlist();

			uint16_t GetFraction(uint16_t general) { return SortedFraction * (general / 100); }

			// percent of lucky population | 0 < SortFraction < 100 < 65535
			uint16_t SortedFraction = 30; 

		};

		class IOs {

		/* up to this point Species class should be defined */

		public:
			std::ofstream Cout, Fout, Statout;
			std::ifstream Cin;

			void WriteBest(Species& Ind)
			{
				Cout.open(output_dir + "/RT/current/best.txt", std::ios_base::out);
				for (auto const& cur : Ind.coefs_values) Cout << cur << std::endl;
				Cout.close();
			}

			void WriteResult(Species& Ind)
			{
				Cout.open(output_dir + "/RT/current/leaders_C.txt", std::ios_base::out | std::ios_base::app);
				for (auto const& cur : Ind.coefs_values) Cout << cur << "\t\t\t";
				Cout << std::endl;
				Cout.close();

				Fout.open(output_dir + "/RT/current/leaders_F.txt", std::ios_base::out | std::ios_base::app);
				Fout << Ind.F_value << std::endl;
				Fout.close();
			}

			void WriteStatData(Species& best_Ind)
			{
				Statout.open(output_dir + "/RT/accumul/stat_N.txt", std::ios_base::out | std::ios_base::app);
				for (auto const& cur : best_Ind.coefs_values)
					Statout << cur << "\t\t\t";
				Statout << best_Ind.F_value << std::endl;
				Statout.close();
			}

			void ReadBest(Species& Ind)
			{
				Ind.coefs_values.clear();
				Cin.open(output_dir + "RT/current/best.txt");
				if (!Cin)
				{
					std::cout << "\n stop this maddness right now, and give me proper direction for best RT-solution data!\n";
					return;
				}
				while (!Cin.eof())
				{
					double_t tmp; Cin >> tmp;
					Ind.coefs_values.emplace_back(tmp);
				}
				Cin.close();
			}


		};


		[[nodiscard]] inline double_t RecombExpr(double_t x, double_t y, double_t alp) noexcept {
			return x + alp * (y - x);
		}

		[[nodiscard]] inline double_t MutExpr(double_t L, double_t R, double_t gamma, double_t mu) noexcept {
			return mu * (L - R) * std::pow(2, -16 * gamma);
		}

		template <typename ST_Method>
		class ISolver: public Parameters{

		public:

			ISolver() {

				// Setting parameters for calc-method: N, \tau, full_amount_of_gaps
				STM.Set(1500, 24.0, 1);

				/* setting up the StraightSask for every BGA iteration:
				* initial data setting and presolved data definition */
				STM.PrepairTheTask();

				// initialising default member by current default state of the ODE System in ST
				for (size_t i = 0; i <= amount_of_attributes; i++) {
					default_member.coefs_values.emplace_back(*CoefsToVariate[i].second);
				}

			
			}

			ISolver(uint16_t iter_amount, uint16_t init_pop_size, uint16_t reg_pop_size, uint16_t Sort_Fract) 
				:Parameters(iter_amount, init_pop_size, reg_pop_size, Sort_Fract)
			{

				// Setting parameters for calc-method: N, \tau, full_amount_of_gaps
				STM.Mthd.Set(1500, 24.0F, 4);

				/* setting up the StraightSask for every BGA iteration:
				* initial data setting and presolved data definition */
				STM.PrepairTheTask();

				// initialising default member by current default state of the ODE System in ST
				for (size_t i = 0; i <= amount_of_attributes; i++) {
					default_member.coefs_values.emplace_back(*CoefsToVariate[i].second);
				}

				// setting initial population size
				Population.assign(p_0, default_member); 
				Population.shrink_to_fit();

				
				F.GatherData();

				srand((unsigned int)(time(0))); // запускаем генератор случайных чисел

				// randomising coefs in population generated above
				for (size_t cind = 0; cind < p_0; cind++) {
					for (size_t j = 0; j <= amount_of_attributes; j++) {
						Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
					}
				}

				// this object is ready for the algorithm, kind of...


				ios.Cout.open(output_dir + "RT/current/leaders_C.txt", std::ios_base::out);
				if (!ios.Cout) {
					std::cout << "\n !!!stream for coeficients was not allocated for some reason.\n Care to attend!";
					getchar();
					return;
				}

				ios.Cout << '#'
					<< "mu = " << mu << "\t"
					<< "d = " << rc << "\t"
					<< "N_gen" << amount_of_iterations << "\t"
					<< "population : " << p_0 << " -> " << p << std::endl;

				ios.Cout << "#";
				for (size_t i = 0; i <= amount_of_attributes; i++)
					ios.Cout << CoefsToVariate[i].first << "\t";
				
				ios.Cout << std::endl;

				ios.Cout.close();



				ios.Fout.open(output_dir + "RT/current/leaders_F.txt", std::ios_base::out);
				if (!ios.Fout) {
					std::cout << "\n !!!stream for F-values was not allocated for some reason.\n Care to attend!";
					getchar();
					return;
				}
				ios.Fout.close();

			}
			
			Species default_member;

			~ISolver() = default;

			
			void SolveForOutput() {
		
				/* let's process the first fraction of the population;
				* by that I mean calculate F_value (cind = current_individual) */
				for (size_t cind = 1; cind < Sp; cind++){

					// вносим текущие значения коэффициентов на конвеер
					for (uint16_t j = 0; j <= amount_of_attributes; j++)
						(*CoefsToVariate[j].second) = Population[cind].coefs_values[j];

					// calculating current Aberration value
					Population[cind].F_value = STM.SolveForBGA(F); //1111111111111111111111

				}

				// let us embrace the BGA
				for (uint16_t cit = 1; cit <= amount_of_iterations; cit++)
				{
					std::cout << "\n" << cit << '/' << amount_of_iterations 
						<< "th iteration of BGA: " << std::endl;

					// 1. Recalculating F_values for current population
					for (size_t cind = Sp; cind < Population.size(); cind++)
					{
						
						// putting current individual "on the conveyor"
						for (uint16_t i = 0; i <= amount_of_attributes; i++)
							(*CoefsToVariate[i].second) = Population[cind].coefs_values[i];

						// calculating current Aberration value
						Population[cind].F_value = STM.SolveForBGA(F);

#ifdef _DEBUG
						// outputs current progress
						CoutProgress(cind, Population.size());
#endif
					}

					// 2. Sorting vector container in ascending order of F_values
					Sort1stFrac();

					// 4. Recombining
					for (size_t cind = Sp; cind < p-p_rec; cind++) {
						Recombine(Population[cind]);
						Mutate(Population[cind]);
					}

					// 5. Accepting new members to the Population
					for (size_t cind = p-p_rec; cind <= p; cind++) {
						for (size_t j = 0; j <= amount_of_attributes; j++) {
							Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
						}
					}

					ios.WriteResult(Population[1]);

					if (cit == 1) {
						Population.erase(Population.begin() + p + 1, Population.end());
					}//*/
				}

				ios.WriteBest(Population[1]);
				ios.WriteStatData(Population[1]);


				std::cout << "\n Reversed task is solved for now \n";

				getchar();
				return;

			}

			void SolveForStatistics();
			
			void OutputBestSolution() {

				std::ifstream in(output_dir + "RT/current/best.txt");
				if (!in) {
					std::cout << "\n GetFucked! \n";
					return;
				}

				for (size_t i = 0; i <= amount_of_attributes; i++) {
					in >> (*CoefsToVariate[i].second);
				}
				in.close();

				STM.SolveAndOutput(1500, 24.0, 4);
			}

		private:

			vector<Species> Population;

			IOs ios;

			ST_Method STM;

			IAggregateControls F;

			bool Sort1stFrac() {

				size_t min_ind;
				double_t tmp = 1000;
				

				bool isBestfound = true; 

				// cycle to sort in ascending order best individials
				for (size_t i = 1; i <= Sp; i++) {

					// cycle to find best individual in the population beginning with (i)-th individual
					for (size_t cind = i; cind < Population.size(); cind++) {

						if (Population[cind] < tmp)
						{
							tmp = Population[cind].F_value;
							min_ind = cind;
						}

					} // min_ind is now an index of the individual whose F_value is to be next in the accending order after (i-1)-th individual

					tmp = 1000;// reseting buffer-variable

					// if min_ind is in the right place... 
					if (min_ind == i) {

						/* and if Population[min_ind] is the leader of the previous iteration, 
						* then new best solution was not found */
						if (i == 1) isBestfound = false;

						// ...then it stays right where it is. And we go on searching next best individual
						continue;
					}

					// otherwise we change its position in the current "leaderboard"
					else SwapInVect(i, min_ind);


					/* outputting current state of the leaderboard with regards to the first-third and SP'th leaders
					* within current population // */
					if (i <= 3) {
						std::cout
							<< "\t" << min_ind << " ===> " << i << " \t\t F = " << Population[i].F_value << std::endl;
						continue;
					}

					if (i == Sp) {
						std::cout << "\t ... \t ... \n"
							<< "\t" << min_ind << " ===> " << i << " \t\t F = " << Population[i].F_value << std::endl;
					}
				}

				return isBestfound;
			}

			inline void SwapInVect(size_t loser, size_t winner) noexcept {
				Species tmp = Population[loser];
				Population[loser] = Population[winner];
				Population[winner] = tmp;
			}

			void Recombine(Species & Ind)
			{
				size_t prns[] = { (size_t)(random(1.0, Sp)), (size_t)(random(1.0, Sp)) };

				for (size_t j = 0; j <= amount_of_attributes; j++)
				{
					double_t alp = random(-rc, 1.0 + rc);

					Ind.coefs_values[j] = RecombExpr(
						Population[prns[0]].coefs_values[j],
						Population[prns[1]].coefs_values[j], alp);
				}
			}

			void Mutate(Species & Ind)
			{
				for (size_t j = 0; j <= amount_of_attributes; j++)
				{
					bool oper = (int32_t)(random(1.0, 100.0)) % 2; // choice of the sign {+,-} randomisation
					double_t gamma = random(0.0, 1.0);
					double_t deviation = MutExpr(bounds[0][j], bounds[1][j], gamma, mu);

					if (oper) Ind.coefs_values[j] += deviation;
					else Ind.coefs_values[j] -= deviation;
				}

			}
		};

	}
}