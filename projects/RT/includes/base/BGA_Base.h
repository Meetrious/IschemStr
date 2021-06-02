#pragma once
#include <string>

#include <array>
#include <cmath>

#include <thread>
#define MAX_AMOUNT_OF_THREADS 8

#include <stdlib.h>     /* srand, rand */
#include <ctime>       /* time */

#include<base/Timer.h>



/* size of the std::array object that contains pointers to parameters that
* are to variate during BGA iterations */
#define CTV_SIZE 30

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

		// returns pseudorandom value from [a,b] interval
		[[nodiscard]] double_t random(double_t a, double_t b) noexcept {
			double_t tmp = (double_t)(rand() % 10001) / 10000.0; // tmp \in [ 0 , 1 ]
			return a + (b - a) * tmp; // return -> \in [ a , b ]
		}

		class Species {

			/* all fields are independent */

		public:

			// attributes of current individual
			vector<double_t> coefs_values;

			// Aberration value
			double_t F_value = 0;

			bool operator< (Species const& other) {
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

		class Parameters {

			friend class IOs;
			/* all fields are independent except SetCVlist member-function.
			* SetCTVlist is to be defined somewhere further after equation.h inclusion. */

		public:

			Parameters(
				// parameters for ReverseTask
				uint16_t iter_amount = 200,
				uint16_t thread_amount = 1,
				uint16_t init_pop_size = 200,
				uint16_t reg_pop_size = 100,
				uint16_t Sort_Fract = 30,
				uint16_t Recr_Fract = 10,
				double_t recombination_parameter = 0.05,
				double_t mutation_parameter = 0.01,

				// parameters for StraightTask 
				//double_t t0 = 0.5,
				uint32_t ST_gridPow = 1500,
				double_t ST_gapWidth = 24.0,
				uint16_t ST_amountOfGaps = 1)
				:
				amount_of_iterations{ iter_amount },
				amount_of_threads{thread_amount},
				p_0{ init_pop_size },
				p{ reg_pop_size },
				SortedFraction{ Sort_Fract },
				RecruitedFraction{Recr_Fract},
				rc{ recombination_parameter },
				mu{ mutation_parameter },

				N{ ST_gridPow },
				gap_width{ ST_gapWidth },
				full_amount_of_gaps{ ST_amountOfGaps }

			{

				bounds.clear();
				bounds.emplace_back();
				bounds.emplace_back();
				bounds.shrink_to_fit();

				//amount_of_attributes = SetCTVlist();

				Sp = SortedFraction * (p / 100);
				p_rec = RecruitedFraction * (p / 100);

			}

			// parameters for StraightTask

			//double_t t_0; // initial time moment
			uint32_t N; // amount of nods in a grid
			double_t gap_width;
			uint16_t full_amount_of_gaps; // crucial parameter in terms of dealing with delayed(retarded) arguments in equations

			~Parameters() = default;

		protected:

			
			size_t amount_of_attributes;

			// list of coefficients varying in the BGA
			//std::array<StraightTask::IConstant*, CTV_SIZE> CoefsToVariate;
			

			/* list of boundaries within which coefficient values are initially set by RNG :
			* dim = (2 x amount_of_attributes) */
			matrix<double_t> bounds;

			// initial population size | 0 < p_0 < 65535
			uint16_t p_0;

			// regular population size | 0 < p <= p_0 < 65535
			uint16_t p;

			// fraction of selected individuals in the population
			uint16_t Sp;

			// percent of selected population | 0 < SortFraction < 100 < 65535
			uint16_t SortedFraction;

			// amount of random recruited population
			uint16_t p_rec;

			// percent of selected population | 0 < SortFraction < 100 < 65535
			uint16_t RecruitedFraction;

			// BGA iteration amount | 0 < iterAmount < 32767
			uint16_t amount_of_iterations;

			// recombination parameter
			double_t rc;

			// mutation parameter
			double_t mu;

			uint16_t amount_of_threads;

			void DisplayParameters() {
				std::cout << "\n ReverseTask parameters: \n"
					<< "\n\t amount of iterations = " << amount_of_iterations
					<< "\n\t amount of threads = " << amount_of_threads
					<< "\n\t initial population size = " << p_0
					<< "\n\t regular population size = " << p
					<< "\n\t percent of selected population = " << SortedFraction
					<< "\n\t percent of recruited population = " << RecruitedFraction
					<< "\n\t recombination parameter = " << rc
					<< "\n\t mutation parameter = " << mu << std::endl;

				std::cout << "\n\n StraightTask parameters: \n"
					<< "\n\t grid-step size = " << N
					<< "\n\t width of a gap in the step method = " << gap_width
					<< "\n\t full amount of gaps in the step method = " << full_amount_of_gaps
					<< std::endl;
			}

		private:

			uint16_t GetFraction(uint16_t general) { return SortedFraction * (general / 100); }

		};

		class IOs {

			/* up to this point Species class should be defined */

		public:
			std::ofstream CSout, Fout, Statout;
			std::ifstream Cin;

			void WriteBest(Species& Ind)
			{
				CSout.open(output_dir + "/RT/current/best.txt", std::ios_base::out);
				for (auto const& cur : Ind.coefs_values) CSout << cur << std::endl;
				CSout.close();
			}

			void WriteResult(Species& Ind)
			{
				CSout.open(output_dir + "/RT/current/leaders_C.txt", std::ios_base::out | std::ios_base::app);
				for (auto const& cur : Ind.coefs_values) CSout << cur << "\t\t\t";
				CSout << std::endl;
				CSout.close();

				Fout.open(output_dir + "/RT/current/leaders_F.txt", std::ios_base::out | std::ios_base::app);
				Fout << Ind.F_value << std::endl;
				Fout.close();
			}

			void RestartCollector(Parameters const & RT_P,
				std::array<StraightTask::IConstant*, CTV_SIZE> & CoefsToVariate) {

				CSout.open(output_dir + "RT/current/leaders_C.txt", std::ios_base::out);
				if (!CSout) {
					std::cout << "\n !!!stream for coeficients was not allocated for some reason.\n Care to attend!";
					getchar();
					return;
				}

				CSout << '#'
					<< "mu = " << RT_P.mu << "\t"
					<< "d = " << RT_P.rc << "\t"
					<< "N_gen" << RT_P.amount_of_iterations << "\t"
					<< "population : " << RT_P.p_0 << " -> " << RT_P.p << std::endl;

				CSout << "#";
				for (size_t i = 0; i <= RT_P.amount_of_attributes; i++)
					CSout << CoefsToVariate[i]->name << "\t";

				CSout << std::endl;

				CSout.close();

				Fout.open(output_dir + "RT/current/leaders_F.txt", std::ios_base::out);
				if (!Fout) {
					std::cout << "\n !!!stream for F-values was not allocated for some reason.\n Care to attend!";
					getchar();
					return;
				}
				Fout.close();

			}

			void WriteCoefDefBorders(matrix<double_t> const & bounds) {
				CSout.open(output_dir + "/RT/accumul/CoefBorders.txt");
				if (!CSout) {
					std::cout << "\n For some reason I couldn't create CoefBorders.txt file" << std::endl;
					getchar();
					return;
				}
				for (size_t i = 0; i < bounds[0].size(); i++) {
					CSout << bounds[0][i] << "  \t" << bounds[1][i] << std::endl;
				}
				CSout.close();
			}

			void WriteOptimisedCoefs(std::array<StraightTask::IConstant*, CTV_SIZE>& CoefsToVariate, size_t AOA) {
				CSout.open(output_dir + "/RT/accumul/OptCoefNames.txt", std::ios_base::out);
				if (!CSout) {
					std::cout << "\n For some reason I couldn't create OptCoefNames.txt file" << std::endl;
					getchar();
					return;
				}
				for (size_t i = 0; i <= AOA; i++) {
					CSout << CoefsToVariate[i]->name << "\n";
				}
				CSout.close();
			}

			void WriteStatData(Species& best_Ind, const char* name)
			{
				Statout.open(output_dir + "/RT/accumul/" + name + ".txt", std::ios_base::out | std::ios_base::app);
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

			void ConstructMultiplotScript() {

				CSout.open
				(output_dir + "RT/plot_scripts/leaders_evo.plt", std::ios_base::out);
				if (!CSout) {
					std::cout << " \n Failed to create the leaders_evo.plt script in the given folder; "
						<< " \n external_file_allocation_error is thrown";
					external_file_allocation_error exception;
					throw(exception);
				}
				CSout << "set terminal wxt size 1900, 1000" << std::endl
					<< "set multiplot " << std::endl
					<< "set size 0.75, 0.96" << std::endl
					<< "set origin 0, 0 " << std::endl
					<< "load \'C_evo.plt\' " << std::endl

					<< "set size 0.27, 0.986" << std::endl
					<< "set origin 0.73, 0 " << std::endl
					<< "#set x2tics(sprintf(\"min(F) = %3.8f\", fmin)fmin)" << std::endl
					<< "set xrange[0.2:0.6] " << std::endl
					<< "load \'F_evo.plt\' " << std::endl

					<< "unset multiplot" << std::flush;
				CSout.close();
			}
			
			void ConstructCoefEvoPlotScript ( matrix<double_t> &bounds,
				std::array<StraightTask::IConstant*,CTV_SIZE> & CoefsToVariate)
			
			{
				double_t R = 0;

				// looking for the biggest Right bound
				for (auto const & cur: bounds[1]) {
					if (cur > R) R = cur;
				}

				CSout.open (output_dir + "RT/plot_scripts/C_evo.plt", std::ios_base::out);
				if (!CSout.is_open()) {
					std::cout << " \n Failed to create the C_evo.plt script in the given folder; "
						<< " \n external_file_allocation_error is thrown";
					external_file_allocation_error exception;
					throw(exception);
				}

				CSout << "set border 9" << std::endl
					<< "set xtics nomirror \n" << std::endl
					<< "set nox2tics" << std::endl
					<< "set xrange[-0.5:" << R << "]" << std::endl
					<< "set grid xtics ytics \n" << std::endl

					<< "plot \"" << output_dir << "RT/current/leaders_C.txt\" \\" << std::endl
					<< "u 1:0 with points ls 1 title \"" << CoefsToVariate[0]->name << "\","
					<< " \"\" u 1:0 with lines ls 1 title \"" << CoefsToVariate[0]->name << "\",\\" << std::endl;
				for (size_t i = 2; i < bounds[0].size() + 1 ; i++) {
					CSout << "\"\" u " << i << ":0 with points ls " << i << " title \"" << CoefsToVariate[i-1]->name << "\", "
						<< "\"\" u " << i << ":0 with lines ls " << i << " title \"" << CoefsToVariate[i-1]->name << "\",\\" << std::endl;
				}
				CSout.close();
			}
			
			void ConstructAberEvoPlotScript() {
				CSout.open
				(output_dir + "RT/plot_scripts/F_evo.plt", std::ios_base::out);
				if (!CSout) {
					std::cout << " \n Failed to create the F_evo.plt script in the given folder; "
						<< " \n external_file_allocation_error is thrown";
					external_file_allocation_error exception;
					throw(exception);
				}
				CSout << "set x2tics scale 0" << std::endl

					<< "set border 3" << std::endl
					<< "set ytics nomirror" << std::endl
					<< "set grid xtics ytics" << std::endl

					<< "plot \"" << output_dir << "RT/current/leaders_F.txt\" \\" << std::endl
					<< "u 1:0 with points lt 8 pt 7 ps 1 t \"F\", \"\" u 1 : 0 with lines lt 1 notitle" << std::endl

					<< "unset grid" << std::flush;
				CSout.close();
			}

		};


		[[nodiscard]] inline double_t RecombExpr(double_t x, double_t y, double_t alp) noexcept {
			return x + alp * (y - x);
		}

		[[nodiscard]] inline double_t MutExpr(double_t L, double_t R, double_t gamma, double_t mu) noexcept {
			return mu * (R - L) * std::pow(2, -16 * gamma);
		}


		template<typename ST_Method>
		class IStraightTaskSolverThread {
		public:
			IStraightTaskSolverThread() {}

			void PrepairToWork(
				size_t IDN,
				uint32_t N,
				double_t gap_width,
				uint16_t full_amount_of_gaps) {

				m_IDN = IDN;

				STM.Mthd.Set(N, gap_width, full_amount_of_gaps);
				STM.PrepairTheTask();

				F.GatherData();

			}
		
		
			void employ(
				vector<Species>& Population,
				size_t first,
				uint16_t AOT,
				size_t last,
				size_t AOA) {

				thr = std::thread([&]() {
					size_t cind = first + m_IDN;
					while (cind < last) {

						// putting current individual "on the conveyor"
						for (size_t i = 0; i <= AOA; i++) {
							CoefsToVariate[i]->value = Population[cind].coefs_values[i];
						}

						// calculating current Aberration value
						Population[cind].F_value = STM.SolveForBGA(F);
#ifdef _DEBUG	
						std::cout << std::this_thread::get_id() << " -- thread finished with " << m_IDN << " - " << cind << " indiv\n" << std::endl;
#endif	
						cind += AOT;
					}
					}
				);

			} // */
		

			~IStraightTaskSolverThread() { if (thr.joinable()) thr.join(); }

			std::thread thr;
			ST_Method STM;
			IAggregateControls F;
		
			// list of coefficients varying in the BGA
			std::array<StraightTask::IConstant*, CTV_SIZE> CoefsToVariate;
		private:
			size_t m_IDN;
		};

		template <typename ST_Method>
		class Task : public Parameters {

		public:

			// custom constructor
			Task(Parameters CurrentPar_s);

			Species default_member;

			~Task() = default;

			/* void ShowFforDefault() {

				// вносим текущие значения коэффициентов на конвеер
				for (uint16_t j = 0; j <= amount_of_attributes; j++)
					Workers[0].CoefsToVariate[j]->value = default_member.coefs_values[j];

				default_member.F_value = STM.SolveForBGA(F);

				std::cout << "\n F(default_member) = " << default_member.F_value << std::endl;
				getchar();

			}

			void LookOverDecentOnes(uint16_t it) {
				
				double_t F_tmp = 1000.0;
				size_t i_tmp = 0;
				
				while (it != 0) {
					for (size_t cind = 1; cind < p_0; cind++) {

						// вносим текущие значения коэффициентов на конвеер
						for (uint16_t j = 0; j <= amount_of_attributes; j++)
							Workers[0].CoefsToVariate[j]->value = Population[cind].coefs_values[j];

						// calculating current Aberration value
						Population[cind].F_value = STM.SolveForBGA(F);
					}

					for (size_t cind = 1; cind < p_0; cind++) {
						if (Population[cind] < F_tmp) {
							F_tmp = Population[cind].F_value;
							i_tmp = cind;
						}
					}
					std::cout << it << " - nd best:\n"
						<< "index = " << i_tmp << std::endl
						<< "F_value = " << F_tmp << std::endl;

					F_tmp = 1000.0;

					ios.WriteResult(Population[i_tmp]);

					if (it != 1) {
						// Accepting new members to the Population
						for (size_t cind = 1; cind < p_0; cind++) {
							for (size_t j = 0; j <= amount_of_attributes; j++) {
								Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
							}
						}
					}
					std::cout << "\n";
					it--;
				}
				std::cout << "\n Decent were looked over and outputed!";
				getchar();
			} // */

			void SolveForOutput();

			void SolveForStatistics(uint16_t amount_of_shots);

			void OutputBestSolution() {

				std::ifstream in(output_dir + "RT/current/best.txt");
				if (!in) {
					std::cout << "\n GetFucked! \n";
					return;
				}

				// вносим текущие значения коэффициентов на конвеер
				for (uint16_t j = 0; j <= amount_of_attributes; j++)
					Workers[0].CoefsToVariate[j]->value = default_member.coefs_values[j];

				// calculating Aberration value of the default_member
				Population[0].F_value = Workers[0].STM.SolveForBGA(Workers[0].F);//

				for (size_t i = 0; i <= amount_of_attributes; i++) {
					in >> (Workers[0].CoefsToVariate[i]->value);
				}
				in.close();
				std::cout << "\n \t initial relative aberration = " << Population[0].F_value << std::endl;

				Population[1].F_value = Workers[0].STM.SolveForBGA(Workers[0].F);
				std::cout << "\n Current best solution: F = " << Population[1].F_value << std::endl;

				Workers[0].STM.SolveAndOutput(1500, 24.0F, 4);
				getchar();
			}


		private:

			std::array<IStraightTaskSolverThread<ST_Method>, MAX_AMOUNT_OF_THREADS> Workers;

			vector<Species> Population;
					
			IOs ios;

			// defined somewhere further after equation.h inclusion 
			size_t SetCVlist();

			bool Sort1stFrac(bool is_time_to_print);

			void SwapInVect(size_t loser, size_t winner) noexcept;

			void Recombine(Species& Ind)
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

			void Mutate(Species& Ind)
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


			void RecombineCarefully(Species& Ind) {
				size_t prns[] = { (size_t)(random(1.0, Sp)), (size_t)(random(1.0, Sp)) };

				for (size_t j = 0; j <= amount_of_attributes; j++) {
					
					double_t a = Population[prns[0]].coefs_values[j];
					double_t b = Population[prns[1]].coefs_values[j];

					if (a > b) {
						double_t tmp = a;
						a = b;
						b = tmp;
					}

					double_t alp;

					if (a < 0.0)
						alp = random(0.9, 1.0 + rc);
					else
						alp = random(-rc, 1.0 + rc);

					Ind.coefs_values[j] = RecombExpr(a, b, alp);

				}
			}

			void MutateCarefully(Species& Ind)
			{
				for (size_t j = 0; j <= amount_of_attributes; j++)
				{
					bool oper = (int32_t)(random(1.0, 100.0)) % 2; // choice of the sign {+,-} randomisation
					double_t gamma = random(0.0, 1.0);
					double_t deviation = MutExpr(bounds[0][j], bounds[1][j], gamma, mu);

					if (Ind.coefs_values[j] < 0) 
						Ind.coefs_values[j] += deviation;
					else
						if (oper) Ind.coefs_values[j] += deviation;
						else Ind.coefs_values[j] -= deviation;
				}

			}
		};


		//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------
		
		template <typename ST_Method>
		Task<ST_Method>::Task(Parameters CurrentPar_s)
			:Parameters{CurrentPar_s}
		{// RT_Task Constructor definition

			/*// Setting parameters for calc-method: N, \tau, full_amount_of_gaps
			STM.Mthd.Set(
				CurrentPar_s.N,
				CurrentPar_s.gap_width,
				CurrentPar_s.full_amount_of_gaps);// */

			/* setting up the StraightSask for every BGA iteration:
			* initial data setting and presolved data definition */
			// STM.PrepairTheTask(); 
			
			//F.GatherData();
			
			amount_of_attributes = SetCVlist();

			


			for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
				Workers[IDN].PrepairToWork (IDN,
					CurrentPar_s.N,
					CurrentPar_s.gap_width,
					CurrentPar_s.full_amount_of_gaps);
			}


			// initialising default member by current default state of the ODE System in ST
			for (size_t i = 0; i <= amount_of_attributes; i++) {
				default_member.coefs_values.emplace_back(Workers[0].CoefsToVariate[i]->value);
			}

			// setting initial population size
			Population.assign(p_0, default_member);
			Population.shrink_to_fit();



			srand((unsigned int)(time(0))); // запускаем генератор случайных чисел

			// randomising coefs in population generated above
			for (size_t cind = 0; cind < p_0; cind++) {
				for (size_t j = 0; j <= amount_of_attributes; j++) {
					Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
				}
			}

			// this object is ready for the algorithm, kind of...

			ios.ConstructCoefEvoPlotScript(bounds, Workers[0].CoefsToVariate);
			ios.ConstructAberEvoPlotScript();
			ios.ConstructMultiplotScript();

			ios.RestartCollector(*this, Workers[0].CoefsToVariate);

			ios.CSout.open(output_dir + "RT/current/leaders_C.txt", std::ios_base::out);
			if (!ios.CSout) {
				std::cout << "\n !!!stream for coeficients was not allocated for some reason.\n Care to attend!";
				getchar();
				return;
			}
	
			ios.CSout << '#'
				<< "mu = " << mu << "\t"
				<< "d = " << rc << "\t"
				<< "N_gen" << amount_of_iterations << "\t"
				<< "population : " << p_0 << " -> " << p << std::endl;

			ios.CSout << "#";
			for (size_t i = 0; i <= amount_of_attributes; i++)
				ios.CSout << Workers[0].CoefsToVariate[i]->name << "\t";

			ios.CSout << std::endl;

			ios.CSout.close();



			ios.Fout.open(output_dir + "RT/current/leaders_F.txt", std::ios_base::out);
			if (!ios.Fout) {
				std::cout << "\n !!!stream for F-values was not allocated for some reason.\n Care to attend!";
				getchar();
				return;
			}
			ios.Fout.close();

		}

		template <typename ST_Method>
		void Task<ST_Method>::SolveForOutput() {

			DisplayParameters();

			Timer timer;

			// вносим текущие значения коэффициентов на конвеер
			for (uint16_t j = 0; j <= amount_of_attributes; j++)
				Workers[0].CoefsToVariate[j]->value = default_member.coefs_values[j];

			// calculating Aberration value of the default_member
			Population[0].F_value = Workers[0].STM.SolveForBGA(Workers[0].F);// 
						
			timer.ClickEnd();

			auto time = timer.CountInterval();

			std::cout << "\n \t initial relative aberration = " << Population[0].F_value << std::endl;
			/*std::cout << "\t time spent on evaluation ~" << std::chrono::duration_cast<std::chrono::microseconds>(time).count() << " ms" << std::endl
				<< "\t approximate computation time:"
				<< "\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(time * (p_0 + p * (amount_of_iterations - 1)) / amount_of_threads).count() << " ms"
				<< "\t\t" << std::chrono::duration_cast<std::chrono::seconds>(time * (p_0 + p * (amount_of_iterations - 1)) / amount_of_threads).count() << " sec"
				<< "\t\t" << std::chrono::duration_cast<std::chrono::minutes>(time * (p_0 + p * (amount_of_iterations - 1)) / amount_of_threads).count() << " min";
			getchar(); // */

			/* let's process the first fraction of the population;
			* by that I mean calculate F_value (cind = current_individual) */
			

			timer.ClickStart();

			for (size_t IDN = 0; IDN < amount_of_threads; IDN++)	{
				Workers[IDN].employ( Population, 1, amount_of_threads, Sp, amount_of_attributes);
			}
			

			for (size_t IDN = 0; IDN < amount_of_threads; IDN++){
				Workers[IDN].thr.join();
			} // */
			
			Population[1] = Population[0];
			
			// let us embrace the BGA
			for (uint16_t cit = 1; cit <= amount_of_iterations; cit++)
			{
				bool is_time_to_print = (cit % (amount_of_threads * 2) == 0);;

				if (is_time_to_print)
					std::cout << "\n" << cit << '/' << amount_of_iterations
					<< "th iteration of BGA: " << std::endl;
				

				for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
					Workers[IDN].employ( Population, Sp, amount_of_threads,	Population.size(), amount_of_attributes);
				}
				
				for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
					Workers[IDN].thr.join();
					
				}

			

				// 2. Sorting vector container in ascending order of F_values
				Sort1stFrac(is_time_to_print);

				// 4. Recombining
				for (size_t cind = Sp; cind < p - p_rec; cind++) {
					RecombineCarefully(Population[cind]);
					MutateCarefully(Population[cind]);
				}

				// 5. Accepting new members to the Population
				for (size_t cind = p - p_rec; cind <= p; cind++) {
					for (size_t j = 0; j <= amount_of_attributes; j++) {
						Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
					}
				}

				ios.WriteResult(Population[1]);

				if (cit == 1) {
					Population.erase(Population.begin() + p + 1, Population.end());
					Population.shrink_to_fit();
				}//*/
			}

			ios.WriteBest(Population[1]);
			ios.WriteStatData(Population[1], "Raw_stat_N");

			timer.ClickEnd();
			time = timer.CountInterval();

			std::cout << "\n Reversed task is solved for now;" << std::endl;
			std::cout << "Final computation time:\n "
				<< "\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " ms" << std::endl
				<< "\t\t" << std::chrono::duration_cast<std::chrono::seconds>(time).count() << " sec" << std::endl
				<< "\t\t" << std::chrono::duration_cast<std::chrono::minutes>(time).count() << " min" << std::endl;
			getchar();
			return;

		}

		template<typename ST_Method>
		void Task<ST_Method>::SolveForStatistics(uint16_t amount_of_shots)	{

			DisplayParameters();
			std::cout << "\n =================== Stat-collector method =================== \n" << std::endl;

			ios.WriteOptimisedCoefs(Workers[0].CoefsToVariate, amount_of_attributes);

			ios.WriteCoefDefBorders(bounds);

			Timer Shot_timer;
			Timer Cit_timer;
			std::chrono::duration<float>shot_durs(0);
			for (size_t csht = 0; csht < amount_of_shots; csht++) {

				ios.RestartCollector(*this, Workers[0].CoefsToVariate);

				Cit_timer.ClickStart();
				std::cout << std::endl << csht << "/" << amount_of_shots << " was processed: \n" << std::endl;

				for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
					Workers[IDN].employ(Population, 1, amount_of_threads, Sp, amount_of_attributes);
				}


				for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
					Workers[IDN].thr.join();
				} // */

				// let us embrace the BGA
				for (uint16_t cit = 1; cit <= amount_of_iterations; cit++)
				{

					if (cit == 1 || cit == amount_of_iterations)
						std::cout << "\n\t\t" << cit << '/' << amount_of_iterations
						<< "th iteration of BGA: " << std::endl;

					for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
						Workers[IDN].employ(Population, Sp, amount_of_threads, Population.size(), amount_of_attributes);
					}

					for (size_t IDN = 0; IDN < amount_of_threads; IDN++) {
						Workers[IDN].thr.join();

					} 

						// 2. Sorting vector container in ascending order of F_values
					Sort1stFrac(cit == 1 || cit == amount_of_iterations);

					// 4. Recombining
					for (size_t cind = Sp; cind < p - p_rec; cind++) {
						Recombine(Population[cind]);
						Mutate(Population[cind]);
					}

					// 5. Accepting new members to the Population
					for (size_t cind = p - p_rec; cind <= p; cind++) {
						for (size_t j = 0; j <= amount_of_attributes; j++) {
							Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
						}
					}

					ios.WriteResult(Population[1]);

					if (cit == 1) {
						Population.erase(Population.begin() + p + 1, Population.end());
						Population.shrink_to_fit();
					}//*/
				}

				Cit_timer.ClickEnd();

				ios.WriteStatData(Population[1], "stat_N");
				std::cout << std::endl << csht + 1 << "'th launch was processed;\n" << std::endl;

				
				auto time_dur = Cit_timer.CountInterval();
				std::cout << std::endl << "computation time under processed launch:\n "
					<< "\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(time_dur).count() << " ms" << std::endl
					<< "\t\t" << std::chrono::duration_cast<std::chrono::seconds>(time_dur).count() << " sec" << std::endl
					<< "\t\t" << std::chrono::duration_cast<std::chrono::minutes>(time_dur).count() << " min" << std::endl;

				
				shot_durs += time_dur ;
				std::cout << std::endl << "approximate time left: \n"
						<< "\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(shot_durs/(csht + 1) * (amount_of_shots - csht - 1)).count() << " ms" << std::endl
						<< "\t\t" << std::chrono::duration_cast<std::chrono::seconds>(shot_durs / (csht + 1) * (amount_of_shots - csht - 1)).count() << " sec" << std::endl
						<< "\t\t" << std::chrono::duration_cast<std::chrono::minutes>(shot_durs / (csht + 1) * (amount_of_shots - csht - 1)).count() << " min" << std::endl;
				

				// setting initial population size
				Population.clear();
				Population.assign(p_0, default_member);
				Population.shrink_to_fit();

				// randomising coefs in population generated above
				for (size_t cind = 0; cind < p_0; cind++) {
					for (size_t j = 0; j <= amount_of_attributes; j++) {
						Population[cind].RandomiseCoef(j, bounds[0][j], bounds[1][j]);
					}
				}
			}

			Shot_timer.ClickEnd();
			auto time_dur = Shot_timer.CountInterval();
			std::cout << "\n -=-=-=-=-=-=-=-=-=-=-=-=-= Statistics on ReverseTask was gathered =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n\n \t "
				<< amount_of_shots << " launches of BGA were processed \n" << std::endl;
			std::cout << "Final computation time:\n "
				<< "\t\t" << std::chrono::duration_cast<std::chrono::milliseconds>(time_dur).count() << " ms" << std::endl
				<< "\t\t" << std::chrono::duration_cast<std::chrono::seconds>(time_dur).count() << " sec" << std::endl
				<< "\t\t" << std::chrono::duration_cast<std::chrono::minutes>(time_dur).count() << " min" << std::endl;
			std::cout << std::flush;
			getchar();
		}

		//====================================================================================================================================================

		template <typename ST_Method>
		bool Task<ST_Method>::Sort1stFrac(bool is_time_to_print) {

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

				/*// if min_ind is in the right place... 
				if (min_ind == i) {

					// and if Population[min_ind] is the leader of the previous iteration,
					// then new best solution was not found 
					if (i == 1) isBestfound = false;

					// ...then it stays right where it is. And we go on searching next best individual
					continue;
				}

				// otherwise we change its position in the current "leaderboard"
				else SwapInVect(i, min_ind); // */

				if(min_ind != i) SwapInVect(i, min_ind);

				/* outputting current state of the leaderboard with regards to the first-third and SP'th leaders
				* within current population // */
				if (is_time_to_print) {
					if (i <= 3) {
						std::cout << "\t" << min_ind << " ===> " << i << " \t\t F = " << Population[i].F_value << std::endl;
						continue;
					}

					if (i == Sp) {
						std::cout << "\t ... \t ... \n"	<< "\t" << min_ind << " ===> " << i << " \t\t F = " << Population[i].F_value << std::endl;
					} // */
				}
			}

			return isBestfound;
		}
		 
		template <typename ST_Method>
		inline void Task<ST_Method>::SwapInVect(size_t loser, size_t winner) noexcept {
			Species tmp = Population[loser];
			Population[loser] = Population[winner];
			Population[winner] = tmp;
		}

	}
}