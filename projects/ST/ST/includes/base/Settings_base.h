#pragma once
#include <base/DirPaths.h>
#include <base/Methods.h>
#include <base/SplRealisation.h>

#include <iomanip>


namespace StraightTask
{

	class IErr{
	
	public:
		matrix<double_t> ErrorArray;
		vector<double_t> Error;
		size_t min_row;

		void SaveError(double_t Tj, double_t current_error) {
			ErrorArray.push_back({ Tj, current_error });
		}

		size_t GetOrdOfMaxError() {
			Error = { 0.0, 0.0 };
			for (size_t i = 0; i < ErrorArray.size(); i++) {
				if (Error[1] < ErrorArray[i][1]){
					Error = ErrorArray[i];
					min_row = i;
				}
			}

			return min_row;
		}


		static double_t GetDif(double_t Calc, double_t Exact, uint16_t power){	return std::pow(abs(Calc - Exact),power); }

		void Sort(size_t n){
			if (n > ErrorArray.size()) n = ErrorArray.size();

			vector<double_t> tmp = { 0.0, 0.0 };
			size_t k = 0;

			for (size_t i = 0; i < ErrorArray.size(); i++)
			{
				tmp[1] = 0.0;
				for (size_t j = i; j < ErrorArray.size(); j++)
				{
					if (tmp[1] < ErrorArray[j][1])
					{
						tmp = ErrorArray[j];
						k = j;
					}					
				}
				ErrorArray[k] = ErrorArray[i];
				ErrorArray[i] = tmp;
				if (i > n) 
					break;
			}
		}

	};

	class IRet{
	public:
		IRet(){}

		void SetParameters(Methods::Parameters ST) {
			N = ST.N;
			t_0 = ST.t_0;
			gap_width = ST.gap_width;
		}

		double_t GetRetValue(double_t tau, uint32_t Nj) {
			return data[CountRetIndex(tau, Nj)]; 
		}

		void ShiftRets(uint32_t Nj, double_t prev_U) { data[Nj - 1] = prev_U; }

		virtual void StateRetArray(uint32_t ST_N, double_t ST_t0, double_t ST_gap) {
			N = ST_N; t_0 = ST_t0; gap_width = ST_gap;
			data.clear();
			for (uint32_t i = 0; i < N; i++) { data.emplace_back(0.0); }
		}

	protected:
		uint32_t N=-1;
		double_t t_0=0;
		double_t gap_width;

		vector<double_t> data;
	private:
		[[nodiscard]] size_t CountRetIndex(double_t tau, size_t Nj) {
			return (Nj + (size_t)( (N * (10.0 * gap_width - (tau - t_0) * 10.0) / (10.0 * gap_width) ))) % N;
		}
		
	};

	class ISpline
	{
	private:
		[[nodiscard]] inline double_t GetValue(double_t t, double_t ti, double_t ai, double_t bi, double_t ci, double_t di) {
			return ai + bi * (t - ti) + ci * std::pow(t - ti, 2) + di * std::pow(t - ti, 3);
		}
		[[nodiscard]] inline double_t GetSplValue(double_t current_time, uint32_t gap_counter) // нерабочий вариант метода!, доработать
		{
			return GetValue(current_time,
				SplineData[gap_counter][0],
				SplineData[gap_counter][1],
				SplineData[gap_counter][2],
				SplineData[gap_counter][3],
				SplineData[gap_counter][4]);
		}
	public:

		bool is_exist;
		bool is_triggered;
		matrix<double_t> ExperimentData;
		matrix<double_t> SplineData;

		size_t max_gap_amount = 0;

		virtual const char* ini_spl_name() { return nullptr; };
		const std::string spl_dir() { return input_dir + "exp/" + ini_spl_name() + ".txt"; }

		void CollectExpData() {
			ExperimentData.clear();

			std::ifstream in(spl_dir());
			if (!in)
			{
				std::cout << "\n something went wrong with exp_data collection for " << ini_spl_name()
					<< ". Care to attend!" << std::endl;
				ExperimentData.push_back({ 0, 0 });
				return;
			}
			double_t tmp[2];
			for (uint16_t i = 0; !in.eof(); i++)
			{
				ExperimentData.emplace_back();
				in >> tmp[0]; in >> tmp[1];
				ExperimentData[i].emplace_back(tmp[0]);
				ExperimentData[i].emplace_back(tmp[1]);
			} // сюда надо добавить определение max_gap_counter
			in.close();
		}

		void CheckShiftInterpGap(double_t current_time, uint32_t& gap_counter) {
			if (!is_triggered) return;
			// Вычисляем промежуток текущей интерполяции для каждой величины

			if (current_time - SplineData[gap_counter][0] > 0)
			{ //  попадение сюда знаменует, что t переменная вышла за пределы промежутка интерполирования текущим кубическим сплайном
				// вывод: надо менять 3-многочлен. 
				// Проверяем, на последнем ли мы промежутке интерполяции для данной i - той величины
				if (gap_counter + 1 == max_gap_amount) // надо k_e + 1
				{ // i-тая величина оказалась на последнем промежутке интерполяции и t вышла в промежуток экстраполяции 
					is_triggered = false; // выключаем передачу значения интерполянта для i-той величины
					return; //continue; // выходим из всех if блоков и переходим на проверку сделующей (i+1)-й величины
				}

				gap_counter++; // обновляем текущий промежуток интерполяции (инкремент счетчика k_e[i]) для текущей i-той величины
			}

			return;
		}

		void ConfigureSpline(bool(*method)(vector<vector<double_t>> const&, matrix<double_t>&)) {
			is_triggered = method(ExperimentData, SplineData);
			this->SplineData.shrink_to_fit();

			if (max_gap_amount < SplineData.size()) max_gap_amount = SplineData.size();
		}



		void OutputSpline(const std::string name, double_t H = 0.016)
		{
			if (!is_triggered)
			{
				std::cout << " \n\t there are no splines for you at " << name
					<< " <- this vatiable. Try another one";
				return;
			}

			std::ofstream out(output_dir + "SPL/SPL_" + name + ".txt");
			size_t k = 1;
			double_t t = SplineData[0][0];
			while (k < max_gap_amount)
			{
				while (t < SplineData[k][0])
				{
					out << t << "\t\t\t"
						<< GetValue(t, SplineData[k][0], SplineData[k][1], SplineData[k][2], SplineData[k][3], SplineData[k][4])
						<< "\n";
					t += H; // диаметр разбиения по умолчанию
				}
				k++;
			}
			out.close();
		}
	public:
		bool QueueRule(double_t& val, double_t tj, uint32_t gap_counter)
		{

			// если сплайн не выключен, вызываем численный метод
			// вопрос: есть ли возможность взаимозаменять вызовы численного метода на вызовы решения?

			// проверить текущий промежуток интерполяции:
			// если кап текущего промежутка достигнут, то сдвигаем gap_counter на +1.
			//		если кап по промежуткам достигнут, то возвращаем false.
			//		выпиливать текущий вызов нет нужды. Мало ли в обратной задаче понадобится
			// возвращаем значение сплайна

			if (!is_triggered) return false;

			// 
			val = GetSplValue(tj, gap_counter); return true;
		}

	};

	class IOs
	{
	protected:
		std::ofstream OutSolutionStream;
		std::ofstream OutBudgetStream;
		std::ofstream OutErrorStream;

		matrix<double_t> PreSavedSolData; // 2-dim table for presaved solution

		// a method that outputs directory of a *.txt file with presaved solution
		const std::string presol_data_dir()noexcept { return input_dir + "preserved_solution/" + sol_name() + ".txt"; }

		// a method that outputs a name of a ODE_System member
		virtual const char* sol_name() = 0;

	public:

		// a method that outputs solution in a proper form 
		void OutputSol(uint32_t Nj, double_t Tj, double_t sol) {
			OutSolutionStream << Nj << "\t\t\t" << std::setprecision(5) << Tj << "\t\t\t";
			OutSolutionStream << std::setprecision(15) << sol << std::endl;
		}

		void OutputError(double_t Tj, double_t Error_value) noexcept {
			OutErrorStream << std::setprecision(5) << Tj << "\t\t\t";
			OutErrorStream << std::setprecision(15) << Error_value << std::endl;
		}

		virtual void OutputBuds(uint32_t Nj, double_t Tj) = 0;

		void AllocateOutputStreams()
		{
			OutSolutionStream.open(output_dir + "ST/SOL/UNI/" + sol_name() + ".txt");
			if (!OutSolutionStream)
			{
				std::cerr << "\n For some reason solution output stream for <" << sol_name() << "> wasn't allocated. \n";
				std::cerr << "\t Imposible to output solution.";
				getchar();
				throw("Output stream for solution was not allocated");
			}
			OutBudgetStream.open(output_dir + "ST/SOL/BUDS/" + sol_name() + ".txt");
			if (!OutBudgetStream)
			{
				std::cerr << "\n For some reason budgets output stream for <" << sol_name() << "> wasn't allocated \n";
				std::cerr << "\t Imposible to output budgets. \n";
				
				// closing opened stream
				OutSolutionStream.close(); // closing opened stream
				
				std::cout << "Opened streams were closed. \n";

				getchar();
				throw("Stream for budgets was not allocated");
			}
			OutErrorStream.open(output_dir + "ST/SOL/ERR/" + sol_name() + ".txt");
			if (!OutErrorStream)
			{
				std::cerr << "\n For some reason Error-output stream for <" << sol_name() << "> wasn't allocated \n";
				std::cerr << "\t Imposible to output calculation errors.";
				
				// closing opened stream
				OutSolutionStream.close();
				OutErrorStream.close();

				std::cout << "Opened streams were closed. \n";

				getchar();
				throw("Stream for Error-output was not allocated");
			}
			// in case of an error this closes only those streams, that related to this object of IOs type.
			// I suggest you to figure out how to close all of those streams that were already allocated higher in the list OutStreamAllocator!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}

		void CollectSolData()
		{
			std::ifstream in(presol_data_dir());
			if (!in)
			{
				std::cerr << "\n !!! For some reason file with presolved_data was not allocated to the in_stream on direction:" << presol_data_dir();
				throw("\n File with presolved_data was not allocated to the in_stream");
			}
			double_t tmp[3];

			while (!in.eof())
			{
				in >> tmp[0]; in >> tmp[1]; in >> tmp[2];
				PreSavedSolData.push_back({ tmp[0], tmp[2] });// */
				/*in >> tmp[0]; in >> tmp[1];
				PreSavedSolData.push_back({ tmp[0], tmp[1] }); // */
			}
			in.close();
			PreSavedSolData.shrink_to_fit();
		}

		double_t GetSolData(uint16_t cur_day, uint32_t Nj, uint32_t N){
			return PreSavedSolData[Nj + cur_day * N][1]; 
		}

		~IOs() {
			if (OutSolutionStream.is_open()) OutSolutionStream.close();
			if (OutBudgetStream.is_open()) OutBudgetStream.close();
			if (OutErrorStream.is_open()) OutErrorStream.close();
		}

	};


	template<typename RightPart>
	class IMember : public IOs, public ISpline, public IRet {
	
	public:

		RightPart RP; // equation of a current member of an ODE system

		// a method that outputs budget contributions
		void OutputBuds(uint32_t Nj, double_t Tj) final {
			OutBudgetStream << Nj << "\t\t\t" << std::setprecision(5) << Tj << "\t\t\t";
			OutBudgetStream << std::setprecision(15);
			for (uint16_t i = 0; i < RP.B.size(); i++) OutBudgetStream << RP.B[i] << "\t\t\t";
			OutBudgetStream << std::endl;
		}

		void GetInitialDataFromOutside(double_t & t0, double_t & val)
		{
			std::ifstream in(input_dir + "list of initials/per_value/" + sol_name() + ".txt");
			if (!in)
			{
				std::cerr
					<< "\n For some reason file with initial data for " << sol_name()
					<< " was not allocated to the in_stream, so initial data will be defined by default \n"
					<< std::endl;
				getchar();
				throw("Output stream for solution was not allocated");
			}

			in >> t0; in >> val;

			in.close();
			return;
		}
		void GetInitialData(double_t& t0, double_t& val) { t0 = this->RP.ini_data[0]; val = this->RP.ini_data[1]; }

	};

}
