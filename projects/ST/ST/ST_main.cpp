#include <Realisation.h>

using namespace StraightTask;

class INI_data {
public:
	INI_data() {}

	INI_data(uint32_t N, double_t t_0, double_t tau, uint16_t Ns) :
		m_N(N), m_t_0(t_0), m_tau(tau), m_Ns(Ns) {}

	~INI_data() = default;

	uint32_t m_N = 200; // размер сетки в итерации метода последовательного интегрирования (МПИ)
	double_t m_t_0 = 1; // начальный момент времени
	double_t m_tau = 0; // величина шага МПИ
	uint16_t m_Ns = 2; // число шагов МПИ
};

// решает диффур для вывода результата
template<typename Method>
void ODE_solver(ISolver<Method> & SYS, INI_data & IDT)
{
	//PhaseTrajOutput PTO("APPROX");
	// здесь определяем параметры численного метода
	SYS.ST.Set(IDT.m_t_0, IDT.m_N, IDT.m_tau, IDT.m_Ns);

	// Инициализируем накопитель запаздывающих аргументов
	for (auto const& cur : SYS.RetInitialiser) { cur(SYS.ST.N, SYS.ST.t_0, SYS.ST.gap_width); }

	// заполняем массивы из уже существующих решений на случай,
	// если надо зафиксировать систему вдоль некоторого решения компоненты системы
	//for (auto const& cur : SYS.DataCollector) { cur(); }

	// и определяем начальные условия из предзаписанного решения
	//SYS.ST.t_0 = SYS.ST.X_init.tim = 0.5;
	//for (auto const& cur : SYS.SolDataInitialiser) { cur(0, 0, SYS.ST.X_init); }


	// и определяем начальные условия
	SYS.ST.X_init.tj = SYS.ST.t_0;
	for (auto const& cur : SYS.IniDataInitialiser) { cur(SYS.ST.X_init); }
	for (auto const& cur : SYS.IniRetInitialiser) { cur(SYS.ST.X_init); }
	//for (auto const& cur : SYS.ExpressSubValues) { cur(SYS.ST.X_init); } // */

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

		if (current_gap == 0) SYS.ApplyPrepStep(Nj, Tj);
	
		// цикл на следующие 24 часа
		for (; Nj <= SYS.ST.N; Nj++)
		{
			// сдвиг на следующий шаг по времени
			Tj = SYS.ST.X_pred.tj += SYS.ST.H;

			//назначаем соответствующие запаздывания
			//for (auto const& cur : SYS.RetUploader) { cur(Nj, SYS.ST.X_pred); }
			SYS.RetUpload(Nj);

			// контролируем промежуток интерполяции
			//ST.X_pred.CheckShiftInterpGap(STpar);

			// предзаписанное решение
			//for (auto const& cur : SYS.SolDataInitialiser) { cur(current_day, Nj, SYS.ST.X_pred); }			
			
			//SYS.SolDataInitialiser[3](current_day, Nj, SYS.ST.X_pred); // цитокины

			SYS.ApplyMethod();

			//вывод, когда предиктор или констатация
			for (auto const& cur : SYS.SolutionOutputter) { cur(Nj, Tj, *SYS.ST.X_sol); }

			// вывод фазового портрета решения
			//PTO.OutputPhaseTraj(Nj, Tj, { (*SYS.ST.X_sol).x, (*SYS.ST.X_sol).y, (*SYS.ST.X_sol).z });
	
			SYS.ErrSaving(Tj);
			SYS.ErrOutput(Tj);
	
			//for (auto const& cur : SYS.BudgetOutputter) { cur(Nj, Tj); }

			// сдвигаем запаздывания (если они есть, лол)
			SYS.RetDataUpdate(Nj);

			// обновляем предшествующий временной ряд для перехода на следующий шаг
			SYS.NodeShift();


		} // конец цикла рассчётов на текущий день for(Nj: 1->N)

		current_gap++; std::cout << " ; ";

		// Проверка, «а надо ли решать дальше?» 
		if (current_gap == SYS.ST.full_amount_of_gaps)
		{
			std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			//TRIG = false; // условие на выход из цикла while
			// else не нужен, ведь если !ID, то внешнее условие while его поймает
		}
		else
		{
			SYS.ST.X_init = *SYS.ST.X_sol;
			// SYS.ST.X_init = SYS.ST.X_pred;
			Nj = 1;
		}
	}


	//SYS.TCOS.GetOrdOfMaxError();
	//SYS.TSIN.GetOrdOfMaxError();
	//SYS.T.GetOrdOfMaxError();

	std::cout << std::endl;
	/*std::cout << "( " << SYS.TCOS.Error[0] <<"\t" <<  SYS.TSIN.Error[0] << "\t" << SYS.T.Error[0] <<  " ) " << "\n"
		<< SYS.TCOS.Error[1] + SYS.TSIN.Error[1] + SYS.T.Error[1] << std::endl;// */
	std::cout << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][0] << "  " << SYS.EXP.ErrorArray[SYS.EXP.GetOrdOfMaxError()][1] << std::endl;
}

//решает диффур для рассчёта погрешности
template<typename MTD>
void ODE_solver(ISolver<MTD>);

void SolutionOutput()
{
	Test::ThreeDim::Ox_ret X;
	Test::ThreeDim::Oy_ret Y;
	Test::ThreeDim::Oz_ret Z;

	PhaseTrajOutput PTO("SOL");

	uint16_t GA = 10;
	uint32_t N = 2000;
	double_t H = pi / N;

	//std::ofstream sol(output_dir + "ST/SOL/PTraj/" + "SOL" + ".txt");
	for (uint32_t i = 0; i < GA * N; i++){
		PTO.OutputPhaseTraj(i, pi + i * H, { X.Solution(pi + i * H), Y.Solution(pi + i * H), Z.Solution(pi + i * H) });
	}
	//sol.close();
	return;
}



int main()
{
	
	Euler SYS; INI_data IDT(200, pi, pi / 2.0, 10);

	//SolutionOutput();
	ODE_solver(SYS, IDT);


	//system("pause");
	
	return EXIT_SUCCESS;
}
