#include "libraries/Methods.h"

using namespace StraightTask;

void ODE_PredCor_solver(int32_t full_amount_of_days, size_t grid_spliting_amount)
{

	Parameters STpar(full_amount_of_days, grid_spliting_amount); 

	// вводим переменные векторы, содержащие всю необходимую информацию о численном решении на каждом шаге
	// и определяем начальные условия
	Methods::PredCor TLayers(STpar);  
	
	Methods STask;
	// открываем потоки для вывода численного решения 
	STask.AllocateOutputStream();

	STask.CollectData("cy_a-spline.dat");
	TLayers.X_init.cy = STask.data[1][0];

	// выводим данные начальных условий во внешний файл
	for (size_t i = 1; i < Methods::N_eq + 5; i++)	STask.Output[i](TLayers.X_init, 0);


	bool TRIG = true;
	uint16_t days_remain = STpar.FAOD;
	Methods::H = STpar.H;

	/*double_t* CFC[12] = {
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
	};*/
	//Methods::setCoefs(CFC, 12, working_directory + "mstat/medX.txt");

//	STask.CollectData("cy_a-spline.dat");


	while (TRIG)
	{
		std::cout << days_remain << " "; // вывод на экран "кол-во дней, которые надо просчитать"
		TLayers.X_prev = TLayers.X_init; // даём начальное условие на предшествующий вектор

		// цикл на следующие 24 часа
		for (size_t Nj = 1; Nj <= STpar.N; Nj++)
		{
			// сдвиг на следующий шаг по времени
			TLayers.X_pred.tim += STpar.H;

			//назначаем соответствующие запаздывания
			TLayers.X_pred.AllocCurRets(Nj, STpar.N);

			// контролируем промежуток интерполяции
			TLayers.X_pred.CheckShiftInterpGap(STpar);
			TLayers.X_cor = TLayers.X_pred;

			// определяем текущее значение цитокинов из данных внешнего файла
			TLayers.X_cor.cy = STask.data[1][(Nj*16) + 24000 * (STpar.FAOD - days_remain)];
			
			launchPredCor(TLayers);

			// вывод во внешние *.txt
			for (size_t i = 1; i < Methods::N_eq + 5; i++) STask.Output[i](TLayers.X_cor, Nj);

			// сдвигаем запаздывания (если они есть, лол)
			rets::ShiftRets(TLayers.X_prev, Nj);

			// обновляем предшествующий временной ряд для перехода на следующий шаг
			TLayers.X_prev = TLayers.X_cor;

		} // конец цикла рассчётов на текущий день for(Nj: 1->N)

		days_remain--; std::cout << " ; ";

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


int main()
{
	// зарание конфигурируем сплайны
	// вне зависимости от того, нужны они в задаче или нет
	Splines::ConfigureSplines();
	
	// отрисовываем построенные выше сплайны, если надо
	//Splines::OutputSplines();

	// вызываем решение диффура (2 дня, по 1500 шагов в сетке)
	ODE_PredCor_solver(2,1500);

	system("pause");
	return EXIT_SUCCESS;
}