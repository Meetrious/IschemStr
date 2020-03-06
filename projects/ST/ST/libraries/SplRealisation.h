#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <vector>


bool QMSmaker(std::string const & SPLway, std::vector<std::vector<double_t>> & SPL_info)
{
	// Примем эксперементальные данные и строим по ним сплайны
	std::vector<double_t> current_t_grid; // буферный массив для содержания узлов на сетке по времени для интерполирования текущих экспериментальных данных
	std::vector<double_t> current_val_grid; // буферный массив для содержания значений текущих экспериментальных данных в соответствующих узлах на сетке

	// буфер для величин промежуточных шагов, между эксперементальными данными текущей k - величины | [current_t_grid[1] - current_t_grid[2]] |
	double_t H_e;
	double_t tmp;

	double_t alp; // первый угол кусочно-линейного приближения
	double_t bet; // второй угол кусочно-линейного приближения

	std::ifstream input;

	input.open(SPLway, std::ifstream::in);
	if (!input)// если вдруг путь к файлам, в которых лежат экспериментальные данные, оказался неверен
	{
		std::cout << "\n Spline-data was not acquired.(( Probably there is some typo in the direction of a *.txt file \n";
		getchar();
		for (int i = 0; i < IND; i++) // инициализируем интерполянты нолями
			SPL_info.push_back(std::vector<double_t>(5, 0));
		system("pause");
		return false; // тогда все сплайны схлопываются в 0 для всех величин
	}

	while (!input.eof())
	{
		input >> tmp;
		current_t_grid.emplace_back(tmp); // заполняем время
		input >> tmp;
		current_val_grid.emplace_back(tmp); // заполняем сырые данные
	}

	input.close();

	
	// узнаём индекс последнего узла интерполяции - кол-во интерполируемых промежутков времени для данной величины
	int N_e = (int)(current_t_grid.size() - 1);

	for (int i = 0; i <= N_e; i++)
		SPL_info.push_back(std::vector<double_t>(5));

	// заполняем начальное условие
	SPL_info[0][0] = current_t_grid[0];
	SPL_info[0][1] = current_val_grid[0];
	SPL_info[0][2] = 0; // первая производная = 0 в начальной точке

	
	for (int i = 1; i <= N_e; i++) // то же самое что и (int i = 1; i < t.size(); i++) 
	{
		SPL_info[i][0] = current_t_grid[i];
		SPL_info[i][1] = current_val_grid[i];
		H_e = current_t_grid[i - 1] - current_t_grid[i];

		if (i == N_e) // если мы на последнем промежутке интерполяции
			SPL_info[i][2] = 0; // то первая производная на правом краю зануляется

		else // 1.1. Проверяем узел на дискретный экстремум
			if ((current_val_grid[i] >= current_val_grid[i - 1] && current_val_grid[i] >= current_val_grid[i + 1]) || (current_val_grid[i] <= current_val_grid[i - 1] && current_val_grid[i] <= current_val_grid[i + 1]))
				SPL_info[i][2] = 0;
		// 1.2. иначе усредним угол между узлами, в которых график имеет монотонное поведение, линейной склейкой
			else
			{
				bet = atan((current_val_grid[i + 1] - current_val_grid[i]) / (current_t_grid[i + 1] - current_t_grid[i]));
				alp = atan((current_val_grid[i] - current_val_grid[i - 1]) / (current_t_grid[i] - current_t_grid[i - 1]));
				// вычисляем первую производную в t[i]
				SPL_info[i][2] = tan(alp - (alp - bet) / 2);
			}

		// 1.3. вычисляем шестую 3й производной в t[i]
		SPL_info[i][4] = (1. / pow(H_e, 2)) * (SPL_info[i][2] + SPL_info[i - 1][2]) + (2. / pow(H_e, 3)) * (SPL_info[i][1] - SPL_info[i - 1][1]);


		// 1.4. вычисляем половину 2й производной в t[i] 
		SPL_info[i][3] = (1. / pow(H_e, 2)) * (SPL_info[i - 1][1] - SPL_info[i][1]) - (1. / H_e) * SPL_info[i][2] - H_e * SPL_info[i][4];
		
		// для отчетности выводим на экран данные сплайнов на данном i-промежутке
		std::cout << "[" << current_t_grid[i - 1] << ", " << SPL_info[i][0] << "]" << "  ==>  "; // интерполируемый промежуток
		
		std::cout << SPL_info[i][1] << " + "; // свободный член 3DEG сплайна 
		std::cout << SPL_info[i][2] << "*x + "; // коэффициент при первой степени переменной 3DEG сплайна
		std::cout << SPL_info[i][3] << "*x^2 + "; // коэффициент при второй степени переменной 3DEG сплайна
		std::cout << SPL_info[i][4] << "*x^3 \n"; // коэффициент при третьей степени переменной 3DEG сплайна
		// Вот всё и готово для построения сплайна на промежутке [x_(i-1), x_i]
	}
	return true;
}
