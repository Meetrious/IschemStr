/* This header constains non-member function that formulates the shape of 
3DEG spline equations depending on initial interpolation data. // */

#pragma once
#include <iostream>
#include <vector>

// my lazy method that doesn't even control the oscillations that may appear 
bool QMSmaker(
	std::vector<std::vector<double_t>> const & exp,
	std::vector<std::vector<double_t>> & SPL_info)
{
	double_t H_e;

	double_t alp; // первый угол кусочно-линейного приближения
	double_t bet; // второй угол кусочно-линейного приближения

	
	// узнаём индекс последнего узла интерполяции - кол-во интерполируемых промежутков времени для данной величины

	size_t N_e = exp.size() - 1;

	for (int i = 0; i <= N_e; i++)
		SPL_info.push_back(std::vector<double_t>(5));

	// заполняем начальное условие
	SPL_info[0][0] = exp[0][0];
	SPL_info[0][1] = exp[0][1]; 
	SPL_info[0][2] = 0; // первая производная = 0 в начальной точке

	
	for (size_t i = 1; i <= N_e; i++) // то же самое что и (int i = 1; i < t.size(); i++) 
	{
		SPL_info[i][0] = exp[i][0]; 
		SPL_info[i][1] = exp[i][1]; 
		H_e = exp[i - 1][0] - exp[i][0]; 

		if (i == N_e) // если мы на последнем промежутке интерполяции
			SPL_info[i][2] = 0; // то первая производная на правом краю зануляется

		else // 1.1. Проверяем узел на дискретный экстремум
			if ((exp[i][1] >= exp[i - 1][1] && exp[i][1] >= exp[i + 1][1]) ||
				(exp[i][1] <= exp[i - 1][1] && exp[i][1] <= exp[i + 1][1])) SPL_info[i][2] = 0;
				
		// 1.2. иначе усредним угол между узлами, в которых график имеет монотонное поведение, линейной склейкой
			else
			{
				bet = atan((exp[i + 1][1] - exp[i][1]) / (exp[i + 1][0] - exp[i][0]));
				alp = atan((exp[i][1] - exp[i - 1][1]) / (exp[i][0] - exp[i - 1][0]));
				// вычисляем первую производную в t[i]
				SPL_info[i][2] = tan(alp - (alp - bet) / 2);
			}

		// 1.3. вычисляем шестую 3й производной в t[i]
		SPL_info[i][4] = (1. / pow(H_e, 2)) * (SPL_info[i][2] + SPL_info[i - 1][2]) + (2. / pow(H_e, 3)) * (SPL_info[i][1] - SPL_info[i - 1][1]);


		// 1.4. вычисляем половину 2й производной в t[i] 
		SPL_info[i][3] = (1. / pow(H_e, 2)) * (SPL_info[i - 1][1] - SPL_info[i][1]) - (1. / H_e) * SPL_info[i][2] - H_e * SPL_info[i][4];
		
		// для отчетности выводим на экран данные сплайнов на данном i-промежутке
		std::cout << "[" << exp[i - 1][0] << ", " << SPL_info[i][0] << "]" << "  ==>  "; // интерполируемый промежуток
		
		std::cout << SPL_info[i][1] << " + "; // свободный член 3DEG сплайна 
		std::cout << SPL_info[i][2] << "*x + "; // коэффициент при первой степени переменной 3DEG сплайна
		std::cout << SPL_info[i][3] << "*x^2 + "; // коэффициент при второй степени переменной 3DEG сплайна
		std::cout << SPL_info[i][4] << "*x^3 \n"; // коэффициент при третьей степени переменной 3DEG сплайна
		// Вот всё и готово для построения сплайна на промежутке [x_(i-1), x_i]
	}
	return true;
}
