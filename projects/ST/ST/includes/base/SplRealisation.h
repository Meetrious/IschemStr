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

	double_t alp; // ������ ���� �������-��������� �����������
	double_t bet; // ������ ���� �������-��������� �����������

	
	// ����� ������ ���������� ���� ������������ - ���-�� ��������������� ����������� ������� ��� ������ ��������

	size_t N_e = exp.size() - 1;

	for (int i = 0; i <= N_e; i++)
		SPL_info.push_back(std::vector<double_t>(5));

	// ��������� ��������� �������
	SPL_info[0][0] = exp[0][0];
	SPL_info[0][1] = exp[0][1]; 
	SPL_info[0][2] = 0; // ������ ����������� = 0 � ��������� �����

	
	for (size_t i = 1; i <= N_e; i++) // �� �� ����� ��� � (int i = 1; i < t.size(); i++) 
	{
		SPL_info[i][0] = exp[i][0]; 
		SPL_info[i][1] = exp[i][1]; 
		H_e = exp[i - 1][0] - exp[i][0]; 

		if (i == N_e) // ���� �� �� ��������� ���������� ������������
			SPL_info[i][2] = 0; // �� ������ ����������� �� ������ ���� ����������

		else // 1.1. ��������� ���� �� ���������� ���������
			if ((exp[i][1] >= exp[i - 1][1] && exp[i][1] >= exp[i + 1][1]) ||
				(exp[i][1] <= exp[i - 1][1] && exp[i][1] <= exp[i + 1][1])) SPL_info[i][2] = 0;
				
		// 1.2. ����� �������� ���� ����� ������, � ������� ������ ����� ���������� ���������, �������� ��������
			else
			{
				bet = atan((exp[i + 1][1] - exp[i][1]) / (exp[i + 1][0] - exp[i][0]));
				alp = atan((exp[i][1] - exp[i - 1][1]) / (exp[i][0] - exp[i - 1][0]));
				// ��������� ������ ����������� � t[i]
				SPL_info[i][2] = tan(alp - (alp - bet) / 2);
			}

		// 1.3. ��������� ������ 3� ����������� � t[i]
		SPL_info[i][4] = (1. / pow(H_e, 2)) * (SPL_info[i][2] + SPL_info[i - 1][2]) + (2. / pow(H_e, 3)) * (SPL_info[i][1] - SPL_info[i - 1][1]);


		// 1.4. ��������� �������� 2� ����������� � t[i] 
		SPL_info[i][3] = (1. / pow(H_e, 2)) * (SPL_info[i - 1][1] - SPL_info[i][1]) - (1. / H_e) * SPL_info[i][2] - H_e * SPL_info[i][4];
		
		// ��� ���������� ������� �� ����� ������ �������� �� ������ i-����������
		std::cout << "[" << exp[i - 1][0] << ", " << SPL_info[i][0] << "]" << "  ==>  "; // ��������������� ����������
		
		std::cout << SPL_info[i][1] << " + "; // ��������� ���� 3DEG ������� 
		std::cout << SPL_info[i][2] << "*x + "; // ����������� ��� ������ ������� ���������� 3DEG �������
		std::cout << SPL_info[i][3] << "*x^2 + "; // ����������� ��� ������ ������� ���������� 3DEG �������
		std::cout << SPL_info[i][4] << "*x^3 \n"; // ����������� ��� ������� ������� ���������� 3DEG �������
		// ��� �� � ������ ��� ���������� ������� �� ���������� [x_(i-1), x_i]
	}
	return true;
}
