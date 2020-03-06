#pragma once
#include <fstream>
#include <string>
#include <iostream>
#include <vector>


bool QMSmaker(std::string const & SPLway, std::vector<std::vector<double_t>> & SPL_info)
{
	// ������ ����������������� ������ � ������ �� ��� �������
	std::vector<double_t> current_t_grid; // �������� ������ ��� ���������� ����� �� ����� �� ������� ��� ���������������� ������� ����������������� ������
	std::vector<double_t> current_val_grid; // �������� ������ ��� ���������� �������� ������� ����������������� ������ � ��������������� ����� �� �����

	// ����� ��� ������� ������������� �����, ����� ������������������ ������� ������� k - �������� | [current_t_grid[1] - current_t_grid[2]] |
	double_t H_e;
	double_t tmp;

	double_t alp; // ������ ���� �������-��������� �����������
	double_t bet; // ������ ���� �������-��������� �����������

	std::ifstream input;

	input.open(SPLway, std::ifstream::in);
	if (!input)// ���� ����� ���� � ������, � ������� ����� ����������������� ������, �������� �������
	{
		std::cout << "\n Spline-data was not acquired.(( Probably there is some typo in the direction of a *.txt file \n";
		getchar();
		for (int i = 0; i < IND; i++) // �������������� ������������ ������
			SPL_info.push_back(std::vector<double_t>(5, 0));
		system("pause");
		return false; // ����� ��� ������� ������������ � 0 ��� ���� �������
	}

	while (!input.eof())
	{
		input >> tmp;
		current_t_grid.emplace_back(tmp); // ��������� �����
		input >> tmp;
		current_val_grid.emplace_back(tmp); // ��������� ����� ������
	}

	input.close();

	
	// ����� ������ ���������� ���� ������������ - ���-�� ��������������� ����������� ������� ��� ������ ��������
	int N_e = (int)(current_t_grid.size() - 1);

	for (int i = 0; i <= N_e; i++)
		SPL_info.push_back(std::vector<double_t>(5));

	// ��������� ��������� �������
	SPL_info[0][0] = current_t_grid[0];
	SPL_info[0][1] = current_val_grid[0];
	SPL_info[0][2] = 0; // ������ ����������� = 0 � ��������� �����

	
	for (int i = 1; i <= N_e; i++) // �� �� ����� ��� � (int i = 1; i < t.size(); i++) 
	{
		SPL_info[i][0] = current_t_grid[i];
		SPL_info[i][1] = current_val_grid[i];
		H_e = current_t_grid[i - 1] - current_t_grid[i];

		if (i == N_e) // ���� �� �� ��������� ���������� ������������
			SPL_info[i][2] = 0; // �� ������ ����������� �� ������ ���� ����������

		else // 1.1. ��������� ���� �� ���������� ���������
			if ((current_val_grid[i] >= current_val_grid[i - 1] && current_val_grid[i] >= current_val_grid[i + 1]) || (current_val_grid[i] <= current_val_grid[i - 1] && current_val_grid[i] <= current_val_grid[i + 1]))
				SPL_info[i][2] = 0;
		// 1.2. ����� �������� ���� ����� ������, � ������� ������ ����� ���������� ���������, �������� ��������
			else
			{
				bet = atan((current_val_grid[i + 1] - current_val_grid[i]) / (current_t_grid[i + 1] - current_t_grid[i]));
				alp = atan((current_val_grid[i] - current_val_grid[i - 1]) / (current_t_grid[i] - current_t_grid[i - 1]));
				// ��������� ������ ����������� � t[i]
				SPL_info[i][2] = tan(alp - (alp - bet) / 2);
			}

		// 1.3. ��������� ������ 3� ����������� � t[i]
		SPL_info[i][4] = (1. / pow(H_e, 2)) * (SPL_info[i][2] + SPL_info[i - 1][2]) + (2. / pow(H_e, 3)) * (SPL_info[i][1] - SPL_info[i - 1][1]);


		// 1.4. ��������� �������� 2� ����������� � t[i] 
		SPL_info[i][3] = (1. / pow(H_e, 2)) * (SPL_info[i - 1][1] - SPL_info[i][1]) - (1. / H_e) * SPL_info[i][2] - H_e * SPL_info[i][4];
		
		// ��� ���������� ������� �� ����� ������ �������� �� ������ i-����������
		std::cout << "[" << current_t_grid[i - 1] << ", " << SPL_info[i][0] << "]" << "  ==>  "; // ��������������� ����������
		
		std::cout << SPL_info[i][1] << " + "; // ��������� ���� 3DEG ������� 
		std::cout << SPL_info[i][2] << "*x + "; // ����������� ��� ������ ������� ���������� 3DEG �������
		std::cout << SPL_info[i][3] << "*x^2 + "; // ����������� ��� ������ ������� ���������� 3DEG �������
		std::cout << SPL_info[i][4] << "*x^3 \n"; // ����������� ��� ������� ������� ���������� 3DEG �������
		// ��� �� � ������ ��� ���������� ������� �� ���������� [x_(i-1), x_i]
	}
	return true;
}
