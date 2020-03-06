#pragma once
#include <iomanip> // нужен для std::setpresicion модификатора в std::fstream потоках
#include "DirPaths.h"
#include "ISCHEM_enum.h"
#include "SplRealisation.h"

namespace StraightTask
{
	struct Splines
	{
		Splines() {}

		// процедура построения массивов информации для сплайнов
		static void ConfigureSplines()
		{
			EST[TIME] =  false;

			EST[NEC] = EST[ACU] = EST[HEL] = EST[LM] = EST[LN] = true;

			EST[CY] = EST[AP_S] = EST[AP_E] = EST[MIA] = EST[MII] = EST[CH] = EST[ADH] = false;

			EST[D_F] = EST[D_INI] = EST[DP_N] = EST[DP_A] = EST[EPS] = EST[ePS] = EST[PSY] = false;

			// буферная строка для наименования текстового файла, в котором лежит информация для текущих экспериментальных данных
			std::string SPLway;

			// "пропушим" каретку последнего элемента в нулевое состояние,
			// чтобы для SP индекс k тоже начинался с единицы
			IPM.push_back(std::vector<std::vector<double_t>>());

			for (int k = 1; k < IND; k++)
			{
				IPM.push_back(std::vector<std::vector<double_t>>());
				if (!EST[k])
				{
					max_gap_amount[k] = 0;
					continue;
				}
				else
				{
					switch (k)
					{
					case NEC:
					{ // Некротические
						SPLway = working_directory + "data/exp/NECR(exact).txt";
					} break;
					case ACU:
					{ // Повреждённые
						SPLway = working_directory + "data/exp/APOP(1).txt";
					} break;
					case HEL:
					{ // Здоровые
						//SPLway = "exp/HEALT(exact).txt";
						SPLway = working_directory + "data/exp/HEALT.txt";
					} break;
					case LM:
					{ // Макрофаги
		//				SPLway += "exp/MF(%).txt";
						SPLway = working_directory + "data/exp/MF(exact).txt"; // 
					} break;
					case LN:
					{// Нейтрофилы
		//				SPLway += "exp/PMN(%).txt";
						SPLway = working_directory + "data/exp/PMN(Par).txt";
					} break;
					case CY:
					{// Цитокины
						SPLway = working_directory + "data/exp/Cyto(exact).txt";
					} break;
					default:
					{
						std::cout << "Some mistake occured in allocating input files with ...uuh dunno what"; system("pause");
					} break;
					}

					// Находим коэф-ты кубического сплайна
					std::cout << "Data for SPL[" << k << "] :\n";
					if (!QMSmaker(SPLway, IPM.at(k)))
					{
						std::cout << "\n An error occured in \"QMSmaker\" prosedure. Care to attend.";
						system("pause");
						return;
					}
					max_gap_amount[k] = IPM[k].size();

					// Вот все и готово для построения сплайна для k-й величины для всех ее промежутков
					IPM.at(k).shrink_to_fit();
				}
				std::cout << "\n\n";
			}
			IPM.shrink_to_fit();
			IsSplinesConfigured = true;
			return;
		}

		inline static double_t GetValue(double_t t, double_t ti, double_t ai, double_t bi, double_t ci, double_t di)
		{
			return ai + bi * (t - ti) + ci * std::pow(t - ti, 2) + di * std::pow(t - ti, 3);
		}

		static void Draw(size_t i);

		static void OutputSplines();

		~Splines() = default;


		static bool EST[IND]; // EACH_SPLINE_TRIGGER
		static std::vector < std::vector < std::vector<double_t >> > IPM; // interpolating polynomial matrix
		static size_t max_gap_amount[IND];
		static bool IsSplinesConfigured;

	};

	bool Splines::IsSplinesConfigured = false;


	size_t Splines::max_gap_amount[IND] = {};
	std::vector < std::vector < std::vector<double_t >> > Splines::IPM = {};
	bool Splines::EST[IND] = {};
	// корректная инициализация последних трёх объектов происходит в самом начале рантайма при вызове ConfigureSplines;
	// без предварительного вызова ConfigureSplines, сплайны не будут работать корректно.

	struct Parameters : Splines
	{
		Parameters()
		{
			std::cout << "Choose the amount of grid spliting: N := ";
			std::cin >> N;
			while (N <= 6 || N % 2 != 0) // защита на случай, когда шаг сетки был введён некорректно
			{
				std::cout << "\nTry again: N must be integer even and more than 6; N := ";
				std::cin >> N;
			}

			H = 24.0 / double_t(N);
			std::cout << " length of a step H = " << H;
			std::cout << "\n\n";

			std::cout << "How many days? \n Amount = ";
			std::cin >> FAOD;
			while (FAOD <= 0)
			{
				std::cout << "Try again: amount has to be more than 0\n";
				std::cin >> FAOD;
			}
			std::cout << "\n\n";

			coef_out_names();

		}

		Parameters(uint32_t full_amount_of_days, size_t grid_split_number)
		{
			FAOD = full_amount_of_days;
			N = grid_split_number;
			H = 24.0 / double_t(grid_split_number);

			//coef_out_names();
			SetCoefs(); // инициализируем коэффициенты из 1го *.txt - шника
		}

		void coef_out_names()
		{
			int32_t i = 0;

			std::string ininum = "1"; // выбор варианта значений коэффициентов
			std::string in1 = "1.Aim.pN_Rep_kN_kA_pR.txt";
			std::string in2 = "2.Micro.p1_cA_cN_cpro__TM1_TM2_cMi1_cMi2__cdMi_cMi_KMi.txt";
			std::string in3 = "3.Leuco.cLm_pdLm_TLm__cLm1_KLm__cLn_pdLn_TLn__cLn1_KLn.txt";
			std::string in4 = "4.Prot.CMa_CLm_ecy__pMach_pLmch_ech.txt";
			std::string in5 = "5.Adhen.pMadhcy1_pMadhcy2_eMadh.txt";
			std::string in6 = "6.Tox.pncy__pnLn_CDLn__pmLm_CDLm__Pnn_D0_pD.txt";
			std::string in7 = "7.Phag.en_enLm_enLn_enmi.txt";
			std::string in8 = "8.Psy.pxcy0_cymax_T'_t0.txt";

			std::string place = "coefs/I/";

			std::ifstream inCoefs;
			inCoefs.open(place + ininum + in1, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in1 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			place = "coefs/II/";
			inCoefs.open(place + ininum + in2, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in2 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			inCoefs.open(place + ininum + in3, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in3 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			inCoefs.open(place + ininum + in4, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in4 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			inCoefs.open(place + ininum + in5, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in5 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			place = "coefs/sub/";
			inCoefs.open(place + ininum + in6, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in6 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			inCoefs.open(place + ininum + in7, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in7 << "> was not created\n";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();

			inCoefs.open(place + ininum + in8, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs for <" << in8 << "> was not created\n ";
				std::cin >> i;
			}
			while (!inCoefs.eof())
			{
				inCoefs >> Coefs.at(i);
				i++;
			}
			inCoefs.close();
			std::cout << " \n all inCoefs were succesfully assigned\n ";
		}

		void SetCoefs()
		{
			std::ifstream inCoefs;
			inCoefs.open(working_directory + coefs_default_dir, 2);
			if (!inCoefs)
			{
				std::cout << "\n inCoefs was not created for some reason \n";
			}
			size_t i = 0;
			while (!inCoefs.eof()) 	inCoefs >> Coefs.at(i++);
			if (i != length_coef_array)
				std::cout << "\n something is wrong with coefs allocation \n";
		}

		~Parameters() = default;


		std::array<double_t, length_coef_array>Coefs;

		const static std::array<std::string, IND> sol_names;

		uint16_t FAOD; // full amount of days  | 0 < FAOD < 32767
		size_t N; // Grid Spliting amount
		double_t H;

	};

	const std::array<std::string, IND> Parameters::sol_names = {
		" ",
		"Necr.txt",
		"Ac_ch.txt",
		"Apop_s.txt",
		"Apop_e.txt",
		"Healt.txt",

		"CY.txt",
		"CH.txt",
		"Adhen.txt",
		"Mi_act.txt",
		"Mi_inact.txt",

		"Macr.txt",
		"Neut.txt",

		"D_Full.txt",

		"D_INI.txt",
		"DN_c.txt",
		"DA_c.txt",
		"Eps_c.txt",
		"PSY_c.txt"
	};

	void Splines::Draw(size_t i)
	{
		if (!EST[i])
		{
			std::cout << " \n\t there are no splines for you at " << i
				<< "th vatiable. Try another one";
			return;
		}

		std::ofstream out(working_directory + "Splines/SPL_" + Parameters::sol_names[i]);

		size_t k = 1;
		double_t t = IPM[i][0][0];
		while (k < max_gap_amount[i])
		{
			while (t < IPM[i][k][0])
			{
				out << t << "\t\t\t"
					<< GetValue(t, IPM[i][k][0], IPM[i][k][1], IPM[i][k][2], IPM[i][k][3], IPM[i][k][4])
					<< "\n";
				t += 0.016; // диаметр разбиения по умолчанию
			}
			k++;
		}
		out.close();
	}


	void Splines::OutputSplines()
	{
		for (size_t i = 0; i < IND; i++)
		{
			if (!Splines::EST[i]) continue;
			else Splines::Draw(i);
		}
	}
}