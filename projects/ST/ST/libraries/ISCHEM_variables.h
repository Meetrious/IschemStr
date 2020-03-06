#pragma once
#include <array>

#include "StraightTaskParameters.h"

namespace StraightTask
{
	struct variables;
	struct rets
	{
		static void ShiftRets(variables & prev_X, size_t Nj);

		static void SetInitialState(size_t N)
		{
			AsRet.clear();
			HelRet.clear();
			AdhRet.clear();
			DRet.clear();

			for (size_t i = 0; i < N; i++)
			{
				AsRet.emplace_back(0.0);
				HelRet.emplace_back(0.0);
				AdhRet.emplace_back(0.0);
				DRet.emplace_back(0.0);
			}
		}


		void AllocCurRets(size_t cur_ind, size_t N)
		{
			a_12 = AsRet[CountRetIndex(12.0, cur_ind, N)];
			hel_12 = HelRet[CountRetIndex(12.0, cur_ind, N)];
			adh_12 = AdhRet[CountRetIndex(12.0, cur_ind, N)];
			d_12 = DRet[CountRetIndex(12.0, cur_ind, N)];

			adh_4 = AdhRet[CountRetIndex(4.0, cur_ind, N)];
			adh_24 = AdhRet[CountRetIndex(24.0, cur_ind, N)];
		}

		double_t a_12;
		double_t hel_12;
		double_t adh_4, adh_12, adh_24;
		double_t d_12;

	private:

		[[nodiscard]] size_t CountRetIndex(double_t ret, size_t cur_ind, size_t N)
		{
			return (cur_ind + (size_t)((N * (240.0 - ret * 10.0) / 240.0))) % N;
		}

		static std::vector <double_t> AsRet;
		static std::vector <double_t> HelRet;
		static std::vector <double_t> AdhRet;
		static std::vector <double_t> DRet;
	};

	std::vector<double_t> rets::AsRet = {};
	std::vector<double_t> rets::HelRet = {};
	std::vector<double_t> rets::AdhRet = {};
	std::vector<double_t> rets::DRet = {};

	struct variables : rets
	{
		variables() {
			SIS[TIME] = SIS[NEC] = SIS[ACU] = SIS[HEL] = false;
			SIS[CY] = SIS[LM] = SIS[LN] = false;
			SIS[AP_S] = SIS[AP_E] = SIS[CH] = SIS[ADH] = SIS[MIA] = SIS[MII] = false;
			SIS[D_F] = SIS[D_INI] = SIS[DP_N] = SIS[DP_A] = SIS[EPS] = SIS[ePS] = SIS[PSY] = false;
		}

		double_t* GetMatchedVar(NumOfDerivative current)
		{
			switch (current)
			{
			case NEC: return &nec;
			case ACU: return &acu_c;
			case AP_S: return &ap_s;
			case AP_E: return &ap_e;
			case HEL: return &hel;


			case CY: return &cy;
			case CH: return &ch;
			case ADH: return &adh;

			case MIA:return &mia;
			case MII:return &mii;
			case LM: return &lm;
			case LN: return &ln;

			case D_F: return &d_F;

			case D_INI: return &d_ini;
			case DP_N: return &dp_N;
			case DP_A: return &dp_A;

			case ePS: return &eps;
			case EPS: return &Eps;
			case PSY: return &psy;
			default: return nullptr;
			}
		}
		double_t* GetMatchedVar(size_t current)
		{
			switch (current)
			{
			case 1: return &nec;
			case 2: return &acu_c;
			case 3: return &ap_s;
			case 4: return &ap_e;
			case 5: return &hel;


			case 6: return &cy;
			case 7: return &ch;
			case 8: return &adh;

			case 9:return &mia;
			case 10:return &mii;
			case 11: return &lm;
			case 12: return &ln;

			case 13: return &d_F;

			case 14: return &d_ini;
			case 15: return &dp_N;
			case 16: return &dp_A;

			case 17: return &eps;
			case 18: return &Eps;
			case 19: return &psy;
			default: return nullptr;
			}
		}

		void SetMatchedVar(NumOfDerivative current, double_t value)
		{
			double_t* var = GetMatchedVar(current);
			(*var) = value;
		}
		void SetMatchedVar(size_t current, double_t value)
		{
			double_t* var = GetMatchedVar(current);
			(*var) = value;
		}


		void CheckShiftInterpGap(Splines const & CurSTPar)
		{
			if (!Splines::IsSplinesConfigured) return;
			// ¬ычисл€ем промежуток текущей интерпол€ции дл€ каждой величины
			for (int i = 1; i < IND; i++)
			{
				if (!CurSTPar.EST[i]) // если i-та€ величина не интерполируетс€ => 
					continue; // ломаем блок if/else, массив for дает следующую величину на проверку
				else
				{
					if (tim - CurSTPar.IPM[i][gap_counter[i]][0] > 0)
					{ //  попадение сюда знаменует, что t переменна€ вышла за пределы промежутка интерполировани€ текущим кубическим сплайном
						// вывод: надо мен€ть 3-многочлен. 
						// ѕровер€ем, на последнем ли мы промежутке интерпол€ции дл€ данной i - той величины
						if (gap_counter[i] + 1 == Splines::max_gap_amount[i]) // надо k_e + 1
						{ // i-та€ величина оказалась на последнем промежутке интерпол€ции и t вышла в промежуток экстрапол€ции 
							SIS[i] = false; // выключаем передачу значени€ интерпол€нта дл€ i-той величины
							continue; // выходим из всех if блоков и переходим на проверку сделующей (i+1)-й величины
						}

						gap_counter[i]++; // обновл€ем текущий промежуток интерпол€ции (инкремент счетчика k_e[i]) дл€ текущей i-той величины

						// переопредел€ем коэффициенты сплайна
						for (int j = 0; j < 5; j++)
							SplineData[i][j] = CurSTPar.IPM[i][gap_counter[i]][j];
					}
				}
			}
			return;
		}

		double_t GetSplineValue(size_t i)
		{
			return Splines::GetValue(tim, SplineData[i][0], SplineData[i][1], SplineData[i][2], SplineData[i][3], SplineData[i][4]);
		}



		// процедура инициализации текущих коэффициентов на полиномы сплайна
		void SetSplineData(Parameters const &CurSTPar)
		{
			for (int i = 1; i < IND; i++)
			{
				int32_t k;
				if (!CurSTPar.EST[i])
				{
					gap_counter[i] = 0;
					SIS[i] = false;
					continue; // дл€ i-той величины у нас нет данных.   следующей величине
				}
				else
				{
					// настройка gap_counter с учЄтом указанного начального момента времени
					k = 1;
					while (tim > CurSTPar.IPM[i][k][0]) k++; // вычисл€ем gap_counter дл€ i-той величины
					gap_counter[i] = k;

					if (gap_counter[i] > CurSTPar.max_gap_amount[i])
					{
						gap_counter[i] = 0;
						SIS[i] = false;
						continue; // t вышло за область интерпол€ции.   следующей величине
					}
					//				SIS[i] = true;
									// заполн€ем данные уравнени€ сплайна дл€ вычисленного выше gap_counter 
					for (int j = 0; j < 5; j++)
						SplineData[i][j] = CurSTPar.IPM[i][gap_counter[i]][j];
				}
			}
		}

		// процедура задани€ начальных условий из внешних файлов/левых краев сплайнов
		void MakeVarsInitial(Parameters const &CurSTPar)
		{

			std::ifstream ini(working_directory + init_dir);
			if (!ini)
			{
				std::cout << "\n Could not create an output.txt in the initial_step. \n";
				getchar();
				return;
			}

			ini >> tim;
			ini >> nec;  ini >> acu_c;  ini >> ap_s;  ini >> ap_e;  ini >> hel;

			ini >> cy;  ini >> ch;  ini >> adh;
			ini >> mia; ini >> mii; ini >> lm; ini >> ln;
			ini.close();

			SetSplineData(CurSTPar);

			if (SplinesAreInitials)
			{
				for (size_t i = 1; i < IND; i++)
				{
					if (gap_counter[i] == 0) continue;
					else SetMatchedVar(static_cast<NumOfDerivative>(i), GetSplineValue(i));
				}
			}

			return;
		}


		~variables() = default;

		double_t tim;

		double_t nec, acu_c, ap_s, ap_e, hel;

		double_t cy, ch, adh;
		double_t eps, Eps;

		double_t mia, mii, lm, ln;

		double_t d_F, d_ini;
		double_t dp_A, dp_N;

		double_t psy;

		uint16_t gap_counter[IND];
		std::array<std::array<double_t, 5>, IND> SplineData;

		bool SIS[IND];

		// триггер:
		// если true, то данные сплайна превалируют в начальных услови€х против тех, что указаны во внешнем *.txt
		static constexpr bool SplinesAreInitials = true;
	};

	void rets::ShiftRets(variables & prev_X, size_t cur_ind)
	{
		AsRet[cur_ind - 1] = prev_X.ap_s;
		HelRet[cur_ind - 1] = prev_X.hel;
		DRet[cur_ind - 1] = prev_X.d_ini;
		AdhRet[cur_ind - 1] = prev_X.adh; // Ui[]; // AdhRet
	}
}