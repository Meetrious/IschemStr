#pragma once
#include <string>
#include <fstream>
#include <iostream>

#include <vector>

#include <base/DirPaths.h>

namespace ReverseTask
{
	template <typename T>
	using vector = std::vector<T>;

	template <typename T>
	using matrix = vector<vector<T>>;


	class Discrete {
	public:

		Discrete(const char* name)
			: m_name{ name }, Nk{ 0 }, m_res{ 0.0 } {}
		~Discrete() = default;
		const char* m_name;

		void GatherData(std::string DaWay) {

			std::ifstream in(DaWay);
			if (!in) {
				std::cout << "\n Get Fucked \n";
				getchar();
				exit(1337);
			}

			for (size_t i = 0; !in.eof(); i++)
			{
				in >> Tjs.emplace_back();
				in >> Proper.emplace_back();
				Calc.emplace_back(-1.0);
			}

			in.close();

			Tjs.shrink_to_fit();
			Proper.shrink_to_fit();
			Calc.shrink_to_fit();
		}

		void CollectCalc(float_t H, float_t Tj, double_t val) {

			if (Tjs[Nk] - Tj > H / 2.0) return;
			else { Calc[Nk] = val;	Nk++; }
		}

		double_t CountResult() {

			for (size_t i = 0; i < Proper.size(); i++) {
				m_res += std::abs(Proper[i] - Calc[i]);
			}

			m_res /= Proper.size();

			return m_res;
		}

		void ResetState() { Nk = 0; m_res = 0.0; }


	private:

		uint32_t Nk;

		vector<float_t> Tjs;
		vector<double_t> Proper;
		vector<double_t> Calc;

		double_t m_res;
	};

	class Continuous {
	public:
		Continuous(const char* name) : m_name{ name } {}
		~Continuous() = default;
		const char* m_name;

		void GatherData(std::string DaWay) {

			std::ifstream in(DaWay);
			if (!in) {
				std::cout << "\n Get Fucked \n";
				getchar();
				exit(1337);
			}

			float_t dummy;
			for (size_t i = 0; !in.eof(); i++)
			{
				in >> dummy;
				in >> dummy;
				in >> Proper.emplace_back();
				Calc.emplace_back();
			}

			in.close();

			//Tjs.shrink_to_fit();
			Proper.shrink_to_fit();
			Calc.shrink_to_fit();
		}

		void CollectCalc(uint32_t Nj, uint16_t gap, uint32_t N, double_t val) {
			Calc[Nj + gap * N] = val;
		}

		double_t CountResult() {

			for (size_t i = 0; i < Proper.size(); i++) {
				m_res += std::abs(Proper[i] - Calc[i]);
			}

			m_res /= Proper.size();

			return m_res;
		}

		void ResetState() { m_res = 0.0; }


	private:
		//vector<float_t> Tjs;
		vector<double_t> Proper;
		vector<double_t> Calc;

		double_t m_res;
	};

	class IAggregateControls;
}