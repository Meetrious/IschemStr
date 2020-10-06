#pragma once
#include <vector>



namespace StraightTask
{
	template<typename T>
	using vector = std::vector<T>;

	template<typename T>
	using matrix = vector<vector<T>>;

	template<typename argType>
	class IRightPart
	{
	protected:
		double_t Sum_B() {
			B[0] = 0.0;
			for (uint16_t i = 1; i < B.size(); i++) { B[0] += B[i]; }
			return B[0];
		}
	
	public:
		
		vector<double_t> B;

		virtual double_t Expression(argType u) = 0;
		bool ret_is;

		virtual double_t Solution(double_t t) { return 0.0; }
	};

} // namespace StraightTask