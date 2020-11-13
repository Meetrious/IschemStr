#pragma once

#include <array>

namespace StraightTask
{

	template<typename argType>
	class IRightPart
	{
	protected:
		//amount of components in the RightPart of equation
		size_t comp_amount = 9;

		double_t Sum_B()
		{
			B[0] = 0.0;
			for (size_t i = 1; i <= comp_amount; i++) { B[0] += B[i]; }
			return B[0];
		}
	
	public:
		IRightPart() {}
		~IRightPart() = default;

		// an array that keeps values of components of the equation 
		std::array<double_t, 10> B;

		virtual double_t Expression(argType u) = 0;
		bool ret_is;

		virtual double_t Solution(double_t t) { return 0.0; }
	};

} // namespace StraightTask