/*This header constains the right parts of equations of ODE_system of the 1st-model //*/


#pragma once
#include <tests/variables.h>
#include <base/Eq_base.h>


namespace StraightTask
{
	static double_t pi = 3.14159265358979323;

	template<typename argType> class IEqMember: public IRightPart<argType>{};

	namespace Test
	{
		

		namespace OneDim 
		{
			static double_t k = 2;

			size_t factorial(size_t n) {
				if (n > 1) return n * factorial(n - 1);	else return 1;
			}

			class Exp_ret :public IEqMember<variables const&>{
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final{
					return u.ret.x_1;
				}
			public:
				Exp_ret() { ret_is = true; }
				const double_t ini_data[2] = { 1.0, 1.0 }; // initial data for t_0 = 0.0;

				const char* name = "Exp_ret";

				double_t Solution(double_t t) noexcept final {
					if (t <= 0) return 1.0;
					double_t value = 1;
					for (size_t n = 1; n < floor(t); n++) {
						value += pow((t - ((double)n - 1.0)), (double)n) / (double)factorial(n);
					}
					return value;
				}
				
				
			};

			class Exp :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return u.x;
				}
			public:
				Exp() { ret_is = false; }
				const double_t ini_data[2] = {0.0, 1.0}; // initial data for t_0 = 0.0;
				const char* name = "Exp";
				double_t Solution(double_t t) noexcept final { return exp(t); }
			};	

			class Sin :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return -k * u.ret.x_pi2 ;
				}
			public:
				Sin() { ret_is = true; }
				const double_t ini_data[2] = { pi, 0 }; // initial data for t_0 = pi;
				const char* name = "RetSin";
				double_t Solution(double_t t) noexcept final { return k * sin(t); }
			};


		}
		namespace ThreeDim
		{

			class Ox_ret :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return (2.0 / pi) * (u.x + u.ret.y_pi2) - u.ret.x_pi2 - (0.5 * pi) * (u.y / u.z);
				}
			public:
				Ox_ret() { ret_is = true; }
				double_t Solution(double_t t)noexcept final { return t * cos(t); }
				const double_t ini_data[2] = { pi,-pi };
				const char* name = "RetSp1"; 
			};
			class Oy_ret :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return (2.0 / pi) * (u.y - u.ret.x_pi2) - u.ret.y_pi2 + (0.5 * pi) * (u.x / u.z);
				}
			public:
				Oy_ret() { ret_is = true; }
				double_t Solution(double_t t)noexcept final { return t * sin(t); }
				const double_t ini_data[2] = { pi, 0 };
				const char* name = "RetSp2"; 
			};
			class Oz_ret :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u)noexcept final {
					return pow(u.ret.x_pi2 * u.ret.x_pi2 + u.ret.y_pi2 * u.ret.y_pi2, 0.5) / u.ret.z_pi2 ;
				}
			public:
				Oz_ret() { ret_is = true; }
				double_t Solution(double_t t)noexcept final { return t; }
				const double_t ini_data[2] = { pi, pi };
				const char* name = "RetSp3"; 
				
			};

			class Ox :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return u.x/u.z - u.y;
				}
			public:
				Ox() { ret_is = false; }
				double_t Solution(double_t t)noexcept final{ return t * cos(t);	}
				const double_t ini_data[2] = { pi,-pi };
				const char* name = "Sp1"; bool deflecting = false;
			};
			class Oy :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return u.y/u.z + u.x;
				}
			public:
				Oy() { ret_is = false; }
				double_t Solution(double_t t)noexcept final { return t * sin(t); }
				const double_t ini_data[2] = { pi,0 };
				const char* name = "Sp2"; bool deflecting = false;
				
			};
			class Oz :public IEqMember<variables const&> {
				[[nodiscard]] inline double_t Expression(variables const& u) final {
					return (u.x * u.x + u.y * u.y) / (u.z * u.z);
				}
			public:
				Oz() { ret_is = false; }
				double_t Solution(double_t t)noexcept final { return t; }
				const double_t ini_data[2] = { pi, pi };
				const char* name = "Sp3"; bool deflecting = false;
			};
		}
	}

} // namespace StraightTask