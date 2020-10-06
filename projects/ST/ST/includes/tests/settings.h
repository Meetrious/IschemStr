//settings for Staight test-task
#pragma once
#include <tests/equations.h>
#include <base/Settings_base.h>


namespace StraightTask
{
	template<typename RightPart>
	class ISysMember : public IMember<RightPart>, public IErr{
		const char* sol_name() final { return this->RP.name; }
	public:
		void GetInitialData(double_t& t0, double_t& val) { t0 = this->RP.ini_data[0]; val = this->RP.ini_data[1]; }
		void StateRetArray(uint32_t ST_N, double_t ST_t0, double_t ST_gap) final {
			this->N = ST_N; this->t_0 = ST_t0; this->gap_width = ST_gap;
			this->data.clear();
			//for (uint32_t i = 0; i < N; i++) { data.emplace_back(1.0); }
			for (uint32_t i = 0; i < this->N; i++) {
				this->data.emplace_back
				(this->RP.Solution(this->t_0 - this->gap_width + i * (this->gap_width / this->N))); 
			}
		}
		bool is_deflecting() { return RP.is_ret; }
	};

//---------------------------------------------------------------------------------
//=================================================================================


	namespace Test
	{
		namespace OneDim {
			class M_ODE : public ISysMember<Exp_ret>{};
		}
		namespace ThreeDim
		{
			class X_m_ODE : public ISysMember<Ox_ret> {};
			class Y_m_ODE : public ISysMember<Oy_ret> {};
			class Z_m_ODE : public ISysMember<Oz_ret> {};
		}
	}
}
