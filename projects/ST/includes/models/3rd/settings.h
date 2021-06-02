/* This header implements the IMember interface for 3rd-model ODE_system members

M_ODE stands for ODE_MEMBER */

#pragma once
#include <models/3rd/equations.h>
#include <base/Settings_base.h>


namespace StraightTask
{

	template<typename RightPart>
	class ISysMember : public IMember<RightPart>, public IErr{};

	namespace Neurons
	{
		namespace NecroticCells
		{
			class M_ODE : public ISysMember<Equation>	{
				const char* sol_name()noexcept final { return "Necr"; };
				const char* ini_spl_name() final { return "necr4SPL"; };
			public:
				M_ODE() { CollectExpData(); }
			};
		}

		namespace AcuteChanges
		{
			class M_ODE : public ISysMember<Equation>	{
				const char* sol_name()noexcept final { return "Ac_ch"; };
				const char* ini_spl_name() final { return "apop4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}
		
		namespace IntactCells
		{
			class M_ODE : public ISysMember<Equation>{
				const char* sol_name()noexcept final { return "Healt"; };
				const char* ini_spl_name() final { return "hel4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}
	}

	namespace Cytokines
	{
		namespace Pro_Inflam
		{
			class M_ODE : public ISysMember<Equation> {
				const char* sol_name()noexcept final { return "Cy"; };
				const char* ini_spl_name() final { return "cyto4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}
	}

	namespace Adhesion
	{
		class M_ODE : public ISysMember<Equation> {
			const char* sol_name()noexcept final { return "Adhes"; };
		};
	}

	namespace LeuMacrophags
	{
		class M_ODE : public ISysMember<Equation> {
			const char* sol_name()noexcept final { return "Lm"; };
			const char* ini_spl_name() final { return "1lm4SPL"; }
		public:
			M_ODE() { CollectExpData(); }
		};
	}

	namespace LeuNeutrophils
	{
		class M_ODE : public ISysMember<Equation> {
			const char* sol_name()noexcept final { return "Ln"; };
			const char* ini_spl_name() final { return "ln4SPL"; }
		public:
			M_ODE() { CollectExpData(); }
		};
	}

	namespace Microglia
	{
		namespace Active
		{
			class M_ODE : public ISysMember<Equation> {
				const char* sol_name()noexcept final { return "Mi_active"; };
			};
		}
		namespace Inactive
		{
			class M_ODE : public ISysMember<Equation> {
				const char* sol_name()noexcept final { return "Mi_inactive"; };
			};
		}
	}

	template<typename RightPart>
	class ISubMember : public IMember<RightPart>{};

	namespace ToxDamage
	{
		namespace Full
		{
			class M_Sub : public ISubMember<SubVal>	{
			public:
				const char* sol_name()noexcept final { return "D_Full"; }
				M_Sub() {}
			};
		}
		namespace Nec_partial
		{
			class M_Sub : public ISubMember<SubVal> {
			public:
				const char* sol_name()noexcept final { return "DN_c"; }
			};
		}
		namespace Apop_partial
		{
			class M_Sub : public ISubMember<SubVal> {
			public:
				const char* sol_name()noexcept final { return "DA_c"; }
			};
		}
	}

	namespace Phagocytosis 
	{
		namespace Strong
		{
			class M_Sub : public ISubMember<SubVal>	{
			public:
				const char* sol_name()noexcept final { return "Eps_s"; }
			};
		}

		namespace Weak
		{
			class M_Sub : public ISubMember<SubVal>	{
			public:
				const char* sol_name()noexcept final { return "Eps_w"; }
			};
		}
	}
}

