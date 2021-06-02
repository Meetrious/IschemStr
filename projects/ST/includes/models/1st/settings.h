// This header constains whole objects for each ODE_system member from the 1st-model

#pragma once
#include <models/1st/equations.h>
#include <base/Settings_base.h>




namespace StraightTask
{
	template<typename RightPart>
	class ISysMember : public IMember<RightPart>, public IErr {};

	namespace Neurons
	{
		namespace NecroticCells
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name() noexcept final { return "Necr"; };
				const char* ini_spl_name() noexcept final { return "necr4SPL"; };
			public:
				M_ODE() { CollectExpData(); }
			};
		}

		namespace AcuteChanges
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name()noexcept final { return "Ac_ch"; };
				const char* ini_spl_name()noexcept final { return "apop4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}

		namespace Apoptosis
		{
			namespace Started
			{
				class M_ODE : public ISysMember<RightPart> {
					const char* sol_name()noexcept final { return "Ap_s"; };
				};
			}
			namespace Ended
			{
				class M_ODE : public ISysMember<RightPart> {
					const char* sol_name()noexcept final { return "Ap_e"; };
					const char* ini_spl_name()noexcept final { return "apop4SPL"; }
				};
			}
		}

		namespace IntactCells
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name()noexcept final { return "Healt"; };
				const char* ini_spl_name()noexcept final { return "hel4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}
	}

	namespace Cytokines
	{
		namespace Pro_Inflam
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name()noexcept final { return "Cy"; };
				const char* ini_spl_name()noexcept final { return "cyto4SPL"; }
			public:
				M_ODE() { CollectExpData(); }
			};
		}
	}

	namespace Adhesion
	{
		class M_ODE : public ISysMember<RightPart> {
			const char* sol_name()noexcept final { return "Adhes"; };
		};
	}

	namespace LeuMacrophags
	{
		class M_ODE : public ISysMember<RightPart> {
			const char* sol_name()noexcept final { return "Lm"; };
			const char* ini_spl_name()noexcept final { return "lm4SPL"; }
		public:
			M_ODE() { CollectExpData(); }
		};
	}

	namespace LeuNeutrophils
	{
		class M_ODE : public ISysMember<RightPart> {
			const char* sol_name()noexcept final { return "Ln"; };
			const char* ini_spl_name()noexcept final { return "ln4SPL"; }
		public:
			M_ODE() { CollectExpData(); }
		};
	}

	namespace Microglia
	{
		namespace Active
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name()noexcept final { return "Mi_active"; };
			};
		}
		namespace Inactive
		{
			class M_ODE : public ISysMember<RightPart> {
				const char* sol_name()noexcept final { return "Mi_inactive"; };
			};
		}
	}

	template<typename RightPart>
	class ISubMember : public IMember<RightPart> {};

	namespace Cytokines
	{
		class M_Sub : public ISubMember<Psy>
		{
			const char* sol_name()noexcept final { return "Psy_c"; }
		};
	}

	namespace ToxDamage
	{
		namespace Full
		{
			class M_Sub : public ISubMember<SubVal>
			{
			public:
				const char* sol_name()noexcept final { return "D_Full"; }
				M_Sub() {}
			};
		}
		namespace Initial
		{
			class M_Sub : public ISubMember<SubVal>
			{
			public:
				const char* sol_name()noexcept final { return "D_ini"; }
				const char* ini_spl_name() noexcept final { return "Orig_d_ini"; }
				M_Sub() {}
			};
		}
	}

	namespace Phagocytosis
	{
		namespace Strong
		{
			class M_Sub : public ISubMember<SubVal>
			{
			public:
				const char* sol_name()noexcept final { return "Eps_s"; }
			};
		}
	}

}
