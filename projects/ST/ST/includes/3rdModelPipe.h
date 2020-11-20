/* This header contains the class IAggregate that lists members of ODE system for 3rd Ischemic Model and 
and ISolver-member-function definitions specifically for the system listed in IAggregate + 
subclasses of ISolver that serve for implementations of various calculation methods.*/

#pragma once
#include <models/3rd/settings.h>

namespace StraightTask
{
	class IAggregate {
	public:

		// list of equations in the processed system
		Neurons::NecroticCells::M_ODE NEC;
		Neurons::AcuteChanges::M_ODE AC;
		Neurons::IntactCells::M_ODE HEL;

		Cytokines::Pro_Inflam::M_ODE CY;
		Adhesion::M_ODE ADH;
		LeuMacrophags::M_ODE LM;
		LeuNeutrophils::M_ODE LN;

		Microglia::Active::M_ODE MIA;
		Microglia::Inactive::M_ODE MII;


		// subbordinate values relevant for a system 
		ToxDamage::Full::M_Sub DF;
		ToxDamage::Apop_partial::M_Sub DP_A;
		ToxDamage::Nec_partial::M_Sub DP_N;
		Phagocytosis::Strong::M_Sub EPS_S;
		Phagocytosis::Weak::M_Sub EPS_W;

		std::array<std::function<void(variables&)>, 5> SubValAssigner = {
			[&](variables& X)->void {X.d_F = DF.RP.Expression(X); },
			[&](variables& X)->void {X.dp_A = DP_A.RP.Expression(X); },
			[&](variables& X)->void {X.dp_N = DP_N.RP.Expression(X); },
			[&](variables& X)->void {X.eps_s = EPS_S.RP.Expression(X); },
			[&](variables& X)->void {X.eps_w = EPS_W.RP.Expression(X); }
		};
	};
}

/* including all common instruments stated in ISolver template class;
it is compulsory that IAggregate class is defined above this include directive */
#pragma once
#include<base/Solver_base.h>

namespace StraightTask
{

#define	T_Func(RETURN_TYPE)\
template<typename Method>\
RETURN_TYPE ISolver<Method>::


	T_Func(void)RetUpload(uint32_t Nj) {
		Mthd.X_pred.ret.adh_4 = ADH.GetRetValue(4.0, Nj);
	}

	T_Func(void)RetDataUpdate(uint32_t Nj) {
		ADH.ShiftRets(Nj, Mthd.X_prev.adh);// */
	}

	T_Func(void)CollectData() {

		/*NEC.CollectSolData();
		AC.CollectSolData();
		HEL.CollectSolData(); //*/

		//CY.CollectSolData();
		/*ADH.CollectSolData();

		LM.CollectSolData();
		LN.CollectSolData();
		MIA.CollectSolData();
		MII.CollectSolData();

		//DF.CollectSolData();
		//EPS_S.CollectSolData(); //*/
	}

	T_Func(void)SetIniData(vector<float_t>& T0s) {

		T0s.clear();

		NEC.SetInitialData(Mthd.X_init.tj, Mthd.X_init.nec); T0s.emplace_back(Mthd.X_init.tj);
		AC.SetInitialData(Mthd.X_init.tj, Mthd.X_init.acu_c);  T0s.emplace_back(Mthd.X_init.tj);
		HEL.SetInitialData(Mthd.X_init.tj, Mthd.X_init.hel);  T0s.emplace_back(Mthd.X_init.tj);

		CY.SetInitialData(Mthd.X_init.tj, Mthd.X_init.cy);  T0s.emplace_back(Mthd.X_init.tj);
		ADH.SetInitialData(Mthd.X_init.tj, Mthd.X_init.adh);  T0s.emplace_back(Mthd.X_init.tj);

		LM.SetInitialData(Mthd.X_init.tj, Mthd.X_init.lm);  T0s.emplace_back(Mthd.X_init.tj);
		LN.SetInitialData(Mthd.X_init.tj, Mthd.X_init.ln);  T0s.emplace_back(Mthd.X_init.tj);

		MIA.SetInitialData(Mthd.X_init.tj, Mthd.X_init.mia);  T0s.emplace_back(Mthd.X_init.tj);
		MII.SetInitialData(Mthd.X_init.tj, Mthd.X_init.mii);  T0s.emplace_back(Mthd.X_init.tj);
	}

	T_Func(void)SetIniDataFromOutside(vector<float_t>& T0s) {

		T0s.clear();

		NEC.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.nec);  T0s.emplace_back(Mthd.X_init.tj);
		AC.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.acu_c); T0s.emplace_back(Mthd.X_init.tj);
		HEL.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.hel); T0s.emplace_back(Mthd.X_init.tj);

		CY.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.cy); T0s.emplace_back(Mthd.X_init.tj);
		ADH.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.adh); T0s.emplace_back(Mthd.X_init.tj);

		LM.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.lm); T0s.emplace_back(Mthd.X_init.tj);
		LN.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.ln); T0s.emplace_back(Mthd.X_init.tj);

		MIA.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.mia); T0s.emplace_back(Mthd.X_init.tj);
		MII.SetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.mii); T0s.emplace_back(Mthd.X_init.tj);
	}

	T_Func(void)InitialiseRetArrays() {
		ADH.StateRetArray(Mthd.N, Mthd.t_0, Mthd.gap_width);
	}

	T_Func(void)InitialiseIniRetValues() {
		Mthd.X_init.ret.adh_4 = ADH.GetRetValue(4.0, 0);
	}

	T_Func(void)AssignSolData(uint16_t gap, uint32_t Nj, variables& X) {

		/*X.nec = NEC.GetSolData(gap, Nj, Mthd.N);
		X.acu_c = AC.GetSolData(gap, Nj, Mthd.N);
		X.hel = HEL.GetSolData(gap, Nj, Mthd.N); //*/

		//X.cy = CY.GetSolData(gap, Nj, Mthd.N);
		/*X.adh = ADH.GetSolData(gap, Nj, Mthd.N);

		X.lm = LM.GetSolData(gap, Nj, Mthd.N);
		X.ln = LN.GetSolData(gap, Nj, Mthd.N);

		X.mia = MIA.GetSolData(gap, Nj, Mthd.N);
		X.mii = MII.GetSolData(gap, Nj, Mthd.N);

		//X.d_F = DF.GetSolData(gap, Nj, Mthd.N);
		//X.eps_s = EPS_S.GetSolData(gap, Nj, Mthd.N); // */
	}

	T_Func(void)ExpressSubValues(variables& X) {

		X.d_F = DF.RP.Expression(X);
		X.dp_A = DP_A.RP.Expression(X);
		X.dp_N = DP_N.RP.Expression(X);
		X.eps_s = EPS_S.RP.Expression(X);
		X.eps_w = EPS_W.RP.Expression(X);

	}

	T_Func(void)AllocateOutputStreams() {

		NEC.AllocateOutputStreams();
		AC.AllocateOutputStreams();
		HEL.AllocateOutputStreams();

		CY.AllocateOutputStreams();
		ADH.AllocateOutputStreams();
		LM.AllocateOutputStreams();
		LN.AllocateOutputStreams();
		MIA.AllocateOutputStreams();
		MII.AllocateOutputStreams();

		DF.AllocateOutputStreams();
		DP_N.AllocateOutputStreams();
		DP_A.AllocateOutputStreams();
		EPS_S.AllocateOutputStreams();
		EPS_W.AllocateOutputStreams(); // 14 */
	}

	T_Func(void)OutputSolution(uint32_t Nj, float_t Tj, variables const& X) {

		NEC.OutputSol(Nj, Tj, X.nec);
		AC.OutputSol(Nj, Tj, X.acu_c);
		HEL.OutputSol(Nj, Tj, X.hel);

		CY.OutputSol(Nj, Tj, X.cy);
		ADH.OutputSol(Nj, Tj, X.adh);

		LM.OutputSol(Nj, Tj, X.lm);
		LN.OutputSol(Nj, Tj, X.ln);

		MIA.OutputSol(Nj, Tj, X.mia);
		MII.OutputSol(Nj, Tj, X.mii);

		DF.OutputSol(Nj, Tj, X.d_F);
		DP_N.OutputSol(Nj, Tj, X.dp_N);
		DP_A.OutputSol(Nj, Tj, X.dp_A);
		EPS_S.OutputSol(Nj, Tj, X.eps_s);
		EPS_W.OutputSol(Nj, Tj, X.eps_w);

	}

	T_Func(void)OutputBudgets(uint32_t Nj, float_t Tj) {

		NEC.OutputBuds(Nj, Tj);
		AC.OutputBuds(Nj, Tj);
		HEL.OutputBuds(Nj, Tj);

		CY.OutputBuds(Nj, Tj);
		ADH.OutputBuds(Nj, Tj);

		LM.OutputBuds(Nj, Tj);
		LN.OutputBuds(Nj, Tj);

		MIA.OutputBuds(Nj, Tj);
		MII.OutputBuds(Nj, Tj);

		DF.OutputBuds(Nj, Tj);
		EPS_S.OutputBuds(Nj, Tj);

	}

	T_Func(void)DeallocateOutputStreams() {
		NEC.DeallocateOutputStreams();
		AC.DeallocateOutputStreams();
		HEL.DeallocateOutputStreams();

		CY.DeallocateOutputStreams();
		ADH.DeallocateOutputStreams();
		LM.DeallocateOutputStreams();
		LN.DeallocateOutputStreams();
		MIA.DeallocateOutputStreams();
		MII.DeallocateOutputStreams();

		DF.DeallocateOutputStreams();
		DP_N.DeallocateOutputStreams();
		DP_A.DeallocateOutputStreams();
		EPS_S.DeallocateOutputStreams();
		EPS_W.DeallocateOutputStreams(); // 14 */
	}

	T_Func(bool)is_SYS_deflecting() {
		bool if_ret = NEC.RP.ret_is + AC.RP.ret_is + HEL.RP.ret_is +
			CY.RP.ret_is + ADH.RP.ret_is + LM.RP.ret_is + LN.RP.ret_is +
			MIA.RP.ret_is + MII.RP.ret_is;
		return if_ret;
	}

#undef T_Func
	//--------------------------------------------------------------------------------------------------------------------------------------------

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { Mthd.X_sol = &Mthd.X_pred; }

		void ApplyMethod() final {

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_prev.hel, HEL.RP);//*/

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_prev.mii, MII.RP);// */

			ExpressSubValues(Mthd.X_pred);// */
		}

	};

	class PredCor : public ISolver<Methods::ModEuler>
	{
	public:
		PredCor() { Mthd.X_sol = &Mthd.X_cor; }
		void ApplyMethod() final
		{
			ApplyEuler();

			Mthd.X_cor.tj = Mthd.X_pred.tj;
			Mthd.X_cor.ret = Mthd.X_pred.ret;

			Mthd.X_cor.nec = Mthd.Corrector(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_cor.acu_c = Mthd.Corrector(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_cor.hel = Mthd.Corrector(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.Corrector(Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.Corrector(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.Corrector(Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.Corrector(Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.Corrector(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.Corrector(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_cor);
		}

		

	private:
		void ApplyEuler() {
			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}


	};

	class RunKut : public ISolver<Methods::RunKut4> {
	public:
		RunKut() { *Mthd.X_sol = Mthd.X_pred; }
		void ApplyMethod() final
		{
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_sub.acu_c, Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_sub.hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_sub.cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_sub.adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_sub.lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_sub.ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_sub.mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_sub.mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred); // */
		}

	};

	template<typename Method>
	class IMultistepMethod : public ISolver<Method>{

	public:
		void ApplyPrepStep(uint32_t& Nj, float_t& Tj) final {
			for (; Nj < this->Mthd.X.size(); Nj++)
			{
				// shifting independent variable on one step further for predicted solution
				Tj = this->Mthd.X_pred.tj += this->Mthd.H;

				// setting ret-values for Tj-time-moment in X_pred
				if (this->is_SYS_deflecting()) this->RetUpload(Nj);

				this->ApplyPrepMethod();

				this->Mthd.X[Nj] = this->Mthd.X_pred;

				// outputting solution in current Tj - time-moment
				this->OutputSolution(Nj, Tj, this->Mthd.X_pred);

				// updating RetArray(s) pushing X_prev.value of solution
				this->RetDataUpdate(Nj);

				// shifting X_prev next step further in one-step-methods and
				this->Mthd.X_prev = this->Mthd.X_pred;
			}
		}
	
	};

	class Gear : public IMultistepMethod<Methods::RKGear> {

	public:
		Gear() { *Mthd.X_sol = Mthd.X_cor; }

		void ApplyMethod() final
		{
			ApplyPred();

			Mthd.X_cor.tj = Mthd.X_pred.tj;
			Mthd.X_cor.ret = Mthd.X_pred.ret;

			Mthd.X_cor.nec = Mthd.Corrector(Mthd.X[1].nec, Mthd.X[2].nec, Mthd.X[3].nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_cor.acu_c = Mthd.Corrector(Mthd.X[1].acu_c, Mthd.X[2].acu_c, Mthd.X[3].acu_c, Mthd.X_prev.acu_c, AC.RP);

			Mthd.X_cor.hel = Mthd.Corrector(Mthd.X[1].hel, Mthd.X[2].hel, Mthd.X[3].hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.Corrector(Mthd.X[1].cy, Mthd.X[2].cy, Mthd.X[3].cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.Corrector(Mthd.X[1].adh, Mthd.X[2].adh, Mthd.X[3].adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.Corrector(Mthd.X[1].lm, Mthd.X[2].lm, Mthd.X[3].lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.Corrector(Mthd.X[1].ln, Mthd.X[2].ln, Mthd.X[3].ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.Corrector(Mthd.X[1].mia, Mthd.X[2].mia, Mthd.X[3].mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.Corrector(Mthd.X[1].mii, Mthd.X[2].mii, Mthd.X[3].mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_cor);
		}

		void NodeShift() final {
			Mthd.X[1] = Mthd.X[2]; Mthd.X[2] = Mthd.X[3];
			Mthd.X[3] = Mthd.X_prev; Mthd.X_prev = *Mthd.X_sol;
		}

	private:
		void ApplyPred() {
			Mthd.X_pred.nec = Mthd.GPred(Mthd.X[1].nec, Mthd.X[2].nec, Mthd.X[3].nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.GPred(Mthd.X[1].acu_c, Mthd.X[2].acu_c, Mthd.X[3].acu_c, Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.GPred(Mthd.X[1].hel, Mthd.X[2].hel, Mthd.X[3].hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.GPred(Mthd.X[1].cy, Mthd.X[2].cy, Mthd.X[3].cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.GPred(Mthd.X[1].adh, Mthd.X[2].adh, Mthd.X[3].adh, Mthd.X_prev.adh, ADH.RP);
			Mthd.X_pred.lm = Mthd.GPred(Mthd.X[1].lm, Mthd.X[2].lm, Mthd.X[3].lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.GPred(Mthd.X[1].ln, Mthd.X[2].ln, Mthd.X[3].ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.GPred(Mthd.X[1].mia, Mthd.X[2].mia, Mthd.X[3].mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.GPred(Mthd.X[1].mii, Mthd.X[2].mii, Mthd.X[3].mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}

		void ApplyPrepMethod() final {
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_sub.acu_c, Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_sub.hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_sub.cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_sub.adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_sub.lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_sub.ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_sub.mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_sub.mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred); // */
		}
	};

	class Adams : public IMultistepMethod<Methods::Adams> {
	public:
		Adams() { Mthd.X_sol = &Mthd.X_pred; }

		void ApplyMethod() {
			Mthd.X_pred.nec = Mthd.A_Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.A_Predictor(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.A_Predictor(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.A_Predictor(Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.A_Predictor(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.A_Predictor(Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.A_Predictor(Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.A_Predictor(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.A_Predictor(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}

		void NodeShift() final {
			Mthd.X[1] = Mthd.X[2]; Mthd.X[2] = Mthd.X[3];
			Mthd.X[3] = Mthd.X_prev; Mthd.X_prev = *Mthd.X_sol;
		}

	private:
		void ApplyPrepMethod() final
		{
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_sub.acu_c, Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_sub.hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_sub.cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_sub.adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_sub.lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_sub.ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_sub.mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_sub.mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}

	};

	class ABM : public IMultistepMethod<Methods::ABM> {
	public:
		ABM() { Mthd.X_sol = &Mthd.X_cor; }

		void ApplyMethod() final
		{
			ApplyAdams();

			Mthd.X_cor.tj = Mthd.X_pred.tj;
			Mthd.X_cor.ret = Mthd.X_pred.ret;

			Mthd.X_cor.nec = Mthd.A_Corrector(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_cor.acu_c = Mthd.A_Corrector(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_cor.hel = Mthd.A_Corrector(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.A_Corrector(Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.A_Corrector(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.A_Corrector(Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.A_Corrector(Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.A_Corrector(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.A_Corrector(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_cor);
		}

		void NodeShift() final {
			Mthd.X[1] = Mthd.X[2]; Mthd.X[2] = Mthd.X[3];
			Mthd.X[3] = Mthd.X_prev; Mthd.X_prev = *Mthd.X_sol;
		}

	private:
		void ApplyAdams() {
			Mthd.X_pred.nec = Mthd.A_Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.A_Predictor(Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.A_Predictor(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.A_Predictor(Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.A_Predictor(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.A_Predictor(Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.A_Predictor(Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.A_Predictor(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.A_Predictor(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}

		void ApplyPrepMethod() final {
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.acu_c = Mthd.Predictor(Mthd.X_sub.acu_c, Mthd.X_prev.acu_c, AC.RP);
			Mthd.X_pred.hel = Mthd.Predictor(Mthd.X_sub.hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.Predictor(Mthd.X_sub.cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.Predictor(Mthd.X_sub.adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.Predictor(Mthd.X_sub.lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.Predictor(Mthd.X_sub.ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.Predictor(Mthd.X_sub.mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.Predictor(Mthd.X_sub.mii, Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}
	};
}
