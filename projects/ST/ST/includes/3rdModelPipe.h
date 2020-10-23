/* This header contains the class IAggregate that lists members of ODE system for 3rd Ischemic Model and 
and ISolver-member-function definitions specifically for the system listed in IAggregate + 
subclasses of ISolver that serve for implementations of various calculation methods*/

#pragma once
#include <models/3rd/settings.h>
//#include <functional>

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

	};
}

/* including all common instruments stated in ISolver template class;
it is compulsory that IAggregate class is defined above this include directive */
#include<base/Solver_base.h>

namespace StraightTask
{

#define	T_Func(RETURN_TYPE)\
template<typename Method>\
RETURN_TYPE ISolver<Method>::


	T_Func(void)RetUpload(uint32_t Nj) {
		ST.X_pred.ret.adh_4 = ADH.GetRetValue(4.0, Nj);
	}

	T_Func(void)RetDataUpdate(uint32_t Nj) {
		ADH.ShiftRets(Nj, ST.X_prev.adh);// */
	}

	T_Func(void)CollectData() {

		NEC.CollectSolData();
		AC.CollectSolData();
		HEL.CollectSolData();

		CY.CollectSolData();
		ADH.CollectSolData();

		LM.CollectSolData();
		LN.CollectSolData();
		MIA.CollectSolData();
		MII.CollectSolData();

		DF.CollectSolData();
		EPS_S.CollectSolData();
	}

	T_Func(void)SetIniData(vector<double_t>& T0s) {

		T0s.clear();

		NEC.GetInitialData(ST.X_init.tj, ST.X_init.nec); T0s.emplace_back(ST.X_init.tj);
		AC.GetInitialData(ST.X_init.tj, ST.X_init.acu_c);  T0s.emplace_back(ST.X_init.tj);
		HEL.GetInitialData(ST.X_init.tj, ST.X_init.hel);  T0s.emplace_back(ST.X_init.tj);

		CY.GetInitialData(ST.X_init.tj, ST.X_init.cy);  T0s.emplace_back(ST.X_init.tj);
		ADH.GetInitialData(ST.X_init.tj, ST.X_init.adh);  T0s.emplace_back(ST.X_init.tj);

		LM.GetInitialData(ST.X_init.tj, ST.X_init.lm);  T0s.emplace_back(ST.X_init.tj);
		LN.GetInitialData(ST.X_init.tj, ST.X_init.ln);  T0s.emplace_back(ST.X_init.tj);

		MIA.GetInitialData(ST.X_init.tj, ST.X_init.mia);  T0s.emplace_back(ST.X_init.tj);
		MII.GetInitialData(ST.X_init.tj, ST.X_init.mii);  T0s.emplace_back(ST.X_init.tj);
	}

	T_Func(void)SetIniDataFromOutside(vector<double_t>& T0s) {

		T0s.clear();

		NEC.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.nec);  T0s.emplace_back(ST.X_init.tj);
		AC.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.acu_c); T0s.emplace_back(ST.X_init.tj);
		HEL.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.hel); T0s.emplace_back(ST.X_init.tj);

		CY.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.cy); T0s.emplace_back(ST.X_init.tj);
		ADH.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.adh); T0s.emplace_back(ST.X_init.tj);

		LM.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.lm); T0s.emplace_back(ST.X_init.tj);
		LN.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.ln); T0s.emplace_back(ST.X_init.tj);

		MIA.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.mia); T0s.emplace_back(ST.X_init.tj);
		MII.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.mii); T0s.emplace_back(ST.X_init.tj);
	}

	T_Func(void)InitialiseRetArrays() {
		ADH.StateRetArray(ST.N, ST.t_0, ST.gap_width);
	}

	T_Func(void)InitialiseIniRetValues() {
		ST.X_init.ret.adh_4 = ADH.GetRetValue(4.0, 0);
	}

	T_Func(void)AssignSolData(uint16_t day, uint32_t Nj, variables& X) {

		X.nec = NEC.GetSolData(day, Nj, ST.N);
		X.acu_c = AC.GetSolData(day, Nj, ST.N);
		X.hel = HEL.GetSolData(day, Nj, ST.N);

		X.cy = CY.GetSolData(day, Nj, ST.N);
		X.adh = ADH.GetSolData(day, Nj, ST.N);

		X.lm = LM.GetSolData(day, Nj, ST.N);
		X.ln = LN.GetSolData(day, Nj, ST.N);

		X.mia = MIA.GetSolData(day, Nj, ST.N);
		X.mii = MII.GetSolData(day, Nj, ST.N);

		X.d_F = DF.GetSolData(day, Nj, ST.N);
		X.eps_s = EPS_S.GetSolData(day, Nj, ST.N);
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

	T_Func(void)OutputSolution(uint32_t Nj, double_t Tj, variables const& X) {

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

	T_Func(void)OutputBudgets(uint32_t Nj, double_t Tj) {

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
		Euler() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}

	};

	class PredCor : public ISolver<Methods::ModEuler>
	{
	public:
		PredCor() { ST.X_sol = &ST.X_cor; }
		void ApplyMethod() final
		{
			ApplyEuler();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.nec = ST.Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.Corrector(ST.X_prev.acu_c, AC.RP);
			ST.X_cor.hel = ST.Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_cor);
		}

	private:
		void ApplyEuler() {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}
	};

	class RunKut : public ISolver<Methods::RunKut4> {
	public:
		RunKut() { *ST.X_sol = ST.X_pred; }
		void ApplyMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred); // */
		}
	};

	template<typename Method>
	class IMultistepMethod : public ISolver<Method>{

	public:
		void ApplyPrepStep(uint32_t& Nj, double_t& Tj) final {
			for (; Nj < this->ST.X.size(); Nj++)
			{
				// shifting independent variable on one step further for predicted solution
				Tj = this->ST.X_pred.tj += this->ST.H;

				// setting ret-values for Tj-time-moment in X_pred
				if (this->is_SYS_deflecting()) this->RetUpload(Nj);

				this->ApplyPrepMethod();

				this->ST.X[Nj] = this->ST.X_pred;

				// outputting solution in current Tj - time-moment
				this->OutputSolution(Nj, Tj, this->ST.X_pred);

				// updating RetArray(s) pushing X_prev.value of solution
				this->RetDataUpdate(Nj);

				// shifting X_prev next step further in one-step-methods and
				this->ST.X_prev = this->ST.X_pred;
			}
		}
	
	};

	class Gear : public IMultistepMethod<Methods::RKGear> {

	public:
		Gear() { *ST.X_sol = ST.X_cor; }

		void ApplyMethod() final
		{
			ApplyPred();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.nec = ST.Corrector(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.Corrector(ST.X[1].acu_c, ST.X[2].acu_c, ST.X[3].acu_c, ST.X_prev.acu_c, AC.RP);

			ST.X_cor.hel = ST.Corrector(ST.X[1].hel, ST.X[2].hel, ST.X[3].hel, ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X[1].cy, ST.X[2].cy, ST.X[3].cy, ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X[1].adh, ST.X[2].adh, ST.X[3].adh, ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X[1].lm, ST.X[2].lm, ST.X[3].lm, ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X[1].ln, ST.X[2].ln, ST.X[3].ln, ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X[1].mia, ST.X[2].mia, ST.X[3].mia, ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X[1].mii, ST.X[2].mii, ST.X[3].mii, ST.X_prev.mii, MII.RP);

			ST.X_cor.ret = ST.X_pred.ret; // in case of emergency to update ret-values on X_cor, but what no earth for

			ExpressSubValues(ST.X_cor);
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPred() {
			ST.X_pred.nec = ST.GPred(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.GPred(ST.X[1].acu_c, ST.X[2].acu_c, ST.X[3].acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.GPred(ST.X[1].hel, ST.X[2].hel, ST.X[3].hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.GPred(ST.X[1].cy, ST.X[2].cy, ST.X[3].cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.GPred(ST.X[1].adh, ST.X[2].adh, ST.X[3].adh, ST.X_prev.adh, ADH.RP);
			ST.X_pred.lm = ST.GPred(ST.X[1].lm, ST.X[2].lm, ST.X[3].lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.GPred(ST.X[1].ln, ST.X[2].ln, ST.X[3].ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.GPred(ST.X[1].mia, ST.X[2].mia, ST.X[3].mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.GPred(ST.X[1].mii, ST.X[2].mii, ST.X[3].mii, ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}

		void ApplyPrepMethod() final {
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}
	};

	class Adams : public IMultistepMethod<Methods::Adams> {
	public:
		Adams() { ST.X_sol = &ST.X_pred; }

		void ApplyMethod() {
			ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.A_Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}

	};

	class ABM : public IMultistepMethod<Methods::ABM> {
	public:
		ABM() { ST.X_sol = &ST.X_cor; }

		void ApplyMethod() final
		{
			ApplyAdams();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.nec = ST.A_Corrector(ST.X_prev.nec, NEC.RP);
			ST.X_cor.acu_c = ST.A_Corrector(ST.X_prev.acu_c, AC.RP);
			ST.X_cor.hel = ST.A_Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.A_Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.A_Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.A_Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.A_Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.A_Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.A_Corrector(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_cor);
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyAdams() {
			ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.A_Predictor(ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.A_Predictor(ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.A_Predictor(ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.A_Predictor(ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.A_Predictor(ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.A_Predictor(ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.A_Predictor(ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.A_Predictor(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}

		void ApplyPrepMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.acu_c = ST.Predictor(ST.X_sub.acu_c, ST.X_prev.acu_c, AC.RP);
			ST.X_pred.hel = ST.Predictor(ST.X_sub.hel, ST.X_prev.hel, HEL.RP);

			ST.X_pred.cy = ST.Predictor(ST.X_sub.cy, ST.X_prev.cy, CY.RP);
			ST.X_pred.adh = ST.Predictor(ST.X_sub.adh, ST.X_prev.adh, ADH.RP);

			ST.X_pred.lm = ST.Predictor(ST.X_sub.lm, ST.X_prev.lm, LM.RP);
			ST.X_pred.ln = ST.Predictor(ST.X_sub.ln, ST.X_prev.ln, LN.RP);

			ST.X_pred.mia = ST.Predictor(ST.X_sub.mia, ST.X_prev.mia, MIA.RP);
			ST.X_pred.mii = ST.Predictor(ST.X_sub.mii, ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_pred);// */
		}
	};
}
