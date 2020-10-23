/* This header contains the class IAggregate that lists members of ODE system for Original Ischemic Model and
and ISolver-member-function definitions specifically for the system listed in IAggregate +
subclasses of ISolver that serve for implementations of various calculation methods*/

#pragma once
#include <models/1st/settings.h>
//#include <functional>

namespace StraightTask
{
	class IAggregate
	{
	public:

		// list of equations in the processed system
		Neurons::NecroticCells::M_ODE NEC;
		Neurons::Apoptosis::Started::M_ODE AS;
		Neurons::Apoptosis::Ended::M_ODE AE;
		Neurons::IntactCells::M_ODE HEL;

		Cytokines::Pro_Inflam::M_ODE CY;
		Adhesion::M_ODE ADH;
		LeuMacrophags::M_ODE LM;
		LeuNeutrophils::M_ODE LN;

		Microglia::Active::M_ODE MIA;
		Microglia::Inactive::M_ODE MII;


		// subbordinate values relevant for the system listed above
		Cytokines::M_Sub PSY;
		ToxDamage::Full::M_Sub DF;
		ToxDamage::Initial::M_Sub D_INI;

		Phagocytosis::Strong::M_Sub EPS_S;

	};

}

/* including all common instruments stated in ISolver template class;
it is compulsory that IAggregate class is defined above this include directive */
#include<base/Solver_base.h>

namespace StraightTask
{

//=================================================================================================================
#define	T_Func(TEMPLATE_NAME, RETURN_TYPE)\
template<typename TEMPLATE_NAME>\
RETURN_TYPE ISolver<TEMPLATE_NAME>::


	T_Func(Method, void)RetUpload(uint32_t Nj) {
		ST.X_pred.ret.hel_12 = HEL.GetRetValue(12.0, Nj);
		ST.X_pred.ret.adh_12 = ADH.GetRetValue(12.0, Nj);
		ST.X_pred.ret.adh_24 = ADH.GetRetValue(24.0, Nj);
		ST.X_pred.ret.d_12 = D_INI.GetRetValue(12.0, Nj);
	}

	T_Func(Method, void)RetDataUpdate(uint32_t Nj) {
		HEL.ShiftRets(Nj, ST.X_prev.hel);
		D_INI.ShiftRets(Nj, ST.X_prev.d_ini);
		ADH.ShiftRets(Nj, ST.X_prev.adh);
	}

	T_Func(Method, void)CollectData() {

		NEC.CollectSolData();
		AE.CollectSolData(); // this container is rather empty
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

	T_Func(Method, void)SetIniData(vector<double_t>& T0s) {

		T0s.clear();

		NEC.GetInitialData(ST.X_init.tj, ST.X_init.nec); T0s.emplace_back(ST.X_init.tj);
		AS.GetInitialData(ST.X_init.tj, ST.X_init.ap_s);  T0s.emplace_back(ST.X_init.tj);
		AE.GetInitialData(ST.X_init.tj, ST.X_init.ap_e);  T0s.emplace_back(ST.X_init.tj);
		HEL.GetInitialData(ST.X_init.tj, ST.X_init.hel);  T0s.emplace_back(ST.X_init.tj);

		CY.GetInitialData(ST.X_init.tj, ST.X_init.cy);  T0s.emplace_back(ST.X_init.tj);
		ADH.GetInitialData(ST.X_init.tj, ST.X_init.adh);  T0s.emplace_back(ST.X_init.tj);

		LM.GetInitialData(ST.X_init.tj, ST.X_init.lm);  T0s.emplace_back(ST.X_init.tj);
		LN.GetInitialData(ST.X_init.tj, ST.X_init.ln);  T0s.emplace_back(ST.X_init.tj);

		MIA.GetInitialData(ST.X_init.tj, ST.X_init.mia);  T0s.emplace_back(ST.X_init.tj);
		MII.GetInitialData(ST.X_init.tj, ST.X_init.mii);  T0s.emplace_back(ST.X_init.tj);
	}

	T_Func(Method, void)SetIniDataFromOutside(vector<double_t>& T0s) {

		T0s.clear();

		NEC.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.nec);  T0s.emplace_back(ST.X_init.tj);
		AS.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.ap_s); T0s.emplace_back(ST.X_init.tj);
		AE.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.ap_e); T0s.emplace_back(ST.X_init.tj);
		HEL.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.hel); T0s.emplace_back(ST.X_init.tj);

		CY.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.cy); T0s.emplace_back(ST.X_init.tj);
		ADH.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.adh); T0s.emplace_back(ST.X_init.tj);

		LM.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.lm); T0s.emplace_back(ST.X_init.tj);
		LN.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.ln); T0s.emplace_back(ST.X_init.tj);

		MIA.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.mia); T0s.emplace_back(ST.X_init.tj);
		MII.GetInitialDataFromOutside(ST.X_init.tj, ST.X_init.mii); T0s.emplace_back(ST.X_init.tj);
	}

	T_Func(Method, void)InitialiseRetArrays() {
		ADH.StateRetArray(ST.N, ST.t_0, ST.gap_width);
		HEL.StateRetArray(ST.N, ST.t_0, ST.gap_width);
		D_INI.StateRetArray(ST.N, ST.t_0, ST.gap_width); // */
	}

	T_Func(Method, void)InitialiseIniRetValues() {
		ST.X_init.ret.hel_12 = HEL.GetRetValue(12.0, 0);
		ST.X_init.ret.adh_24 = ADH.GetRetValue(24.0, 0);
		ST.X_init.ret.d_12 = D_INI.GetRetValue(12.0, 0);
		ST.X_init.ret.adh_12 = ADH.GetRetValue(12.0, 0);
	}

	T_Func(Method, void)AssignSolData(uint16_t day, uint32_t Nj, variables& X) {

		X.nec = NEC.GetSolData(day, Nj, ST.N);
		X.ap_e = AE.GetSolData(day, Nj, ST.N); // this container is rather empty
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

	T_Func(Method, void)ExpressSubValues(variables& X) {

		X.d_F = DF.RP.Expression(X);
		X.d_ini = D_INI.RP.Expression(X);
		X.eps_s = EPS_S.RP.Expression(X);
	}

	T_Func(Method, void)AllocateOutputStreams() {

		NEC.AllocateOutputStreams();
		AS.AllocateOutputStreams();
		AE.AllocateOutputStreams();
		HEL.AllocateOutputStreams();

		CY.AllocateOutputStreams();
		ADH.AllocateOutputStreams();

		LM.AllocateOutputStreams();
		LN.AllocateOutputStreams();
		MIA.AllocateOutputStreams();
		MII.AllocateOutputStreams();

		DF.AllocateOutputStreams();
		D_INI.AllocateOutputStreams();
		EPS_S.AllocateOutputStreams();
	}

	T_Func(Method, void)OutputSolution(uint32_t Nj, double_t Tj, variables const& X) {

		NEC.OutputSol(Nj, Tj, X.nec);
		AS.OutputSol(Nj, Tj, X.ap_s);
		AE.OutputSol(Nj, Tj, X.ap_e);
		HEL.OutputSol(Nj, Tj, X.hel);

		CY.OutputSol(Nj, Tj, X.cy);
		ADH.OutputSol(Nj, Tj, X.adh);

		LM.OutputSol(Nj, Tj, X.lm);
		LN.OutputSol(Nj, Tj, X.ln);

		MIA.OutputSol(Nj, Tj, X.mia);
		MII.OutputSol(Nj, Tj, X.mii);

		DF.OutputSol(Nj, Tj, X.d_F);
		D_INI.OutputSol(Nj, Tj, X.d_ini);
		EPS_S.OutputSol(Nj, Tj, X.eps_s);
		
	}

	T_Func(Method, void)OutputBudgets(uint32_t Nj, double_t Tj) {

		NEC.OutputBuds(Nj, Tj);
		AS.OutputBuds(Nj, Tj);
		AE.OutputBuds(Nj, Tj);
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

	T_Func(Method, bool)is_SYS_deflecting() {
		bool if_ret = NEC.RP.ret_is + AS.RP.ret_is + AE.RP.ret_is + HEL.RP.ret_is +
			CY.RP.ret_is + ADH.RP.ret_is + LM.RP.ret_is + LN.RP.ret_is +
			MIA.RP.ret_is + MII.RP.ret_is;
		return if_ret;
	}

#undef T_Func


//----------------------------------------------------------------------------------------------------------------------------------

	

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_prev.ap_e, AE.RP);
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
			ST.X_cor.ap_s = ST.Corrector(ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.Corrector(ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_cor);// */
		}

	private:
		void ApplyEuler() {
			ST.X_pred.nec = ST.Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_prev.ap_e, AE.RP);
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
		RunKut() { ST.X_sol = &ST.X_pred; }
		void ApplyMethod() final
		{
			ST.X_sub = ST.X_prev;

			ST.X_pred.nec = ST.Predictor(ST.X_sub.nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
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
	class IMultistepMethod : public ISolver<Method> {

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
		Gear() { ST.X_sol = &ST.X_cor; }
		

		void ApplyMethod() final
		{
			ApplyPred();

			ST.X_cor.tj = ST.X_pred.tj;
			ST.X_cor.ret = ST.X_pred.ret;

			ST.X_cor.nec = ST.Corrector(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_cor.ap_s = ST.Corrector(ST.X[1].ap_s, ST.X[2].ap_s, ST.X[3].ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.Corrector(ST.X[1].ap_e, ST.X[2].ap_e, ST.X[3].ap_e, ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.Corrector(ST.X[1].hel, ST.X[2].hel, ST.X[3].hel, ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.Corrector(ST.X[1].cy, ST.X[2].cy, ST.X[3].cy, ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.Corrector(ST.X[1].adh, ST.X[2].adh, ST.X[3].adh, ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.Corrector(ST.X[1].lm, ST.X[2].lm, ST.X[3].lm, ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.Corrector(ST.X[1].ln, ST.X[2].ln, ST.X[3].ln, ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.Corrector(ST.X[1].mia, ST.X[2].mia, ST.X[3].mia, ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.Corrector(ST.X[1].mii, ST.X[2].mii, ST.X[3].mii, ST.X_prev.mii, MII.RP);

			ST.X_cor.ret = ST.X_pred.ret; // in case of emergency to update ret-values on X_cor, but what no earth for

			ExpressSubValues(ST.X_cor); // */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyPred() {
			ST.X_pred.nec = ST.GPred(ST.X[1].nec, ST.X[2].nec, ST.X[3].nec, ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.GPred(ST.X[1].ap_s, ST.X[2].ap_s, ST.X[3].ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.GPred(ST.X[1].ap_e, ST.X[2].ap_e, ST.X[3].ap_e, ST.X_prev.ap_e, AE.RP);
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
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
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
			ST.X_pred.ap_s = ST.A_Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.A_Predictor(ST.X_prev.ap_e, AE.RP);
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
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
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
			ST.X_cor.ap_s = ST.A_Corrector(ST.X_prev.ap_s, AS.RP);
			ST.X_cor.ap_e = ST.A_Corrector(ST.X_prev.ap_e, AE.RP);
			ST.X_cor.hel = ST.A_Corrector(ST.X_prev.hel, HEL.RP);

			ST.X_cor.cy = ST.A_Corrector(ST.X_prev.cy, CY.RP);
			ST.X_cor.adh = ST.A_Corrector(ST.X_prev.adh, ADH.RP);

			ST.X_cor.lm = ST.A_Corrector(ST.X_prev.lm, LM.RP);
			ST.X_cor.ln = ST.A_Corrector(ST.X_prev.ln, LN.RP);

			ST.X_cor.mia = ST.A_Corrector(ST.X_prev.mia, MIA.RP);
			ST.X_cor.mii = ST.A_Corrector(ST.X_prev.mii, MII.RP);

			ExpressSubValues(ST.X_cor);// */
		}

		void NodeShift() final {
			ST.X[1] = ST.X[2]; ST.X[2] = ST.X[3];
			ST.X[3] = ST.X_prev; ST.X_prev = *ST.X_sol;
		}

	private:
		void ApplyAdams() {
			ST.X_pred.nec = ST.A_Predictor(ST.X_prev.nec, NEC.RP);
			ST.X_pred.ap_s = ST.A_Predictor(ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.A_Predictor(ST.X_prev.ap_e, AE.RP);
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
			ST.X_pred.ap_s = ST.Predictor(ST.X_sub.ap_s, ST.X_prev.ap_s, AS.RP);
			ST.X_pred.ap_e = ST.Predictor(ST.X_sub.ap_e, ST.X_prev.ap_e, AE.RP);
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