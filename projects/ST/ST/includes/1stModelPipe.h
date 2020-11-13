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
#define	T_Func(RETURN_TYPE)\
template<typename Method>\
RETURN_TYPE ISolver<Method>::


	T_Func(void)RetUpload(uint32_t Nj) {
		Mthd.X_pred.ret.hel_12 = HEL.GetRetValue(12.0, Nj);
		Mthd.X_pred.ret.adh_12 = ADH.GetRetValue(12.0, Nj);
		Mthd.X_pred.ret.adh_24 = ADH.GetRetValue(24.0, Nj);
		Mthd.X_pred.ret.d_12 = D_INI.GetRetValue(12.0, Nj);
	}

	T_Func(void)RetDataUpdate(uint32_t Nj) {
		HEL.ShiftRets(Nj, Mthd.X_prev.hel);
		D_INI.ShiftRets(Nj, Mthd.X_prev.d_ini);
		ADH.ShiftRets(Nj, Mthd.X_prev.adh);
	}

	T_Func(void)CollectData() {

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

	T_Func(void)SetIniData(vector<float_t>& T0s) {

		T0s.clear();

		NEC.GetInitialData(Mthd.X_init.tj, Mthd.X_init.nec); T0s.emplace_back(Mthd.X_init.tj);
		AS.GetInitialData(Mthd.X_init.tj, Mthd.X_init.ap_s);  T0s.emplace_back(Mthd.X_init.tj);
		AE.GetInitialData(Mthd.X_init.tj, Mthd.X_init.ap_e);  T0s.emplace_back(Mthd.X_init.tj);
		HEL.GetInitialData(Mthd.X_init.tj, Mthd.X_init.hel);  T0s.emplace_back(Mthd.X_init.tj);

		CY.GetInitialData(Mthd.X_init.tj, Mthd.X_init.cy);  T0s.emplace_back(Mthd.X_init.tj);
		ADH.GetInitialData(Mthd.X_init.tj, Mthd.X_init.adh);  T0s.emplace_back(Mthd.X_init.tj);

		LM.GetInitialData(Mthd.X_init.tj, Mthd.X_init.lm);  T0s.emplace_back(Mthd.X_init.tj);
		LN.GetInitialData(Mthd.X_init.tj, Mthd.X_init.ln);  T0s.emplace_back(Mthd.X_init.tj);

		MIA.GetInitialData(Mthd.X_init.tj, Mthd.X_init.mia);  T0s.emplace_back(Mthd.X_init.tj);
		MII.GetInitialData(Mthd.X_init.tj, Mthd.X_init.mii);  T0s.emplace_back(Mthd.X_init.tj);
	}

	T_Func(void)SetIniDataFromOutside(vector<float_t>& T0s) {

		T0s.clear();

		NEC.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.nec);  T0s.emplace_back(Mthd.X_init.tj);
		AS.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.ap_s); T0s.emplace_back(Mthd.X_init.tj);
		AE.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.ap_e); T0s.emplace_back(Mthd.X_init.tj);
		HEL.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.hel); T0s.emplace_back(Mthd.X_init.tj);

		CY.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.cy); T0s.emplace_back(Mthd.X_init.tj);
		ADH.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.adh); T0s.emplace_back(Mthd.X_init.tj);

		LM.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.lm); T0s.emplace_back(Mthd.X_init.tj);
		LN.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.ln); T0s.emplace_back(Mthd.X_init.tj);

		MIA.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.mia); T0s.emplace_back(Mthd.X_init.tj);
		MII.GetInitialDataFromOutside(Mthd.X_init.tj, Mthd.X_init.mii); T0s.emplace_back(Mthd.X_init.tj);
	}

	T_Func(void)InitialiseRetArrays() {
		ADH.StateRetArray(Mthd.N, Mthd.t_0, Mthd.gap_width);
		HEL.StateRetArray(Mthd.N, Mthd.t_0, Mthd.gap_width);
		D_INI.StateRetArray(Mthd.N, Mthd.t_0, Mthd.gap_width); // */
	}

	T_Func(void)InitialiseIniRetValues() {
		Mthd.X_init.ret.hel_12 = HEL.GetRetValue(12.0, 0);
		Mthd.X_init.ret.adh_24 = ADH.GetRetValue(24.0, 0);
		Mthd.X_init.ret.d_12 = D_INI.GetRetValue(12.0, 0);
		Mthd.X_init.ret.adh_12 = ADH.GetRetValue(12.0, 0);
	}

	T_Func(void)AssignSolData(uint16_t day, uint32_t Nj, variables& X) {

		X.nec = NEC.GetSolData(day, Nj, Mthd.N);
		X.ap_e = AE.GetSolData(day, Nj, Mthd.N); // this container is rather empty
		X.hel = HEL.GetSolData(day, Nj, Mthd.N);

		X.cy = CY.GetSolData(day, Nj, Mthd.N);
		X.adh = ADH.GetSolData(day, Nj, Mthd.N);

		X.lm = LM.GetSolData(day, Nj, Mthd.N);
		X.ln = LN.GetSolData(day, Nj, Mthd.N);

		X.mia = MIA.GetSolData(day, Nj, Mthd.N);
		X.mii = MII.GetSolData(day, Nj, Mthd.N);

		X.d_F = DF.GetSolData(day, Nj, Mthd.N);
		X.eps_s = EPS_S.GetSolData(day, Nj, Mthd.N);
	}

	T_Func(void)ExpressSubValues(variables& X) {

		X.d_F = DF.RP.Expression(X);
		X.d_ini = D_INI.RP.Expression(X);
		X.eps_s = EPS_S.RP.Expression(X);
	}

	T_Func(void)AllocateOutputStreams() {

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

	T_Func(void)OutputSolution(uint32_t Nj, float_t Tj, variables const& X) {

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

	T_Func(void)OutputBudgets(uint32_t Nj, float_t Tj) {

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

	T_Func(void)DeallocateOutputStreams(){
		NEC.DeallocateOutputStreams();
		AS.DeallocateOutputStreams();
		AE.DeallocateOutputStreams();
		HEL.DeallocateOutputStreams();

		CY.DeallocateOutputStreams();
		ADH.DeallocateOutputStreams();

		LM.DeallocateOutputStreams();
		LN.DeallocateOutputStreams();
		MIA.DeallocateOutputStreams();
		MII.DeallocateOutputStreams();

		DF.DeallocateOutputStreams();
		D_INI.DeallocateOutputStreams();
		EPS_S.DeallocateOutputStreams();
	}

	T_Func(bool)is_SYS_deflecting() {
		bool if_ret = NEC.RP.ret_is + AS.RP.ret_is + AE.RP.ret_is + HEL.RP.ret_is +
			CY.RP.ret_is + ADH.RP.ret_is + LM.RP.ret_is + LN.RP.ret_is +
			MIA.RP.ret_is + MII.RP.ret_is;
		return if_ret;
	}

#undef T_Func


//----------------------------------------------------------------------------------------------------------------------------------

	

	class Euler : public ISolver<Methods::Euler> {
	public:
		Euler() { Mthd.X_sol = &Mthd.X_pred; }
		void ApplyMethod() final {
			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_prev.ap_e, AE.RP);
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
			Mthd.X_cor.ap_s = Mthd.Corrector(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_cor.ap_e = Mthd.Corrector(Mthd.X_prev.ap_e, AE.RP);
			Mthd.X_cor.hel = Mthd.Corrector(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.Corrector(Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.Corrector(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.Corrector(Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.Corrector(Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.Corrector(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.Corrector(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_cor);// */
		}

	private:
		void ApplyEuler() {
			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_prev.ap_e, AE.RP);
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
		RunKut() { Mthd.X_sol = &Mthd.X_pred; }
		void ApplyMethod() final
		{
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_sub.ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_sub.ap_e, Mthd.X_prev.ap_e, AE.RP);
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
	class IMultistepMethod : public ISolver<Method> {

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
		Gear() { Mthd.X_sol = &Mthd.X_cor; }
		

		void ApplyMethod() final
		{
			ApplyPred();

			Mthd.X_cor.tj = Mthd.X_pred.tj;
			Mthd.X_cor.ret = Mthd.X_pred.ret;

			Mthd.X_cor.nec = Mthd.Corrector(Mthd.X[1].nec, Mthd.X[2].nec, Mthd.X[3].nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_cor.ap_s = Mthd.Corrector(Mthd.X[1].ap_s, Mthd.X[2].ap_s, Mthd.X[3].ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_cor.ap_e = Mthd.Corrector(Mthd.X[1].ap_e, Mthd.X[2].ap_e, Mthd.X[3].ap_e, Mthd.X_prev.ap_e, AE.RP);
			Mthd.X_cor.hel = Mthd.Corrector(Mthd.X[1].hel, Mthd.X[2].hel, Mthd.X[3].hel, Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.Corrector(Mthd.X[1].cy, Mthd.X[2].cy, Mthd.X[3].cy, Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.Corrector(Mthd.X[1].adh, Mthd.X[2].adh, Mthd.X[3].adh, Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.Corrector(Mthd.X[1].lm, Mthd.X[2].lm, Mthd.X[3].lm, Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.Corrector(Mthd.X[1].ln, Mthd.X[2].ln, Mthd.X[3].ln, Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.Corrector(Mthd.X[1].mia, Mthd.X[2].mia, Mthd.X[3].mia, Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.Corrector(Mthd.X[1].mii, Mthd.X[2].mii, Mthd.X[3].mii, Mthd.X_prev.mii, MII.RP);

			Mthd.X_cor.ret = Mthd.X_pred.ret; // in case of emergency to update ret-values on X_cor, but what no earth for

			ExpressSubValues(Mthd.X_cor); // */
		}

		void NodeShift() final {
			Mthd.X[1] = Mthd.X[2]; Mthd.X[2] = Mthd.X[3];
			Mthd.X[3] = Mthd.X_prev; Mthd.X_prev = *Mthd.X_sol;
		}

	private:
		void ApplyPred() {
			Mthd.X_pred.nec = Mthd.GPred(Mthd.X[1].nec, Mthd.X[2].nec, Mthd.X[3].nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.GPred(Mthd.X[1].ap_s, Mthd.X[2].ap_s, Mthd.X[3].ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.GPred(Mthd.X[1].ap_e, Mthd.X[2].ap_e, Mthd.X[3].ap_e, Mthd.X_prev.ap_e, AE.RP);
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
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_sub.ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_sub.ap_e, Mthd.X_prev.ap_e, AE.RP);
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

	class Adams : public IMultistepMethod<Methods::Adams> {
	public:
		Adams() { Mthd.X_sol = &Mthd.X_pred; }

		void ApplyMethod() {
			Mthd.X_pred.nec = Mthd.A_Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.A_Predictor(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.A_Predictor(Mthd.X_prev.ap_e, AE.RP);
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
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_sub.ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_sub.ap_e, Mthd.X_prev.ap_e, AE.RP);
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
			Mthd.X_cor.ap_s = Mthd.A_Corrector(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_cor.ap_e = Mthd.A_Corrector(Mthd.X_prev.ap_e, AE.RP);
			Mthd.X_cor.hel = Mthd.A_Corrector(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_cor.cy = Mthd.A_Corrector(Mthd.X_prev.cy, CY.RP);
			Mthd.X_cor.adh = Mthd.A_Corrector(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_cor.lm = Mthd.A_Corrector(Mthd.X_prev.lm, LM.RP);
			Mthd.X_cor.ln = Mthd.A_Corrector(Mthd.X_prev.ln, LN.RP);

			Mthd.X_cor.mia = Mthd.A_Corrector(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_cor.mii = Mthd.A_Corrector(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_cor);// */
		}

		void NodeShift() final {
			Mthd.X[1] = Mthd.X[2]; Mthd.X[2] = Mthd.X[3];
			Mthd.X[3] = Mthd.X_prev; Mthd.X_prev = *Mthd.X_sol;
		}

	private:
		void ApplyAdams() {
			Mthd.X_pred.nec = Mthd.A_Predictor(Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.A_Predictor(Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.A_Predictor(Mthd.X_prev.ap_e, AE.RP);
			Mthd.X_pred.hel = Mthd.A_Predictor(Mthd.X_prev.hel, HEL.RP);

			Mthd.X_pred.cy = Mthd.A_Predictor(Mthd.X_prev.cy, CY.RP);
			Mthd.X_pred.adh = Mthd.A_Predictor(Mthd.X_prev.adh, ADH.RP);

			Mthd.X_pred.lm = Mthd.A_Predictor(Mthd.X_prev.lm, LM.RP);
			Mthd.X_pred.ln = Mthd.A_Predictor(Mthd.X_prev.ln, LN.RP);

			Mthd.X_pred.mia = Mthd.A_Predictor(Mthd.X_prev.mia, MIA.RP);
			Mthd.X_pred.mii = Mthd.A_Predictor(Mthd.X_prev.mii, MII.RP);

			ExpressSubValues(Mthd.X_pred);// */
		}

		void ApplyPrepMethod() final
		{
			Mthd.X_sub = Mthd.X_prev;

			Mthd.X_pred.nec = Mthd.Predictor(Mthd.X_sub.nec, Mthd.X_prev.nec, NEC.RP);
			Mthd.X_pred.ap_s = Mthd.Predictor(Mthd.X_sub.ap_s, Mthd.X_prev.ap_s, AS.RP);
			Mthd.X_pred.ap_e = Mthd.Predictor(Mthd.X_sub.ap_e, Mthd.X_prev.ap_e, AE.RP);
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