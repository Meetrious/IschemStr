enum NumOfDerivative
{
	// tim, nec, /*a_st, a_en,*/ apop, hel, mia, mii, chem, madh, mac, neut, cyto, dpN, ep, psy, ind
	TIME,

	NEC,
	ACU,
	AP_S,
	AP_E,
	HEL,

	CY,
	CH,
	ADH,

	MIA,
	MII,
	LM,
	LN,


	D_F,

	D_INI,
	DP_N,
	DP_A,

	ePS,
	EPS,
	PSY,

	IND

};
/*
#define	TIME 0

#define	NEC	1
#define	ACU 2
#define	AP_S 3
#define	AP_E 4
#define	HEL 5

#define	CY 6
#define	CH 7
#define	ADH 8

#define	MIA 9
#define	MII 10
#define	LM 11
#define	LN 12
	

#define	D_F 13

#define	D_INI 14
#define	DP_N 15
#define	DP_A 16

#define	EPS 17
#define	PSY 18

#define	IND 19*/


enum CoefName
{
	/*1*/
	pN, rep, kN, kA, pR,

	/*2*/
	p1, cA, cN, cpro,
	TM1, TM2, cMi1, cMi2,
	cdMi, cMi, KMi,

	/*3*/
	cLm, pdLm, TLm,
	cdLm, KLm,
	cLn, pdLn, TLn,
	cdLn, KLn,

	/*4*/
	CMa, CLm, CLn, ecy,
	pMach, pLmch, ech,

	/*5*/
	pMadhcy1, pMadhcy2, eMadh,

	/*6*/
	pncy, CDcy,
	pLn, CDLn, // CDLm = CLm
	pLm, CDLm, CD,
	Pnn, D0, pD,

	/*7*/
	en, enLm, enLn, enmi,

	/*8*/
	pxcy0, cymax, T__, t0,

	length_coef_array
};
//
//#define pN 1 
//#define	REP 2
//#define kN 3
//#define kA 4
//#define pR 5
//
//
//#define	p1 6
//#define cA 7
//#define cN 8
//#define cpro 9
//#define TM1 10
//#define TM2 11
//#define cMi1 12
//#define cMi2 13
//
//#define cdMi 14
//#define cMi 15
//#define KMi 16
//
//
//#define	cLm 17
//#define pdLm 18
//#define TLm 19
//#define cLm1 20
//#define KLm 21
//
//
//#define cLn 22
//#define pdLn 23
//#define TLn 24
//
//#define cLn1 25
//#define KLn 26
//
//
//#define CMa 27
//#define CLm 28
//#define ecy 29
//#define pMach 30
//#define pLmch 31
//#define ech 32
//
//#define pMadhcy1 33
//#define pMadhcy2 34
//#define eMadh 35
//
//#define pncy 36
//#define pnLn 37
//#define CDLn 38 // CDLm = CLm
//#define pmLm 39
//#define CDLm 40
//
//#define Pnn 41
//#define D0 42
//#define pD 43
//
//#define en 44
//#define enLm 45
//#define enLn 46
//#define enmi 47
//
//#define pxcy0 48
//#define cymax 49
//#define _T 50
//#define t0 51
//
//#define LCA 52


std::string DerEnumToString(NumOfDerivative const & current)
{
	switch (current)
	{
	case TIME: return "tim";
	case NEC: return "nec";
	case ACU: return "acu_c";
	case AP_S: return "a_st";
	case AP_E: return "a_en";
	case HEL: return "hel";
	case CY: return "cyto";
	case CH: return "chem";

	
	case ADH: return "madh";

	case MIA: return "mia";
	case MII: return "mii";
	case LM: return "mac";
	case LN: return "neut";
	
	
	case D_F: return "dpF";
	case D_INI: return "d_in";

	case DP_N: return "dpN";
	case DP_A: return "dpA";
	case EPS: return "eps";
	case PSY: return "psy";
	default: return "DerivatNUM";
	}
}
