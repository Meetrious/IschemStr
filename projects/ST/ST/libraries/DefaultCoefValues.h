#define	DST(STRUCT, COEF, VAL) double_t StraightTask::STRUCT::COEF = VAL;

DST(Neurons, p_N, 0.9)
DST(Neurons, Rep, 1.0e-4)
DST(Neurons, k_N, 0.56366) // DST(Neurons, k_N, 2.15)
DST(Neurons, k_A, 3.0895) // DST(Neurons, k_A, 3.0)
DST(Neurons, p_R, 0.267460) // DST(Neurons, p_R, 0.25894) //DST(Neurons, p_R, 0.27)
double* StraightTask::Neurons::D_0 = &StraightTask::ToxDamage::D_0;

//  по модели 3
DST(Neurons,q_A, 1)
DST(Neurons, q_A_h, 0.006)
DST(Neurons, q_H, 8.13)
DST(Neurons, p_R3, 0.1)


DST(Neurons, q_N, 1.5)
DST(Neurons, q_A_n, 0.406)
DST(Neurons, c_H, 0.006)
DST(Neurons, c_A, 1.2e-2)
DST(Neurons, q_epN, 0.15)
DST(Neurons, q_epA, 0.001)



DST (Microglia, p_1, 1.0)
DST (Microglia, c_A, 0.06)
DST (Microglia, c_N, 2.1)
DST (Microglia, c_pro, 1.0)

DST(Microglia, T_M1, 60.0)
DST(Microglia, T_M2, 60.0)
DST(Microglia, c_Mi1, 0.3)
DST(Microglia, c_Mi2, 0.3)

DST(Microglia, c_dMi, 0.2)
DST(Microglia, c_Mi, 1.0)
DST(Microglia, K_Mi, 0.1)




DST(Cytokines, p_xcy0, 0.5)
DST(Cytokines, cy_max, 1.0)
DST(Cytokines, T_, 1.0)
DST(Cytokines, t_0, 8.0)


DST(Cytokines, p_Mach, 4.5)
DST(Cytokines, C_Ma, 0.4)
DST(Cytokines, p_Lmch, 4.5)
DST(Cytokines, C_Lm, 0.18)
DST(Cytokines, e_ch, 0.18)

DST(Cytokines, e_cy, 0.1)
DST(Cytokines, C_Ln, 0.18)


DST(Adhension, o_cy_1, 2.62e-2) //DST(Adhension, o_cy_1, 2.0)
DST(Adhension, o_cy_2, 3.14e-2) // DST(Adhension, o_cy_2, 1.0)
DST(Adhension, e_adh, 1.31e-2) //DST(Adhension, e_adh, 1.0)



DST(LeuMacrophags, c_Lm, 9.63e-2)
DST(LeuMacrophags, p_dLm, 3.6)
DST(LeuMacrophags, T_Lm, 1.0)

DST(LeuMacrophags, c_dLm, 0.08)
DST(LeuMacrophags, K_Lm, 43.9)
DST(LeuMacrophags, d_Lm, 0.139)


DST(LeuNeutrophils, c_Ln, 1.47)
DST(LeuNeutrophils, p_dLn, 3.2)
DST(LeuNeutrophils, T_Ln, 1.0)


DST(LeuNeutrophils, K_Ln, 1.0)
DST(LeuNeutrophils, c_dLn, 1.0)
DST(LeuNeutrophils, c_dLn1, 16.2)
DST(LeuNeutrophils, K_Ln1, 3.96)
DST(LeuNeutrophils, c_dLn2, 4.43)
DST(LeuNeutrophils, K_Ln2, 9.899)
DST(LeuNeutrophils, d_Ln, 84)


DST(ToxDamage, p_ncy, 17.83860)  // DST(ToxDamage, p_ncy, 16.652) // DST(ToxDamage, p_ncy, 0.5)
DST(ToxDamage, p_Ln, 27.078500)// DST(ToxDamage, p_Ln, 23.241) // DST(ToxDamage, p_Ln, 0.4)
DST(ToxDamage, C_DLn, 40.5670) // DST(ToxDamage, C_DLn, 44.165) // DST(ToxDamage, C_DLn, 0.6)
DST(ToxDamage, P_nn, 4.820900)// DST(ToxDamage, P_nn, 4.3495850) //DST(ToxDamage, P_nn, 5.4762)// DST(ToxDamage, P_nn, 0.05)

DST(ToxDamage, p_Lm, 34.70) // DST(ToxDamage, p_Lm, 33.394) // DST(ToxDamage, p_Lm, 0.4)
DST(ToxDamage, C_DLm, 8.6732) // DST(ToxDamage, C_DLm, 0.6)
DST(ToxDamage, C_D, 1.0)
DST(ToxDamage, C_Dcy, 9.7557)  // DST(ToxDamage, C_Dcy, 0.5) 

DST(ToxDamage, D_0, 0.001) //DST(ToxDamage, D_0, 1.643800) // DST(ToxDamage, D_0, 1.6789) // DST(ToxDamage, D_0, 0.25)
DST(ToxDamage, p_D, 0.28817) //  DST(ToxDamage, p_D, 0.25)

DST(ToxDamage, p_q1, 0.1)	DST(ToxDamage, c_q1, 13.206)
DST(ToxDamage, p_q2, 0.4)	DST(ToxDamage, c_q2, 1.32)
DST(ToxDamage, p_q3, 0.4)	DST(ToxDamage, c_q3, 0.1)
DST(ToxDamage, p_q4, 0.1)	DST(ToxDamage, c_q4, 0.1)

DST(Phagocytosis, e_Ma, 0.05) // DST(Phagocytosis, e_Ma, 0.005)
DST(Phagocytosis, e_Ln, 5.68e-3) //DST(Phagocytosis, e_Ln, 0.005)
DST(Phagocytosis, e_Lm, 0.05) // DST(Phagocytosis, e_Lm, 0.003)
DST(Phagocytosis, e_Mi, 1.0e-6) //DST(Phagocytosis, e_Mi, 0.001)


#undef DST