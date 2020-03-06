#include "libraries/Methods.h"

using namespace StraightTask;

void ODE_PredCor_solver(int32_t full_amount_of_days, size_t grid_spliting_amount)
{

	Parameters STpar(full_amount_of_days, grid_spliting_amount); 

	// ������ ���������� �������, ���������� ��� ����������� ���������� � ��������� ������� �� ������ ����
	// � ���������� ��������� �������
	Methods::PredCor TLayers(STpar);  
	
	Methods STask;
	// ��������� ������ ��� ������ ���������� ������� 
	STask.AllocateOutputStream();

	STask.CollectData("cy_a-spline.dat");
	TLayers.X_init.cy = STask.data[1][0];

	// ������� ������ ��������� ������� �� ������� ����
	for (size_t i = 1; i < Methods::N_eq + 5; i++)	STask.Output[i](TLayers.X_init, 0);


	bool TRIG = true;
	uint16_t days_remain = STpar.FAOD;
	Methods::H = STpar.H;

	/*double_t* CFC[12] = {
	&StraightTask::Neurons::k_N,
	&StraightTask::Neurons::k_A,
	&StraightTask::Neurons::p_R,

	&StraightTask::ToxDamage::p_ncy,
	&StraightTask::ToxDamage::C_Dcy,
	&StraightTask::ToxDamage::C_DLn,
	&StraightTask::ToxDamage::C_DLm,
	&StraightTask::ToxDamage::P_nn,
	&StraightTask::ToxDamage::p_Lm,
	&StraightTask::ToxDamage::p_Ln,
	&StraightTask::ToxDamage::D_0,
	&StraightTask::ToxDamage::p_D
	};*/
	//Methods::setCoefs(CFC, 12, working_directory + "mstat/medX.txt");

//	STask.CollectData("cy_a-spline.dat");


	while (TRIG)
	{
		std::cout << days_remain << " "; // ����� �� ����� "���-�� ����, ������� ���� ����������"
		TLayers.X_prev = TLayers.X_init; // ��� ��������� ������� �� �������������� ������

		// ���� �� ��������� 24 ����
		for (size_t Nj = 1; Nj <= STpar.N; Nj++)
		{
			// ����� �� ��������� ��� �� �������
			TLayers.X_pred.tim += STpar.H;

			//��������� ��������������� ������������
			TLayers.X_pred.AllocCurRets(Nj, STpar.N);

			// ������������ ���������� ������������
			TLayers.X_pred.CheckShiftInterpGap(STpar);
			TLayers.X_cor = TLayers.X_pred;

			// ���������� ������� �������� ��������� �� ������ �������� �����
			TLayers.X_cor.cy = STask.data[1][(Nj*16) + 24000 * (STpar.FAOD - days_remain)];
			
			launchPredCor(TLayers);

			// ����� �� ������� *.txt
			for (size_t i = 1; i < Methods::N_eq + 5; i++) STask.Output[i](TLayers.X_cor, Nj);

			// �������� ������������ (���� ��� ����, ���)
			rets::ShiftRets(TLayers.X_prev, Nj);

			// ��������� �������������� ��������� ��� ��� �������� �� ��������� ���
			TLayers.X_prev = TLayers.X_cor;

		} // ����� ����� ��������� �� ������� ���� for(Nj: 1->N)

		days_remain--; std::cout << " ; ";

		// ��������, �� ���� �� ������ ������?� 
		if (days_remain == 0)
		{
			std::cout << "\t\a Finita! \n\n All assigned days were rendered;\n";
			TRIG = false; // ������� �� ����� �� ����� while
			// else �� �����, ���� ���� !ID, �� ������� ������� while ��� �������
		}
		else TLayers.X_init = TLayers.X_cor;
	}
}


int main()
{
	// ������� ������������� �������
	// ��� ����������� �� ����, ����� ��� � ������ ��� ���
	Splines::ConfigureSplines();
	
	// ������������ ����������� ���� �������, ���� ����
	//Splines::OutputSplines();

	// �������� ������� ������� (2 ���, �� 1500 ����� � �����)
	ODE_PredCor_solver(2,1500);

	system("pause");
	return EXIT_SUCCESS;
}