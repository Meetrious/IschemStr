
#include <4thModelOptTaskPipe.h>

using namespace ReverseTask;

int main() {
	
	BGA::Task<StraightTask::Euler> RT (
		BGA::Parameters(
			250, 8, 3000, 500,
			30,	10, 0.1, 0.01,
			1500, 24.0, 4)
	);
	
	//RT.ShowFforDefault();

	//RT.LookOverDecentOnes(10);

	//RT.SolveForStatistics(2);

	//RT.SolveForOutput();

	RT.OutputBestSolution();

	return EXIT_SUCCESS;
}