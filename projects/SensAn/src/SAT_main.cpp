// SAT is for Sensitivity Analysis Task
#include <SensAn_base.h>
using namespace SensAnalysisTask;

int main() {

	SensAnTask<Euler> Task(
		Parameters(0.1f, 1500, 24.0, 1));

	Task.SolveForOutput();

	//Task.STM.SolveAndOutput(1500, 24.0, 4);

	return EXIT_SUCCESS;
}