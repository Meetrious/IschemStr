
#include <3rdModelOptTaskPipe.h>


int main()
{
	//using namespace ReverseTask;
		
	ReverseTask::BGA::ISolver<StraightTask::Euler> Task (100, 5000, 1000, 20);
	
	//Task.SolveForOutput();


	Task.OutputBestSolution();

	return EXIT_SUCCESS;
}