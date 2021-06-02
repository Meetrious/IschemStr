/* Here we include *Pipe.h with the task we want to solve. */
#include <4thModelPipe.h>

using namespace StraightTask;

// 'cause I want to
int main(){

	//calling for calculations
	try {

		Euler Task;
		//Task.AC.ConfigureSpline(QMSmaker);
		//Task.AC.OutputSpline("2AC");// */

		Task.SolveAndOutput(1500, 24.0, 4);
	}

	catch (MyException const & bzz) {
		std::cout << "\n Calculation process was interrupted due to the problem occured:\n\t"
			<< bzz.what() << "\n\t" << bzz.what_exactly();
		getchar();
		return EXIT_FAILURE;
	} 

	catch (std::exception const & bzz){
		std::cout << "\n Calculation process was interrupted due to the problem occured:\n\t"
			<< bzz.what();
		getchar();
		return EXIT_FAILURE;
	}

	catch (...){
		std::cout << "I've caught something I can't explain.\n "
			<< "Even std::exception can't cover that.";
		getchar();
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}