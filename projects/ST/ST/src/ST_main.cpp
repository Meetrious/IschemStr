#include <1stModelPipe.h>
using namespace StraightTask;

// 'cause I want to
auto main()->int{

	//calling for calculations
	try {
		PredCor Task;
		Task.SolveAndOutput(6000, 24.0, 4);
	}

	catch (MyException const & bzz) {
		std::cout << "\n Calculation process was interrupted due to the problem occured:\n\t"
			<< bzz.what() << "\n\t" << bzz.what_exactly();
		return EXIT_FAILURE;
	} 

	catch (std::exception const & bzz){
		std::cout << "\n Calculation process was interrupted due to the problem occured:\n\t"
			<< bzz.what();
		getchar();
		return EXIT_FAILURE;
	}

	catch (...){
		std::cout << "I've caught something I can't even explain.\n "
			<< "Even std::exception can't cover that.";
		getchar();
		return EXIT_FAILURE;
	}//*/

	return EXIT_SUCCESS;
}