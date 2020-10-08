#include <3rdModel_Pipe.h>

using namespace StraightTask;


int main()
{
	//calling for calculatings
	try { ODE_solver_output();}
	catch (const char* exception) {
		std::cout << "\n Calculation process was interrupted due to the problem occured:\n\t" << exception << std::endl;
	} //*/

	getchar();

	return EXIT_SUCCESS;
}