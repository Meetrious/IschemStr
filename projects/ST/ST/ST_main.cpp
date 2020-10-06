<<<<<<< HEAD
#include <3rdModel_Pipe.h>
=======
//âûáèðàåì çàäà÷ó, êîòîðóþ áóäåì ðåøàòü
#include <Test1D_Pipe.h>
>>>>>>> 83b0c79880d264b11dc44cc5f219e5a8f06899ba

using namespace StraightTask;


int main()
{
<<<<<<< HEAD
	//calling for calculatings
	try {ODE_solver();}
	catch (const char* exception) {
		std::cout << "\n Calculation process were interrupted due to the problem occured:\n\t" << exception << std::endl;
	} //*/

	getchar();
=======
	//SolutionOutput();
	
	//calling for calculations
	ODE_solver();
>>>>>>> 83b0c79880d264b11dc44cc5f219e5a8f06899ba

	return EXIT_SUCCESS;
}
