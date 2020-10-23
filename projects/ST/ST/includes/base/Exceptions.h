#pragma once
#include <exception>
#include <stdio.h>


class MyException : public std::exception {
public:
	MyException() 
	{
		oflog = fopen("Error_log.txt", "w");
		if (oflog == NULL) { printf("Cannot create error_log. Watch the console output! \n"); }
	}
	~MyException() { fclose(oflog); }

	FILE* oflog;

	const virtual char* what()const noexcept final {
		return "Some local shenanigan is happening:";
	}
	const virtual char* what_exactly() const noexcept;

};

class external_file_allocation_error : public MyException {
public:
	const virtual char* what_exactly() const noexcept final {
		return "Cannot allocate given IO-stream to external file.";
	}
};

class inconsistent_initial_data : public MyException
{
public:
	const virtual char* what_exactly() const noexcept final {
		return "Initial data for the ODE system is not consistent";
	}
	const virtual char* problem() const noexcept;
};
