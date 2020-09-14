#pragma once
#include<iostream>
#include<string>


class Gnuplot
{
public:
	Gnuplot();
	void wxt(int32_t term_ID, uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size);
	void gif(int32_t term_ID, uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size, uint16_t delay);
	void png(uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size, std::string dir, std::string name);
	void cd(const std::string & directory);
	~Gnuplot(); 
	void operator ()(const std::string & command);
	
//protected:
	FILE *gnuplotpipe; 
};

Gnuplot::Gnuplot()
{
	// with -persist option you will see the windows as your program ends
	//gnuplotpipe=_popen("gnuplot -persist","w");

	//without that option you will not see the window
	 // because I choose the terminal to output files so I don't want to see the window
	gnuplotpipe = _popen("gnuplot", "w"); // ф-ция открывающая консоль gnuplot

	if (!gnuplotpipe)
		std::cerr << ("Gnuplot not found !");

}

Gnuplot::~Gnuplot()
{
	fprintf(gnuplotpipe, "exit\n");
	_pclose(gnuplotpipe);
}

void Gnuplot::operator()(const std::string & command)
{
	fprintf(gnuplotpipe, "%s\n", command.c_str());
	fflush(gnuplotpipe);
	// flush is necessary, nothing gets plotted else
}

void Gnuplot::wxt
(int32_t term_ID, uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size)
{
	this->operator()
		("set terminal wxt " + std::to_string(term_ID) + " size " + std::to_string(Width) + ',' + std::to_string(Height) 
			+ " enhanced font \"Times Roman, " + std::to_string(xtick_inscript_size) + "\" ");
}

void Gnuplot::gif
(int32_t term_ID, uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size, uint16_t delay)
{
	this->operator()("set terminal gif " + std::to_string(term_ID) + " size " + std::to_string(Width) + ',' + std::to_string(Height)
		+ " enhanced font \"Times Roman, " + std::to_string(xtick_inscript_size) + "\" animate delay " + std::to_string(delay));
}

void Gnuplot::png(uint16_t Width, uint16_t Height, uint16_t xtick_inscript_size, std::string dir, std::string name)
{
	this->operator()
		("set terminal pngcairo size " + std::to_string(Width) + ", " + std::to_string(Height) +
		" enhanced font \"times roman, " + std::to_string(xtick_inscript_size) + "\" ");
	this->operator()
		("set output \"" + dir + name + ".png\"");
}

void Gnuplot::cd(const std::string & directory) { this->operator()("cd \'" + directory + "\'"); }

