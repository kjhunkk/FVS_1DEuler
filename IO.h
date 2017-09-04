#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define gamma 1.4

typedef int Type;

class IO
{
public:
	IO();
	~IO();

public:
	// functions
	inline Type getInitial() const { return m_initial; };

	inline Type getScheme() const { return m_scheme; };

	inline int getNFX() const { return m_NFX; };

	inline double getTargetTime() const { return m_TargetTime; };

	inline double getCFL() const { return m_CFL; };

	// read input & scheme type / r.t. True/False
	bool initialize();

	// export solutions / p.m. index, x coordinate, density, velocity, pressure
	void output(std::string, std::vector<double>, std::vector<double>, std::vector<double>, std::vector<double>);
	
protected:
	// variables
	Type m_initial;
	Type m_scheme;
	int m_NFX;
	double m_TargetTime;
	double m_CFL;

protected:
	// functions
	// read input file / r.t. T : success, F : fail
	bool read_file();

	// type scheme / r.t. T : success, F : fail
	bool type_input();

	// print inputs
	void print();
};