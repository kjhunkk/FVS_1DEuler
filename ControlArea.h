#pragma once

#include <cmath>
#include <iostream>
#include <array>
#include <vector>

#define gamma 1.4

typedef int Type;

class ControlArea
{
public:
	// constructor / p.m. initial condition, NFX
	ControlArea(Type, int);
	~ControlArea();

public:
	// functions
	inline std::vector<double> getX() const { return m_x; };

	inline std::vector<double> getDensity() const { return m_rho; };

	inline std::vector<double> getVelocity() const { return m_u; };

	inline std::vector<double> getPressure() const { return m_p; };

	inline std::array<std::vector<double>, 3> getConserv() const { return m_U; };

	inline std::array<std::vector<double>, 3> getEigen() const { return m_Eigen; };

	inline double getDX() const { return m_dx; };

	inline int getSideGhost() const { return m_nSideGhost; }

	inline void setDensity(std::vector<double> rho) { m_rho = rho; };

	inline void setVelocity(std::vector<double> u) { m_u = u; };

	inline void setPressure(std::vector<double> p) { m_p = p; };

	inline void setConserv(std::array<std::vector<double>, 3> U) { m_U = U; };

	// generate control area (primitive variables) / p.m. number of side ghost cell / r.t. number of total cell
	int generate_primitive(int);

	// generate control area (conservative variables) / r.t. number of total cell
	void generate_conserv();

	// update conservative variables from primitive
	void update_conserv();

	// update primitive variables from conservative
	void update_primitive();

	// update eigen values
	void update_eigen();

	// calculate maximum absolute eigenvalue / r.t. maximum absolute value of eigenvalue
	double cal_max_eigen() const;

protected:
	// variables
	Type m_initial;
	int m_NFX;
	int m_nSideGhost;
	double m_dx;
	std::vector<double> m_x;
	std::vector<double> m_rho;
	std::vector<double> m_u;
	std::vector<double> m_p;
	std::vector<double> m_a;
	std::array<std::vector<double>, 3> m_U;
	std::array<std::vector<double>, 3> m_Eigen;

protected:
	// functions
	// initial conditions / r.t. primitive variables
	std::array<std::vector<double>, 3> initial_condition_1();

	std::array<std::vector<double>, 3> initial_condition_2();

	std::array<std::vector<double>, 3> initial_condition_3();

	std::array<std::vector<double>, 3> initial_condition_4();

	std::array<std::vector<double>, 3> initial_condition_5();
	
	// clear variables
	void clear_primitive();

	void clear_conserv();
};