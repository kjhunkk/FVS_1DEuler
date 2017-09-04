#pragma once

#include <vector>
#include <array>

class FVS_Flux
{
public:
	// constructor / p.m. sign(+ : 1, - : 0), flux type, gamma
	FVS_Flux(bool, int, double);

	// destructor
	~FVS_Flux();

public:
	// functions
	inline static void setPrimitive(std::vector<double> rho, std::vector<double> u, std::vector<double> p, std::array<std::vector<double>, 3> eigen) { m_rho = rho; m_u = u; m_p = p; m_eigen = eigen; };

	// calculate flux / r.t. flux
	std::array<std::vector<double>, 3> solve();

protected:
	// variables
	bool m_sign;
	int m_scheme;
	double m_gamma;

	// static variables
	static std::vector<double> m_rho;
	static std::vector<double> m_u;
	static std::vector<double> m_p;
	static std::array<std::vector<double>, 3> m_eigen;

protected:
	// functions
	// Steger-Warming Splitting
	std::array<std::vector<double>, 3> Steger_Warming();

	// van Leer Splitting
	std::array<std::vector<double>, 3> van_Leer();

	// AUSM
	std::array<std::vector<double>, 3> AUSM();
};