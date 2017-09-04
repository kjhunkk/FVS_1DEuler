#include "FVS_Flux.h"

std::vector<double> FVS_Flux::m_rho;
std::vector<double> FVS_Flux::m_u;
std::vector<double> FVS_Flux::m_p;
std::array<std::vector<double>, 3> FVS_Flux::m_eigen;

FVS_Flux::FVS_Flux(bool sign, int flux, double gamma)
{
	m_sign = sign;
	m_scheme = flux;
	m_gamma = gamma;
}

FVS_Flux::~FVS_Flux()
{

}

std::array<std::vector<double>, 3> FVS_Flux::solve()
{
	std::array<std::vector<double>, 3> flux;
	switch (m_scheme)
	{
	case 1:
		flux = Steger_Warming();
		break;
	case 2:
		flux = van_Leer();
		break;
	case 3:
		flux = AUSM();
		break;
	default:
		for (int i = 0; i < 3; ++i)
		{
			flux[i].clear();
		}
	}

	return flux;
}

std::array<std::vector<double>, 3> FVS_Flux::Steger_Warming()
{
	std::array<std::vector<double>, 3> flux;
	std::array<double, 3> signed_eigen;
	double a;
	double H;

	for (int i = 0; i < m_eigen[0].size(); ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			signed_eigen[j] = 0.5*(m_eigen[j][i] + (2.0*m_sign - 1.0)*abs(m_eigen[j][i]));
		}
		a = sqrt(m_gamma*m_p[i] / m_rho[i]);
		H = 0.5*m_u[i] * m_u[i] + a*a / (m_gamma - 1.0);
		flux[0].push_back((signed_eigen[0] + 2.0*(m_gamma - 1.0)*signed_eigen[1] + signed_eigen[2])*m_rho[i] / (2.0*m_gamma));
		flux[1].push_back(((m_u[i] - a)*signed_eigen[0] + 2.0*(m_gamma - 1.0)*m_u[i] * signed_eigen[1] + (m_u[i] + a)*signed_eigen[2])*m_rho[i] / (2.0*m_gamma));
		flux[2].push_back(((H - m_u[i] * a)*signed_eigen[0] + (m_gamma - 1.0)*m_u[i] * m_u[i] * signed_eigen[1] + (H + m_u[i] * a)*signed_eigen[2])*m_rho[i] / (2.0*m_gamma));
	}

	return flux;
}

std::array<std::vector<double>, 3> FVS_Flux::van_Leer()
{
	std::array<std::vector<double>, 3> flux;
	double sign = (m_sign == 1 ? 1 : -1);
	double a;
	double M;
	double C;

	for (int i = 0; i < m_eigen[0].size(); ++i)
	{
		a = sqrt(m_gamma*m_p[i] / m_rho[i]);
		M = m_u[i] / a;
		if (M >= 1.0)
		{
			if (m_sign == true)
			{
				flux[0].push_back(m_rho[i] * a*M);
				flux[1].push_back(m_rho[i] * a*a*(M*M + (1.0 / m_gamma)));
				flux[2].push_back(m_rho[i] * pow(a, 3.0)*M*(0.5*M*M + 1.0 / (m_gamma - 1.0)));
			}
			else
			{
				flux[0].push_back(0.0);
				flux[1].push_back(0.0);
				flux[2].push_back(0.0);
			}
		}
		else if (M <= -1.0)
		{
			if (m_sign == false)
			{
				flux[0].push_back(m_rho[i] * a*M);
				flux[1].push_back(m_rho[i] * a*a*(M*M + (1.0 / m_gamma)));
				flux[2].push_back(m_rho[i] * pow(a, 3.0)*M*(0.5*M*M + 1.0 / (m_gamma - 1.0)));
			}
			else
			{
				flux[0].push_back(0.0);
				flux[1].push_back(0.0);
				flux[2].push_back(0.0);
			}
		}
		else
		{
			C = sign*0.25*m_rho[i] * a*pow((1.0 + sign*M), 2.0);
			flux[0].push_back(C);
			flux[1].push_back(C*2.0*a*(0.5*(m_gamma - 1.0)*M + sign) / m_gamma);
			flux[2].push_back(C*2.0*pow(a*(0.5*(m_gamma - 1.0)*M + sign), 2.0) / (m_gamma*m_gamma - 1.0));
		}
	}

	return flux;
}

std::array<std::vector<double>, 3> FVS_Flux::AUSM()
{
	std::array<std::vector<double>, 3> flux;

	return flux;
}