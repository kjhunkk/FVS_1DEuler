#include "ControlArea.h"

ControlArea::ControlArea(Type initial, int NFX)
{
	// initialize member variables
	m_initial = initial;
	m_NFX = NFX;
	m_nSideGhost = 0;
	m_dx = 1.0 / double(NFX);

	// clear vector container
	m_x.clear();
	m_rho.clear();
	m_u.clear();
	m_p.clear();
	for (int i = 0; i < 3; ++i)
	{
		m_U[i].clear();
		m_Eigen[i].clear();
	}
}

ControlArea::~ControlArea()
{

}

int ControlArea::generate_primitive(int nSideGhost)
{
	m_nSideGhost = nSideGhost;

	clear_primitive();
	m_x.resize(m_NFX + 2 * m_nSideGhost);
	m_rho.resize(m_NFX + 2 * m_nSideGhost);
	m_u.resize(m_NFX + 2 * m_nSideGhost);
	m_p.resize(m_NFX + 2 * m_nSideGhost);

	for (int i = 0; i < m_x.size(); ++i)
	{
		m_x[i] = m_dx*((double)(-m_nSideGhost + i) + 0.5);
	}

	std::array<std::vector<double>, 3> initial;
	
	switch (m_initial)
	{
	case 1:
		initial = initial_condition_1();
		break;
	case 2:
		initial = initial_condition_2();
		break;
	case 3:
		initial = initial_condition_3();
		break;
	case 4:
		initial = initial_condition_4();
		break;
	case 5:
		initial = initial_condition_5();
		break;
	default:
		std::cout << "initial condition error";
	}
	m_rho = initial[0];
	m_u = initial[1];
	m_p = initial[2];

	return m_x.size();
}

void ControlArea::generate_conserv()
{
	clear_conserv();

	for (int i = 0; i < 3; ++i)
	{
		m_U[i].resize(m_NFX + 2 * m_nSideGhost);
		m_Eigen[i].resize(m_NFX + 2 * m_nSideGhost);
		m_a.resize(m_NFX + 2 * m_nSideGhost);
	}
	for (int i = 0; i < m_x.size(); ++i)
	{
		m_U[0][i] = m_rho[i];
		m_U[1][i] = m_rho[i] * m_u[i];
		m_U[2][i] = 0.5*m_rho[i] * m_u[i] * m_u[i] + m_p[i] / (gamma - 1.0);
		m_a[i] = sqrt(gamma*m_p[i] / m_rho[i]);
		m_Eigen[0][i] = m_u[i] - m_a[i];
		m_Eigen[1][i] = m_u[i];
		m_Eigen[2][i] = m_u[i] + m_a[i];
	}
}

void ControlArea::update_conserv()
{
	for (int i = 0; i < m_x.size(); ++i)
	{
		m_U[0][i] = m_rho[i];
		m_U[1][i] = m_rho[i] * m_u[i];
		m_U[2][i] = (0.5*m_rho[i] * m_u[i] * m_u[i]) + (m_p[i] / (gamma - 1.0));
	}
}

void ControlArea::update_primitive()
{
	for (int i = 0; i < m_x.size(); ++i)
	{
		m_rho[i] = m_U[0][i];
		m_u[i] = m_U[1][i] / m_U[0][i];
		m_p[i] = (gamma - 1.0)*(m_U[2][i] - 0.5*m_U[1][i] * m_U[1][i] / m_U[0][i]);
	}
}

void ControlArea::update_eigen()
{
	for (int i = 0; i < m_x.size(); ++i)
	{
		m_a[i] = sqrt(gamma*m_p[i] / m_rho[i]);
		m_Eigen[0][i] = m_u[i] - m_a[i];
		m_Eigen[1][i] = m_u[i];
		m_Eigen[2][i] = m_u[i] + m_a[i];
	}
}

double ControlArea::cal_max_eigen() const
{
	std::vector<double> eigen;
	for (int i = 0; i < m_Eigen[0].size(); ++i)
	{
		eigen.push_back(abs(m_Eigen[1][i]) + m_a[i]);
	}

	return *std::max_element(eigen.begin(), eigen.end());
}

std::array<std::vector<double>, 3> ControlArea::initial_condition_1()
{
	std::array<std::vector<double>, 3> W;
	for (int i = 0; i < 3; ++i)
	{
		W[i].resize(m_x.size());
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		if (m_x[i] <= 0.3)
		{
			W[0][i] = 1.0;
			W[1][i] = 0.75;
			W[2][i] = 1.0;
		}
		else
		{
			W[0][i] = 0.125;
			W[1][i] = 0.0;
			W[2][i] = 0.1;
		}
	}

	return W;
}

std::array<std::vector<double>, 3> ControlArea::initial_condition_2()
{
	std::array<std::vector<double>, 3> W;
	for (int i = 0; i < 3; ++i)
	{
		W[i].resize(m_x.size());
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		if (m_x[i] <= 0.5)
		{
			W[0][i] = 1.0;
			W[1][i] = -2.0;
			W[2][i] = 3.0;
		}
		else
		{
			W[0][i] = 1.0;
			W[1][i] = 2.0;
			W[2][i] = 3.0;
		}
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		W[2][i] = (W[2][i] - 0.5*W[0][i] * W[1][i] * W[1][i]) / (gamma - 1);
	}

	return W;
}

std::array<std::vector<double>, 3> ControlArea::initial_condition_3()
{
	std::array<std::vector<double>, 3> W;
	for (int i = 0; i < 3; ++i)
	{
		W[i].resize(m_x.size());
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		if (m_x[i] <= 0.5)
		{
			W[0][i] = 1.0;
			W[1][i] = 0.0;
			W[2][i] = 1000.0;
		}
		else
		{
			W[0][i] = 1.0;
			W[1][i] = 0.0;
			W[2][i] = 0.01;
		}
	}

	return W;
}

std::array<std::vector<double>, 3> ControlArea::initial_condition_4()
{
	std::array<std::vector<double>, 3> W;
	for (int i = 0; i < 3; ++i)
	{
		W[i].resize(m_x.size());
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		if (m_x[i] <= 0.4)
		{
			W[0][i] = 5.99924;
			W[1][i] = 19.5975;
			W[2][i] = 460.894;
		}
		else
		{
			W[0][i] = 5.99242;
			W[1][i] = -6.19633;
			W[2][i] = 46.0950;
		}
	}

	return W;
}

std::array<std::vector<double>, 3> ControlArea::initial_condition_5()
{
	std::array<std::vector<double>, 3> W;
	for (int i = 0; i < 3; ++i)
	{
		W[i].resize(m_x.size());
	}

	for (int i = 0; i < m_x.size(); ++i)
	{
		if (m_x[i] <= 0.8)
		{
			W[0][i] = 1.0;
			W[1][i] = -19.59745;
			W[2][i] = 1000.0;
		}
		else
		{
			W[0][i] = 1.0;
			W[1][i] = -19.59745;
			W[2][i] = 0.01;
		}
	}

	return W;
}

void ControlArea::clear_primitive()
{
	m_x.clear();
	m_rho.clear();
	m_u.clear();
	m_p.clear();
}

void ControlArea::clear_conserv()
{
	for (int i = 0; i < 3; ++i)
	{
		m_U[i].clear();
		m_Eigen[i].clear();
	}
	m_a.clear();
}
