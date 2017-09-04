#include "IO.h"

IO::IO()
{
	m_initial = 0;
	m_scheme = 0;
	m_NFX = 0;
	m_TargetTime = 0.0;
	m_CFL = 0;
}

IO::~IO()
{

}

bool IO::initialize()
{
	bool input_file = read_file();
	bool scheme_type = type_input();

	print();

	return input_file*scheme_type;
}

void IO::output(std::string index, std::vector<double> x, std::vector<double> rho, std::vector<double> u, std::vector<double> p)
{
	std::ofstream file;
	std::string FileName = "./output/";
	FileName += "1DEuler_";
	std::string Scheme = "";
	switch (m_initial)
	{
	case 1:
		FileName += "ShockTube_";
		break;
	case 2:
		FileName += "Expansion_";
		break;
	case 3:
		FileName += "Blast_";
		break;
	case 4:
		FileName += "DoubleShock_";
		break;
	case 5:
		FileName += "Contact_";
		break;
	default:
		FileName += "Error_";
	}
	switch (m_scheme)
	{
	case 1:
		Scheme = "Steger-Warming";
		break;
	case 2:
		Scheme = "van Leer";
		break;
	case 3:
		Scheme = "AUSM";
		break;
	default:
		Scheme = "Scheme error";
	}
	FileName += "CFL";
	FileName += std::to_string(m_CFL);
	FileName += "_";
	FileName += index;
	FileName += ".plt";
	file.open(FileName, std::ios::trunc);
	if (file.is_open())
	{
		std::cout << "output file open\n";
		file << "Variables = \"X\", \"density\", \"velocity\", \"pressure\", \"total energy\"\n";
		file << "Zone t = \"" << Scheme << "\", i=" << x.size() << ", f=point\n";
		for (int i = 0; i < x.size(); ++i)
		{
			file << x[i] << "\t\t" << rho[i] << "\t\t" << u[i] << "\t\t" << p[i] << "\t\t" << (0.5*rho[i]*u[i]*u[i] + p[i]/(gamma - 1)) << "\n";
		}
	}
	else std::cout << "output file error\n";

	file.close();
}

bool IO::read_file()
{
	bool ind = 0;
	std::ifstream file;
	std::string buff;
	std::string FileName = "./input.inp";
	file.open(FileName);

	if (file.is_open())
	{
		std::cout << "----------Input File Open----------\n";
		file >> buff >> buff >> buff >> buff >> m_initial >> m_NFX >> m_TargetTime >> m_CFL;
		ind = 1;
	}
	file.close();
	return ind;
}

bool IO::type_input()
{
	std::cout << "Scheme type = ";
	std::cin >> m_scheme;
	if ((m_scheme <= 3) && (m_scheme > 0)) return 1;
	else return 0;
}

void IO::print()
{
	std::cout << "----------Conditions----------\n";
	std::cout << "$ Initial condition : " << m_initial << "\n";
	std::cout << "$ Scheme type       : " << m_scheme << "\n";
	std::cout << "$ Number of cells   : " << m_NFX << "\n";
	std::cout << "$ Target Time       : " << m_TargetTime << "\n";
	std::cout << "$ CFL number        : " << m_CFL << "\n";
	std::cout << "------------------------------\n";
}