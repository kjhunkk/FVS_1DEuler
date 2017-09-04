#include <memory>
#include <vector>
#include <array>

#include "ControlArea.h"
#include "IO.h"
#include "FVS_Flux.h"

#define gamma 1.4
#define forward 1
#define backward 0

int main()
{
	// Input Output object
	std::shared_ptr<IO> Euler_IO = std::make_shared<IO>();
	if(Euler_IO->initialize() == 0) return 0;
	int NFX = Euler_IO->getNFX();

	// Control area object
	std::shared_ptr<ControlArea> Euler_area = std::make_shared<ControlArea>(Euler_IO->getInitial(), NFX);
	int Total_NFX = Euler_area->generate_primitive(2);
	Euler_area->generate_conserv();

	Euler_IO->output("initial", Euler_area->getX(), Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure());

	// Flux object
	std::shared_ptr<FVS_Flux> Euler_flux_p = std::make_shared<FVS_Flux>(forward, Euler_IO->getScheme(), gamma);
	std::shared_ptr<FVS_Flux> Euler_flux_m = std::make_shared<FVS_Flux>(backward, Euler_IO->getScheme(), gamma);

	// Time
	double TargetTime = Euler_IO->getTargetTime();
	double dt = 0.0;
	double Time;

	// CFL
	double CFL = Euler_IO->getCFL();
	double dx = Euler_area->getDX();

	// Solution variables
	std::array<std::vector<double>, 3> U;
	std::array<std::vector<double>, 3> Flux_p;
	std::array<std::vector<double>, 3> Flux_m;

	// test
	//int iter = 0;

	Euler_area->update_primitive();
	Euler_area->update_eigen();

	for (Time = 0.0; (Time + dt) < TargetTime; Time += dt)
	{
		// test
		//iter++;

		U = Euler_area->getConserv();
		FVS_Flux::setPrimitive(Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure(), Euler_area->getEigen());
		Flux_p = Euler_flux_p->solve();
		Flux_m = Euler_flux_m->solve();
		dt = CFL*dx / Euler_area->cal_max_eigen();

		for (int i = 2; i < Total_NFX - 2; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				U[j][i] = U[j][i] - (dt / dx)*((Flux_p[j][i] + Flux_m[j][i + 1]) - (Flux_p[j][i - 1] + Flux_m[j][i]));
			}
		}
		//if ((iter % 10) == 0)
		//{
		//	Euler_IO->output(std::to_string(iter), Euler_area->getX(), Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure());
		//}
		Euler_area->setConserv(U);
		Euler_area->update_primitive();
		Euler_area->update_eigen();
		//Euler_IO->output("test", Euler_area->getX(), Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure());
	}
	if ((TargetTime - Time) > 10.0e-8)
	{
		dt = TargetTime - Time;
		U = Euler_area->getConserv();
		FVS_Flux::setPrimitive(Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure(), Euler_area->getEigen());
		Flux_p = Euler_flux_p->solve();
		Flux_m = Euler_flux_m->solve();

		for (int i = 2; i < Total_NFX - 2; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				U[j][i] = U[j][i] - (dt / dx)*((Flux_p[j][i] - Flux_m[j][i + 1]) - (Flux_p[j][i - 1] - Flux_m[j][i]));
			}
		}

		Euler_area->setConserv(U);
		Euler_area->update_primitive();
		Euler_area->update_eigen();
		Time += dt;
	}

	std::cout << "Finished in " << Time << "<sec>\n";
	std::cout << "Last time step = " << dt << "<sec>\n";

	Euler_IO->output("result", Euler_area->getX(), Euler_area->getDensity(), Euler_area->getVelocity(), Euler_area->getPressure());

	return 0;
}