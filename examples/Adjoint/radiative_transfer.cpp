#include<iostream>
#include<aad/App>

#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[])
{
	if(argc != 2)
	{
		std::cout << "Usage: ./radiative_transfer <configuration filename>" << std::endl;

		return 1;
	}

	aad::RadiativeTransfer rtm(argv[1]);

	const auto& geo = rtm.geometry();
	std::vector<aad::core::RadiativeLayer> adj_source(1);
	adj_source[0].resize(geo.Ntheta, geo.M);
	adj_source[0].reflectance_m_top_cos[0](0, 0) = 1.0;

	rtm.run_adjoint(adj_source);

	const auto& result = rtm.result_adjoint();

	std::ofstream output("vertical_sensitivity.dat");
	output << "#altitude, dJ/dn, dJ/dq, dJ/d(ln q), dJ/dn, dJ/dq, dJ/d(ln q), ..." << std::endl;
	for(int i = 0; i < result.layers.size(); ++i)
	{
		output << std::setprecision(8) << std::scientific << result.altitude[i];

		for(int j = 0; j < result.layers[i].species.size(); ++j)
		{
			output << std::setprecision(8) << std::scientific << ", " << result.layers[i].species[j].number_density << ", " << result.layers[i].species[j].mixing_ratio << ", " << result.layers[i].species[j].log_mixing_ratio;
		}

		output << std::setprecision(8) << std::scientific << std::endl;
	}

	return 0;
}