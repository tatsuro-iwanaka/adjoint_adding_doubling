#include<filesystem>
#include<iostream>
#include<sstream>
#include<chrono>
#include<vector>
#include<string>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<algorithm>
#include<netcdf>

int main(int argc, char* argv[])
{
	std::string filename = std::string(argv[1]);
	std::string filename_out;
	for(int i = filename.size() - 1; i >= 0; --i)
	{
		if(filename[i] == '.')
		{
			filename_out = filename.substr(0, i) + ".dat";
		}
	}
	
	if(argc == 3)
	{
		filename_out = std::string(argv[2]);
	}
	
	netCDF::NcFile input(filename, netCDF::NcFile::read);

	netCDF::NcDim dim_theta_e = input.getDim("theta_e");
	int num_theta_e = dim_theta_e.getSize();
	netCDF::NcDim dim_theta_i = input.getDim("theta_i");
	int num_theta_i = dim_theta_i.getSize();
	netCDF::NcDim dim_phi = input.getDim("delta_phi");
	int num_phi = dim_phi.getSize();
	netCDF::NcDim dim_layer = input.getDim("layer");
	int num_layer = dim_layer.getSize();
	netCDF::NcDim dim_wavelength = input.getDim("wavelength");
	int num_wavelength = dim_wavelength.getSize();

	std::cout << filename << " -> " << filename_out << std::endl;
	std::cout << num_theta_e << ", " << num_layer << ", " << num_wavelength << std::endl;

	std::vector<double> zenith_angle(num_theta_e);
	netCDF::NcVar nc_za = input.getVar("theta_e");
	nc_za.getVar(zenith_angle.data());

	std::vector<double> reflectance_flat(num_wavelength * num_theta_e * num_theta_i * num_phi);
	netCDF::NcVar nc_r = input.getVar("reflectance");
	nc_r.getVar(reflectance_flat.data());

	std::vector<double> thermal_emission_flat(num_wavelength * num_theta_e);
	netCDF::NcVar nc_b = input.getVar("thermal_emission");
	nc_b.getVar(thermal_emission_flat.data());

	std::vector<std::vector<std::vector<std::vector<double>>>> reflectance(num_wavelength, std::vector<std::vector<std::vector<double>>>(num_theta_e, std::vector<std::vector<double>>(num_theta_i, std::vector<double>(num_phi))));
	std::vector<std::vector<double>> thermal_emission(num_wavelength, std::vector<double>(num_theta_e));

	// export emission angle dependencies of reflecntace and thermal emission
	for(int i = 0; i < num_wavelength; ++i)
	{
		for(int j = 0; j < num_theta_e; ++j)
		{
			thermal_emission[i][j] = thermal_emission_flat[i * num_theta_e + j];

			for(int k = 0; k < num_theta_i; ++k)
			{
				for(int l = 0; l < num_phi; ++l)
				{
					reflectance[i][j][k][l] = reflectance_flat[i * num_theta_e * num_theta_i * num_phi + j * num_theta_i * num_phi + k * num_phi + l];
				}
			}
		}
	}

	std::ofstream output(filename_out);

	for(int e = 0; e < num_theta_e; ++e)
	{
		output << std::scientific << std::setprecision(12) << zenith_angle[e] * 180.0 / M_PI << ", " << reflectance[0][e][0][0] << ", " << thermal_emission[0][e] << std::endl;
	}

	return 0;
}
