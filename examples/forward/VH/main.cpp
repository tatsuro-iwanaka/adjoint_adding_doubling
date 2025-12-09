#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <aad/Core>
#include <aad/utils/Spline>

aad::core::RadiativeLayer initializeLayer(double tau, double omega, const std::vector<double>& theta_grid, const std::vector<double>& phase_function, const aad::core::Geometry& geo, double initial_tau = 1.0E-8)
{
	aad::core::RadiativeLayer layer;
	layer.optical_thickness = tau;
	layer.is_surface = false;
	layer.enable_atmospheric_emission = false;
	layer.n_doubling = 0;

	while(layer.optical_thickness > initial_tau)
	{
		layer.optical_thickness *= 0.5;
		layer.n_doubling++;
	}

	auto R_data = std::vector<std::vector<std::vector<double>>>(geo.Ntheta, std::vector<std::vector<double>>(geo.Ntheta, std::vector<double>(geo.Nphi, 0.0)));
	auto T_data = R_data;

	double base_coeff = omega * layer.optical_thickness / 4.0;

	for (int i = 0; i < geo.Ntheta; ++i)
	{
		for (int j = 0; j < geo.Ntheta; ++j)
		{
			double mu_i = geo.mu(i);
			double mu_j = geo.mu(j);
			double coef = base_coeff / (mu_i * mu_j);

			double sin_i = std::sqrt(std::max(0.0, 1.0 - mu_i * mu_i));
			double sin_j = std::sqrt(std::max(0.0, 1.0 - mu_j * mu_j));

			for (int k = 0; k < geo.Nphi; ++k)
			{
				double phi = geo.phi(k);
				double cos_phi = std::cos(phi);

				double cos_theta_refl = -mu_i * mu_j + sin_i * sin_j * cos_phi;
				double theta_refl = std::acos(std::clamp(cos_theta_refl, -1.0, 1.0));
				
				double cos_theta_tran = mu_i * mu_j + sin_i * sin_j * cos_phi;
				double theta_tran = std::acos(std::clamp(cos_theta_tran, -1.0, 1.0));

				double P_refl = 0.0; 
				double P_tran = 0.0;
				
				{
					// 反射
					auto it = std::lower_bound(theta_grid.begin(), theta_grid.end(), theta_refl);
					size_t idx = std::distance(theta_grid.begin(), it);
					if(idx == 0)
					{
						P_refl = phase_function[0];
					}
					else if(idx >= theta_grid.size())
					{
						P_refl = phase_function.back();
					}
					else
					{
						double r = (theta_refl - theta_grid[idx-1]) / (theta_grid[idx] - theta_grid[idx-1]);
						P_refl = phase_function[idx-1] * (1.0 - r) + phase_function[idx] * r;
					}

					// 透過
					it = std::lower_bound(theta_grid.begin(), theta_grid.end(), theta_tran);
					idx = std::distance(theta_grid.begin(), it);
					if(idx == 0)
					{
						P_tran = phase_function[0];
					}
					else if(idx >= theta_grid.size())
					{
						P_tran = phase_function.back();
					}
					else
					{
						double r = (theta_tran - theta_grid[idx-1]) / (theta_grid[idx] - theta_grid[idx-1]);
						P_tran = phase_function[idx-1] * (1.0 - r) + phase_function[idx] * r;
					}
				}

				R_data[i][j][k] = coef * P_refl;
				T_data[i][j][k] = coef * P_tran;
			}
		}
	}

	auto R_coeffs = aad::core::computeFourierSeriesCoefficients(R_data, geo);
	auto T_coeffs = aad::core::computeFourierSeriesCoefficients(T_data, geo);

	layer.reflectance_m_top_cos = R_coeffs[0];
	layer.reflectance_m_top_sin = R_coeffs[1];

	layer.reflectance_m_bottom_cos = R_coeffs[0];
	layer.reflectance_m_bottom_sin = R_coeffs[1];

	layer.transmittance_m_top_cos = T_coeffs[0];
	layer.transmittance_m_top_sin = T_coeffs[1];
	layer.transmittance_m_bottom_cos = T_coeffs[0];
	layer.transmittance_m_bottom_sin = T_coeffs[1];

	layer.source_up = Eigen::VectorXd::Zero(geo.Ntheta);
	layer.source_down = Eigen::VectorXd::Zero(geo.Ntheta);

	return layer;
}

int main(void)
{
	int n_layer = 1;
	int n_scattering_angle = 1201;
	

	std::ofstream output_12("vh_table12.dat");
	{
		std::vector<double> theta_grid(n_scattering_angle);
		std::vector<double> isotropical_scattering_phase_angle(n_scattering_angle, 1.0);
		for(int i = 0; i < n_scattering_angle; ++i)
		{
			theta_grid[i] = std::numbers::pi / double(n_scattering_angle - 1) * double(i);
		}
		
		std::vector<double> mu = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
		std::vector<double> mu0 = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};

		std::vector<double> a = {0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0};
		std::vector<double> b = {0.03125, 0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 1.0E3};

		int n_theta = 12;
		auto geo = aad::core::generateGeometryGaussRadau(n_theta, n_theta * 4 - 1, (n_theta * 4 - 2) / 2);
		std::vector<double> mu_grid(geo.mu.size());
		for(int i = 0; i < geo.mu.size(); ++i)
		{
			mu_grid[geo.mu.size() - 1 - i] = geo.mu(i);
		}

		for(int i = 0; i < b.size(); ++i)
		{
			output_12 << "===========================================================================================================" << std::endl;
			output_12 << "VH. Table 12. b (total optical thickness) = " << b[i] << ". Intensities out at Top. Interpolated by spline." << std::endl; 
			output_12 << "===========================================================================================================" << std::endl;
			output_12 << "a, mu=0.1, mu=0.3, mu=0.5, mu=0.7, mu=0.9, mu=1.0" << std::endl;

			std::vector<std::vector<std::vector<double>>> flux(a.size(), std::vector<std::vector<double>>(mu.size(), std::vector<double>(mu0.size(), 0.0)));
			std::vector<aad::core::RadiativeLayer> layers(n_layer);

			for(int j = 0; j < a.size(); ++j)
			{
				std::vector<std::vector<double>> flux_raw(geo.mu.size(), std::vector<double>(geo.mu.size(), 0.0));

				for(int k = 0; k < layers.size(); ++k)
				{
					layers[k] = initializeLayer(b[i], a[j], theta_grid, isotropical_scattering_phase_angle, geo);
				}

				auto result = aad::core::computeAtmosphere(layers, geo);

				for(int k = 0; k < geo.mu.size(); ++k)
				{
					for(int l = 0; l < geo.mu.size(); ++l)
					{
						flux_raw[k][l] = result.reflectance_m_bottom_cos[0](geo.mu.size() - 1 - k, geo.mu.size() - 1 - l);
					}
				}

				for(int k = 0; k < mu.size(); ++k)
				{
					for(int l = 0; l < mu0.size(); ++l)
					{
						flux[j][k][l] = aad::utils::interpolate2DSpline(mu[k], mu0[l], mu_grid, mu_grid, flux_raw);
					}
				}
			}

			for(int j = 0; j < mu0.size(); ++j)
			{
				output_12 << "mu0 = " << mu0[j] << std::endl;

				for(int k = 0; k < a.size(); ++k)
				{
					output_12 << a[k];
					// output_12 << mu[j];

					for(int l = 0; l < mu.size(); ++l)
					{
						output_12 << ", " << flux[k][l][j];	
					}

					output_12 << std::endl;
				}

				output_12 << std::endl;
			}

			output_12 << std::endl;
		}		
	}

	std::ofstream output_35("vh_table35.dat");
	{
		std::vector<double> theta_grid(n_scattering_angle);
		for(int i = 0; i < n_scattering_angle; ++i)
		{
			theta_grid[i] = std::numbers::pi / double(n_scattering_angle - 1) * double(i);
		}

		std::vector<double> mu = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};
		std::vector<double> mu0 = {0.1, 0.3, 0.5, 0.7, 0.9, 1.0};

		std::vector<double> a = {0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 0.99, 1.0};
		std::vector<double> b = {0.5, 1.0, 2.0, 4.0, 8.0};
		std::vector<double> g = {0.5, 0.75};

		int n_theta = 12;
		auto geo = aad::core::generateGeometryGaussRadau(n_theta, n_theta * 4 - 1, (n_theta * 4 - 2) / 2);
		std::vector<double> mu_grid(geo.mu.size());
		for(int i = 0; i < geo.mu.size(); ++i)
		{
			mu_grid[geo.mu.size() - 1 - i] = geo.mu(i);
		}

		for(int i = 0; i < b.size(); ++i)
		{
			for(int m = 0; m < g.size(); ++m)
			{
				output_35 << "=========================================================================================" << std::endl;
				output_35 << "VH. Table 35. g = " << g[m] << ", b = " << b[i] << ". Reflection. Interpolated by spline." << std::endl; 
				output_35 << "=========================================================================================" << std::endl;
				output_35 << "a, mu=0.1, mu=0.3, mu=0.5, mu=0.7, mu=0.9, mu=1.0" << std::endl;

				std::vector<aad::core::RadiativeLayer> layers(n_layer);
				std::vector<double> HG_scattering_phase_angle(n_scattering_angle, 1.0);

				for(int k = 0; k < n_scattering_angle; ++k)
				{
					HG_scattering_phase_angle[k] = (1.0 - g[m] * g[m]) / std::pow(1.0 + g[m] * g[m] - 2.0 * g[m] * cos(theta_grid[k]), 1.5);
				}

				std::vector<std::vector<std::vector<double>>> flux(a.size(), std::vector<std::vector<double>>(mu.size(), std::vector<double>(mu0.size(), 0.0)));

				for(int j = 0; j < a.size(); ++j)
				{
					std::vector<std::vector<double>> flux_raw(geo.mu.size(), std::vector<double>(geo.mu.size(), 0.0));

					for(int k = 0; k < layers.size(); ++k)
					{
						layers[k] = initializeLayer(b[i], a[j], theta_grid, HG_scattering_phase_angle, geo);
					}

					auto result = aad::core::computeAtmosphere(layers, geo);

					for(int k = 0; k < geo.mu.size(); ++k)
					{
						for(int l = 0; l < geo.mu.size(); ++l)
						{
							flux_raw[k][l] = result.reflectance_m_bottom_cos[0](geo.mu.size() - 1 - k, geo.mu.size() - 1 - l);
						}
					}

					for(int k = 0; k < mu.size(); ++k)
					{
						for(int l = 0; l < mu0.size(); ++l)
						{
							flux[j][k][l] = aad::utils::interpolate2DSpline(mu[k], mu0[l], mu_grid, mu_grid, flux_raw);
						}
					}
				}

				for(int j = 0; j < mu0.size(); ++j)
				{
					output_35 << "mu0 = " << mu0[j] << std::endl;

					for(int k = 0; k < a.size(); ++k)
					{
						output_35 << a[k];
						// output << mu[j];

						for(int l = 0; l < mu.size(); ++l)
						{
							output_35 << ", " << flux[k][l][j];	
						}

						output_35 << std::endl;
					}

					output_35 << std::endl;
				}

				output_35 << std::endl;
			}
		}		
	}

	return 0;
}