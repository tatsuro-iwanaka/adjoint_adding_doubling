#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

#include <aad/Core>

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

aad::core::RadiativeLayer initializeIsotropicLayer(double target_tau, double omega, const aad::core::Geometry& geo, double initial_tau = 1.0E-8)
{
	aad::core::RadiativeLayer layer;
	layer.optical_thickness = target_tau;
	layer.enable_atmospheric_emission = false;
	layer.is_surface = false;
	layer.n_doubling = 0;

	while(layer.optical_thickness > initial_tau)
	{
		layer.optical_thickness *= 0.5;
		layer.n_doubling++;
	}

	double tau = layer.optical_thickness;

	layer.source_up = Eigen::VectorXd::Zero(geo.Ntheta);
	layer.source_down = Eigen::VectorXd::Zero(geo.Ntheta);
	layer.reflectance_m_top_cos.resize(geo.M + 1);
	layer.reflectance_m_bottom_cos.resize(geo.M + 1);
	layer.transmittance_m_top_cos.resize(geo.M + 1);
	layer.transmittance_m_bottom_cos.resize(geo.M + 1);
	
	layer.reflectance_m_top_sin.resize(geo.M + 1);
	layer.reflectance_m_bottom_sin.resize(geo.M + 1);
	layer.transmittance_m_top_sin.resize(geo.M + 1);
	layer.transmittance_m_bottom_sin.resize(geo.M + 1);

	for(int m=0; m <= geo.M; ++m)
	{
		layer.reflectance_m_top_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.reflectance_m_bottom_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_top_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_bottom_cos[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);

		layer.reflectance_m_top_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.reflectance_m_bottom_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_top_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
		layer.transmittance_m_bottom_sin[m] = Eigen::MatrixXd::Zero(geo.Ntheta, geo.Ntheta);
	}

	for (int i = 0; i < geo.Ntheta; ++i)
	{
		for (int j = 0; j < geo.Ntheta; ++j)
		{
			double mu_i = geo.mu(i);
			double mu_j = geo.mu(j);
			
			// Reflection
			double r_val = (omega * tau) / (4.0 * (mu_i * mu_j));
			layer.reflectance_m_top_cos[0](i, j) = r_val;
			layer.reflectance_m_bottom_cos[0](i, j) = r_val;

			// Transmission
			double t_val = (omega * tau) / (4.0 * (mu_i * mu_j));
			layer.transmittance_m_top_cos[0](i, j) = t_val;
			layer.transmittance_m_bottom_cos[0](i, j) = t_val;
		}
	}

	return layer;
}

int main(void)
{
	int n_layer = 21;
	int n_theta = 10;
	std::vector<double> layer_taus(n_layer, 0.2);
	std::vector<double> layer_omegas(n_layer);

	int n_scattering_angle = 1201;
	std::vector<double> theta_grid(n_scattering_angle);
	for(int i = 0; i < n_scattering_angle; ++i)
	{
		theta_grid[i] = std::numbers::pi / double(n_scattering_angle - 1) * double(i);
	}

	for(int i = 0; i < n_layer; ++i)
	{
		if (i < n_layer / 2)
		{
			layer_omegas[i] = 0.9;
		}
		else
		{
			layer_omegas[i] = 0.1;
		}
	}

	auto geo = aad::core::generateGeometryGaussRadau(n_theta, n_theta * 4 - 1, (n_theta * 4 - 2) / 2);

	std::vector<aad::core::RadiativeLayer> layers(n_layer);
		
	for(int i = 0; i < layers.size(); ++i)
	{
		layers[i] = initializeIsotropicLayer(layer_taus[i], layer_omegas[i], geo);
	}

	auto result = aad::core::computeAtmosphere(layers, geo);

	aad::core::RadiativeLayer adj_result = result;
	
	adj_result.optical_thickness = 0.0;
	adj_result.source_up.setZero();
	adj_result.source_down.setZero();
	for(int m=0; m<=geo.M; ++m)
	{
		adj_result.reflectance_m_top_cos[m].setZero();
		adj_result.reflectance_m_bottom_cos[m].setZero();
		adj_result.transmittance_m_top_cos[m].setZero();
		adj_result.transmittance_m_bottom_cos[m].setZero();
		adj_result.reflectance_m_top_sin[m].setZero();
		adj_result.reflectance_m_bottom_sin[m].setZero();
		adj_result.transmittance_m_top_sin[m].setZero();
		adj_result.transmittance_m_bottom_sin[m].setZero();
	}

	adj_result.reflectance_m_top_cos[0](0, 0) = 1.0;
	auto adj_layers = aad::core::computeAtmosphere_adjoint(layers, geo, adj_result);

	/*
	{
		std::cout << "Sensitivity on optical thickness" << std::endl;
		std::cout << "Central Finite Difference" << std::endl;
		for(int i = 0; i < n_layer; ++i)
		{
			std::vector<aad::core::RadiativeLayer> layers_p(n_layer);
			std::vector<aad::core::RadiativeLayer> layers_m(n_layer);

			double dtau = layers[i].optical_thickness * 1.0E-4;
			double domega = layer_omegas[i] * 1.0E-4;

			for(int j = 0; j < layers.size(); ++j)
			{
				if(i == j)
				{
					layers_p[j] = initializeIsotropicLayer(layer_taus[j] + dtau, layer_omegas[j], geo);
					layers_m[j] = initializeIsotropicLayer(layer_taus[j] - dtau, layer_omegas[j], geo);
				}
				else
				{
					layers_p[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
					layers_m[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
				}
			}

			auto result_p = aad::core::computeAtmosphere(layers_p, geo);
			auto result_m = aad::core::computeAtmosphere(layers_m, geo);

			double dj = (result_p.reflectance_m_top_cos[0](0, 0) - result_m.reflectance_m_top_cos[0](0, 0)) / (2.0 * dtau);

			std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << dj << std::endl;
		}

		std::cout << "Adjoint" << std::endl;
		for(int i = 0; i < n_layer; ++i)
		{
			double grad_tau = 0.0;
			double grad_omega = 0.0;
			std::vector<double> grad_P(n_scattering_angle, 0.0);
			double scaling_factor = 1.0 / std::pow(2.0, layers[i].n_doubling);
			aad::core::computeInitializationSensitivities(adj_layers[i], layers[i], layer_omegas[i], geo, theta_grid, grad_tau, grad_omega, grad_P);
			grad_tau *= scaling_factor;
			std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << grad_tau << std::endl;
		}
	}
	
	{
		std::cout << "Sensitivity on single scattering albedo" << std::endl;
		std::cout << "Central Finite Difference" << std::endl;
		for(int i = 0; i < n_layer; ++i)
		{
			std::vector<aad::core::RadiativeLayer> layers_p(n_layer);
			std::vector<aad::core::RadiativeLayer> layers_m(n_layer);

			double dtau = layers[i].optical_thickness * 1.0E-4;
			double domega = layer_omegas[i] * 1.0E-4;

			for(int j = 0; j < layers.size(); ++j)
			{
				if(i == j)
				{
					// layers_p[j] = initializeIsotropicLayer(layer_taus[j] + dtau, layer_omegas[j], geo);
					// layers_m[j] = initializeIsotropicLayer(layer_taus[j] - dtau, layer_omegas[j], geo);
					layers_p[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j] + domega, geo);
					layers_m[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j] - domega, geo);
				}
				else
				{
					layers_p[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
					layers_m[j] = initializeIsotropicLayer(layer_taus[j], layer_omegas[j], geo);
				}
			}

			auto result_p = aad::core::computeAtmosphere(layers_p, geo);
			auto result_m = aad::core::computeAtmosphere(layers_m, geo);

			// double dj = (result_p.reflectance_m_top_cos[0](0, 0) - result_m.reflectance_m_top_cos[0](0, 0)) / (2.0 * dtau);
			double dj = (result_p.reflectance_m_top_cos[0](0, 0) - result_m.reflectance_m_top_cos[0](0, 0)) / (2.0 * domega);

			std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << dj << std::endl;
		}

		std::cout << "Adjoint" << std::endl;
		for(int i = 0; i < n_layer; ++i)
		{
			double grad_tau = 0.0;
			double grad_omega = 0.0;
			std::vector<double> grad_P(n_scattering_angle, 0.0);
			double scaling_factor = 1.0 / std::pow(2.0, layers[i].n_doubling);
			aad::core::computeInitializationSensitivities(adj_layers[i], layers[i], layer_omegas[i], geo, theta_grid, grad_tau, grad_omega, grad_P);
			grad_tau *= scaling_factor;
			// std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << grad_tau << std::endl;
			std::cout << "Layer " << i << " (Top=" << (n_layer-1) << "): " << std::scientific << std::setprecision(8) << grad_omega << std::endl;
		}
	}

	{
		for(int i = 0; i < n_layer; i ++)
		{
			int target_layer = i; 

			std::cout << "\n--- Phase Function Sensitivity Verification (Layer " << target_layer << ") ---" << std::endl;

			std::vector<double> base_P(n_scattering_angle, 1.0);

			double grad_tau = 0.0;
			double grad_omega = 0.0;
			std::vector<double> grad_P_adj(n_scattering_angle, 0.0);
			
			aad::core::computeInitializationSensitivities(adj_layers[target_layer], layers[target_layer], layer_omegas[target_layer], geo, theta_grid, grad_tau, grad_omega, grad_P_adj);

			std::cout << "Angle_Index, Theta(deg), Adjoint, FiniteDifference, Error" << std::endl;
			
			std::vector<int> check_indices = {0, n_scattering_angle/4, n_scattering_angle/2, n_scattering_angle*3/4, n_scattering_angle-1};

			for (int k : check_indices)
			{
				double epsilon = 1.0e-4;
				double dp = epsilon * base_P[k];
				
				std::vector<double> P_plus = base_P;
				P_plus[k] += dp;
				
				std::vector<double> P_minus = base_P;
				P_minus[k] -= dp;

				std::vector<aad::core::RadiativeLayer> layers_p(n_layer);
				std::vector<aad::core::RadiativeLayer> layers_m(n_layer);

				for(int j = 0; j < layers.size(); ++j)
				{
					if(j == target_layer)
					{
						layers_p[j] = initializeLayer(layer_taus[j], layer_omegas[j], theta_grid, P_plus, geo);
						layers_m[j] = initializeLayer(layer_taus[j], layer_omegas[j], theta_grid, P_minus, geo);
					}
					else
					{
						layers_p[j] = initializeLayer(layer_taus[j], layer_omegas[j], theta_grid, base_P, geo);
						layers_m[j] = initializeLayer(layer_taus[j], layer_omegas[j], theta_grid, base_P, geo);
					}
				}

				auto res_p = aad::core::computeAtmosphere(layers_p, geo);
				auto res_m = aad::core::computeAtmosphere(layers_m, geo);

				double J_p = res_p.reflectance_m_top_cos[0](0, 0);
				double J_m = res_m.reflectance_m_top_cos[0](0, 0);

				double grad_fd = (J_p - J_m) / (2.0 * dp);

				std::cout << k << ", " << theta_grid[k] * 180.0 / std::numbers::pi << ", "<< grad_P_adj[k] << ", " << grad_fd << ", " << (grad_P_adj[k] - grad_fd) / grad_P_adj[k] << std::endl;
			}
		}
	}
	*/

	std::ofstream output("layer_sensitivity.dat");

	{
		for(int i = 0; i < n_layer; i ++)
		{
			int target_layer = i;

			double grad_tau = 0.0;
			double grad_omega = 0.0;
			std::vector<double> grad_P_adj(n_scattering_angle, 0.0);

			aad::core::computeInitializationSensitivities(adj_layers[target_layer], layers[target_layer], layer_omegas[target_layer], geo, theta_grid, grad_tau, grad_omega, grad_P_adj);
			
			std::vector<int> check_indices = {0, n_scattering_angle/6*1, n_scattering_angle/6*2, n_scattering_angle/6*3, n_scattering_angle/6*4, n_scattering_angle/6*5, n_scattering_angle - 1};

			double scaling_factor = 1.0 / std::pow(2.0, layers[target_layer].n_doubling);

			output << i << ", " << grad_tau * scaling_factor << ", " << grad_omega;

			for (int k : check_indices)
			{
				output << ", " << grad_P_adj[k];
			}

			output << std::endl;
		}
	}

	return 0;
}