#ifndef ISING_H
#define ISING_H

#include <vector>
#include <algorithm>
#include <random>

#include <iostream>
#include <iomanip>

#include "crystal.h"
#include "site.h"
#include "lattice.h"
#include "GSLpp/matrix.h"

template<size_t dim>
class Ising_t{
	private:
		std::array<size_t, dim> size_m;
		Crystal_t<dim> cr_m;
		std::vector<int8_t> field_m;
		std::vector<std::vector<Neighbours<dim>>> nn_shells_m;
		std::vector<double> J_m;
		double H_m;
		double beta_m;

		size_t calc_length()
		{
			size_t length = 1;
			for(auto val : size_m){
				length *= val;
			}
			return length;
		}

		void setup_field()
		{
			size_t length = calc_length();
			field_m = std::vector<int8_t>(length);
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<int8_t> dist(0,1);
			for(auto& val : field_m){
				val = static_cast<int8_t>(2*dist(gen) - 1);
			}
		}

		void setup_crystal()
		{
			cr_m.set_size(size_m);
		}

		void setup_crystal_sites()
		{
			size_t length = calc_length();
			std::vector<GSL::Vector> positions(length);
			size_t offset = 1, tmp_index;
			std::array<size_t, dim> coord;
			for(size_t idx = 0; idx < length; idx++){
				tmp_index = idx;
				for(size_t i = dim; i > 0; i--)
				{
					offset = 1;
					for(size_t j = 0; j < i - 1; j++ ){
						offset *= size_m[j];
					}
					coord[i-1] = tmp_index / offset;
					tmp_index -= coord[i-1]*offset;
				}
				positions[idx] = GSL::Vector(dim);
				for(size_t i = 0; i < coord.size(); i++){
					positions[idx][i] = static_cast<double>(coord[i]);
				}
			}
			cr_m.add_sites(positions);
		}

		void setup_crystal_lattice_vectors(bool periodic)
		{
			if(periodic){
				cr_m.set_Rn(1);
			}else{
				cr_m.set_Rn(0);
			}
		}

		void setup_nearest_neighbour_shells(const size_t n_steps = 1)
		{
			std::vector<Neighbours<dim>> nn = cr_m.calc_nearest_neighbours(n_steps);
			nn_shells_m = cr_m.determine_nn_shells(nn);
		}

		double site_energy(size_t index)
		{
			double energy = 0;
			int other_spins = 0;
			// Loop over all interaction constants provided
			for(size_t i = 0; i < J_m.size(); i++){
				for(const auto neighbour : nn_shells_m[index][i]){
					other_spins += field_m[neighbour.index()];
				}
				energy += J_m[i]/2 * other_spins;
				other_spins = 0;
			}
			energy += H_m;
			energy *= -field_m[index];

			return energy;
		}

	public:
		Ising_t() : size_m(), cr_m(), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m() {}
		Ising_t(const Lattice_t<dim>& l, const std::array<size_t, dim> & s, bool periodic = false)
			: size_m(s), cr_m(l), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m()
		{
			setup_field();
			setup_crystal();
			setup_crystal_sites();
			setup_crystal_lattice_vectors(periodic);
			setup_nearest_neighbour_shells(1);
		}

		void set_interaction_parameters(const std::vector<double>& J)
		{
			J_m = J;
			if(J_m.size() > 2){
				std::cout << "Recalculating nearest neighbours to enable inclusion of (at least) " << J_m.size() << " nearest neighbour shells\n";
				setup_nearest_neighbour_shells(J_m.size()/2 + 1);
			}
		}

		void set_H(const double H){H_m = H;}
		void set_beta(const double beta){beta_m = beta;}

		void update()
		{
			double e_site, e_trial;

			size_t length = calc_length();
			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<size_t> dist_i(0, length - 1);
			std::uniform_real_distribution<double> dist_d(0., 1.);
			size_t index;

			index = dist_i(gen);
			e_site = site_energy(index);
			field_m[index] = static_cast<int8_t>(-field_m[index]);
			e_trial = site_energy(index);

			if( e_trial > e_site && dist_d(gen) > GSL::exp(-beta_m*(e_trial - e_site)).val){
				field_m[index] = static_cast<int8_t>(-field_m[index]);
			}

		}

		double total_energy()
		{
			double res = 0;
			for(size_t i = 0; i < field_m.size(); i++){
				res += site_energy(i);
			}
			return res;
		}

		double average_site_energy()
		{
			size_t n_sites = calc_length();
			return total_energy()/static_cast<double>(n_sites);
		}


		double magnetization()
		{
			double length = static_cast<double>(calc_length());
			int tot_spin = 0;
			for(auto spin : field_m){
				tot_spin += spin;
			}
			return static_cast<double>(tot_spin)/length;
		}

		std::vector<int8_t>& field(){

			return field_m;
		}	
};

#endif // ISING_H
