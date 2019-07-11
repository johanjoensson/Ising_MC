#ifndef ISING_H
#define ISING_H

#include <vector>
#include <tuple>
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
	using Site = Site_t<dim>;
	private:
		std::array<size_t, dim> size_m;
		Crystal_t<dim> cr_m;
		std::vector<int8_t> field_m;
		std::vector<std::vector<Neighbours<dim>>> nn_shells_m;
		std::vector<double> J_m;
		double H_m;
		double beta_m;
		std::vector<std::pair<size_t, size_t>> correlators_m;

		size_t calc_length() const 
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

		GSL::Vector calc_coords(const size_t index)
		{
			GSL::Vector res(dim);
			size_t offset = 1;
			size_t tmp_index = index;
			for(size_t i = dim; i > 0; i--)
			{
				offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size_m[j];
				}
				res[i-1] = static_cast<double>(tmp_index / offset);
				tmp_index -= (tmp_index/offset)*offset;
			}
			return res;
		}

		size_t calc_index(const std::array<size_t, dim>& coords)
		{
			size_t res = 0, offset;
			for(size_t i = dim; i > 0; i--)
			{
				offset = 1;
				for(size_t j = 0; j < i - 1; j++ ){
					offset *= size_m[j];
				}
				res += coords[i]*offset;
			}
			return res;

		}

		void setup_crystal_sites()
		{
			size_t length = calc_length();
			std::vector<GSL::Vector> positions(length);
			GSL::Vector coord(dim);
			for(size_t idx = 0; idx < length; idx++){
				coord = calc_coords(idx);
				positions[idx] = 1./cr_m.lat().scale()*cr_m.lat().lat()*coord;
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

		double site_energy(const size_t index) const
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
		Ising_t() : size_m(), cr_m(), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m(), correlators_m() {}
		Ising_t(const Lattice_t<dim>& l, const std::array<size_t, dim> & s, bool periodic = false)
			: size_m(s), cr_m(l), field_m(), nn_shells_m(), J_m(), H_m(0), beta_m(), correlators_m()
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

		void add_spin_spin_correlator(const size_t i, const size_t j)
		{
			size_t dist = i > j ? i - j : j - i;
			size_t dcol = dist/size_m[0];
			size_t drow = dist % size_m[0];
			if(nn_shells_m[i].size() <= std::max(dcol, drow)){
				setup_nearest_neighbour_shells(std::max(drow, dcol)/2 + 1);
			}
			correlators_m.push_back(std::make_pair(i, j));
		}

		void add_spin_correlator(const size_t index, const double r_max)
		{
			double r = 0;
			size_t j = 0;
			while(r <= r_max && j < nn_shells_m[index].size() && nn_shells_m[index][j].size() > 0){
				r = nn_shells_m[index][j][0].pos(). template norm<double>();
				for(auto site : nn_shells_m[index][j]){
					correlators_m.push_back(std::make_pair(index, site.index()));
				}
				j++;
			}

		}

		void add_spin_correlator(const std::array<size_t, dim>& i, const double r_max)
		{
			add_spin_correlator(calc_index(i), r_max);
		}

		void add_spin_correlator(const std::array<size_t, dim>& i, const std::array<size_t, dim>& j)
		{
			size_t i_index = 0;
			size_t j_index = 0;

			add_spin_correlator(i_index, j_index);
		}

		std::vector<int8_t>& field(){return field_m;}	

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

		int spin_spin(const size_t i, const size_t j) const
		{
			return static_cast<int>(field_m[i]*field_m[j]);
		}
		
		std::vector<std::tuple<double, int, double, double>> measure_spin_correlators()
		{
			std::vector<std::tuple<double, int, double, double>> res(correlators_m.size());
			size_t i, j;
			int si, ss, sj;
			double r = 0;
			std::pair<size_t, size_t> corr;
			for(auto index = 0; index < correlators_m.size(); index++){
				corr = correlators_m[index];
				i = corr.first;
				j = corr.second;
				si = field_m[i];
				sj = field_m[j];
				ss = si*sj;
				r = (1./cr_m.lat().scale()*cr_m.lat().lat()*
					(calc_coords(i) - calc_coords(j))).
					 template norm<double>();
				res[index] = std::make_tuple(r, ss, si, sj);
			}
			return res;
		}                                                     

		double total_energy() const
		{
			double res = 0;
			for(size_t i = 0; i < field_m.size(); i++){
				res += site_energy(i);
			}
			return res;
		}

		double average_site_energy() const
		{
			size_t n_sites = calc_length();
			return total_energy()/static_cast<double>(n_sites);
		}


		double magnetization() const
		{
			double length = static_cast<double>(calc_length());
			int tot_spin = 0;
			for(auto spin : field_m){
				tot_spin += spin;
			}
			return static_cast<double>(tot_spin)/length;
		}

};

#endif // ISING_H
