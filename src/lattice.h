#ifndef LATTICE_H
#define LATTICE_H

#include "GSLpp/matrix.h"
#include "GSLpp/vector.h"
#include "GSLpp/linalg.h"

template<size_t dim, class T = double, class M = GSL::Matrix, class V = GSL::Vector>
class Lattice_t
{
private:
    M lat_m, recip_lat_m;
    T scale_m;
public:
    Lattice_t() : lat_m{}, recip_lat_m{}, scale_m{} {};
    Lattice_t(const std::initializer_list<V>& l)
	    : lat_m(l), recip_lat_m(), scale_m(0)
    {
	    for(auto elem : lat_m[0]){
		    scale_m += elem*elem;
	    }
	    scale_m = std::sqrt(scale_m);
	    lat_m *= 1./scale_m;
	    recip_lat_m = 2*M_PI*lat_m.inverse();
    }


    T scale() {return scale_m;}
    M lat() {return scale_m*lat_m;}
    M recip_lat() {return 1/scale_m*recip_lat_m;}
};

#endif // LATTICE_H
