#include <bitset>
#include <iostream> 
#include <iomanip>
#include <fstream> 
#include <random>
#include "lattice.h"
#include "crystal.h"
#include "ising.h"
#include "GSLpp/error.h"

struct Bitmap_header{	
	uint16_t bm;
	uint32_t size, reserved, offset;
	Bitmap_header() 
		: bm(0x4D42), size(0), reserved(0), offset(14 + 56) 
	{} 
};

std::ostream& operator<<(std::ostream& os, const Bitmap_header& bmh)
{
	os << bmh.bm << " " << bmh.size << " " << bmh.reserved << " ";
        os << bmh.offset;
	return os;
}

struct DIB_header{	
	uint32_t size;
	int32_t width, height;
	uint16_t planes, bit_per_pixel;
	uint32_t compression, bitmap_size;
	int32_t h_resolution, v_resolution;
	uint32_t colors, important_colors;
	uint16_t resolution_units, padding, halftoning;
	uint32_t halftoning_param1, halftoning_param2,
		 color_encoding, reserved,
		red_bitmask, green_bitmask, blue_bitmask, 
		alpha_bitmask;//, color_space, reserved[16]{0} ; 
	DIB_header()
		: size(124), width(0), height(0), planes(1), bit_per_pixel(32), 
		  compression(3), bitmap_size(0), h_resolution(2835), 
		  v_resolution(2835), colors(0), important_colors(0),
		  resolution_units(0), padding(0),  halftoning(0),
		  halftoning_param1(0), halftoning_param2(0),
		  color_encoding(0), reserved(0),
		  red_bitmask(0x00FF0000), green_bitmask(0x0000FF00), 
		  blue_bitmask(0x000000FF), alpha_bitmask(0xFF000000)
	{} 
};

std::ostream& operator<<(std::ostream& os, const DIB_header& dibh)
{
	os << dibh.size << " " << dibh.width << " " << dibh.height << " ";
	os << dibh.planes << " " << dibh.bit_per_pixel << " ";
	os << dibh.compression << " " << dibh.bitmap_size << " ";
	os << dibh.h_resolution << " " << dibh.v_resolution << " "; 
	os << dibh.colors << " " << dibh.important_colors << " ";
	os << dibh.resolution_units << " " << dibh.padding << " ";
	os << dibh.halftoning << " " << dibh.halftoning_param1 << " ";
	os << dibh.halftoning_param2 << " " << dibh.color_encoding << " ";
	os << dibh.reserved << " " << dibh.red_bitmask << " ";
	os << dibh.green_bitmask << " " << dibh.blue_bitmask << " ";
	os << dibh.alpha_bitmask;
	return os;
}

void bitmap_print(const std::vector<int> field, const std::string filename, const std::array<size_t, 2>& size) 
{
	Bitmap_header bmh; 
	DIB_header dibh;
	dibh.width = static_cast<int32_t>(size[0]);
	dibh.height = static_cast<int32_t>(size[1]);
	dibh.bitmap_size = static_cast<uint32_t>(dibh.width * dibh.height * dibh.bit_per_pixel/8);
	bmh.size = bmh.offset  + dibh.bitmap_size;

	std::ofstream file(filename, std::ios::out | std::ios::binary);
	file.write(reinterpret_cast<char*>(&bmh), 14);
	file.write(reinterpret_cast<char*>(&dibh), 124);

	int8_t red[4] {0x00, 0x7F, -1, -0x0F};
	//int8_t red[4] {-1, -1, -1, -1};
	int8_t blue[4] {-1, 0x7F, 0x00, -0x0F};
	//int8_t blue[4] {0, 0, 0, -1};
	for(auto row = size[0]; row > 0; row--){
		for(size_t column = 0; column < size[1]; column++){
			if(field[(row-1)*size[1] + column] > 0){
			    file.write(reinterpret_cast<char*>(&red), 4);
			}else{
			    file.write(reinterpret_cast<char*>(&blue), 4);
			}
		}
	}
	file.close();
}


int main() 
{ 
	GSL::Error_handler e_handler;
	e_handler.off();

	size_t width = 500, height = 500;
	Lattice_t<2> lat{{static_cast<double>(width), 0}, 
		{0, static_cast<double>(height)}}; 
	Ising_t<2> ising(lat, {width, height}, true);
	ising.set_beta(10);

	// One interaction parameter per nearest neighbour shell to consider
	ising.set_interaction_parameters({1.0, 0, -0.2});
	ising.set_H(0);

	bitmap_print(ising.field(), "initial.bmp", {width, height});
	size_t num_iterations = 100*width*height;
	for(size_t it = 0; it < num_iterations; it++){
		ising.update();
	}

	std::cout << "Total energy = " << ising.total_energy() << "\n";
	std::cout << "Magnetization = " << ising.magnetization() << "\n";
	bitmap_print(ising.field(), "final.bmp", {width, height});

	return 0; 
}
