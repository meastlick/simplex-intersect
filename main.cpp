#include "intersect.hpp"

#include <iostream>
#include <cstdlib>

template <typename value_type>
bool intersect(const simplex<2, 3, value_type> &triangle_a, const simplex<2, 3, value_type> &triangle_b);

template <typename value_type>
value_type random(const value_type &maximum)
{
    return maximum * rand() / RAND_MAX ;
}

template <typename value_type>
point<3, value_type> make_random_point3(const value_type &maximum)
{
    return make_point(random(maximum), random(maximum), random(maximum)) ;
}

// generates a triangle in 3-dimensions whose vertices are uniformly spread within cube (0, 0, 0).(maximum, maximum, maximum)

template <typename value_type>
simplex<2, 3, value_type> make_random_simplex23(const value_type &maximum)
{
    return make_simplex(make_random_point3(maximum), make_random_point3(maximum), make_random_point3(maximum)) ;
}

int main(int argc, const char *argv[])
{
    typedef double coord_type ;
    typedef point<3, coord_type> point3_type ;
    typedef simplex<2, 3, coord_type> simplex23_type ;

	srand (13); //initialise random seed = 13

	const int num_trials=1000000;
	int num_intersects=0;

	for(int i=0;i<num_trials;i++)
	{
		const simplex23_type triangle_a = make_random_simplex23(100.0),
                         triangle_b = make_random_simplex23(100.0) ;
		if(intersect(triangle_a, triangle_b))
			num_intersects++;
	}

	std::cout << "Probability of intersection: : " << (coord_type) num_intersects/num_trials << std::endl ;

	std::cin.ignore();

	return 0 ;
}
