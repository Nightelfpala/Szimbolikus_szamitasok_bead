
#include <iostream>
#include "arbpp.hpp"

//TODO - causes segfault on specialization
#include <boost/math/constants/constants.hpp>
/*
namespace boost{
namespace math{
namespace constants{
template<> arbpp::arb half<arbpp::arb>() { return arbpp::arb("5.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-01"); }	// ezzel mukodik, de valoszinuleg ertelmetlenne teszi
}}}
*/
void arbpp_demo_const()
{
	std::cout << "half:\t" << boost::math::constants::half<arbpp::arb>() << std::endl;
	//std::cout << "pi:\t" << boost::math::constants::pi<arbpp::arb>() << std::endl;
}


int main()
{
	arbpp_demo_const();
	
	return 0;
}
