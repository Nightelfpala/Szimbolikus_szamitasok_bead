
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
#include "arbpp.hpp"

void arbpp_test_eq()
{
	// operator==
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	arbpp::arb D(3);
	std::cout << A << " == " << B << " : " << (A == B) << std::endl;
	std::cout << A << " == " << C << " : " << (A == C) << std::endl;
	std::cout << A << " == " << 1 << " : " << (A == 1) << std::endl;
	std::cout << 1.0f << " == " << B << " : " << (1 == B) << std::endl;
	std::cout << A << " == " << 2 << " : " << (A == 2) << std::endl;
	std::cout << 2.0f << " == " << B << " : " << (2.0 == B) << std::endl;
}

void arbpp_add_test()
{
	arbpp_test_eq();
}

int main()
{
	arbpp_add_test();
	
	
	return 0;
}