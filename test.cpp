
#include <iostream>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include "arbpp.hpp"

// http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/real_concepts.html megvalositasok tesztelese
	// https://github.com/bluescarni/arbpp kiegeszitesebol

// round.hpp
void arbpp_demo_boostround()
{
	arbpp::arb A(1.4d);
	arbpp::arb B(2.5d);
	// megjegyzes double-rol konvertalva arb tipusra belekerulhet kerekitesi hiba
	
	// kerekito fuggvenyek
	std::cout << "round(" << A << ") : " << boost::math::round(A) << std::endl;
	std::cout << "round(" << B << ") : " << boost::math::round(B) << std::endl;
	// a megfelelo boostos iround, lround, llround fuggvenyek nincsenek implementalva
		// ehhez az arb tipusnak static_cast-olhatonak kellene lennie a megfelelo tipusokra
		// long long eseteben pedig konstruktor is kellene hozza
	
	std::cout << std::endl;
}

// trunc.hpp
void arbpp_demo_boosttrunc()
{
	arbpp::arb A(1.4d);
	arbpp::arb B(-2.5d);
	
	// trunc fuggvenyek
	std::cout << "trunc(" << A << ") : " << boost::math::trunc(A) << std::endl;
	std::cout << "trunc(" << B << ") : " << boost::math::trunc(B) << std::endl;	// negativ esetben felfele kerekit
	// itrunc, ltrunc, lltrunc - ugyanaz mint a kerekitesnel
}

int main()
{
	arbpp_demo_boostround();
	arbpp_demo_boosttrunc();	
	
	return 0;
}