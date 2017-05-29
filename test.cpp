
#include <iostream>
#include "arbpp.hpp"

// boost 1.55 verzioval dolgoztam, ami ezert hianyozhat:
	// polinom aritmetika (1.61)
	// hyperexponential eloszlas (1.57)
	// bernoulli_numbers (1.56)
	// eloszlasok max ertek helyett overflow_error -t dobhatnak (1.56)	-	https://svn.boost.org/trac/boost/ticket/10111a 

#include <boost/math/special_functions/round.hpp>
void arbpp_demo_boostround()
{
	arbpp::arb A(1.4d);
	arbpp::arb B(-2.5d);
	// megjegyzes double-rol konvertalva arb tipusra belekerulhet kerekitesi hiba (a tobbi fv-nel is)
	
	// kerekito fuggvenyek
	std::cout << "round(" << A << ") : " << boost::math::round(A) << std::endl;
	std::cout << "round(" << B << ") : " << boost::math::round(B) << std::endl;	// -2.5-ot -3-ra kerekit - "round away from zero" vagy nem teljesen pontos kiindulo ertek
	// a megfelelo boostos iround, lround, llround fuggvenyek nincsenek implementalva
		// ehhez az arb tipusnak static_cast-olhatonak kellene lennie a megfelelo tipusokra
		// long long eseteben pedig konstruktor is kellene hozza
	std::cout << std::endl;
}

#include <boost/math/special_functions/trunc.hpp>
void arbpp_demo_boosttrunc()
{
	arbpp::arb A(1.4d);
	arbpp::arb B(-2.5d);
	
	// trunc fuggvenyek
	std::cout << "trunc(" << A << ") : " << boost::math::trunc(A) << std::endl;
	std::cout << "trunc(" << B << ") : " << boost::math::trunc(B) << std::endl;	// -2.5-ot -2-re vag
	// itrunc, ltrunc, lltrunc - ugyanaz mint a kerekitesnel
	std::cout << std::endl;
}

#include <boost/math/special_functions/modf.hpp>
void arbpp_demo_modf()
{
	arbpp::arb A(1.4d);
	arbpp::arb B(-2.5d);
	arbpp::arb C;
	arbpp::arb D;
	
	// modf fuggvenyek - egeszresz-tortresz bontas, mindket eredmeny elojele megegyezik az eredetivel
	C = boost::math::modf(A, &D);	// D = int(A), C = A - D
	std::cout << "modf(" << A << ") i: " << D << " frac: " << C << std::endl;
	C = boost::math::modf(B, &D);
	std::cout << "modf(" << B << ") i: " << D << " frac: " << C << std::endl;
	std::cout << std::endl;
}

//#include <boost/math/special_functions/fpclassify.hpp>	// nem tudtam beincludeolni, feltehetoleg std::numeric_limits kell
	// vegesseg, NaN-sag, denormalizaltsag ellenorzesek, ezek alapjan osztalyozas

#include <boost/math/special_functions/sign.hpp>
void arbpp_demo_sign()
{
	arbpp::arb A(3);
	arbpp::arb B(0);
	arbpp::arb C(-2);
	
	// elojelbit
	std::cout << "signbit(" << A << ") : " << boost::math::signbit(A) << std::endl;
	std::cout << "signbit(" << C << ") : " << boost::math::signbit(C) << std::endl;
	// elojelfv
	std::cout << "sign(" << A << ") : " << boost::math::sign(A) << std::endl;
	std::cout << "sign(" << B << ") : " << boost::math::sign(B) << std::endl;
	std::cout << "sign(" << C << ") : " << boost::math::sign(C) << std::endl;
	// elojelmasolas
	std::cout << "copysign(" << A << ", " << C <<") : " << boost::math::copysign(A, C) << std::endl;
}

//#include <boost/math/special_functions/nonfinite_num_facets.hpp>	// nem tudtam beincludeolni, feltehetoleg std::numeric_limits kell
	// string konverzio es IO a vegtelen es NaN ertekekhez (kiiras es beolvasas) - a wrapper mar ezeket tudja operatorokkal
void arbpp_demo_io()
{
	arbpp::arb A(2.0);
	arbpp::arb inf = arbpp::arb::pos_inf();
	arbpp::arb minf = arbpp::arb::neg_inf();
	
	std::cout << A << "\t" << inf << "\t" << minf << std::endl;
	
	std::cout << "Irj be egy szamot!" << std::endl;
	arbpp::arb B;
	std::cin >> B;	// tesztelve: +-inf, NaN, egesz, tizedes tort, tizedes tort e10 hatvannyal
		// TODO egyelore std::invalid_argument-et dob, ha a konverzio nem mukodik, kivant mukodesre le lehetne cserelni
	std::cout << B << std::endl;
}

/* //TODO - causes segfault on specialization
#include <boost/math/constants/constants.hpp>
namespace boost{
namespace math{
namespace constants{
template<> arbpp::arb half<arbpp::arb>() { return arbpp::arb("0.5"); }	// this makes it work
}}}
void arbpp_demo_const()
{
	std::cout << "half:\t" << boost::math::constants::half<arbpp::arb>() << std::endl;
	//std::cout << "pi:\t" << boost::math::constants::pi<arbpp::arb>() << std::endl;
}
*/

//#include <boost/math/complex/asin.hpp>
//#include <boost/math/complex/acos.hpp>
//#include <boost/math/complex/atan.hpp>
//#include <boost/math/complex/asinh.hpp>
//#include <boost/math/complex/acosh.hpp>
//#include <boost/math/complex/atanh.hpp>
	// a headerek includeolasa nelkul is mukodnek
void arbpp_demo_complex()
{
	arbpp::arb r(2.3);
	arbpp::arb i(-0.9);
	std::complex<arbpp::arb> C(r, i);
	
	// complex
	std::cout << "complex: " << C << std::endl;
	std::cout << "asin(" << C << ") : " << asin(C) << std::endl;
	//std::cout << "acos(" << C << ") : " << acos(C) << std::endl;	// TODO long double arb konstruktor
	std::cout << "atan(" << C << ") : " << atan(C) << std::endl;
	std::cout << "asinh(" << C << ") : " << asinh(C) << std::endl;
	std::cout << "acosh(" << C << ") : " << acosh(C) << std::endl;
	std::cout << "atanh(" << C << ") : " << atanh(C) << std::endl;
	std::cout << std::endl;
}

/*
#include <boost/math/common_factor_rt.hpp>
void arbpp_demo_gcd_lcm()
{
	arbpp::arb A(12);
	arbpp::arb B(8);
	
	// legkisebb kozos oszto, legnagyobb kozos tobbszoros
		// TODO mukodeshez operator%= implementalasa
	std::cout << "gcd(" << A << ", " << B << ") : " << boost::math::gcd(A, B) << std::endl;
	std::cout << "lcm(" << A << ", " << B << ") : " << boost::math::lcm(A, B) << std::endl;
	std::cout << std::endl;
}
*/

int main()
{
	//arbpp_demo_boostround();
	//arbpp_demo_boosttrunc();
	//arbpp_demo_modf();
	//arbpp_demo_sign();
	//arbpp_demo_io();
	//arbpp_demo_const();	// TODO
	arbpp_demo_complex();
	//arbpp_demo_gcd_lcm();	// TODO
	
	return 0;
}
