
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
template<> arbpp::arb half<arbpp::arb>() { return arbpp::arb("5.00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000e-01"); }	// ezzel mukodik, de valoszinuleg ertelmetlenne teszi
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
	//std::cout << "acos(" << C << ") : " << acos(C) << std::endl;	// TODO long double arb konstruktor implementalasa
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
		// vagy kell-e egyaltalan? vegulis az arb az lebegopontosan szamol, ez meg egeszekkel dolgozik
	std::cout << "gcd(" << A << ", " << B << ") : " << boost::math::gcd(A, B) << std::endl;
	std::cout << "lcm(" << A << ", " << B << ") : " << boost::math::lcm(A, B) << std::endl;
	std::cout << std::endl;
}
*/

#include <boost/math/tools/roots.hpp>
	// fuggveny gyokenek keresese
namespace deriv_roots
{
struct roots_func_D1	// f(x) = x^2 - 1 -> gyokei -1, 1
{
	boost::math::tuple<arbpp::arb, arbpp::arb> operator()(const arbpp::arb &a)
	{
		return boost::math::make_tuple( a * a, 2 * a );	// f, f'
	}
};
struct roots_func_D2	// f(x) = x^2 - 1 -> gyokei -1, 1
{
	boost::math::tuple<arbpp::arb, arbpp::arb, arbpp::arb> operator()(const arbpp::arb &a)
	{
		return boost::math::make_tuple( a * a, 2 * a, arbpp::arb(2));	// f, f', f''
	}
};
void arbpp_demo_roots()
{
	roots_func_D1 fD1;
	roots_func_D2 fD2;
	arbpp::arb startv(0.6);
	arbpp::arb minv(0.25);
	arbpp::arb maxv(5);
	int prec = 20;
	// opcionalis: max_iter, kulonben gyakorlatilag vegtelensegig futhat ha nem eri el a kivant pontossagot
	
	// derivaltat hasznalo iterativ gyokkereses
	//arbpp::arb ret = boost::math::tools::newton_raphson_iterate(fD1, startv, minv, maxv, prec);	// elso derivalt kell iteraciohoz
	//arbpp::arb ret = boost::math::tools::halley_iterate(fD2, startv, minv, maxv, prec);	// masodik derivalt kell
	arbpp::arb ret = boost::math::tools::schroeder_iterate(fD2, startv, minv, maxv, prec);
	
	std::cout << "x^2 - 1 gyoke [-1, 1]-ben: " << ret << std::endl;
	std::cout << ret << "^2 + 1 = " << (pow(ret, arbpp::arb(2)) + 1) << std::endl;
	std::cout << std::endl;
}
}	// namespace deriv_roots

namespace noderiv_roots
{
struct roots_func	// f(x) = x^3 - 1 -> gyoke 1
{
	arbpp::arb operator()(const arbpp::arb &a)
	{
		return a * a * a - 1;
	}
};
struct tolerance
{
	bool operator()(const arbpp::arb &a, const arbpp::arb &b)	// ha az intervallumszelekkel meghivva igazat ad vissza, akkor terminal az algoritmus
	{
		return (abs(b - a) <= arbpp::arb(1e-5));
	}
};
void arbpp_demo_roots()
{
	roots_func f;
	arbpp::arb minv(-0.5);
	arbpp::arb maxv(3.25);
	tolerance tol;
	// opcionalis: max_iter
	
	//std::pair<arbpp::arb, arbpp::arb> ret = boost::math::tools::bisect(f, minv, maxv, tol);
		// feltetel: f(min) * f(max) <= 0 && min < max, kulonben hibat ad
	
	arbpp::arb guess(0.5);
	arbpp::arb factor(2);	// ennyival oszt/szoroz lepesenkent
	bool rising = true;	// mon novo a fv
	boost::uintmax_t max_iter = 50;	// itt mar nem opcionalis
	//std::pair<arbpp::arb, arbpp::arb> ret = boost::math::tools::bracket_and_solve_root(f, guess, factor, rising, tol, max_iter);
		// csak monoton fv-re mukodik
		// erdekes modon, ha itt a guess rossz volt (-0.5), akkor nem tudta megtalalni
		
	std::pair<arbpp::arb, arbpp::arb> ret = boost::math::tools::toms748_solve(f, minv, maxv, tol, max_iter);
		// ennek elvileg nem kell monoton fv
	
	std::cout << "x^3 - 1 gyoke a kovetkezo intervallumban van: [" << ret.first << ", " << ret.second << "]" << std::endl;
	std::cout << std::endl;
}
}	// namespace noderiv_roots

#include <boost/math/tools/minima.hpp>
	// fuggveny minimumhely kereses
namespace minimize_func
{
struct min_func	// f(x) = x^2 minimalizalasa [-1, 1]-en
{
	arbpp::arb operator()(const arbpp::arb &a)
	{
		return a * a;
	}
};
void arbpp_demo_minimize()
{
	min_func f;
	arbpp::arb minv(-1);
	arbpp::arb maxv(1);
	int bits = 40;	// pontossag
	// opcionalis: max_iter
	
	std::pair<arbpp::arb, arbpp::arb> ret = boost::math::tools::brent_find_minima(f, minv, maxv, bits);
		// az intervallumban nem lehet maximum (?)
	std::cout << "x^2 minimuma a kovetkezo intervallumban van: [" << ret.first << ", " << ret.second << "]" << std::endl;
	std::cout << std::endl;
}
}	// namespace minimize_func

int main()
{
	//arbpp_demo_boostround();
	//arbpp_demo_boosttrunc();
	//arbpp_demo_modf();
	//arbpp_demo_sign();
	//arbpp_demo_io();
	//arbpp_demo_const();	// TODO
	//arbpp_demo_complex();
	//arbpp_demo_gcd_lcm();	// TODO
	//deriv_roots::arbpp_demo_roots();
	//noderiv_roots::arbpp_demo_roots();
	minimize_func::arbpp_demo_minimize();
	
	return 0;
}
