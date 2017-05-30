
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
	//std::cout << "half:\t" << boost::math::constants::half<arbpp::arb>() << std::endl;
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
		// vagy kell-e egyaltalan? vegulis az arb az lebegopontosan szamol, ez meg egeszekkel dolgozik - megbeszelesunk alapjan nem kell, mindenesetre benn hagyom
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

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/laguerre.hpp>
#include <boost/math/special_functions/hermite.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
namespace special_functions
{
	void arbpp_demo_gamma()
	{
		arbpp::arb A(4.2);
		arbpp::arb B(5);
		
		//std::cout << boost::math::tgamma(A) << std::endl;	// TODO invalid static_cast from arbpp::arb to int
		//std::cout << boost::math::tgamma(B) << std::endl;
		std::cout << std::endl;
	}
	
	void arbpp_demo_factorial()
	{
		arbpp::arb A(5);
		
		//std::cout << A << "! : " << boost::math::factorial<arbpp::arb>(5) << std::endl;	// TODO invalid static_cast to int
		std::cout << std::endl;
	}
	
	void arbpp_demo_legendre()
	{
		arbpp::arb A(0.5);	// abs(A) <= 1, kulonben domain_error
		int n = 2;
		int m = 3;
		
		// legendre_p(n, x) = 1 / (2^n * n!) * (d^l/dx^l)((x^2 - 1)^n)
		// legendre_p(n, m, x) = (-1)^m * (1 - x^2)^(m/2) * (d^m/dx^m)legendre_p(n,x)
		
		// legendre_q(0, x) = 1/2 ln((1 + x) / (1 - x))
		// legendre_q(1, x) = x/2 ln((1 + x) / (1 - x)) - 1
			// a dokumentacio nem ad rendes kepletet, csak ezt a ket peldat
		
		std::cout << "legendre_p(" << n << ", " << A << ") : " << boost::math::legendre_p(n, A) << std::endl;
		//std::cout << "legendre_p(" << n << ", " << m << ", " << A << ") : " << boost::math::legendre_p(n, m, A) << std::endl;	// TODO invalid static_cast to int
		std::cout << "legendre_q(" << 0 << ", " << A << ") : " << boost::math::legendre_q(0, A) << std::endl;
		std::cout << "legendre_q(" << n << ", " << A << ") : " << boost::math::legendre_q(n, A) << std::endl;
		std::cout << std::endl;
	}
	
	void arbpp_demo_laguerre()
	{
		arbpp::arb A(5.3);
		int n = 3;
		int m = 2;
		
		// laguerre(n, x) = e^x / n! * (d^n / dx^n) (x^n * e^-x)
		// laguerre(n, m, x) = (-1)^m * (d^m / dx^m) (laguerre(n + m, x))
		
		std::cout << "laguerre(" << n << ", " << A << ") : " << boost::math::laguerre(n, A) << std::endl;
		std::cout << "laguerre(" << n << ", " << m << ", " << A << ") : " << boost::math::laguerre(n, m, A) << std::endl;
		std::cout << std::endl;
	}
	
	void arbpp_demo_hermite()
	{
		arbpp::arb A(0.8);
		int n = 3;
		
		// hermite(n, A) = (-1)^n * e^(x^2) * (d^2 / dx^2) (e^(-x^2))
		
		std::cout << "hermite(" << n << ", " << A << ") : " << boost::math::hermite(n, A) << std::endl;
		std::cout << std::endl;
	}
	
	void arbpp_demo_spherical_harmonic()
	{
		arbpp::arb A(1.2);	// [0, pi]
		arbpp::arb B(-0.6);	// [0, 2pi)
		int n = 3;
		int m = 2;
		
		// spherical_harmonic(n, m, theta, phi) = sqrt(((2 * n + 1) * (n - m)!) / ((4 * pi) * (n + m)!)) * Pmn(cos theta) * e^(i * m * phi)
			// ez std::complex<type>-ot ad vissza, _r es _i valosakat
		
		//std::cout << "spherical_harmonic_r(" << n << ", " << m << ", " << A << ", " << B << ") : " << boost::math::spherical_harmonic_r(n, m, A, B) << std::endl;	// TODO static cast
		std::cout << std::endl;
	}
}

int main()
{
	//arbpp_demo_boostround();
	//arbpp_demo_boosttrunc();
	//arbpp_demo_modf();
	//arbpp_demo_sign();
	//arbpp_demo_io();
	//arbpp_demo_const();	// unfinished
	//arbpp_demo_complex();	// partially unfinished
	//arbpp_demo_gcd_lcm();	// unfinished
	//deriv_roots::arbpp_demo_roots();
	//noderiv_roots::arbpp_demo_roots();
	//minimize_func::arbpp_demo_minimize();
	//special_functions::arbpp_demo_gamma();	// unfinished
	//special_functions::arbpp_demo_factorial();	// unfinished
	//special_functions::arbpp_demo_legendre();	// partially unfinished
	//special_functions::arbpp_demo_laguerre();
	//special_functions::arbpp_demo_hermite();
	//special_functions::arbpp_demo_spherical_harmonic();	// unfinished
	
	return 0;
}
