
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

#include <boost/math/special_functions/fpclassify.hpp>	// std::numeric_limits kell hogy jol mukodjon
void arbpp_demo_fpclass()
{
	arbpp::arb A(2);
	arbpp::arb B(arbpp::arb::pos_inf());
	
	int cl = boost::math::fpclassify(A);
	std::string str = (cl == FP_ZERO ? "FP_ZERO"
		: cl == FP_NORMAL ? "FP_NORMAL"
		: cl == FP_INFINITE ? "FP_INFINITE"
		: cl == FP_NAN ? "FP_NAN"
		: cl == FP_SUBNORMAL ? "FP_SUBNORMAL"
		: "not found");
	
	// vegesseg, NaN-sag, denormalizaltsag ellenorzesek, ezek alapjan osztalyozas
	std::cout << "fpclassify(" << A << ") : " << str << std::endl;
	std::cout << "isfinite(" << A << ") : " << boost::math::isfinite(A) << std::endl;
	std::cout << "isinf(" << A << ") : " << boost::math::isinf(A) << std::endl;
	std::cout << "isinf(" << B << ") : " << boost::math::isinf(B) << std::endl;	// nem jo, TODO numeric_limits
	std::cout << "isnan(" << A << ") : " << boost::math::isnan(A) << std::endl;
	std::cout << "isnormal(" << A << ") : " << boost::math::isnormal(A) << std::endl;
	std::cout << std::endl;
}

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

//#include <boost/math/special_functions/nonfinite_num_facets.hpp>	// feltehetoleg std::numeric_limits kell
	// string konverzio es IO a vegtelen es NaN ertekekhez (kiiras es beolvasas) - a wrapper mar ezeket tudja operatorokkal, bar ez nem standardizalt, de a wrapper kiegeszitheto
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

#include <boost/math/constants/constants.hpp>
	// TODO - causes segfault on specialization
//#include "arbppconsts.h"
	// this explicitly specializes constants with the values taken from the boost header, avoids segmentation faults
void arbpp_demo_const()
{
	//std::cout << "half:\t" << boost::math::constants::half<arbpp::arb>() << std::endl;
	//std::cout << "pi:\t" << boost::math::constants::pi<arbpp::arb>() << std::endl;
}


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
		// mindharom fuggveny mukodik
	
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
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/hankel.hpp>	// a header nev nincs a dokumentacioban benne
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/ellint_rf.hpp>
#include <boost/math/special_functions/ellint_rd.hpp>
#include <boost/math/special_functions/ellint_rj.hpp>
#include <boost/math/special_functions/ellint_rc.hpp>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/ellint_3.hpp>
#include <boost/math/special_functions/jacobi_elliptic.hpp>
#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/special_functions/expint.hpp>
//#include <boost/math/special_functions/sin_pi.hpp>
//#include <boost/math/special_functions/cos_pi.hpp>
//#include <boost/math/special_functions/log1p.hpp>
//#include <boost/math/special_functions/expm1.hpp>
//#include <boost/math/special_functions/cbrt.hpp>
//#include <boost/math/special_functions/sqrt1pm1.hpp>
//#include <boost/math/special_functions/powm1.hpp>
//#include <boost/math/special_functions/hypot.hpp>
//#include <boost/math/special_functions/sinc.hpp>
#include <boost/math/special_functions/sinhc.hpp>
#include <boost/math/special_functions/asinh.hpp>
#include <boost/math/special_functions/acosh.hpp>
#include <boost/math/special_functions/atanh.hpp>
#include <boost/math/special_functions/owens_t.hpp>
namespace special_functions
{
	void gamma()
	{
		arbpp::arb A(4.2);
		arbpp::arb B(5);
		
		//std::cout << boost::math::tgamma(A) << std::endl;	// TODO invalid static_cast from arbpp::arb to int
		//std::cout << boost::math::tgamma(B) << std::endl;
		std::cout << std::endl;
	}
	
	void factorial()
	{
		arbpp::arb A(5);
		
		//std::cout << A << "! : " << boost::math::factorial<arbpp::arb>(5) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
	void legendre()
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
		//std::cout << "legendre_p(" << n << ", " << m << ", " << A << ") : " << boost::math::legendre_p(n, m, A) << std::endl;	// TODO static_cast int
		std::cout << "legendre_q(" << 0 << ", " << A << ") : " << boost::math::legendre_q(0, A) << std::endl;
		std::cout << "legendre_q(" << n << ", " << A << ") : " << boost::math::legendre_q(n, A) << std::endl;
		std::cout << std::endl;
	}
	
	void laguerre()
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
	
	void hermite()
	{
		arbpp::arb A(0.8);
		int n = 3;
		
		// hermite(n, A) = (-1)^n * e^(x^2) * (d^2 / dx^2) (e^(-x^2))
		
		std::cout << "hermite(" << n << ", " << A << ") : " << boost::math::hermite(n, A) << std::endl;
		std::cout << std::endl;
	}
	
	void spherical_harmonic()
	{
		arbpp::arb A(1.2);	// [0, pi]
		arbpp::arb B(-0.6);	// [0, 2pi)
		int n = 3;
		int m = 2;
		
		// spherical_harmonic(n, m, theta, phi) = sqrt(((2 * n + 1) * (n - m)!) / ((4 * pi) * (n + m)!)) * Pmn(cos theta) * e^(i * m * phi)
			// ez std::complex<type>-ot ad vissza, _r es _i valosakat
		
		//std::cout << "spherical_harmonic_r(" << n << ", " << m << ", " << A << ", " << B << ") : " << boost::math::spherical_harmonic_r(n, m, A, B) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
	void bessel()	// tgamma()-ra epul
	{
		arbpp::arb A(2);
		arbpp::arb B(3.1);
		
		// cyl_bessel_j(v, x) = (x / 2)^v * sum(k=0..inf, (-x^2 / 4)^k / (k! * gamma(v + k + 1)))
		// cyl_neumann(v, x) = (cyl_bessel_j(v, x) * cos(v * pi) - cyl_bessel_j(-v, x)) / sin(v * pi)
		
		//std::cout << "cyl_bessel_j(" << A << ", " << B << ") : " << boost::math::cyl_bessel_j(A, B) << std::endl;	// TODO static_cast int
		//std::cout << "cyl_neumann(" << A << ", " << B << ") : " << boost::math::cyl_neumann(A, B) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
	void hankel()	// cyl_bessel_j()-re epul
	{
		arbpp::arb A(2);
		arbpp::arb B(3.1);
		
		// cyl_hankel_1(v, x) = cyl_bessel_j(v, x) + i * cyl_neumann(v, x)
		// cyl_hankel_2(v, x) = cyl_bessel_j(v, x) - i * cyl_neumann(v, x)
		// sph_hankel_1(v, x) = sqrt(pi/2) * 1/sqrt(x) * cyl_hankel_1(v + 1/2, x)
		// sph_hankel_2(v, x) = sqrt(pi/2) * 1/sqrt(x) * cyl_hankel_2(v + 1/2, x)
			// std::complex<type>-ot adnak vissza
		
		//std::cout << "cyl_hankel_1(" << A << ", " << B << ") : " << boost::math::cyl_hankel_1(A, B) << std::endl;	// TODO static_cast int
		//std::cout << "cyl_hankel_2(" << A << ", " << B << ") : " << boost::math::cyl_hankel_2(A, B) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
	void airy()	// cyl_bessel_j()-re es cyl_bessel_k()-ra epul
	{
		arbpp::arb A(0);
		
		// (d^2 / dz^2)w = z * w diffegyenlet megoldasai es a derivaltjaik
			// a harmadik megoldas airy_ai(z * exp(+- 2 * pi * i / 3))
		
		//std::cout << "airy_ai(" << A << ") : " << boost::math::airy_ai(A) << std::endl;	// TODO numeric_limits
		//std::cout << "airy_bi(" << A << ") : " << boost::math::airy_bi(A) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "airy_ai_prime(" << A << ") : " << boost::math::airy_ai_prime(A) << std::endl;	// TODO numeric_limits
		//std::cout << "airy_bi_prime(" << A << ") : " << boost::math::airy_bi_prime(A) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void carlson_elliptic()
	{
		arbpp::arb A(1);
		arbpp::arb B(0.5);
		arbpp::arb C(1.5);
		arbpp::arb D(2.6);
		
		// ellint_rf(x, y, z) = 1/2 * integral(0..+inf) ((t + x) * (t + y) * (t + z))^(-1/2) dt
			// mind nemnegativ, legfeljebb az egyik 0	->	kulonben domain_error
		// ellint_rd(x, y, z) = 3/2 * integral(0..+inf) ((t + x) * (t + y))^(-1/2) * (t + z)^(-3/2) dt
			// x, y nemnegativ, legfeljebb az egyik 0; z nemnegativ
		// ellint_rj(x, y, z, p) = 3/2 * integral(0..+inf) (t + p)^(-1) ((t + x) * (t + y) * (t + z))^(-1/2) dt
			// x, y, z nemnegativ, legf az egyik 0; p nem 0
		// ellint_rc(x, y) = 1/2 * integral(0..+inf) (t + x)^(-1/2) * (t + y)^(-1) dt
			// x pozitiv; y nem 0
		
		std::cout << "ellint_rf(" << A << ", " << B << ", " << C << ") : " << boost::math::ellint_rf(A, B, C) << std::endl;
		std::cout << "ellint_rd(" << A << ", " << B << ", " << C << ") : " << boost::math::ellint_rd(A, B, C) << std::endl;
		std::cout << "ellint_rj(" << A << ", " << B << ", " << C << ", " << D <<  ") : " << boost::math::ellint_rj(A, B, C, D) << std::endl;
		std::cout << "ellint_rc(" << A << ", " << B << ") : " << boost::math::ellint_rc(A, B) << std::endl;
		std::cout << std::endl;
	}
	
	void legendre_elliptic()
	{
		arbpp::arb A(1);
		arbpp::arb B(0.5);
		arbpp::arb C(2);
		
		// ellint_1(phi, k) = integral(0..phi) 1 / sqrt(1 - k^2 * sin^2(t)) dt
		// ellint_1(k) = ellint_1(pi/2, k)
			// abs(k) <= 1 	-	elvileg, de k = 1-el tesztelve elszallt overflow_error-ral
		// ellint_2(phi, k) = integral(0..phi) sqrt(1 - k^2 * sin^2(t)) dt
		// ellint_2(k) = ellint_2(pi/2, k)
			// abs(k) <= 1 	-	ez nem szallt el k = 1-el
		// ellint_3(n, phi, k) = integral(0..phi) 1 / ((1 - n *sin^2(t)) * sqrt(1 - k^2 * sin^2(t))) dt
		// ellint_3(n, k) = ellint_3(n, pi/2, k)
			// abs(k) <= 1; n < 1
		
		//std::cout << "ellint_1(" << A << ", " << B << ") : " << boost::math::ellint_1(A, B) << std::endl;	// TODO fmod(arb, arb) not implemented	-	arb_t doesn't have the required function
		std::cout << "ellint_1(" << B << ") : " << boost::math::ellint_1(B) << std::endl;
		//std::cout << "ellint_2(" << A << ", " << B << ") : " << boost::math::ellint_2(A, B) << std::endl;	// TODO fmod
		std::cout << "ellint_2(" << B << ") : " << boost::math::ellint_2(B) << std::endl;
		//std::cout << "ellint_3(" << A << ", " << C << ", " << B << ") : " << boost::math::ellint_3(A, C, B) << std::endl;	// TODO fmod
		//std::cout << "ellint_3(" << B << ", " << A << ") : " << boost::math::ellint_3(B, A) << std::endl;	// TODO segfault - valoszinuleg constants.hpp-bol akar valamit arb-al lepeldanyositani (a main() elott szall el)
			// itt valoszinuleg a parameterezesem nem jo es azert dob hibat, de a segfault megszunt a konstansok behackelesevel
		std::cout << std::endl;
	}
	
	void jacobi_elliptic()
	{
		arbpp::arb A(1);
		arbpp::arb B(0.5);
		arbpp::arb C;
		arbpp::arb D;
		
		// u = integral(0..phi) 1 / sqrt(1 - k^2 *sin^2(t)) dt
			// jacobi_sn(k, u) = sn(u, k) := sin(phi)
			// jacobi_cn(k, u) = cn(u, k) := cos(phi)
			// jacobi_dn(k, u) = dn(u, k) := sqrt(1 - k^2 *sin^2(t))
		// jacobi_elliptic(k, u, &C, &D) = sn(u, k)
			// C = cn(u, k), D = dn(u, k)
		// jacobi_[A][B] = jacobi_[A] / jacobi_[B]
			// pl. jacobi_cd = jacobi_cn / jacobi_dn
		// jacobi_n[A] = 1 / jacobi_[A]
			// pl. jacobi_nd = 1 / jacobi_dn
		
		std::cout << "jacobi_elliptic(" << A << ", " << B << ", " << C << ", " << D << ") : " << boost::math::jacobi_elliptic(A, B, &C, &D) << std::endl;
		std::cout << "\tC: " << C << std::endl << "\tD: " << D << std::endl;
		std::cout << "jacobi_cd(" << A << ", " << B << ") : " << boost::math::jacobi_cd(A, B) << std::endl;
		std::cout << "jacobi_nd(" << A << ", " << B << ") : " << boost::math::jacobi_nd(A, B) << std::endl;
		std::cout << std::endl;
	}
	
	void zeta()
	{
		arbpp::arb A(4);
		
		// zeta(s) = sum(k=1..+inf) (1 / (k^s))
		
		//std::cout << "zeta(" << A << ") : " << boost::math::zeta(A) << std::endl;	// TODO numeric_limits
		std::cout << std::endl;
	}
	
	void exponential_integral()
	{
		arbpp::arb A(1.5);
		int n = 2;
		
		// expint(n, x) = integral(1..+inf) e^(-x * t) / t^n dt
		// expint(x) = -expint(1, -x)
		
		//std::cout << "expint(" << n << ", " << A << ") : " << boost::math::expint(n, A) << std::endl;	// TODO segfault, konstans hack nem eleg
		//std::cout << "expint(" << A << ") : " << boost::math::expint(A) << std::endl;	// TODO segfault, konstans hack nem eleg
		std::cout << std::endl;
	}
	
	void basics()
	{
		arbpp::arb A(0.5);
		arbpp::arb B(1);
		
		// sin|cos (pi * x), ln(1+x), exp(x) - 1, kobgyok, sqrt(1+x) - 1, x^y - 1, sqrt(x^2 + y^2)
		
		//std::cout << "sin_pi(" << A << ")" << boost::math::sin_pi(A) << std::endl;	// TODO static_cast int
		//std::cout << "cos_pi(" << A << ")" << boost::math::cos_pi(A) << std::endl;	// TODO static_cast int
		std::cout << "log1p(" << A << ")" << boost::math::log1p(A) << std::endl;
		//std::cout << "expm1(" << A << ")" << boost::math::expm1(A) << std::endl;	// TODO numeric_limits
		//std::cout << "cbrt(" << A << ")" << boost::math::cbrt(A) << std::endl;	// TODO frexp(arb, arb) not implemented	-	arb doesn't have the required function
		//std::cout << "sqrt1pm1(" << A << ")" << boost::math::sqrt1pm1(A) << std::endl;	// TODO numeric_limits
		//std::cout << "powm1(" << A << ", " << B << ")" << boost::math::powm1(A, B) << std::endl;	// TODO numeric_limits
		std::cout << "hypot(" << A << ", " << B << ")" << boost::math::hypot(A, B) << std::endl;
		std::cout << std::endl;
	}
	
	void sinc_sinhc()
	{
		arbpp::arb A(0);
		
		// sinc_pi(x) = sin(x) / x
		// sinhc_pi(x) = sinh(x) / x
		std::cout << "sinc_pi(" << A << ") : " << boost::math::sinc_pi(A) << std::endl;
		std::cout << "sinhc_pi(" << A << ") : " << boost::math::sinhc_pi(A) << std::endl;
		std::cout << std::endl;
	}
	
	void inverse_hyperbolic()
	{
		arbpp::arb A(2);
		arbpp::arb B(-0.6);
		
		//std::cout << "asinh(" << A << ") : " << boost::math::asinh(A) << std::endl;	// TODO numeric_limits
		//std::cout << "acosh(" << A << ") : " << boost::math::acosh(A) << std::endl;	// TODO segfault
		std::cout << "atanh(" << B << ") : " << boost::math::atanh(B) << std::endl;
		std::cout << std::endl;
	}
	
	void owenst()
	{
		arbpp::arb A(2);
		arbpp::arb B(-0.5);
		
		// owens_t(h, a) = 1 / (2*pi) * integral(0..a) exp(-1/2 * h^2 * (1 + x^2)) / (1 + x^2) dx
		
		//std::cout << "owens_t(" << A << ", " << B << ") : " << boost::math::owens_t(A, B) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
}

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/bernoulli.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/exponential.hpp>
#include <boost/math/distributions/extreme_value.hpp>	// a dokumentacioban extreme.hpp neven van, de az nem letezik
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/gamma.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/inverse_gamma.hpp>
#include <boost/math/distributions/inverse_gaussian.hpp>
#include <boost/math/distributions/laplace.hpp>
#include <boost/math/distributions/logistic.hpp>
#include <boost/math/distributions/lognormal.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/non_central_beta.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/non_central_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <boost/math/distributions/pareto.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/rayleigh.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/weibull.hpp>
namespace distributions
{
	// itt a kiirasok tukrozni fogjak az eloszlas paramereterit, de a fvhivasok nem egyeznek meg azzal, amit kiirunk
	// az elso peldat leszamitva csak az eloszlasfv, a surusegfv, a varhato ertek, a szoras es a kvantilis lesz tesztelve (illetve amelyik ertelmes az adott eloszlasra)
	void normal()
	{
		boost::math::normal_distribution<arbpp::arb> A(1.5, 2);
		arbpp::arb x(1.2);
		
		// normal / Gauss-eloszlas
		
		// cdf - eloszlasfv
		// pdf - surusegfv
		// range - ertekkeszlet
		// support - tarto
		// variance - szorasnegyzet
		// skewness - ferdeseg
		// kurtosis - lapultsag
			// kurtosis_excess = kurtosis - 3
		
		//std::cout << "cdf(normal_distribution<arbpp::arb>(1.5, 2), " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(normal_distribution<arbpp::arb>(1.5, 2), " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "mean(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::mean(A) << std::endl;
		std::cout << "median(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::median(A) << std::endl;
		std::cout << "mode(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::mode(A) << std::endl;
		std::pair<arbpp::arb, arbpp::arb> range = boost::math::range(A);
		std::cout << "range(normal_distribution<arbpp::arb>(1.5, 2)) : " << range.first << " : " << range.second << std::endl;
		//std::cout << "quantile(normal_distribution<arbpp::arb>(1.5, 2), 0.4) : " << boost::math::quantile(A, 0.4) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "standard_deviation(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::standard_deviation(A) << std::endl;
		std::pair<arbpp::arb, arbpp::arb> support = boost::math::support(A);
		std::cout << "support(normal_distribution<arbpp::arb>(1.5, 2)) : " << support.first << " : " << support.second << std::endl;
		std::cout << "variance(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::variance(A) << std::endl;
		std::cout << "skewness(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::skewness(A) << std::endl;
		std::cout << "kurtosis(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::kurtosis(A) << std::endl;
		std::cout << "kurtosis_excess(normal_distribution<arbpp::arb>(1.5, 2)) : " << boost::math::kurtosis_excess(A) << std::endl;
		std::cout << std::endl;
	}
	
	void bernoulli()
	{
		arbpp::arb val(0.4);
		boost::math::bernoulli_distribution<arbpp::arb> A(val);
		
		// Bernoulli-eloszlas: egy kiserlet, p valoszinuseggel sikeres
		
		std::cout << "cdf(bernoulli_distribution<arbpp::arb>(" << val << ", 0) : " << boost::math::cdf(A, 0) << std::endl;
		std::cout << "pdf(bernoulli_distribution<arbpp::arb>(" << val << ", 1) : " << boost::math::pdf(A, 1) << std::endl;
		std::cout << "mean(bernoulli_distribution<arbpp::arb>(" << val << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(bernoulli_distribution<arbpp::arb>(" << val << ") : " << boost::math::standard_deviation(A) << std::endl;
		std::cout << "quantile(standard_deviation<arbpp::arb>(" << val << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void beta_distr()
	{
		arbpp::arb a1(2);
		arbpp::arb a2(3);
		boost::math::beta_distribution<arbpp::arb> A(a1, a2);	// mindketto pozitiv, kulonben domain_error
		arbpp::arb x(0.4);
		
		// surusegfv: f(x, a1, a2) = x^(a1-1) * (1 - x)^(a2-1) / B(a1, a2)
		
		//std::cout << "cdf(beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void binomial_distr()
	{
		arbpp::arb n(6);	// n >= 0
		arbpp::arb p(0.3);	// 0 <= p <= 1
		boost::math::binomial_distribution<arbpp::arb> A(n, p);
		arbpp::arb x(4);
		
		//std::cout << "cdf(binomial_distribution<arbpp::arb>(" << n << ", " << p << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(binomial_distribution<arbpp::arb>(" << n << ", " << p << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(binomial_distribution<arbpp::arb>(" << n << ", " << p << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(binomial_distribution<arbpp::arb>(" << n << ", " << p << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(binomial_distribution<arbpp::arb>(" << n << ", " << p << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void cauchy_lorentz()
	{
		arbpp::arb a1(1);
		arbpp::arb a2(1.4);
		boost::math::cauchy_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(4);
		
		// surusegfv: f(x, a1, a2) = 1/pi * (a2 / (x - a1)^2 + a2^2)
			// NINCSEN: mean, std, variance, "stb"
			// ha ilyet kerdezunk le, BOOST_STATIC_ASSERTION_FAILURE forditaskor
		
		//std::cout << "cdf(cauchy_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "pdf(cauchy_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "quantile(cauchy_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << std::endl;
	}
	
	void chi_sqr()
	{
		arbpp::arb v(4);	// v pozitiv, kulonben domain_error
		boost::math::chi_squared_distribution<arbpp::arb> A(v);
		arbpp::arb x(2.7);
		
		//std::cout << "cdf(chi_squared_distribution<arbpp::arb>(" << v << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(chi_squared_distribution<arbpp::arb>(" << v << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "mean(chi_squared_distribution<arbpp::arb>(" << v << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(chi_squared_distribution<arbpp::arb>(" << v << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(chi_squared_distribution<arbpp::arb>(" << v << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void exponential()
	{
		arbpp::arb lambda(1.5);	// lambda pozitiv, kulonben domain_error
		boost::math::exponential_distribution<arbpp::arb> A(lambda);
		arbpp::arb x(2.2);
		
		// surusegfv: f(x) = lambda * e^(-lambda * x)
		
		//std::cout << "cdf(exponential_distribution<arbpp::arb>(" << lambda << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "pdf(exponential_distribution<arbpp::arb>(" << lambda << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(exponential_distribution<arbpp::arb>(" << lambda << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(exponential_distribution<arbpp::arb>(" << lambda << ") : " << boost::math::standard_deviation(A) << std::endl;
		std::cout << "quantile(exponential_distribution<arbpp::arb>(" << lambda << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void extreme()
	{
		arbpp::arb a1(2);
		arbpp::arb a2(1.3);	// a2 pozitiv, kulonben domain_error
		boost::math::extreme_value_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(2.2);
		
		// f(x) = 1/a2 * e^(-(x - a1) / a2) * e^(-e^((x - a1) / a2))
		
		std::cout << "cdf(extreme_value_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;
		std::cout << "pdf(extreme_value_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		//std::cout << "mean(extreme_value_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "standard_deviation(extreme_value_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "quantile(extreme_value_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void fisher_f()
	{
		arbpp::arb a1(4);	// pozitiv
		arbpp::arb a2(10);	// pozitiv
		boost::math::fisher_f_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1.2);
		
		// surusegfv: f(x, n, m) = m^(m/2) * n^(n/2) * x^(n/2 - 1) / ((m + n * x)^((n + m) / 2) * B(n / 2, m / 2))
		
		//std::cout << "cdf(fisher_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(fisher_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(fisher_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(fisher_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(fisher_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void gamma_distr()
	{
		arbpp::arb a1(4);	// pozitiv
		arbpp::arb a2(1.5);	// pozitiv
		boost::math::gamma_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(3);
		
		// surusegfv: f(x, a1, a2) = x^(a1 - 1) * e^(-x/a2) / (a2^a1 * gamma(a1))
		
		//std::cout << "cdf(gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "mean(gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void geometric()
	{
		arbpp::arb p(0.4);
		boost::math::geometric_distribution<arbpp::arb> A(p);
		arbpp::arb x(4);
		
		//std::cout << "cdf(geometric_distribution<arbpp::arb>(" << p << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "pdf(geometric_distribution<arbpp::arb>(" << p << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(geometric_distribution<arbpp::arb>(" << p << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(geometric_distribution<arbpp::arb>(" << p << ") : " << boost::math::standard_deviation(A) << std::endl;
		std::cout << "quantile(geometric_distribution<arbpp::arb>(" << p << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void hypergeometric()
	{
		// ez az eloszlas a template parameteret a publikus value_type tipusan kivul sehol sem hasznalja
		unsigned int r(20);
		unsigned int n(10);
		unsigned int N(200);
		boost::math::hypergeometric_distribution<arbpp::arb> A(r, n, N);
		arbpp::arb x(5);
		
		//std::cout << "cdf(hypergeometric_distribution<arbpp::arb>(" << r << ", " << n << ", " << N << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(hypergeometric_distribution<arbpp::arb>(" << r << ", " << n << ", " << N << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(hypergeometric_distribution<arbpp::arb>(" << r << ", " << n << ", " << N << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(hypergeometric_distribution<arbpp::arb>(" << r << ", " << n << ", " << N << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(hypergeometric_distribution<arbpp::arb>(" << r << ", " << n << ", " << N << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << std::endl;
	}
	
	void inverse_chi_sq()
	{
		arbpp::arb a1(5);	// pozitiv, bizonyos kuszobok alatt egyes jellemzok (pl atlag) nem leteznek
		arbpp::arb a2(0.8);
		boost::math::inverse_chi_squared_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(5);
		
		//std::cout << "cdf(inverse_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(inverse_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(inverse_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(inverse_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(inverse_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void inverse_gamma()
	{
		arbpp::arb a1(4);	// pozitiv
		arbpp::arb a2(1.5);	// pozitiv
		boost::math::inverse_gamma_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(3);
		
		// surusegfv: f(x, a1, a2) = a2^a1 * (1/x)^(a1 + 1) * e^(-a2/x) / gamma(a1)
		
		//std::cout << "cdf(inverse_gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(inverse_gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "mean(inverse_gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(inverse_gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(inverse_gamma_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void inverse_normal()
	{
		arbpp::arb a1(3);
		arbpp::arb a2(5);
		boost::math::inverse_gaussian_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1);
		
		// surusegfv: f(x, a1, a2) = sqrt(a2 / (2 * pi * x^3)) * e^(-a2 * (x - a1)^2 / (2 * a1^2 *x))
		
		//std::cout << "cdf(inverse_gaussian_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(inverse_gaussian_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "mean(inverse_gaussian_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(inverse_gaussian_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(inverse_gaussian_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO static_cast int
		std::cout << std::endl;
	}
	
	void laplace()
	{
		arbpp::arb a1(3);
		arbpp::arb a2(5);	// pozitiv, kulonben domain_error
		boost::math::laplace_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1);
		
		// surusegfv: f(x, a1, a2) = 1 / (2 * a2) * e^(-abs(x - a1) / a2)
		
		std::cout << "cdf(laplace_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;
		std::cout << "pdf(laplace_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(laplace_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		//std::cout << "standard_deviation(laplace_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "quantile(laplace_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void logistic()
	{
		arbpp::arb a1(3);
		arbpp::arb a2(5);	// pozitiv, kulonben domain_error
		boost::math::logistic_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1);
		
		//std::cout << "cdf(logistic_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits
		//std::cout << "pdf(logistic_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "mean(logistic_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		//std::cout << "standard_deviation(logistic_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "quantile(logistic_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void lognormal()
	{
		arbpp::arb a1(3);	// pozitiv
		arbpp::arb a2(0.4);	// 0 <= a2 <= 1
		boost::math::lognormal_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1);	// diszkret elo
		
		// surusegfv: f(x, a1, a2) = gamma(a1 + x) / (x! * gamma(a1)) * a2^a1 * (1 - a2)^x
		
		//std::cout << "cdf(lognormal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(lognormal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "mean(lognormal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		//std::cout << "standard_deviation(lognormal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO numeric_limits
		//std::cout << "quantile(lognormal_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void negbinom()
	{
		arbpp::arb a1(3);	// pozitiv
		arbpp::arb a2(0.4);	// 0 <= a2 <= 1
		boost::math::negative_binomial_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(1);	// diszkret elo
		
		// surusegfv: f(x, a1, a2) = gamma(a1 + x) / (x! * gamma(a1)) * a2^a1 * (1 - a2)^x
		
		//std::cout << "cdf(negative_binomial_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(negative_binomial_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(negative_binomial_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(negative_binomial_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO numeric_limits
		//std::cout << "quantile(negative_binomial_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void noncentral_beta()
	{
		arbpp::arb a1(8);	// pozitiv
		arbpp::arb a2(12);	// pozitiv
		arbpp::arb a3(30);	// nemnegativ
		boost::math::non_central_beta_distribution<arbpp::arb> A(a1, a2, a3);
		arbpp::arb x(0.8);
		
		//std::cout << "cdf(non_central_beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(non_central_beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(non_central_beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(non_central_beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(non_central_beta_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void noncentral_chi_sq()
	{
		arbpp::arb a1(4);	// pozitiv
		arbpp::arb a2(2.1);	// nemnegativ
		boost::math::non_central_chi_squared_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(6.3);
		
		//std::cout << "cdf(non_central_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(non_central_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(non_central_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(non_central_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(non_central_chi_squared_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void noncentral_f()
	{
		arbpp::arb a1(8);	// pozitiv
		arbpp::arb a2(12);	// pozitiv
		arbpp::arb a3(30);	// nemnegativ
		boost::math::non_central_f_distribution<arbpp::arb> A(a1, a2, a3);
		arbpp::arb x(0.8);
		
		//std::cout << "cdf(non_central_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(non_central_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "mean(non_central_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(non_central_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(non_central_f_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits
		std::cout << std::endl;
	}
	
	void noncentral_t()
	{
		arbpp::arb a1(4);	// pozitiv
		arbpp::arb a2(2.1);	// veges
		boost::math::non_central_t_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(6.3);
		
		//std::cout << "cdf(non_central_t_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(non_central_t_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "mean(non_central_t_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "standard_deviation(non_central_t_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "quantile(non_central_t_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void pareto()
	{
		arbpp::arb a1(1.5);	// pozitiv
		arbpp::arb a2(4);	// pozitiv
		boost::math::pareto_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(2.3);
		
		// surusegfv: f(x, a2, a1) = a2 * a1^a2 / x^(a1 + 1)
		
		//std::cout << "cdf(pareto_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "pdf(pareto_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(pareto_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "standard_deviation(pareto_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << "quantile(pareto_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void poisson()
	{
		arbpp::arb a1(4);
		boost::math::poisson_distribution<arbpp::arb> A(a1);
		arbpp::arb x(3);
		
		// surusegfv: f(k, a1) = e^(-a1) * a1^k / k!
		
		//std::cout << "cdf(poisson_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(poisson_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "mean(poisson_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(poisson_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(poisson_distribution<arbpp::arb>(" << a1 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void rayleigh()
	{
		arbpp::arb a1(3.6);	// pozitiv
		boost::math::rayleigh_distribution<arbpp::arb> A(a1);
		arbpp::arb x(3.4);
		
		// surusegfv: f(x, a1) = x * e^(-x^2 / a1^2) / a1^2
		
		//std::cout << "cdf(rayleigh_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "pdf(rayleigh_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		//std::cout << "mean(rayleigh_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::mean(A) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "standard_deviation(rayleigh_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO segfault, konstans hack megoldja
		std::cout << "quantile(rayleigh_distribution<arbpp::arb>(" << a1 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void skew_normal()
	{
		arbpp::arb a1(1);
		arbpp::arb a2(1.2);	// pozitiv
		arbpp::arb a3(4);
		boost::math::skew_normal_distribution<arbpp::arb> A(a1, a2, a3);
		arbpp::arb x(0.8);
		
		//std::cout << "cdf(skew_normal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(skew_normal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "mean(skew_normal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::mean(A) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "standard_deviation(skew_normal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO segfault, konstans hack megoldja
		//std::cout << "quantile(skew_normal_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void students_t()
	{
		arbpp::arb a1(4);	// pozitiv
		boost::math::students_t_distribution<arbpp::arb> A(a1);
		arbpp::arb x(1.6);
		
		// surusegfv: f(x, a1) = gamma((a1 + 1) / 2) / sqrt(a1 * pi) * gamma(a1/2) * (1 + x^2/a1)^((a1 + 1) / 2)
		
		//std::cout << "cdf(students_t_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(students_t_distribution<arbpp::arb>(" << a1 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO static_cast int
		std::cout << "mean(students_t_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(students_t_distribution<arbpp::arb>(" << a1 << ") : " << boost::math::standard_deviation(A) << std::endl;
		//std::cout << "quantile(students_t_distribution<arbpp::arb>(" << a1 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;	// TODO numeric_limits & static_cast int
		std::cout << std::endl;
	}
	
	void triangular()
	{
		arbpp::arb a1(-1);	// veges
		arbpp::arb a2(0.8);	// veges
		arbpp::arb a3(4);	// veges
		boost::math::triangular_distribution<arbpp::arb> A(a1, a2, a3);
		arbpp::arb x(0.8);
		
		std::cout << "cdf(triangular_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;
		std::cout << "pdf(triangular_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(triangular_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(triangular_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ") : " << boost::math::standard_deviation(A) << std::endl;
		std::cout << "quantile(triangular_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << a3 << ", 0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void uniform()
	{
		arbpp::arb a1(-1.5);	// veges
		arbpp::arb a2(4);	// veges
		boost::math::uniform_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(2.3);
		
		std::cout << "cdf(uniform_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;
		std::cout << "pdf(uniform_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		std::cout << "mean(uniform_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;
		std::cout << "standard_deviation(uniform_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;
		std::cout << "quantile(uniform_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
	
	void weibull()
	{
		arbpp::arb a1(3);	// pozitiv
		arbpp::arb a2(2);	// pozitiv
		boost::math::weibull_distribution<arbpp::arb> A(a1, a2);
		arbpp::arb x(2.3);
		
		//std::cout << "cdf(weibull_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits
		std::cout << "pdf(weibull_distribution<arbpp::arb>(" << a1 << ", " << a2 << ", " << x << ") : " << boost::math::pdf(A, x) << std::endl;
		//std::cout << "mean(weibull_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::mean(A) << std::endl;	// TODO static_cast int
		//std::cout << "standard_deviation(weibull_distribution<arbpp::arb>(" << a1 << ", " << a2 << ") : " << boost::math::standard_deviation(A) << std::endl;	// TODO static_cast int
		std::cout << "quantile(weibull_distribution<arbpp::arb>(" << a1 << ", " << a2 << "0.5) : " << boost::math::quantile(A, 0.5) << std::endl;
		std::cout << std::endl;
	}
}

int main()
{
	//arbpp_demo_boostround();
	//arbpp_demo_boosttrunc();
	//arbpp_demo_modf();
	//arbpp_demo_fpclass();	// unfinished, compiles & runs
	//arbpp_demo_sign();
	//arbpp_demo_io();
	//arbpp_demo_const();	// unfinished, "hacked"
	//arbpp_demo_complex();	// partially unfinished
	//arbpp_demo_gcd_lcm();	// unfinished, probably unneeded
	
	//deriv_roots::arbpp_demo_roots();
	//noderiv_roots::arbpp_demo_roots();
	//minimize_func::arbpp_demo_minimize();
	
	// unfinished - egyik fuggveny sem mukodik (forditasideju hiba vagy segfault)
	// partially unfinished - legalabb egy fuggveny nem mukodik
	
	//special_functions::gamma();	// unfinished
	//special_functions::factorial();	// unfinished
	//special_functions::legendre();	// partially unfinished
	//special_functions::laguerre();
	//special_functions::hermite();
	//special_functions::spherical_harmonic();	// unfinished
	//special_functions::bessel();	// unfinished
	//special_functions::hankel();	// unfinished
	//special_functions::airy();	// unfinished
	//special_functions::carlson_elliptic();
	//special_functions::legendre_elliptic();	// partially unfinished
	//special_functions::jacobi_elliptic();
	//special_functions::zeta();	// unfinished
	//special_functions::exponential_integral();	// unfinished
	//special_functions::basics();	// partially unfinished
	//special_functions::sinc_sinhc();
	//special_functions::inverse_hyperbolic();	// partially unfinished
	//special_functions::owenst();	// unfinished

	//distributions::normal();	// partially unfinished
	//distributions::bernoulli();
	//distributions::beta_distr();	// partially unfinished
	//distributions::binomial_distr();	// partially unfinished
	//distributions::cauchy_lorentz();	// unfinished, hackable
	//distributions::chi_sqr();	// partially unfinished
	//distributions::exponential();	// partially unfinished
	//distributions::extreme();	// partially unfinished
	//distributions::fisher_f();	// partially unfinished
	//distributions::gamma_distr();	// partially unfinished
	//distributions::geometric();	// partially unfinished
	//distributions::hypergeometric();	// partially unfinished, partially hackable
	//distributions::inverse_chi_sq();	// partially unfinished
	//distributions::inverse_gamma();	// partially unfinished
	//distributions::inverse_normal();	// partially unfinished, partially hackable
	//distributions::laplace();	// partially unfinished, hackable
	//distributions::logistic();	// partially unfinished, partially hackable
	//distributions::lognormal();	// partiall unfinished, partially hackable
	//distributions::negbinom();	// partially unfinished
	//distributions::noncentral_beta();	// partially unfinished
	//distributions::noncentral_f();	// partially unfinished
	//distributions::noncentral_t();	// unfinished
	//distributions::pareto();	// partially unfinished
	//distributions::poisson();	// partially unfinished
	//distributions::rayleigh();	// partially unfinished, partially hackable
	//distributions::skew_normal();	// unfinished, partially hackable
	//distributions::students_t();	// partially unfinished
	//distributions::triangular();
	//distributions::uniform();
	//distributions::weibull();	// partially unfinished
	
	
	return 0;
}
