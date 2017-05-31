
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
		
		//std::cout << "expint(" << n << ", " << A << ") : " << boost::math::expint(n, A) << std::endl;	// TODO segfault
		//std::cout << "expint(" << A << ") : " << boost::math::expint(A) << std::endl;	// TODO segfault
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
namespace distributions
{
	// itt a kiirasok tukrozni fogjak az eloszlas paramereterit, de a fvhivasok nem egyeznek meg azzal, amit kiirunk
	void normal()
	{
		boost::math::normal_distribution<arbpp::arb> A(1.5, 2);	// normal / Gauss-eloszlas
		arbpp::arb x(1.2);
		
		// cdf - eloszlasfv
		// pdf - surusegfv
		// range - ertekkeszlet
		// support - tarto
		// variance - szorasnegyzet
		// skewness - ferdeseg
		// kurtosis - lapultsag
			// kurtosis_excess = kurtosis - 3
		
		//std::cout << "cdf(normal_distribution<arbpp::arb>(1.5, 2), " << x << ") : " << boost::math::cdf(A, x) << std::endl;	// TODO numeric_limits & static_cast int
		//std::cout << "pdf(normal_distribution<arbpp::arb>(1.5, 2), " << x << ") : " << boost::math::pdf(A, x) << std::endl;	// TODO segfault
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
		boost::math::bernoulli_distribution<arbpp::arb> A(val);	// Bernoulli-eloszlas: egy kiserlet, p valoszinuseggel sikeres
		
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
}

int main()
{
	//arbpp_demo_boostround();
	//arbpp_demo_boosttrunc();
	//arbpp_demo_modf();
	//arbpp_demo_fpclass();	// unfinished, compiles & runs
	//arbpp_demo_sign();
	//arbpp_demo_io();
	//arbpp_demo_const();	// unfinished
	//arbpp_demo_complex();	// partially unfinished
	//arbpp_demo_gcd_lcm();	// unfinished
	//deriv_roots::arbpp_demo_roots();
	//noderiv_roots::arbpp_demo_roots();
	//minimize_func::arbpp_demo_minimize();
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
	distributions::beta_distr();	// partially unfinished
	
	return 0;
}
