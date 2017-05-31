
#include <iostream>
#include "arbpp.hpp"

// http://www.boost.org/doc/libs/1_63_0/libs/math/doc/html/math_toolkit/real_concepts.html megvalositasok tesztelese
	// https://github.com/bluescarni/arbpp kiegeszitesebol

void arbpp_test_eq()
{
	// operator==
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " == " << B << " : " << (A == B) << std::endl;
	std::cout << A << " == " << C << " : " << (A == C) << std::endl;
	std::cout << C << " == " << A << " : " << (C == A) << std::endl;
	std::cout << A << " == " << 1 << " : " << (A == 1) << std::endl;
	std::cout << 1 << " == " << B << " : " << (1 == B) << std::endl;
	std::cout << A << " == " << 2 << " : " << (A == 2) << std::endl;
	std::cout << 2 << " == " << B << " : " << (2 == B) << std::endl;
	
	std::cout << A << " == " << 1.0f << " : " << (A == 1.0f) << std::endl;
	std::cout << 1.0f << " == " << B << " : " << (1.0f == B) << std::endl;
	std::cout << A << " == " << 2.0f << " : " << (A == 2.0f) << std::endl;
	std::cout << 2.0f << " == " << B << " : " << (2.0f == B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_ineq()
{
	// operator!=
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " != " << B << " : " << (A != B) << std::endl;
	std::cout << A << " != " << C << " : " << (A != C) << std::endl;
	std::cout << C << " != " << A << " : " << (C != A) << std::endl;
	std::cout << A << " != " << 1 << " : " << (A != 1) << std::endl;
	std::cout << 1 << " != " << B << " : " << (1 != B) << std::endl;
	std::cout << A << " != " << 2 << " : " << (A != 2) << std::endl;
	std::cout << 2 << " != " << B << " : " << (2 != B) << std::endl;
	
	std::cout << A << " != " << 1.0f << " : " << (A != 1.0f) << std::endl;
	std::cout << 1.0f << " != " << B << " : " << (1.0f != B) << std::endl;
	std::cout << A << " != " << 2.0f << " : " << (A != 2.0f) << std::endl;
	std::cout << 2.0f << " != " << B << " : " << (2.0f != B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_le()
{
	// operator<=
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " <= " << B << " : " << (A <= B) << std::endl;
	std::cout << A << " <= " << C << " : " << (A <= C) << std::endl;
	std::cout << C << " <= " << A << " : " << (C <= A) << std::endl;
	std::cout << A << " <= " << 1 << " : " << (A <= 1) << std::endl;
	std::cout << 1 << " <= " << B << " : " << (1 <= B) << std::endl;
	std::cout << A << " <= " << 2 << " : " << (A <= 2) << std::endl;
	std::cout << 2 << " <= " << B << " : " << (2 <= B) << std::endl;
	
	std::cout << A << " <= " << 1.0f << " : " << (A <= 1.0f) << std::endl;
	std::cout << 1.0f << " <= " << B << " : " << (1.0f <= B) << std::endl;
	std::cout << A << " <= " << 2.0f << " : " << (A <= 2.0f) << std::endl;
	std::cout << 2.0f << " <= " << B << " : " << (2.0f <= B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_ge()
{
	// operator>=
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " >= " << B << " : " << (A >= B) << std::endl;
	std::cout << A << " >= " << C << " : " << (A >= C) << std::endl;
	std::cout << C << " >= " << A << " : " << (C >= A) << std::endl;
	std::cout << A << " >= " << 1 << " : " << (A >= 1) << std::endl;
	std::cout << 1 << " >= " << B << " : " << (1 >= B) << std::endl;
	std::cout << A << " >= " << 2 << " : " << (A >= 2) << std::endl;
	std::cout << 2 << " >= " << B << " : " << (2 >= B) << std::endl;
	
	std::cout << A << " >= " << 1.0f << " : " << (A >= 1.0f) << std::endl;
	std::cout << 1.0f << " >= " << B << " : " << (1.0f >= B) << std::endl;
	std::cout << A << " >= " << 2.0f << " : " << (A >= 2.0f) << std::endl;
	std::cout << 2.0f << " >= " << B << " : " << (2.0f >= B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_lt()
{
	// operator<
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " < " << B << " : " << (A < B) << std::endl;
	std::cout << A << " < " << C << " : " << (A < C) << std::endl;
	std::cout << C << " < " << A << " : " << (C < A) << std::endl;
	std::cout << A << " < " << 1 << " : " << (A < 1) << std::endl;
	std::cout << 1 << " < " << B << " : " << (1 < B) << std::endl;
	std::cout << A << " < " << 2 << " : " << (A < 2) << std::endl;
	std::cout << 2 << " < " << B << " : " << (2 < B) << std::endl;
	
	std::cout << A << " < " << 1.0f << " : " << (A < 1.0f) << std::endl;
	std::cout << 1.0f << " < " << B << " : " << (1.0f < B) << std::endl;
	std::cout << A << " < " << 2.0f << " : " << (A < 2.0f) << std::endl;
	std::cout << 2.0f << " < " << B << " : " << (2.0f < B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_gt()
{
	// operator>
	arbpp::arb A(1);
	arbpp::arb B(1);
	arbpp::arb C(2);
	
	std::cout << A << " > " << B << " : " << (A > B) << std::endl;
	std::cout << A << " > " << C << " : " << (A > C) << std::endl;
	std::cout << C << " > " << A << " : " << (C > A) << std::endl;
	std::cout << A << " > " << 1 << " : " << (A > 1) << std::endl;
	std::cout << 1 << " > " << B << " : " << (1 > B) << std::endl;
	std::cout << A << " > " << 2 << " : " << (A > 2) << std::endl;
	std::cout << 2 << " > " << B << " : " << (2 > B) << std::endl;
	
	std::cout << A << " > " << 1.0f << " : " << (A > 1.0f) << std::endl;
	std::cout << 1.0f << " > " << B << " : " << (1.0f > B) << std::endl;
	std::cout << A << " > " << 2.0f << " : " << (A > 2.0f) << std::endl;
	std::cout << 2.0f << " > " << B << " : " << (2.0f > B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_abs()
{
	arbpp::arb A(1);
	arbpp::arb B(-1);
	std::cout << "abs(" << A << ") : " << abs(A) << std::endl;
	std::cout << "abs(" << B << ") : " << abs(B) << std::endl;
	
	A = 2.5d;
	B = -2.5d;
	std::cout << "abs(" << A << ") : " << abs(A) << std::endl;
	std::cout << "abs(" << B << ") : " << abs(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_ceil()
{
	arbpp::arb A(1);
	arbpp::arb B(-1);
	std::cout << "ceil(" << A << ") : " << ceil(A) << std::endl;
	std::cout << "ceil(" << B << ") : " << ceil(B) << std::endl;
	
	A = 2.5d;
	B = -2.5d;
	std::cout << "ceil(" << A << ") : " << ceil(A) << std::endl;
	std::cout << "ceil(" << B << ") : " << ceil(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_floor()
{
	arbpp::arb A(1);
	arbpp::arb B(-1);
	std::cout << "floor(" << A << ") : " << floor(A) << std::endl;
	std::cout << "floor(" << B << ") : " << floor(B) << std::endl;
	
	A = 2.5d;
	B = -2.5d;
	std::cout << "floor(" << A << ") : " << floor(A) << std::endl;
	std::cout << "floor(" << B << ") : " << floor(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_exp()
{
	arbpp::arb A(1);
	arbpp::arb B(0);
	arbpp::arb C(-1);
	std::cout << "exp(" << A << ") : " << exp(A) << std::endl;
	std::cout << "exp(" << B << ") : " << exp(B) << std::endl;
	std::cout << "exp(" << C << ") : " << exp(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_pow()
{
	arbpp::arb A(1);
	arbpp::arb B(2);
	arbpp::arb C(3);
	arbpp::arb D(-1);
	std::cout << "pow(" << A << ", " << B << ") : " << pow(A, B) << std::endl;
	std::cout << "pow(" << B << ", " << C << ") : " << pow(B, C) << std::endl;
	std::cout << "pow(" << B << ", " << D << ") : " << pow(B, D) << std::endl;
	
	// nem mukodik, de elvileg csak pow(RealType, RealType) kell, alap tipussal nem kell egyutt mukodnie
	//std::cout << "pow(" << B << ", " << 5 << ") : " << pow(B, 5) << std::endl;
	//std::cout << "pow(" << B << ", " << 5.0d << ") : " << pow(B, 5.0d) << std::endl;
	
	std::cout << std::endl;
}

void arbpp_test_sqrt()
{
	arbpp::arb A(1);
	arbpp::arb B(4);
	arbpp::arb C(10);
	std::cout << "sqrt(" << A << ") : " << sqrt(A) << std::endl;
	std::cout << "sqrt(" << B << ") : " << sqrt(B) << std::endl;
	std::cout << "sqrt(" << C << ") : " << sqrt(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_log()
{
	arbpp::arb A(1);
	arbpp::arb B(exp(arbpp::arb(1)));
	arbpp::arb C(10);
	std::cout << "log(" << A << ") : " << log(A) << std::endl;
	std::cout << "log(" << B << ") : " << log(B) << std::endl;
	std::cout << "log(" << C << ") : " << log(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_sin()	// cos mar meg volt irva
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	arbpp::arb C(asin(arbpp::arb(1)));
	std::cout << "sin(" << A << ") : " << sin(A) << std::endl;
	std::cout << "sin(" << B << ") : " << sin(B) << std::endl;
	std::cout << "sin(" << C << ") : " << sin(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_asin()
{
	arbpp::arb A(0);
	arbpp::arb B(0.5);
	arbpp::arb C(1);
	std::cout << "asin(" << A << ") : " << asin(A) << std::endl;
	std::cout << "asin(" << B << ") : " << asin(B) << std::endl;
	std::cout << "asin(" << C << ") : " << asin(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_acos()	// boost nem irja hogy kene, de azert implementalom
{
	arbpp::arb A(0);
	arbpp::arb B(0.5);
	arbpp::arb C(1);
	std::cout << "acos(" << A << ") : " << acos(A) << std::endl;
	std::cout << "acos(" << B << ") : " << acos(B) << std::endl;
	std::cout << "acos(" << C << ") : " << acos(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_tan()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	std::cout << "tan(" << A << ") : " << tan(A) << std::endl;
	std::cout << "tan(" << B << ") : " << tan(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_atan()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	arbpp::arb C(tan(arbpp::arb(0.5d)));
	std::cout << "atan(" << A << ") : " << atan(A) << std::endl;
	std::cout << "atan(" << B << ") : " << atan(B) << std::endl;
	std::cout << "atan(" << C << ") : " << atan(C) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_atan2()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	std::cout << "atan2(" << A << ", " << B << ") : " << atan2(A, B) << std::endl;
	std::cout << "atan2(" << B << ", " << A << ") : " << atan2(B, A) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_sinh()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	
	std::cout << "sinh(" << A << ") : " << sinh(A) << std::endl;
	std::cout << "sinh(" << B << ") : " << sinh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_cosh()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	
	std::cout << "cosh(" << A << ") : " << cosh(A) << std::endl;
	std::cout << "cosh(" << B << ") : " << cosh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_tanh()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	
	std::cout << "tanh(" << A << ") : " << tanh(A) << std::endl;
	std::cout << "tanh(" << B << ") : " << tanh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_coth()
{
	arbpp::arb A(0);
	arbpp::arb B(1);
	
	std::cout << "coth(" << A << ") : " << coth(A) << std::endl;
	std::cout << "coth(" << B << ") : " << coth(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_asinh()
{
	arbpp::arb A(0.5);
	arbpp::arb B(-4.3);
	
	std::cout << "asinh(" << A << ") : " << asinh(A) << std::endl;
	std::cout << "asinh(" << B << ") : " << asinh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_acosh()
{
	arbpp::arb A(2);
	arbpp::arb B(4.3);
	
	std::cout << "acosh(" << A << ") : " << acosh(A) << std::endl;
	std::cout << "acosh(" << B << ") : " << acosh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_atanh()
{
	arbpp::arb A(0.5);
	arbpp::arb B(-0.8);
	
	std::cout << "atanh(" << A << ") : " << atanh(A) << std::endl;
	std::cout << "atanh(" << B << ") : " << atanh(B) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_ldexp()
{
	// ldexp(x, ex) = x * 2^ex
	arbpp::arb A(2);
	arbpp::arb B(3);
	std::cout << "ldexp(" << A << ", 3) : " << ldexp(A, 3) << std::endl;
	std::cout << "ldexp(" << B << ", -2) : " << ldexp(B, -2) << std::endl;
	std::cout << std::endl;
}

void arbpp_test_boostmathtools()
{
	std::cout << "digits():\t" << boost::math::tools::digits<arbpp::arb>() << std::endl;
	std::cout << "min_value():\t" << boost::math::tools::min_value<arbpp::arb>() << std::endl;
	std::cout << "max_value():\t" << boost::math::tools::max_value<arbpp::arb>() << std::endl;
	std::cout << "epsilon():\t" << boost::math::tools::epsilon<arbpp::arb>() << std::endl;
	std::cout << std::endl;
}

int main()
{
	//arbpp_test_eq();	// op==
	//arbpp_test_ineq();	// op!=
	//arbpp_test_le();	// op<=
	//arbpp_test_ge();	// op>=
	//arbpp_test_lt();	// op<
	//arbpp_test_gt();	// op>
	//arbpp_test_abs();	// abs()
	//arbpp_test_ceil();	// ceil()
	//arbpp_test_floor()	// floor()
	//arbpp_test_exp();	// exp()
	//arbpp_test_pow();	// pow()
	//arbpp_test_sqrt();	// sqrt()
	//arbpp_test_log();	// log()
	//arbpp_test_sin();	// sin()
	//arbpp_test_asin();	// asin()
	//arbpp_test_acos();	// acos()
	//arbpp_test_tan();	// tan()
	//arbpp_test_atan();	// atan()
	//arbpp_test_sinh();	// sinh()
	//arbpp_test_cosh();	// cosh()
	//arbpp_test_tanh();	// tanh()
	//arbpp_test_coth();	// coth()
	//arbpp_test_asinh();
	//arbpp_test_acosh();
	//arbpp_test_atanh();
	//arbpp_test_atan2();	// atan2()
	//arbpp_test_ldexp();	// ldexp()
	//arbpp_test_boostmathtools();
	
	return 0;
}