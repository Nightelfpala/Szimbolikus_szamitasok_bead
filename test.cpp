
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>
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

void arbpp_test_boostmathtools()
{
	std::cout << "digits():\t" << boost::math::tools::digits<arbpp::arb>() << std::endl;
	std::cout << "min_value():\t" << boost::math::tools::min_value<arbpp::arb>() << std::endl;
	std::cout << "max_value():\t" << boost::math::tools::max_value<arbpp::arb>() << std::endl;
}

void arbpp_add_test()
{
	//arbpp_test_eq();
	//arbpp_test_ineq();
	//arbpp_test_le();
	//arbpp_test_ge();
	//arbpp_test_lt();
	//arbpp_test_gt();
	//arbpp_test_abs();
	//arbpp_test_ceil();
	//arbpp_test_floor()
	//arbpp_test_exp();
	//arbpp_test_pow();
	//arbpp_test_sqrt();
	//arbpp_test_log();
	//arbpp_test_sin();
	//arbpp_test_asin();
	//arbpp_test_acos();
	//arbpp_test_tan();
	//arbpp_test_atan();
	arbpp_test_boostmathtools();
}

int main()
{
	arbpp_add_test();
	
	
	return 0;
}