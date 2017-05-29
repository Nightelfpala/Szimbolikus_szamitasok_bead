
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

void arbpp_add_test()
{
	//arbpp_test_eq();
	//arbpp_test_ineq();
	//arbpp_test_le();
	//arbpp_test_ge();
	//arbpp_test_lt();
	//arbpp_test_gt();
	//arbpp_test_abs();
	arbpp_test_ceil();
	arbpp_test_floor();
}

int main()
{
	arbpp_add_test();
	
	
	return 0;
}