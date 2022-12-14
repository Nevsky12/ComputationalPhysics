#include <gtest/gtest.h>
#include <interpolators/newton.h>
#include <interpolators/cubic_spline.h>
#include <iostream>
#include <fstream>

// 0.326019 - Newton in KR
// 0.327426 - Natural spline in KR

TEST(INTERPOLATION, NEWTON)
{
    ComputationalPhysics::Types::BasicTypes::vecSTL f_x = { 92228496. , 106021537., 123202624., 132164569., 151325798.
                                                          , 179323175., 203211926., 226545805., 248709873., 281421906.
                                                          };
    ComputationalPhysics::Types::BasicTypes::vecSTL   x = { 1910., 1920., 1930., 1940., 1950.
                                                          , 1960., 1970., 1980., 1990., 2000.
                                                          };

    ComputationalPhysics::Types::BasicTypes::scalar const x0 = 2010.;

    auto const result = ComputationalPhysics::Interpolators::newton(x, f_x, x0);
    std::cout << "Newton interpolation: " << result << std::endl;
}

TEST(INTERPOLATION, CUBIC_SPLINE)
{
    ComputationalPhysics::Types::BasicTypes::vecSTL f_x = { 92228496. , 106021537., 123202624., 132164569., 151325798.
                                                          , 179323175., 203211926., 226545805., 248709873., 281421906.
                                                          };
    ComputationalPhysics::Types::BasicTypes::vecSTL   x = { 1910., 1920., 1930., 1940., 1950.
                                                          , 1960., 1970., 1980., 1990., 2000.
                                                          };

    ComputationalPhysics::Types::BasicTypes::scalar const h = 10.;
    ComputationalPhysics::Types::BasicTypes::scalar const x0 = 1990.;

    auto const resSpline = ComputationalPhysics::Interpolators::naturalSpline(x, f_x, h);
    auto const result        = ComputationalPhysics::Interpolators::naturalSplineValue(resSpline, x0);
    std::cout << "Cubic spline interpolation: " << result << std::endl;
}