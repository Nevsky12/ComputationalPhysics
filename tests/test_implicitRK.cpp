#include <gtest/gtest.h>
#include <integrators/robust_tasks/implicitRK.h>
#include <iostream>
#include <fstream>

TEST(INTEGRATORS, IMPLICIT_RUNGE_KUTTA)
{
    ComputationalPhysics::Types::BasicTypes::scalar const mu = 1000.;
    ComputationalPhysics::Types::BasicTypes::scalar const h = 0.001;
    ComputationalPhysics::Types::BasicTypes::scalar const T = 1000.;
    ComputationalPhysics::Types::BasicTypes::scalar const t0 = 0.;
    ComputationalPhysics::Types::BasicTypes::scalar const p0 = 0.001;
    ComputationalPhysics::Types::BasicTypes::scalar const x0 = 0.;

    ComputationalPhysics::Types::BasicTypes::vec px {{p0, x0}};
    ComputationalPhysics::Types::CoreTypes::State const init =
    {
        .state = px,
        .t = t0,
    };

    // Implicit Euler
//    ComputationalPhysics::Types::BasicTypes::vec const b {{1.}};
//    ComputationalPhysics::Types::BasicTypes::vec const c {{1.}};
//    ComputationalPhysics::Types::BasicTypes::mat const A {{1.}};

    // Rado method
//    ComputationalPhysics::Types::BasicTypes::vec const c { {(4. - std::sqrt(6.)) / 10.
//                                                          , (4 + std::sqrt(6)) / 10.
//                                                          , 1.}
//                                                         };
//    ComputationalPhysics::Types::BasicTypes::vec const b { {(16. - std::sqrt(6)) / 36.
//                                                          , (16. + std::sqrt(6)) / 36.
//                                                          , 1. / 9.}
//                                                         };
//    ComputationalPhysics::Types::BasicTypes::mat const A { {(88. - 7. * std::sqrt(6.)) / 360., (296. - 169. * std::sqrt(6)) / 1800., (-2. + 3. * std::sqrt(6)) / 255.}
//                                                         , {(296. + 169. * std::sqrt(6)) / 1800., (88. + 7. * std::sqrt(6.)) / 360., (-2. - 3. * std::sqrt(6)) / 255.}
//                                                         , {(16. - std::sqrt(6)) / 36., (16. + std::sqrt(6)) / 36., 1. / 9.}
//                                                         };

    // Gauss method
    ComputationalPhysics::Types::BasicTypes::vec const c { {0.5 - std::sqrt(15) / 10.
                                                          , 0.5
                                                          , 0.5 + std::sqrt(15) / 10.}
                                                         };
    ComputationalPhysics::Types::BasicTypes::vec const b { {5. / 18.
                                                          , 4. / 9.
                                                          , 5. / 18.}
                                                         };
    ComputationalPhysics::Types::BasicTypes::mat const A { {5. / 36., 2. / 9. - std::sqrt(15.) / 15., 5. / 36. - std::sqrt(15.) / 30.}
                                                         , {5. / 36. + std::sqrt(15.) / 24., 2. / 9., 5. / 36. - std::sqrt(15.) / 24.}
                                                         , {5. / 36. + std::sqrt(15.) / 30., 2. / 9. + std::sqrt(15.) / 15., 5. / 36.}
                                                         };
    ComputationalPhysics::Types::CoreTypes::ButcherTable<3u> const &butcherTable =
    {
           .b = b,
           .c = c,
           .A = A,
    };

    auto const &rightPart = [=]( ComputationalPhysics::Types::BasicTypes::scalar const  t
                                                                 , ComputationalPhysics::Types::BasicTypes::   vec const &st
                                                                 ) noexcept
                                                                -> ComputationalPhysics::Types::BasicTypes::vec
    {
        ComputationalPhysics::Types::BasicTypes::scalar const p = st.x();
        ComputationalPhysics::Types::BasicTypes::scalar const x = st.y();
        return ComputationalPhysics::Types::BasicTypes::vec
        {{
            p * mu * (1 - p * p) * p - x,
            p,
        }};
    };
    auto const &jacobi = [=]( ComputationalPhysics::Types::BasicTypes::scalar const t
                                                              , ComputationalPhysics::Types::BasicTypes::   vec const &st
                                                              ) noexcept
                                                            -> ComputationalPhysics::Types::BasicTypes::mat
    {
        ComputationalPhysics::Types::BasicTypes::scalar const p = st.x();
        return  ComputationalPhysics::Types::BasicTypes::mat
       {
               {mu * (1. - 3 * p * p),  -1.},
               {1.,                      0.},
       };
    };

    std::ofstream out("implicitRK.csv");

    out << "Time" << '\t';
    out << "Coord" << '\t';
    out << "Velocity" << '\n';
    auto const &result = ComputationalPhysics::Integrators::Implicit::implicitRKFinite<3u>( init
                                                                                                          , h
                                                                                                          , T
                                                                                                          , butcherTable
                                                                                                          , rightPart
                                                                                                          , jacobi
                                                                                                          );
    for (auto const stuff: result)
    {
        out << stuff.t << '\t';
        out << stuff.state.y() << '\t';
        out << stuff.state.x() << '\n';
    }
    out.close();
}