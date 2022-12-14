#include <gtest/gtest.h>
#include <integrators/embedded_methods/embeddedMethods.h>
#include <iostream>
#include <fstream>

TEST(INTEGRATORS, DORMAN_PRINCE)
{
    ComputationalPhysics::Types::BasicTypes::scalar const mu = 0.012277471;
    ComputationalPhysics::Types::BasicTypes::scalar const eta = 1. - mu;
    ComputationalPhysics::Types::BasicTypes::scalar const T = 17.0652165601579625588917206249;
    ComputationalPhysics::Types::BasicTypes::scalar const broad = T;
    ComputationalPhysics::Types::BasicTypes::scalar const tolerance = 1e-12;
    ComputationalPhysics::Types::BasicTypes::scalar const h0 = 1e-4;
    ComputationalPhysics::Types::BasicTypes::scalar const x0 = 0.994;
    ComputationalPhysics::Types::BasicTypes::scalar const y0 = 0.;
    ComputationalPhysics::Types::BasicTypes::scalar const u0 = 0.;
    ComputationalPhysics::Types::BasicTypes::scalar const v0 = -2.00158510637908252240537862224;

    ComputationalPhysics::Types::BasicTypes::vec px {{x0, u0, y0, v0}};
    ComputationalPhysics::Types::CoreTypes::State const init =
            {
                    .state = px,
                    .t = 0.,
            };

    // Dormann-Prince
    ComputationalPhysics::Types::BasicTypes::vec const c { {0.
                                                           , 1. / 5.
                                                           , 3. / 10.
                                                           , 4. / 5.
                                                           , 8. / 9.
                                                           , 1.
                                                           , 1.
                                                           }
                                                         };
    ComputationalPhysics::Types::BasicTypes::vec const b { { 35. / 384.
                                                           , 0.
                                                           , 500. / 1113.
                                                           , 125. / 192.
                                                           , -2187. / 6784.
                                                           , 11. / 84.
                                                           , 0.
                                                           }
                                                         };
    ComputationalPhysics::Types::BasicTypes::vec const bp { { 5179. / 57600.
                                                            , 0.
                                                            , 7571. / 16695.
                                                            , 393. / 640.
                                                            , -92097. / 339200.
                                                            , 187. / 2100.
                                                            , 1. / 40.
                                                            }
    };
    ComputationalPhysics::Types::BasicTypes::mat const A { {0.             , 0.             , 0.             , 0.            , 0., 0.         , 0.}
                                                         , {1. / 5.        , 0.             , 0.             , 0.            , 0., 0.         , 0.}
                                                         , {3. / 40.       , 9. / 40.       , 0.             , 0.            , 0.             , 0., 0.}
                                                         , {44. / 45.      , -56. / 15.     , 32. / 9., 0.   , 0.            , 0.             , 0.}
                                                         , {19372. / 6561. , -25360. / 2187., 64448. / 6561. , -212. / 729.  , 0.             , 0., 0.}
                                                         , {9017. / 3168.  , -355. / 33.    , -46732. / 5247., 49. / 176.    , -5103. / 18656., 0., 0.}
                                                         , { 35. / 384., 0., 500. / 1113.   , 125. / 192.    , -2187. / 6784., 11. / 84.      , 0.}
                                                         };
    ComputationalPhysics::Types::CoreTypes::EmbeddedTable<7u> const &embeddedTable =
            {
                    ComputationalPhysics::Types::CoreTypes::ButcherTable<7u>
                    {
                            .b = b,
                            .c = c,
                            .A = A,
                    },
                    bp,
            };

    auto const &rightPart = [=]( ComputationalPhysics::Types::BasicTypes::scalar const  t
                                                                 , ComputationalPhysics::Types::BasicTypes::   vec const &st
                                                                 ) noexcept -> ComputationalPhysics::Types::BasicTypes::vec
    {
        auto const &A = [=]( ComputationalPhysics::Types::BasicTypes::scalar const x
                                                   , ComputationalPhysics::Types::BasicTypes::scalar const y
                                                   ) noexcept -> ComputationalPhysics::Types::BasicTypes::scalar
        {
            return std::sqrt(std::pow( (x + mu)
                                           * (x + mu)
                                           + y * y, 3
                                        )
                            );
        };
        auto const &B = [=]( ComputationalPhysics::Types::BasicTypes::scalar const x
                                                   , ComputationalPhysics::Types::BasicTypes::scalar const y
                                                   ) noexcept -> ComputationalPhysics::Types::BasicTypes::scalar
        {
            return std::sqrt(std::pow( (x - eta)
                                           * (x - eta)
                                           + y * y, 3
                                        )
                            );
        };
        ComputationalPhysics::Types::BasicTypes::scalar const x = st.x();
        ComputationalPhysics::Types::BasicTypes::scalar const u = st.y();
        ComputationalPhysics::Types::BasicTypes::scalar const y = st.z();
        ComputationalPhysics::Types::BasicTypes::scalar const v = st.w();
        return ComputationalPhysics::Types::BasicTypes::vec
                {{
                    u,
                    x + 2 * v
                      - eta * ( (x + mu)  / A(x, y) )
                      - mu  * ( (x - eta) / B(x, y) ),
                    v,
                    y - 2 * u
                      - eta * ( y / A(x, y) )
                      - mu  * ( y / B(x, y) ),
                 }};
    };
    auto const &errFunc = +[]( ComputationalPhysics::Types::BasicTypes::vec const &v1
                                                           , ComputationalPhysics::Types::BasicTypes::vec const &v2
                                                           ) noexcept -> ComputationalPhysics::Types::BasicTypes::scalar
    {
        return (v1 - v2).squaredNorm();
    };

    std::ofstream out("aristofen.csv");

    out << "X" << '\t';
    out << "Y" << '\t' << '\n';
    auto const &result = ComputationalPhysics::Integrators::Embedded::embeddedMethod<7u>( init
                                                                                        , h0
                                                                                        , broad
                                                                                        , tolerance
                                                                                        , embeddedTable
                                                                                        , rightPart
                                                                                        , errFunc
                                                                                        );
    for (unsigned i = 0u; i < result.size(); i += 100u)
    {
//        out << stuff.t << '\t';
        out << result[i].state.x() << '\t';
        out << result[i].state.z() << '\n';
    }
    out.close();
}