#pragma once
#include <types/core_types.h>
#include <solvers/three_diagonal_solver.h>
namespace ComputationalPhysics::Interpolators
{

Types::CoreTypes::splineSeries naturalSpline( Types::BasicTypes::vecSTL const &x
                                            , Types::BasicTypes::vecSTL const &y
                                            , Types::BasicTypes::scalar const h
                                            ) noexcept
{
    Types::BasicTypes::vecSTL a;
    Types::BasicTypes::vecSTL b;
    Types::BasicTypes::vecSTL c;
    Types::BasicTypes::vecSTL d;
    Types::BasicTypes::vecSTL m;
    Types::BasicTypes::vecSTL l;
    Types::BasicTypes::vecSTL z;
    Types::BasicTypes::vecSTL f;

    a.insert(a.begin(), y.begin(), y.end());
    unsigned const n = a.size() - 1;

            b.resize(n);
    d.resize(n);
    c.resize(n + 1);

    l.resize(n + 1);
    m.resize(n + 1);
    z.resize(n + 1);

    l[0] = 1.;
    z[0] = 0.;
    m[0] = 0.;
    f.push_back(0.);

    for (unsigned i = 1u; i < n; ++i)
        f.push_back( 3. * (a[i + 1] - a[i]) / h
                   - 3. * (a[i] - a[i - 1]) / h
                   );

    for (unsigned i = 1u; i < n; ++i)
    {
        l[i] = 2 * h - h * m[i - 1];
        m[i] = h / l[i];
        z[i] = (f[i] - h * z[i - 1]) / l[i];
    }

    l[n] = 1.;
    z[n] = 0.;
    c[n] = 0.;

    for (unsigned i = n - 2u; i > 0u; --i)
    {
        c[i] = z[i] - m[i] * c[i + 1];
        b[i] = (a[i + 1] -     a[i]) / h
             - (c[i + 1] + 2 * c[i]) * h / 3.;
        d[i] = (c[i + 1] - c[i]) / (3. * h);
    }

    Types::CoreTypes::splineSeries result(n);
    for (unsigned i = 0u; i < n; ++i)
    {
        result[i].a = a[i];
        result[i].b = b[i];
        result[i].c = c[i];
        result[i].d = d[i];
        result[i].x = x[i];
    }

    return result;
}

Types::BasicTypes::scalar naturalSplineValue( Types::CoreTypes::splineSeries const &data
                                            , Types::BasicTypes::scalar const x
                                            ) noexcept
{
    unsigned const n = data.size();
    Types::BasicTypes::vecSTL intervals;

    for (auto const item: data)
        intervals.push_back(item.x);

    if (x >= intervals.back())
    {
        auto const [a, b, c, d, x0] = data.back();
        return a + b * (x - x0) + c * (x - x0) * (x - x0) + d * (x - x0) * (x - x0) * (x - x0);
    } else if (x <= intervals.front())
    {
        auto const [a, b, c, d, x0] = data.front();
        return a + b * (x - x0) + c * (x - x0) * (x - x0) + d * (x - x0) * (x - x0) * (x - x0);
    }

    for (unsigned i = 0u; i < n - 1u; ++i)
    {
        if (intervals[i] <= x <= intervals[i + 1])
        {
            auto const [a, b, c, d, x0] = data[i];
            return a + b * (x - x0) + c * (x - x0) * (x - x0) + d * (x - x0) * (x - x0) * (x - x0);
        }
    }
}

} // ComputationalPhysics::Interpolators
