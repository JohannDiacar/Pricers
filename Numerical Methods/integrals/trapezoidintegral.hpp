/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Roman Gitlin
 Copyright (C) 2003 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file trapezoidintegral.hpp
    \brief integral of a one-dimensional function using the trapezoid formula
*/

#ifndef quantlib_trapezoid_integral_hpp
#define quantlib_trapezoid_integral_hpp

#include "integral.hpp"
#include "null.hpp"
#include "errors.hpp"

namespace QuantLib {

    //! Integral of a one-dimensional function
    /*! Given a target accuracy \f$ \epsilon \f$, the integral of
        a function \f$ f \f$ between \f$ a \f$ and \f$ b \f$ is
        calculated by means of the trapezoid formula
        \f[
        \int_{a}^{b} f \mathrm{d}x =
        \frac{1}{2} f(x_{0}) + f(x_{1}) + f(x_{2}) + \dots
        + f(x_{N-1}) + \frac{1}{2} f(x_{N})
        \f]
        where \f$ x_0 = a \f$, \f$ x_N = b \f$, and
        \f$ x_i = a+i \Delta x \f$ with
        \f$ \Delta x = (b-a)/N \f$. The number \f$ N \f$ of intervals
        is repeatedly increased until the target accuracy is reached.

        \test the correctness of the result is tested by checking it
              against known good values.
    */
    template <class IntegrationPolicy>
    class TrapezoidIntegral : public Integrator {
      public:
        TrapezoidIntegral(double accuracy,
                          size_t maxIterations)
        : Integrator(accuracy, maxIterations){}

      protected:
        double integrate(const ext::function<double(double)>& f, double a, double b) const override {

            // start from the coarsest trapezoid...
            size_t N = 1;
            double I = (f(a)+f(b))*(b-a)/2.0, newI;
            // ...and refine it
            size_t i = 1;
            do {
                newI = IntegrationPolicy::integrate(f,a,b,I,N);
                N *= IntegrationPolicy::nbEvalutions();
                // good enough? Also, don't run away immediately
                if (std::fabs(I-newI) <= absoluteAccuracy() && i > 5)
                    // ok, exit
                    return newI;
                // oh well. Another step.
                I = newI;
                i++;
            } while (i < maxEvaluations());
            QL_FAIL("max number of iterations reached");
        }
    };

    // Integration policies
    struct Default {
        inline static double integrate(const ext::function<double (double)>& f, 
                                     double a, 
                                     double b, 
                                     double I, 
                                     size_t N)
        {
            double sum = 0.0;
            double dx = (b-a)/N;
            double x = a + dx/2.0;
            for (size_t i=0; i<N; x += dx, ++i)
                sum += f(x);
            return (I + dx*sum)/2.0;
        }
        inline static size_t nbEvalutions(){ return 2;}
    };

    struct MidPoint {
        inline static double integrate(const ext::function<double (double)>& f,
                                     double a, 
                                     double b, 
                                     double I, 
                                     size_t N)
        {
            double sum = 0.0;
            double dx = (b-a)/N;
            double x = a + dx/6.0;
            double D = 2.0*dx/3.0;
            for (size_t i=0; i<N; x += dx, ++i)
                sum += f(x) + f(x+D);
            return (I + dx*sum)/3.0;
        }
        inline static size_t nbEvalutions(){ return 3;}
    };

}

#endif
