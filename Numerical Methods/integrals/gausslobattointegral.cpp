/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2008 Klaus Spanderen

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

/*! \file gausslabottointegral.cpp
    \brief integral of a one-dimensional function using the adaptive 
    Gauss-Lobatto integral
*/

#include "gausslobattointegral.hpp"
#include <algorithm>
#include <vector>
namespace QuantLib {

    const double GaussLobattoIntegral::alpha_ = std::sqrt(2.0/3.0); 
    const double GaussLobattoIntegral::beta_  = 1.0/std::sqrt(5.0);
    const double GaussLobattoIntegral::x1_    = 0.94288241569547971906; 
    const double GaussLobattoIntegral::x2_    = 0.64185334234578130578;
    const double GaussLobattoIntegral::x3_    = 0.23638319966214988028;

    GaussLobattoIntegral::GaussLobattoIntegral(size_t maxIterations,
                                               double absAccuracy,
                                               double relAccuracy,
                                               bool useConvergenceEstimate)
    : Integrator(absAccuracy, maxIterations),
      relAccuracy_(relAccuracy),
      useConvergenceEstimate_(useConvergenceEstimate) {
    }

    double GaussLobattoIntegral::integrate(
                                     const ext::function<double (double)>& f, 
                                     double a, double b) const {

        setNumberOfEvaluations(0);
        const double calcAbsTolerance = calculateAbsTolerance(f, a, b);

        increaseNumberOfEvaluations(2);
        return adaptivGaussLobattoStep(f, a, b, f(a), f(b), calcAbsTolerance);
    }
    double GaussLobattoIntegral::calculateAbsTolerance(
                                     const ext::function<double (double)>& f, 
                                     double a, double b) const {
        

        double relTol = std::max(relAccuracy_, QL_EPSILON);
        
        const double m = (a+b)/2; 
        const double h = (b-a)/2;
        const double y1 = f(a);
        const double y3 = f(m-alpha_*h);
        const double y5 = f(m-beta_*h);
        const double y7 = f(m);
        const double y9 = f(m+beta_*h);
        const double y11= f(m+alpha_*h);
        const double y13= f(b);

        const double f1 = f(m-x1_*h);
        const double f2 = f(m+x1_*h);
        const double f3 = f(m-x2_*h);
        const double f4 = f(m+x2_*h);
        const double f5 = f(m-x3_*h);
        const double f6 = f(m+x3_*h);

        double acc=h*(0.0158271919734801831*(y1+y13)
                  +0.0942738402188500455*(f1+f2)
                  +0.1550719873365853963*(y3+y11)
                  +0.1888215739601824544*(f3+f4)
                  +0.1997734052268585268*(y5+y9) 
                  +0.2249264653333395270*(f5+f6)
                  +0.2426110719014077338*y7);  
        
        increaseNumberOfEvaluations(13);
        if (acc == 0.0 && (   f1 != 0.0 || f2 != 0.0 || f3 != 0.0
                           || f4 != 0.0 || f5 != 0.0 || f6 != 0.0)) {
            QL_FAIL("can not calculate absolute accuracy "
                    "from relative accuracy");
        }

        double r = 1.0;
        if (useConvergenceEstimate_) {
            const double integral2 = (h/6)*(y1+y13+5*(y5+y9));
            const double integral1 = (h/1470)*(77*(y1+y13)+432*(y3+y11)+
                                             625*(y5+y9)+672*y7);
        
            if (std::fabs(integral2-acc) != 0.0) 
                r = std::fabs(integral1-acc)/std::fabs(integral2-acc);
            if (r == 0.0 || r > 1.0)
                r = 1.0;
        }

        if (relAccuracy_ != Null<double>())
            return std::min(absoluteAccuracy(), acc*relTol)/(r*QL_EPSILON);
        else {
            return absoluteAccuracy()/(r*QL_EPSILON);
        }
    }
    
    double GaussLobattoIntegral::adaptivGaussLobattoStep(
                                     const ext::function<double (double)>& f,
                                     double a, double b, double fa, double fb,
                                     double acc) const {
        QL_REQUIRE(numberOfEvaluations() < maxEvaluations(),
                   "max number of iterations reached");
        
        const double h=(b-a)/2; 
        const double m=(a+b)/2;
        
        const double mll=m-alpha_*h; 
        const double ml =m-beta_*h; 
        const double mr =m+beta_*h; 
        const double mrr=m+alpha_*h;
        
        const double fmll= f(mll);
        const double fml = f(ml);
        const double fm  = f(m);
        const double fmr = f(mr);
        const double fmrr= f(mrr);
        increaseNumberOfEvaluations(5);
        
        const double integral2=(h/6)*(fa+fb+5*(fml+fmr));
        const double integral1=(h/1470)*(77*(fa+fb)
                                       +432*(fmll+fmrr)+625*(fml+fmr)+672*fm);
        
        // avoid 80 bit logic on x86 cpu
        volatile double dist = acc + (integral1-integral2);
        if(const_cast<double&>(dist)==acc || mll<=a || b<=mrr) {
            QL_REQUIRE(m>a && b>m,"Interval contains no more machine number");
            return integral1;
        }
        else {
            return  adaptivGaussLobattoStep(f,a,mll,fa,fmll,acc)  
                  + adaptivGaussLobattoStep(f,mll,ml,fmll,fml,acc)
                  + adaptivGaussLobattoStep(f,ml,m,fml,fm,acc)
                  + adaptivGaussLobattoStep(f,m,mr,fm,fmr,acc)
                  + adaptivGaussLobattoStep(f,mr,mrr,fmr,fmrr,acc)
                  + adaptivGaussLobattoStep(f,mrr,b,fmrr,fb,acc);
        }
    }
}
