/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2020 Klaus Spanderen

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

/*! \file exponentialintegrals.hpp
*/

#include "types.hpp"

#include <complex>

namespace QuantLib {
    /*! References:

        B. Rowe et al: GALSIM: The modular galaxy image simulation toolkit
        https://arxiv.org/abs/1407.7676
    */
    namespace ExponentialIntegral {
        double Si(double x);
        double Ci(double x);

        std::complex<double> Ci(std::complex<double> z);
        std::complex<double> Si(std::complex<double> z);
        std::complex<double> E1(std::complex<double> z);
        std::complex<double> Ei(std::complex<double> z);
    }
}
