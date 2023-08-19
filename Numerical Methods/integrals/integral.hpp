/*
 Copyright (C) 2007 François du Vignaud

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

/*! \file integral.hpp
\brief Integrators base class definition
*/

#ifndef quantlib_math_integrator_hpp
#define quantlib_math_integrator_hpp

#include "functional.hpp"

namespace QuantLib {

    class Integrator{
      public:
        Integrator(double absoluteAccuracy,
                   size_t maxEvaluations);
        virtual ~Integrator() = default;

        double operator()(const ext::function<double (double)>& f,
                        double a,
                        double b) const;

        //! \name Modifiers
        //@{
        void setAbsoluteAccuracy(double);
        void setMaxEvaluations(size_t);
        //@}

        //! \name Inspectors
        //@{
        double absoluteAccuracy() const;
        size_t maxEvaluations() const;
        //@}

        double absoluteError() const ;

        size_t numberOfEvaluations() const;

        virtual bool integrationSuccess() const;

      protected:
        virtual double integrate(const ext::function<double (double)>& f,
                               double a,
                               double b) const = 0;
        void setAbsoluteError(double error) const;
        void setNumberOfEvaluations(size_t evaluations) const;
        void increaseNumberOfEvaluations(size_t increase) const;
      private:
        double absoluteAccuracy_;
        mutable double absoluteError_;
        size_t maxEvaluations_;
        mutable size_t evaluations_;
    };

}


#endif
