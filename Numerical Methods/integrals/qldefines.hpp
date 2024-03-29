/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2005, 2006, 2007 StatPro Italia srl
 Copyright (C) 2015 CompatibL

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

/*! \file qldefines.hpp
    \brief Global definitions and compiler switches.
*/

#ifndef quantlib_defines_hpp
/* install-hook */
#define quantlib_defines_hpp

#ifdef _MSC_VER
/* Microsoft-specific, but needs to be defined before
   including <boost/config.hpp> which somehow includes
   <math.h> under VC++10
*/
#define _USE_MATH_DEFINES
#endif
#include <boost/config.hpp>
#include <boost/version.hpp>
#if BOOST_VERSION < 104800
    #error using an old version of Boost, please update.
#endif
#if !defined(BOOST_ENABLE_ASSERT_HANDLER)
    #define BOOST_ENABLE_ASSERT_HANDLER
#endif

/* This allows one to include a given file at this point by
   passing it as a compiler define (e.g., -DQL_INCLUDE_FIRST=foo.hpp).

   The idea is to provide a hook for defining QL_REAL and at the
   same time including any necessary headers for the new type.
*/
#define INCLUDE_FILE(F) INCLUDE_FILE_(F)
#define INCLUDE_FILE_(F) #F
#ifdef QL_INCLUDE_FIRST
#    include INCLUDE_FILE(QL_INCLUDE_FIRST)
#endif
#undef INCLUDE_FILE_
#undef INCLUDE_FILE

/* Eventually these might go into userconfig.hpp.
   For the time being, we hard code them here.
   They can be overridden by passing the #define to the compiler.
*/
#ifndef QL_INTEGER
#    define QL_INTEGER int
#endif

#ifndef QL_BIG_INTEGER
#    define QL_BIG_INTEGER long
#endif

#ifndef QL_REAL
#   define QL_REAL double
#endif


/*! \defgroup macros QuantLib macros

    Global definitions and a few macros which help porting the
    code to different compilers.

    @{
*/

#if (defined(_DEBUG) || defined(DEBUG))
    #define QL_DEBUG
#endif

#if   defined(HAVE_CONFIG_H)    // Dynamically created by configure
   #include <ql/config.hpp>
/* Use BOOST_MSVC instead of _MSC_VER since some other vendors (Metrowerks,
   for example) also #define _MSC_VER
*/
#elif defined(BOOST_MSVC)       // Microsoft Visual C++
   #include "config.msvc.hpp"
#elif defined(__MINGW32__)      // Minimalistic GNU for Windows
   #include <ql/config.mingw.hpp>
#elif defined(__SUNPRO_CC)      // Sun Studio
   #include <ql/config.sun.hpp>
#else                           // We hope that the compiler follows ANSI
   #include <ql/config.ansi.hpp>
#endif


// extra debug checks
#ifdef QL_DEBUG
    #ifndef QL_EXTRA_SAFETY_CHECKS
        #define QL_EXTRA_SAFETY_CHECKS
    #endif
#endif

#ifdef QL_ENABLE_THREAD_SAFE_OBSERVER_PATTERN
    #if BOOST_VERSION < 105800
        #error Boost version 1.58 or higher is required for the thread-safe observer pattern
    #endif
#endif

#ifdef QL_ENABLE_PARALLEL_UNIT_TEST_RUNNER
    #if BOOST_VERSION < 105900
        #error Boost version 1.59 or higher is required for the parallel unit test runner
    #endif
#endif

// ensure that needed math constants are defined
#ifndef quantlib_math_constants_hpp
#define quantlib_math_constants_hpp

#include <cmath>

#ifndef M_E
#define M_E         2.71828182845904523536
#endif

#ifndef M_LOG2E
#define M_LOG2E     1.44269504088896340736
#endif

#ifndef M_LOG10E
#define M_LOG10E    0.434294481903251827651
#endif

#ifndef M_IVLN10
#define M_IVLN10    0.434294481903251827651
#endif

#ifndef M_LN2
#define M_LN2       0.693147180559945309417
#endif

#ifndef M_LOG2_E
#define M_LOG2_E    0.693147180559945309417
#endif

#ifndef M_LN10
#define M_LN10      2.30258509299404568402
#endif

#ifndef M_PI
#    define M_PI 3.141592653589793238462643383280
#endif

#ifndef M_TWOPI
#define M_TWOPI     (M_PI * 2.0)
#endif

#ifndef M_PI_2
#define M_PI_2      1.57079632679489661923
#endif

#ifndef M_PI_4
#define M_PI_4      0.785398163397448309616
#endif

#ifndef M_3PI_4
#define M_3PI_4     2.3561944901923448370E0
#endif

#ifndef M_SQRTPI
#define M_SQRTPI    1.77245385090551602792981
#endif

#ifndef M_1_PI
#define M_1_PI      0.318309886183790671538
#endif

#ifndef M_2_PI
#define M_2_PI      0.636619772367581343076
#endif

#ifndef M_1_SQRTPI
#define M_1_SQRTPI  0.564189583547756286948
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI  1.12837916709551257390
#endif

#ifndef M_SQRT2
#define M_SQRT2     1.41421356237309504880
#endif

#ifndef M_SQRT_2
#define M_SQRT_2    0.7071067811865475244008443621048490392848359376887
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2   0.7071067811865475244008443621048490392848359376887
#endif

#ifndef M_LN2LO
#define M_LN2LO     1.9082149292705877000E-10
#endif

#ifndef M_LN2HI
#define M_LN2HI     6.9314718036912381649E-1
#endif

#ifndef M_SQRT3
#define M_SQRT3     1.73205080756887719000
#endif

#ifndef M_INVLN2
#define M_INVLN2    1.4426950408889633870E0
#endif

/* This should ensure that no macro are redefined if we happen to
   include <math.h> again, whether or not we're using our macros
   or theirs. We can't know in advance, since it depends on the
   order of inclusion of headers in client code. */
#ifdef _MSC_VER
#undef _USE_MATH_DEFINES
#endif

#endif




// import global functions into std namespace
#if defined(BOOST_NO_STDC_NAMESPACE)
    #include <cmath>
    namespace std {
        using ::sqrt; using ::abs; using ::fabs;
        using ::exp; using ::log; using ::pow;
        using ::sin; using ::cos; using ::asin; using ::acos;
        using ::sinh; using ::cosh;
        using ::floor; using ::fmod; using ::modf;
    }
#endif


/*! \defgroup limitMacros Numeric limits

    Some compilers do not give an implementation of
    <code>\<limits\></code> yet.  For the code to be portable
    these macros should be used instead of the corresponding method of
    <code>std::numeric_limits</code> or the corresponding macro
    defined in <code><limits.h></code>.

    @{
*/
/*! \def QL_MIN_INTEGER
    Defines the value of the largest representable negative integer value
*/
/*! \def QL_MAX_INTEGER
    Defines the value of the largest representable integer value
*/
/*! \def QL_MIN_REAL
    Defines the value of the largest representable negative
    floating-point value
*/
/*! \def QL_MIN_POSITIVE_REAL
    Defines the value of the smallest representable positive double value
*/
/*! \def QL_MAX_REAL
    Defines the value of the largest representable floating-point value
*/
/*! \def QL_EPSILON
    Defines the machine precision for operations over doubles
*/
#include <limits>
// limits used as such
#define QL_MIN_INTEGER         ((std::numeric_limits<QL_INTEGER>::min)())
#define QL_MAX_INTEGER         ((std::numeric_limits<QL_INTEGER>::max)())
#define QL_MIN_REAL           -((std::numeric_limits<QL_REAL>::max)())
#define QL_MAX_REAL            ((std::numeric_limits<QL_REAL>::max)())
#define QL_MIN_POSITIVE_REAL   ((std::numeric_limits<QL_REAL>::min)())
#define QL_EPSILON             ((std::numeric_limits<QL_REAL>::epsilon)())
/*! \def QL_NULL_INTEGER
    \deprecated Don't use this macro.
                Deprecated in version 1.27.
*/
#define QL_NULL_INTEGER        ((std::numeric_limits<int>::max)())
/*! \def QL_NULL_REAL
    \deprecated Don't use this macro.
                Deprecated in version 1.27.
*/
#define QL_NULL_REAL           ((std::numeric_limits<float>::max)())
/*! @} */

/*! @}  */


// For the time being we're keeping a QL_DEPRECATED macro because
// of <https://stackoverflow.com/questions/38378693/>.  We need to
// use it to deprecate constructors until we drop support for VC++2015.
// Other features (methods, typedefs etc.) can use [[deprecated]] and
// possibly add a message.

// emit warning when using deprecated features
// clang-format off
#if defined(BOOST_MSVC)       // Microsoft Visual C++
#    define QL_DEPRECATED __declspec(deprecated)
#    define QL_DEPRECATED_DISABLE_WARNING \
        __pragma(warning(push))           \
        __pragma(warning(disable : 4996))
#    define QL_DEPRECATED_ENABLE_WARNING \
        __pragma(warning(pop))
#elif defined(__GNUC__)
#    define QL_DEPRECATED __attribute__((deprecated))
#    define QL_DEPRECATED_DISABLE_WARNING                               \
        _Pragma("GCC diagnostic push")                                  \
        _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")
#    define QL_DEPRECATED_ENABLE_WARNING \
        _Pragma("GCC diagnostic pop")
#elif defined(__clang__)
#    define QL_DEPRECATED __attribute__((deprecated))
#    define QL_DEPRECATED_DISABLE_WARNING                                 \
        _Pragma("clang diagnostic push")                                  \
        _Pragma("clang diagnostic ignored \"-Wdeprecated-declarations\"")
#    define QL_DEPRECATED_ENABLE_WARNING \
        _Pragma("clang diagnostic pop")
#else
// we don't know how to enable it, just define the macros away
#    define QL_DEPRECATED
#    define QL_DEPRECATED_DISABLE_WARNING
#    define QL_DEPRECATED_ENABLE_WARNING
#endif
// clang-format on

/*! \deprecated Use the noexcept keyword instead.
                Deprecated in version 1.27.
*/
#define QL_NOEXCEPT noexcept

/*! \deprecated Use the constexpr keyword instead.
                Deprecated in version 1.27.
*/
#define QL_CONSTEXPR constexpr

/*! \deprecated Do not check; always use std::unique_ptr.
                Deprecated in version 1.27
*/
#ifndef QL_USE_STD_UNIQUE_PTR
#define QL_USE_STD_UNIQUE_PTR 1
#endif

#endif
