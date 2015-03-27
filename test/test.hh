/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
 * Copyright (c) 2014 Frederik Beaujean
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#ifndef WILSON_FITTER_GUARD_TEST_TEST_HH
#define WILSON_FITTER_GUARD_TEST_TEST_HH 1

#include <config.h>
#include <eos/utils/stringify.hh>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace test
{
    class TestCase
    {
        private:
            std::string _name;

        public:
            TestCase(const std::string & name);

            virtual ~TestCase();

            std::string name() const;

            virtual void run() const = 0;
    };

    class TestCaseFailedException
    {
        private:
            int _line;

            std::string _file;

            std::string _reason;

        public:
            TestCaseFailedException(int line, const std::string & file, const std::string & reason);

            const std::string & reason() const;

            std::string where() const;
    };

#define TEST_SECTION(name, body) \
    do \
    { \
        std::cout << name << "> begins" << std::endl; \
        body \
        std::cout << name << "> ends" << std::endl; \
    } \
    while (false)

#define TEST_CHECK(a) \
    do \
    { \
        if (! (a)) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is false"); \
    } \
    while (false)

#define TEST_CHECK_MSG(a, msg) \
    do \
    { \
        if (! (a)) \
            throw TestCaseFailedException(__LINE__, __FILE__, msg); \
    } \
    while (false)

#define TEST_CHECK_EQUAL(a, b) \
    do \
    { \
        if (! ((a) == (b))) \
            throw TestCaseFailedException(__LINE__, __FILE__, "'" #a "' is not equal to '" #b "'"); \
    } \
    while (false)

#define TEST_CHECK_FAILED(s) \
    do \
    { \
        throw TestCaseFailedException(__LINE__, __FILE__, s); \
    } \
    while (false)

#define TEST_CHECK_NEARLY_EQUAL(a, b, eps) \
    do \
    { \
        auto a_val = (a); \
        auto b_val = (b); \
        if (std::abs((a_val - b_val)) <= eps) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' = " + stringify(a_val, 16u) + " is not nearly-equal to '" #b "' = " + stringify(b_val, 16u) + " within '" + stringify(eps, 16u) + "'" \
                    + ", difference is '" + stringify(a_val - b_val, 16u) + "'"); \
    } \
    while (false)

#define TEST_CHECK_RELATIVE_ERROR(a, b, eps) \
    do \
    { \
        auto a_val = (a); \
        auto b_val = (b); \
        if (std::sqrt(std::fabs(a_val)) < std::numeric_limits<double>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #a "' has been evaluated to the zero within computational accuracy, result = " + stringify(a_val)); \
         \
        if (std::sqrt(std::fabs(b_val)) < std::numeric_limits<double>::epsilon()) \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "'" #b "' has been evaluated to the zero within computational accuracy, result = " + stringify(b_val)); \
         \
        if (((std::abs((a_val - b_val) / a_val)) <= eps) && ((std::abs((a_val - b_val) / b_val)) <= eps)) \
            break; \
        else \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "One relative error of '" #a "' = '" + stringify(a_val, 16u) + "' and '" #b "' = '" + stringify(b_val, 16u) + "' is greater than " + stringify(eps, 16u) + ". The results are " + \
                    stringify(std::abs((a_val - b_val) / a_val)) + " and " + stringify(std::abs((a_val - b_val) / b_val), 16u)); \
    } \
    while (false)

#define TEST_CHECK_NO_THROW(expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception in '" #expression "'"); \
        } \
    } \
    while (false)

#define TEST_CHECK_THROWS(exception, expression) \
    do \
    { \
        try \
        { \
            expression; \
        } \
        catch (exception & e) \
        { \
            break; \
        } \
        catch (...) \
        { \
            throw TestCaseFailedException(__LINE__, __FILE__, \
                    "Caught unexpected exception when expecting " #exception " in '" #expression "'"); \
        } \
        \
        throw TestCaseFailedException(__LINE__, __FILE__, \
                "Caught no exception in " #expression " when expecting '" #exception "'"); \
    } \
    while (false)

#define TEST_CHECK_DIAGNOSTICS(diagnostics, reference) \
    do \
    { \
        if (diagnostics.size() != reference.size()) \
            throw TestCaseFailedException(__LINE__, __FILE__, "Number of diagnostics and reference entries differ!"); \
        auto i = diagnostics.begin(), i_end = diagnostics.end(); \
        auto j = reference.begin(); \
        for ( ; i != i_end ; ++i, ++j) \
        { \
            if (std::abs(j->first - i->value) > j->second) \
                throw TestCaseFailedException(__LINE__, __FILE__, "Diagnostic error: " + i->description + "\n\tevaluates to " + stringify(i->value, 7u) + "\n\tdelta to reference value " + stringify(j->first, 7u) + " is " + stringify(j->first - i->value, 7u)); \
        } \
    } \
    while (false)

}
#endif
