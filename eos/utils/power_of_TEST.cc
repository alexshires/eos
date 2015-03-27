/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010 Danny van Dyk
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

#include <test/test.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class PowerOfTest :
    public TestCase
{
    public:
        PowerOfTest() :
            TestCase("power_of_test")
        {
        }

        virtual void run() const
        {
            {
                // Test power_of<n, double>
                static const double eps = 1e-14;

                TEST_CHECK_NEARLY_EQUAL(1.0,      power_of<0>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(1.2,      power_of<1>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(1.44,     power_of<2>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(1.728,    power_of<3>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(2.0736,   power_of<4>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(2.48832,  power_of<5>(1.2), eps);
                TEST_CHECK_NEARLY_EQUAL(2.985984, power_of<6>(1.2), eps);

                TEST_CHECK_NEARLY_EQUAL(1.0,      power_of<0>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.4,      power_of<1>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.16,     power_of<2>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.064,    power_of<3>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.0256,   power_of<4>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.01024,  power_of<5>(0.4), eps);
                TEST_CHECK_NEARLY_EQUAL(0.004096, power_of<6>(0.4), eps);
            }

            {
                // Test power_of<n, std::complex<double>>
                static const double eps = 1e-14;

                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(1.0, 0.0),       power_of<0>(std::complex<double>(1.2, 0.9)), eps);
                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(0.0, -1.0),      power_of<3>(std::complex<double>(0.0, 1.0)), eps);
                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(1.0, 0.0),       power_of<4>(std::complex<double>(0.0, 1.0)), eps);
                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(-9.6719, -1.56), power_of<4>(std::complex<double>(1.2, 1.3)), eps);
                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(2.985984, 0.0),  power_of<6>(std::complex<double>(1.2, 0.0)), eps);
                TEST_CHECK_NEARLY_EQUAL(std::complex<double>(0.004096, 0.0),  power_of<6>(std::complex<double>(0.4, 0.0)), eps);
            }
        }
} power_of_test;
