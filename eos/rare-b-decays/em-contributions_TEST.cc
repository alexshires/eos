/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/em-contributions.hh>

#include <cmath>
#include <iostream>

using namespace test;
using namespace eos;

class EMContributionsTest :
    public TestCase
{
    public:
        EMContributionsTest() :
            TestCase("em_contributions_test")
        {
        }

        virtual void run() const
        {
            /* Comparison with Christophs results from October 2010 */

            static const double eps = 1.0e-5;

            static const double s_hat = 0.2, m_b = 4.6, m_l = 0.105658, mu = 4.2, log_m_l_hat = std::log(m_l / m_b);

            TEST_CHECK_NEARLY_EQUAL(+8.41414,   EMContributions::omegaem_22(s_hat, log_m_l_hat, mu), eps);
            TEST_CHECK_NEARLY_EQUAL(+2.72822,   real(EMContributions::omegaem_27(s_hat, log_m_l_hat, mu)), eps);
            TEST_CHECK_NEARLY_EQUAL(+1.07985,   imag(EMContributions::omegaem_27(s_hat, log_m_l_hat, mu)), eps);
            TEST_CHECK_NEARLY_EQUAL(+9.41651,   real(EMContributions::omegaem_29(s_hat, log_m_l_hat, mu)), eps);
            TEST_CHECK_NEARLY_EQUAL(+3.26257,   imag(EMContributions::omegaem_29(s_hat, log_m_l_hat, mu)), eps);
            TEST_CHECK_NEARLY_EQUAL(-3.77437,   EMContributions::omegaem_77(s_hat, log_m_l_hat), eps);
            TEST_CHECK_NEARLY_EQUAL(+0.0518519, EMContributions::omegaem_79(s_hat, log_m_l_hat), eps);
            TEST_CHECK_NEARLY_EQUAL(+2.11315,   EMContributions::omegaem_99(s_hat, log_m_l_hat), eps);
            TEST_CHECK_NEARLY_EQUAL(+2.02214,   EMContributions::omegaem_1010(s_hat, log_m_l_hat), eps);
        }
} em_contributions_test;
