/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2013 Danny van Dyk
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
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>

#include <cmath>
#include <limits>
#include <vector>

using namespace test;
using namespace eos;

class BToXsGammaNLOTest :
    public TestCase
{
    public:
        BToXsGammaNLOTest() :
            TestCase("b_to_x_s_gamma_nlo_test")
        {
        }

        virtual void run() const
        {
            /* NLO */

            // Standard Model
            {
                Parameters p = Parameters::Defaults();
                // Taken from [CMM1996], p. 6, Eq. (28)
                p["c1"] = -0.480;
                p["c2"] = +1.023;
                p["c3"] = -0.0045;
                p["c4"] = -0.0640;
                p["c5"] = +0.0004;
                p["c6"] = +0.0009;
                p["Abs{c7}"] = 0.32372;
                p["Arg{c7}"] = M_PI;
                p["c8"] = -0.159167;
                p["mu"] = 5.0;
                // PDG 2010 CKM parameters
                p["CKM::A"] = 0.812;
                p["CKM::lambda"] = 0.2243;
                p["CKM::rhobar"] = 0.144;
                p["CKM::etabar"] = 0.342;
                // QED coupling as used in [CMM1996], Sec 5.(iii), p. 11
                p["QED::alpha_e(m_b)"] = 1/130.3;
                // b quark mass
                p["mass::b(MSbar)"] = 4.19;
                // HQE parameters
                p["B->B::mu_pi^2@1GeV"] = 0.4;
                p["B->B::mu_G^2@1GeV"] = 0.36;

                Options oo;
                oo.set("model", "WilsonScan");

                const double eps = 1e-9;
                BToXsGamma<NLO> decay(p, oo);

                /* Diagnostics */
                {
                    Diagnostics diagnostics = decay.diagnostics();
                    static const std::vector<std::pair<double, double>> reference
                    {
                        /* f_ij */
                        std::make_pair(+0.0151009, 1e-5), std::make_pair( 0.000000,  1e-5), // f_22(z = 0.22^2, delta = 1/6)
                        std::make_pair(-0.0143716, 1e-5), std::make_pair( 0.000000,  1e-5), // f_27(z = 0.22^2, delta = 1/6)
                        std::make_pair(+0.147744,  1e-5), std::make_pair( 0.000000,  1e-5), // f_78(z = 0.22^2, delta = 1/6)

                        /* f_ij */
                        std::make_pair(+0.0111111, 1e-5), std::make_pair( 0.000000,  1e-5), // f_22(z = 0.29^2, delta = 1/6)
                        std::make_pair(-0.0005781, 1e-5), std::make_pair( 0.000000,  1e-5), // f_27(z = 0.29^2, delta = 1/6)
                        std::make_pair(+0.147744,  1e-5), std::make_pair( 0.000000,  1e-5), // f_78(z = 0.29^2, delta = 1/6)

                        /* f_ij */
                        std::make_pair(+0.0101370, 1e-5), std::make_pair( 0.000000,  1e-5), // f_22(z = 0.22^2, delta = 1/9)
                        std::make_pair(-0.0105441, 1e-5), std::make_pair( 0.000000,  1e-5), // f_27(z = 0.22^2, delta = 1/9)
                        std::make_pair(+0.10490,   1e-5), std::make_pair( 0.000000,  1e-5), // f_78(z = 0.22^2, delta = 1/9)

                        /* f_ij */
                        std::make_pair(+0.0074849, 1e-5), std::make_pair( 0.000000,  1e-5), // f_22(z = 0.29^2, delta = 1/9)
                        std::make_pair(-0.0009897, 1e-5), std::make_pair( 0.000000,  1e-5), // f_27(z = 0.29^2, delta = 1/9)
                        std::make_pair(+0.10490,   1e-5), std::make_pair( 0.000000,  1e-5), // f_78(z = 0.29^2, delta = 1/9)

                        /* r_i */
                        std::make_pair(+0.830895,  1e-4), std::make_pair(+0.149794,  1e-4), // r_1(0.22^2)
                        std::make_pair(-4.98537,   1e-4), std::make_pair(-0.898762,  1e-4), // r_2(0.22^2)
                        std::make_pair(+1.00589e1, 1e-4), std::make_pair(+1.08598,   1e-4), // r_3(0.22^2)
                        std::make_pair(-1.08907,   1e-4), std::make_pair(-1.2409,    1e-4), // r_4(0.22^2)
                        std::make_pair(+1.858395e2,1e-4), std::make_pair(+1.73757e1, 1e-4), // r_5(0.22^2)
                        std::make_pair(+2.78472,   1e-4), std::make_pair(-1.81132e1, 1e-4), // r_6(0.22^2)
                        std::make_pair(-1.21063e1, 1e-4), std::make_pair( 0.0,       1e-4), // r_7(0.22^2)
                        std::make_pair(+1.96456,   1e-4), std::make_pair(+2.79253,   1e-4), // r_8(0.22^2)

                        /* r_i */
                        std::make_pair(+0.682146,  1e-4), std::make_pair(+0.075000,  1e-4), // r_1(0.29^2)
                        std::make_pair(-4.09288,   1e-4), std::make_pair(-0.450002,  1e-4), // r_2(0.29^2)
                        std::make_pair(+1.00589e1, 1e-4), std::make_pair(+1.08598,   1e-4), // r_3(0.29^2)
                        std::make_pair(-1.06555,   1e-4), std::make_pair(-1.17315,   1e-4), // r_4(0.29^2)
                        std::make_pair(+1.858395e2,1e-4), std::make_pair(+1.73757e1, 1e-4), // r_5(0.29^2)
                        std::make_pair(+8.23379,   1e-4), std::make_pair(-1.51493e1, 1e-4), // r_6(0.29^2)
                        std::make_pair(-1.21063e1, 1e-4), std::make_pair( 0.0,       1e-4), // r_7(0.29^2)
                        std::make_pair(+1.96456,   1e-4), std::make_pair(+2.79253,   1e-4), // r_8(0.29^2)

                        /* g(z) */
                        std::make_pair(0.698828,   1e-6), // g(0.22^2)
                        std::make_pair(0.542035,   1e-6), // g(0.29^2)

                        /* h(z) */
                        std::make_pair(1.93647,    1e-5), // h(0.22^2)
                        std::make_pair(1.37502,    1e-5), // h(0.22^2)

                        /* R_quark(z, delta, ckm) */
                        std::make_pair(2.88e-3,    1e-5), // R_quark(0.29^2, 0.29^2, 0.976)
                        std::make_pair(3.12e-3,    1e-5), // R_quark(0.29^2, 1/3,    0.976)
                    };

                    TEST_CHECK_DIAGNOSTICS(diagnostics, reference);
                }

                TEST_CHECK_NEARLY_EQUAL(3.13833e-4, decay.integrated_branching_ratio(1.6), eps);
                TEST_CHECK_NEARLY_EQUAL(3.09344e-4, decay.integrated_branching_ratio(1.8), eps);
                TEST_CHECK_NEARLY_EQUAL(3.00152e-4, decay.integrated_branching_ratio(2.0), eps);

                TEST_CHECK_NEARLY_EQUAL(2.21890,    decay.photon_energy_moment_1(1.6), 1e-5);
                TEST_CHECK_NEARLY_EQUAL(2.22731,    decay.photon_energy_moment_1(1.8), 1e-5);
                TEST_CHECK_NEARLY_EQUAL(2.24520,    decay.photon_energy_moment_1(2.0), 1e-5);

                TEST_CHECK_NEARLY_EQUAL(0.0423421,  decay.photon_energy_moment_2(1.6), 1e-7);
                TEST_CHECK_NEARLY_EQUAL(0.0379255,  decay.photon_energy_moment_2(1.8), 1e-7);
                TEST_CHECK_NEARLY_EQUAL(0.0323984,  decay.photon_energy_moment_2(2.0), 1e-7);
            }
        }
} b_to_x_s_gamma_nlo_test;
