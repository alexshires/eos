/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#include <eos/form-factors/form-factors.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/hard-scattering.hh>
#include <eos/rare-b-decays/long-distance.hh>
#include <eos/rare-b-decays/qcdf_integrals.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/save.hh>

#include <cmath>
#include <functional>

#include <gsl/gsl_sf.h>

#include <iostream>

namespace eos
{
    using namespace eos::btovll;
    using std::norm;

    /*
     * Decay: B -> K^* l lbar at Large Recoil, cf. [BHP2008]
     */
    template <>
    struct Implementation<BToKstarDilepton<LargeRecoil>>
    {
        std::shared_ptr<Model> model;

        Parameters parameters;

        UsedParameter hbar;

        UsedParameter m_b_MSbar;

        UsedParameter m_c;

        UsedParameter m_B;

        UsedParameter m_Kstar;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter f_B;

        UsedParameter f_Kstar_par;

        UsedParameter f_Kstar_perp;

        UsedParameter lambda_B_p;

        UsedParameter a_1_par;

        UsedParameter a_2_par;

        UsedParameter a_1_perp;

        UsedParameter a_2_perp;

        UsedParameter uncertainty_par_left;

        UsedParameter uncertainty_par_right;

        UsedParameter uncertainty_perp_left;

        UsedParameter uncertainty_perp_right;

        UsedParameter uncertainty_long_left;

        UsedParameter uncertainty_long_right;

        UsedParameter uncertainty_xi_perp;

        UsedParameter uncertainty_xi_par;

        UsedParameter tau;

        double e_q;

        char q;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToV>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WilsonScan"), p, o)),
            parameters(p),
            hbar(p["hbar"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_Kstar(p["mass::K^*_d"], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            f_B(p["decay-constant::B_" + o.get("q", "d")], u),
            f_Kstar_par(p["B->K^*::f_Kstar_par"], u),
            f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], u),
            lambda_B_p(p["lambda_B_p"], u),
            a_1_par(p["B->K^*::a_1_par"], u),
            a_2_par(p["B->K^*::a_2_par"], u),
            a_1_perp(p["B->K^*::a_1_perp"], u),
            a_2_perp(p["B->K^*::a_2_perp"], u),
            uncertainty_par_left(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_par^L") + "_uncertainty@LargeRecoil"], u),
            uncertainty_par_right(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_par^R") + "_uncertainty@LargeRecoil"], u),
            uncertainty_perp_left(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_perp^L") + "_uncertainty@LargeRecoil"], u),
            uncertainty_perp_right(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_perp^R") + "_uncertainty@LargeRecoil"], u),
            uncertainty_long_left(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_0^L") + "_uncertainty@LargeRecoil"], u),
            uncertainty_long_right(p["B->K^*ll::" + std::string(destringify<bool>(o.get("simple-sl")) ? "sl" : "A_0^R") + "_uncertainty@LargeRecoil"], u),
            uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"], u),
            uncertainty_xi_par(p["formfactors::xi_par_uncertainty"], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            e_q(-1.0/3.0),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToV>::create("B->K^*@" + o.get("form-factors", "KMPW2010"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            u.uses(*form_factors);
            u.uses(*model);

            std::string spectator_quark = o.get("q", "d");
            if (spectator_quark.size() != 1)
                throw InternalError("Option q should only be one character!");

            q = spectator_quark[0];
            if (q == 'd')
            {
                e_q = -1.0 / 3.0;
            }
            else if (q == 'u')
            {
                e_q = 2.0 / 3.0;
            }
            else
                throw InternalError("Unsupported spectator quark");
#if 0
            p["Abs{c7}"]  =  0.33670; // c7eff = 0.30726
            p["Abs{c9}"]  =  4.27305;
            p["Abs{c10}"] =  4.17166;
            p["c1"]       = -0.32196;
            p["c2"]       =  1.00925;
            p["c3"]       = -0.00519;
            p["c4"]       = -0.08787;
            p["c5"]       =  0.00036;
            p["c6"]       =  0.00101;
            p["c8"]       = -0.18262; // c8eff = -0.16925
            p["decay-constant::B_d"] = 0.200;
            p["decay-constant::B_u"] = 0.200;
            p["B->K^*::f_Kstar_perp@2GeV"] = 0.165138862;

            Amplitudes a = amplitudes(3.0);
            std::cout << "A_pp^L(q^2 = 3) = " << a.a_perp_left << std::endl;
            std::cout << "A_pp^R(q^2 = 3) = " << a.a_perp_right << std::endl;
            std::cout << "A_pa^L(q^2 = 3) = " << a.a_par_left << std::endl;
            std::cout << "A_pa^R(q^2 = 3) = " << a.a_par_right << std::endl;
            std::cout << "A_00^L(q^2 = 3) = " << a.a_long_left << std::endl;
            std::cout << "A_00^R(q^2 = 3) = " << a.a_long_right << std::endl;

            std::cout << "lam_up    = " << (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts())) << std::endl;
            std::cout << "xi_par(3) = " << xi_par(3.0) << std::endl;
            std::cout << "xi_perp(3)= " << xi_perp(3.0) << std::endl;
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s();
            std::cout << "C9       = " << wc.c9() << std::endl;
            std::cout << "C10      = " << wc.c10() << std::endl;
            std::cout << "C9 - C10 = " << wc.c9() - wc.c10() << std::endl;
            std::cout << "C9 + c10 = " << wc.c9() + wc.c10() << std::endl;
            DipoleFormFactors dff = calT(3.0);
            std::cout << "s = 3" << std::endl;
            std::cout << "calT_perp_left  = " << dff.calT_perp_left << std::endl;
            std::cout << "calT_perp_right = " << dff.calT_perp_right << std::endl;
            std::cout << "calT_parallel   = " << dff.calT_parallel << std::endl;
            dff = calT(0.0);
            std::cout << "s = 0" << std::endl;
            std::cout << "calT_perp_left  = " << dff.calT_perp_left << std::endl;
            std::cout << "calT_perp_right = " << dff.calT_perp_right << std::endl;
            std::cout << "calT_parallel   = " << dff.calT_parallel << std::endl;
            throw std::string("foo");
#endif
        }

        struct DipoleFormFactors
        {
            complex<double> calT_perp_left;
            complex<double> calT_perp_right;
            complex<double> calT_parallel;
        };

        DipoleFormFactors calT(const double & s, const WilsonCoefficients<BToS> & wc) const
        {
            // charges of down- and up-type quarks
            static const double e_d = -1.0/3.0;
            static const double e_u = +2.0/3.0;

            // spectator contributions
            double delta_qu = (q == 'u' ? 1.0 : 0.0);

            // kinematics
            double m_c_pole = model->m_c_pole();
            double m_b_PS = this->m_b_PS(), m_b_PS2 = m_b_PS * m_b_PS;
            double energy = this->energy(s);
            double L = -1.0 * (m_b_PS2 - s) / s * std::log(1.0 - s / m_b_PS2);

            // couplings
            double alpha_s_mu = model->alpha_s(mu()); // alpha_s at the hard scale
            double a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI;
            double alpha_s_mu_f = model->alpha_s(std::sqrt(mu() * 0.5)); // alpha_s at the factorization scale
            double a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;
            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
            if (cp_conjugate)
                lambda_hat_u = std::conj(lambda_hat_u);

            // Compute the QCDF Integrals
            double invm1_par = 3.0 * (1.0 + a_1_par + a_2_par); // <ubar^-1>_par
            double invm1_perp = 3.0 * (1.0 + a_1_perp + a_2_perp); // <ubar^-1>_perp
            QCDFIntegrals::Results qcdf_0 = QCDFIntegrals::dilepton_massless_case(s, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
            QCDFIntegrals::Results qcdf_c = QCDFIntegrals::dilepton_charm_case(s, m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);
            QCDFIntegrals::Results qcdf_b = QCDFIntegrals::dilepton_bottom_case(s, m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);

            // inverse of the "negative" moment of the B meson LCDA
            // cf. [BFS2001], Eq. (54), p. 15
            double omega_0 = lambda_B_p, lambda_B_p_inv = 1.0 / lambda_B_p;
            complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

            /* Y(s) for the up and the top sector */
            // cf. [BFS2001], Eq. (10), p. 4
            complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
            complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6());
            complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

            // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
                 + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
                 + Y_top_0 * CharmLoops::h(mu, s)
                 + Y_top_;
            // cf. [BFS2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            /* Effective wilson coefficients */
            // cf. [BFS2001], below Eq. (9), p. 4
            complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
            // cf. [BFS2001], below Eq. (26), p. 8
            complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            /* perpendicular, top sector */
            // cf. [BFS2001], Eqs. (12), (15), p. 5, in comparison with \delta_1 = 1
            complex<double> C0_top_perp_left  = (c7eff - wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
            complex<double> C0_top_perp_right = (c7eff + wc.c7prime()) + s / (2.0 * m_b_PS * m_B) * Y_top;
            // cf. [BFS2004], Eq. (44), p. 24
            complex<double> C1f_top_perp_left  = (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
            complex<double> C1f_top_perp_right = (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu()) - L - 4.0 * (1.0 - mu_f() / m_b_PS));
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            complex<double> C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

            /* perpendicular, up sector */
            // cf. [BFS2004], comment before Eq. (43), p. 24
            complex<double> C0_up_perp = s / (2.0 * m_b_PS * m_B) * Y_up;
            // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
            // cf. [BFS2001], Eqs. (34), (37), p. 9
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            complex<double> C1nf_up_perp = (-1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (s / (2.0 * m_b_PS * m_B)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

            /* parallel, top sector */
            // cf. [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
            complex<double> C0_top_par = -1.0 * (c7eff - wc.c7prime() + m_B / (2.0 * m_b_PS) * Y_top);
            // cf. [BFS2004], Eq. (45), p. 24
            complex<double> C1f_top_par = -1.0 * (c7eff - wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
            // cf. [BFS2001], Eqs. (38), p. 9
            complex<double> C1nf_top_par = (+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

            /* parallel, up sector */
            // cf. [BFS2004], comment before Eq. (43), p. 24
            complex<double> C0_up_par = -1.0 * m_B / (2.0 * m_b_PS) * Y_up;
            // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
            // cf. [BFS2004], last paragraph in Sec A.1, p. 24
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            complex<double> C1nf_up_par = (+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

            // compute the factorizing contributions
            complex<double> C_perp_left  = C0_top_perp_left  + lambda_hat_u * C0_up_perp
                + a_mu * (C1f_top_perp_left  + C1nf_top_perp + lambda_hat_u * C1nf_up_perp);
            complex<double> C_perp_right = C0_top_perp_right + lambda_hat_u * C0_up_perp
                + a_mu * (C1f_top_perp_right + C1nf_top_perp + lambda_hat_u * C1nf_up_perp);
            complex<double> C_par = C0_top_par + lambda_hat_u * C0_up_par
                + a_mu * (C1f_top_par + C1nf_top_par + lambda_hat_u * C1nf_up_par);


            /* perpendicular, top sector */
            // T0_top_perp_{p,m} = 0, cf. [BFS2001], Eq. (17), p. 6
            // cf. [BFS2004], Eq. (49)
            complex<double> T1f_top_perp_p_left  = (c7eff - wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
            complex<double> T1f_top_perp_p_right = (c7eff + wc.c7prime()) * (2.0 * m_B / energy) * invm1_perp * lambda_B_p_inv;
            // T1f_top_perp_m = 0, cf. [BFS2001], Eq. (22), p. 7
            // cf. [BFS2001], Eq. (23), p. 7
            // [Christoph] Use c8 instead of c8eff
            complex<double> T1nf_top_perp_p = (-4.0 * e_d * c8eff * qcdf_0.j0bar_perp
                + m_B / (2.0 * m_b_PS) * (
                        e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                        + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv;
            // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

            /* perpendicular, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eq. (50), p. 25
            complex<double> T1nf_up_perp_p = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv;

            /* parallel, top sector */
            // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
            // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
            complex<double> T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
            // cf. [BFS2004], Eq. (49), p. 25
            complex<double> T1f_top_par_p  = (c7eff - wc.c7prime()) * (4.0 * m_B / energy) * invm1_par * lambda_B_p_inv;
            // T1f_top_par_m = 0, cf. [BFS2001], Eq. (22), p. 7
            // cf. [BFS2001], Eq. (25), p. 7
            complex<double> T1nf_top_par_p = m_B / m_b_PS * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2001], Eq. (26), pp. 7-8
            complex<double> T1nf_top_par_m = e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                    + 6.0 * m_B / m_b_PS * (
                        (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                        -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

            /* parallel, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
            complex<double> T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
            // cf. [BFS2004], Eq. (50), p. 25
            complex<double> T1nf_up_par_p = +e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
            complex<double> T1nf_up_par_m = +e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


            // Compute the nonfactorizing contributions
            complex<double> T_perp_left  = a_mu_f * (T1f_top_perp_p_left + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
            complex<double> T_perp_right = a_mu_f * (T1f_top_perp_p_right + T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p);
            complex<double> T_par = a_mu_f * (T1f_top_par_p + T1nf_top_par_p + lambda_hat_u * T1nf_up_par_p)
                + (T0_top_par_m + lambda_hat_u * T0_up_par_m + a_mu_f * (T1nf_top_par_m + lambda_hat_u * T1nf_up_par_m));

            // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
            // cf. [BFS2004], Eq. (51)
            complex<double> Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                    -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                    + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                        (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p));
            complex<double> Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B *
                (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p) * 3.0 * delta_qu * wc.c2();
            // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
            // cf. [BFS2004], Eqs. (52), (53)
            complex<double> Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                        (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                        - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (
                        (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                        - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6())));
            complex<double> Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                    + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                    - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                        * (qcdf_c.j6_perp - qcdf_0.j6_perp));

            // Compute the sum of the numerically leading power-suppressed contributions
            complex<double> Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp;
            complex<double> Delta_T_up_perp = Delta_T_ann_up_perp + Delta_T_hsa_up_perp;
            complex<double> Delta_T_perp = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;


            // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
            DipoleFormFactors result;
            result.calT_perp_left  = xi_perp(s) * C_perp_left
                + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_left
                + Delta_T_perp;
            result.calT_perp_right = xi_perp(s) * C_perp_right
                + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp_right
                + Delta_T_perp;
            result.calT_parallel = xi_par(s) * C_par
                + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_par * m_Kstar) / (m_B * energy) * T_par;

            return result;
        }

        /* Form factors */
        //  cf. [BHP2008], Eq. (E.4), p. 23
        double xi_perp(const double & s) const
        {
            const double factor = m_B() / (m_B() + m_Kstar());
            double result = uncertainty_xi_perp * factor * form_factors->v(s);

            return result;
        }

        double xi_par(const double & s) const
        {
            const double factor1 = (m_B() + m_Kstar()) / (2.0 * energy(s));
            const double factor2 = (1.0 - m_Kstar() / m_B());
            double result = uncertainty_xi_par * (factor1 * form_factors->a_1(s) - factor2 * form_factors->a_2(s));

            return result;
        }

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double lam(const double & s) const
        {
            return lambda(m_B() * m_B(), m_Kstar() * m_Kstar(), s);
        }

        double norm(const double & s) const
        {
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_ts()));

            return std::sqrt(power_of<2>(g_fermi() * alpha_e()) / 3.0 / 1024 / power_of<5>(M_PI) / m_B()
                    * lambda_t * lambda_t * s_hat(s)
                    * std::sqrt(lam(s)) * beta_l(s)); // cf. [BHP2008], Eq. (C.6), p. 21
        }

        inline double s_hat(double s) const
        {
            return s / m_B() / m_B();
        }

        inline double mu_f() const
        {
            return 1.5;
        }

        inline double m_b_PS() const
        {
            // Actually use the PS mass at mu_f = 1.5 GeV
            return model->m_b_ps(mu_f());
        }

        inline double energy(const double & s) const
        {
            return (m_B() * m_B() + m_Kstar() * m_Kstar() - s) / (2.0 * m_B());
        }

        /* Amplitudes */
        // cf. [BHP2008], p. 20
        Amplitudes amplitudes(const double & s) const
        {
            Amplitudes result;

            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            double shat = s_hat(s);
            double mbhat = m_b_PS() / m_B;
            double mKhat = m_Kstar / m_B;
            double norm_s = norm(s);

            DipoleFormFactors dff = calT(s, wc);

            // longitudinal amplitude
            complex<double> wilson_long_right = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime());
            complex<double> wilson_long_left  = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime());
            double prefactor_long = -1.0 / (2.0 * m_Kstar() * std::sqrt(s));

            complex<double> a = ((m_B() * m_B() - m_Kstar() * m_Kstar() - s) * 2.0 * energy(s) * xi_perp(s)
                               - lam(s) * m_B() / (m_B() * m_B() - m_Kstar() * m_Kstar()) * (xi_perp(s) - xi_par(s)));
            complex<double> b = 2.0 * m_b_PS() * (((m_B() * m_B() + 3.0 * m_Kstar() * m_Kstar() - s) * 2.0 * energy(s) / m_B()
                               - lam(s) / (m_B() * m_B() - m_Kstar() * m_Kstar())) * dff.calT_perp_left
                               - lam(s) / (m_B() * m_B() - m_Kstar() * m_Kstar()) * dff.calT_parallel);

            result.a_long_right = norm_s * uncertainty_long_right * prefactor_long * (wilson_long_right * a + b);
            result.a_long_left  = norm_s * uncertainty_long_left  * prefactor_long * (wilson_long_left  * a + b);

            // perpendicular amplitude
            double prefactor_perp = +std::sqrt(2.0) * m_B() * std::sqrt(lambda(1.0, mKhat * mKhat, shat));
            complex<double> wilson_perp_right = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime());
            complex<double> wilson_perp_left  = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime());

            result.a_perp_right = norm_s * uncertainty_perp_right * prefactor_perp * (wilson_perp_right * xi_perp(s) + (2.0 * mbhat / shat) * dff.calT_perp_right);
            result.a_perp_left  = norm_s * uncertainty_perp_left  * prefactor_perp * (wilson_perp_left  * xi_perp(s) + (2.0 * mbhat / shat) * dff.calT_perp_right);

            // parallel amplitude
            double prefactor_par = -std::sqrt(2.0) * m_B() * (1.0 - shat);
            complex<double> wilson_par_right = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime());
            complex<double> wilson_par_left  = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime());

            result.a_par_right = norm_s * uncertainty_par_right * prefactor_par * (wilson_par_right * xi_perp(s) +
                                    (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * dff.calT_perp_left);
            result.a_par_left  = norm_s * uncertainty_par_left  * prefactor_par * (wilson_par_left  * xi_perp(s) +
                                    (2.0 * mbhat / shat) * (1.0 - mKhat * mKhat) * dff.calT_perp_left);

            // timelike amplitude
            double m_Kstarhat = m_Kstar / m_B;

            result.a_timelike = norm_s * m_B * sqrt(lambda(1.0, power_of<2>(m_Kstarhat), shat) / shat) *
                complex<double>(0.0, 2.0) * (wc.c10() - wc.c10prime()) * form_factors->a_0(s);

            return result;
        }

        std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitudes(s), s, m_l());
        }

        AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return array_to_angular_coefficients(angular_coefficients_array(amplitudes(s), s, m_l()));
        }

        AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<BToKstarDilepton<LargeRecoil>>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate(integrand, 64, s_min, s_max);

            return array_to_angular_coefficients(integrated_angular_coefficients_array);
        }

        double a_fb_zero_crossing() const
        {
            // We trust QCDF results in a validity range from 0.5 GeV^2 < s < 6.0 GeV^2
            static const double min_result = 0.5;
            static const double max_result = 7.0;

            // use calT_perp / xi_perp = C_7 as start point
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);
            const double start = -2.0 * model->m_b_msbar(mu()) * m_B() * real(wc.c7() / wc.c9());

            double result = start;
            // clamp result to QCDF validity region
            result = std::max(min_result, result);
            result = std::min(max_result, result);

            // perform a couple of Newton-Raphson steps
            for (unsigned i = 0 ; i < 100 ; ++i)
            {
                double xplus = result * 1.03;
                double xminus = result * 0.97;

                AngularCoefficients a_c_central = differential_angular_coefficients(result);
                double f = a_c_central.j6s;
                AngularCoefficients a_c_minus   = differential_angular_coefficients(xminus);
                double f_xminus = a_c_minus.j6s;
                AngularCoefficients a_c_plus    = differential_angular_coefficients(xplus);
                double f_xplus = a_c_plus.j6s;

                double fprime = (f_xplus - f_xminus) / (xplus - xminus);

                if (std::abs(f / fprime) < 1e-8)
                    break;

                result = result - f / fprime;
                // clamp result to QCDF validity region
                result = std::max(min_result, result);
                result = std::min(max_result, result);
            }

            return result;
        }
    };

    BToKstarDilepton<LargeRecoil>::BToKstarDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarDilepton<LargeRecoil>>(new Implementation<BToKstarDilepton<LargeRecoil>>(parameters, options, *this))
    {
    }

    BToKstarDilepton<LargeRecoil>::~BToKstarDilepton()
    {
    }

    double
    BToKstarDilepton<LargeRecoil>::xi_perp(const double & s) const
    {
        return _imp->xi_perp(s);
    }

    double
    BToKstarDilepton<LargeRecoil>::xi_para(const double & s) const
    {
        return _imp->xi_par(s);
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_long(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_long_left;
        else
            return amp.a_long_right;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_perp(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_perp_left;
        else
            return amp.a_perp_right;
    }

    complex<double>
    BToKstarDilepton<LargeRecoil>::a_par(const Helicity & h, const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        if (h == -1)
            return amp.a_par_left;
        else
            return amp.a_par_right;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_decay_width(const double & s) const
    {
        return decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_isospin_asymmetry(const double & s) const
    {
        Save<char> save_q(_imp->q, 'd');
        Save<double> save_e_q(_imp->e_q, -1.0/3.0);

        double gamma_zero = differential_decay_width(s);
        _imp->q = 'u';
        _imp->e_q = +2.0/3.0;
        double gamma_minus = differential_decay_width(s);

        return (gamma_zero - gamma_minus) / (gamma_zero + gamma_minus);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.8)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.11)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.12)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((power_of<2>(_imp->beta_l(s) * a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_re(const double & s) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.25 * _imp->beta_l(s) * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_im(const double & s) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_4(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_5(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_p_prime_6(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_transversal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_1(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_2(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_3(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_h_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // differential angular coefficients
    double
    BToKstarDilepton<LargeRecoil>::differential_j_1c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_1s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_2c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_2s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3_normalized(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_3_normalized_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return (a_c.j3 + a_c_bar.j3) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j4;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j5;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_6c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6c;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_6s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_7(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j7;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_8(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j8;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9;
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9_normalized(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::differential_j_9_normalized_cp_averaged(const double & s) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        return (a_c.j9 + a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double br = integrated_branching_ratio(s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrated_branching_ratio(s_min, s_max);

        return 0.5 * (br + br_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double gamma = integrated_decay_width(s_min, s_max);
        _imp->cp_conjugate = true;
        double gamma_bar = integrated_decay_width(s_min, s_max);

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_isospin_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<char> save_q(_imp->q, 'd');
        Save<double> save_e_q(_imp->e_q, -1.0/3.0);

        double gamma_zero = integrated_decay_width(s_min, s_max);
        _imp->q = 'u';
        _imp->e_q = +2.0/3.0;
        double gamma_minus = integrated_decay_width(s_min, s_max);

        return (gamma_zero - gamma_minus) / (gamma_zero + gamma_minus);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.8), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double a_fb = integrated_forward_backward_asymmetry(s_min, s_max);
        _imp->cp_conjugate = true;
        double a_fb_bar = integrated_forward_backward_asymmetry(s_min, s_max);

        return 0.5 * (a_fb + a_fb_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j1c - a_c.j2c / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double f_l = integrated_longitudinal_polarisation(s_min, s_max);
        _imp->cp_conjugate = true;
        double f_l_bar = integrated_longitudinal_polarisation(s_min, s_max);

        return 0.5 * (f_l + f_l_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], p. 5, eq. (3.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        double f_t = integrated_transversal_polarisation(s_min, s_max);
        _imp->cp_conjugate = true;
        double f_t_bar = integrated_transversal_polarisation(s_min, s_max);

        return 0.5 * (f_t + f_t_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2_cp_averaged(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        Save<bool> save(_imp->cp_conjugate, false);

        double a_t_2 = integrated_transverse_asymmetry_2(s_min, s_max);
        _imp->cp_conjugate = true;
        double a_t_2_bar = integrated_transverse_asymmetry_2(s_min, s_max);

        return 0.5 * (a_t_2 + a_t_2_bar);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.11), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.12), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((power_of<2>(a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.25 * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_4(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_5(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_p_prime_6(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_1(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return  a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_h_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    double
    BToKstarDilepton<LargeRecoil>::a_fb_zero_crossing() const
    {
        return _imp->a_fb_zero_crossing();
    }

    // integrated angular coefficients
    double
    BToKstarDilepton<LargeRecoil>::integrated_j_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_1s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_2s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j3 + a_c_bar.j3) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j4;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j5;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_6c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6c;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_7(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j7;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_8(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j8;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9 / decay_width(a_c);
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j9 + a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_a_9(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return (a_c.j9 - a_c_bar.j9) / (decay_width(a_c) + decay_width(a_c_bar));
    }

    double
    BToKstarDilepton<LargeRecoil>::four_differential_decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const
    {
        // compute d^4 Gamma, cf. [BHvD2010], p. 5, Eq. (2.6)
        // Cosine squared of the angles
        double c_theta_k_2 = c_theta_k * c_theta_k;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = cos(phi);
        // Sine squared of the angles
        double s_theta_k_2 = 1.0 - c_theta_k_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_k = sqrt(s_theta_k_2);
        double s_theta_l = sqrt(s_theta_l_2);
        double s_phi = sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_k = 2.0 * c_theta_k_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_k = 2.0 * s_theta_k * c_theta_k;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = sin(2.0 * phi);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        return 3.0 / 8.0 / M_PI * (
                 a_c.j1s + (a_c.j1c - a_c.j1s) * c_theta_k_2
                +  (a_c.j2s + (a_c.j2c - a_c.j2s) * c_theta_k_2) * c_2_theta_l
                +  a_c.j3 * s_theta_k_2 * s_theta_l_2 * c_2_phi
                +  a_c.j4 * s_2_theta_k * s_2_theta_l * c_phi
                +  a_c.j5 * s_2_theta_k * s_theta_l * c_phi
                +  (a_c.j6s * s_theta_k_2 + a_c.j6c * c_theta_k_2) * c_theta_l
                +  a_c.j7 * s_2_theta_k * s_theta_l * s_phi
                +  a_c.j8 * s_2_theta_k * s_2_theta_l * s_phi
                +  a_c.j9 * s_theta_k_2 * s_theta_l_2 * s_2_phi
                );
    }

    //lepton universality
    //
    double
    BToKstarDilepton<LargeRecoil>::differential_ratio_muons_electrons(const double & s) const
    {
        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            br_electrons = BToKstarDilepton<LargeRecoil>::differential_branching_ratio(s);
        }
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            br_muons = BToKstarDilepton<LargeRecoil>::differential_branching_ratio(s);
        }
        return br_muons / br_electrons;
    }

    double
    BToKstarDilepton<LargeRecoil>::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToKstarDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1) ;
        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            br_electrons = integrate(integrand, 64, s_min, s_max);
        }
        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            br_muons = integrate(integrand, 64, s_min, s_max);
        }
        // cf. [BHP2007], Eq. (4.10), p. 6
        return br_muons / br_electrons;
    }

    //magnitudes of amplitudes
    double 
    BToKstarDilepton<LargeRecoil>::mag_amp_perp(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return std::abs( amp.a_perp_left ) + std::abs( amp.a_perp_right ) ;
    }
    double 
    BToKstarDilepton<LargeRecoil>::mag_amp_para(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return std::abs( amp.a_par_left ) + std::abs( amp.a_par_right ) ;
    }
    double
    BToKstarDilepton<LargeRecoil>::mag_amp_long(const double & s) const
    {
        Amplitudes amp = _imp->amplitudes(s);
        return std::abs( amp.a_long_left ) + std::abs( amp.a_long_right ) ;
    }

}
