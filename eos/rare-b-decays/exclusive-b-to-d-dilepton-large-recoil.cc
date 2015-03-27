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
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-large-recoil.hh>
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
    //using namespace eos::btovll;
    using std::norm;

    /*
     * Decay: B -> pi l lbar at Large Recoil
     */
    template <>
    struct Implementation<BToPiDilepton<LargeRecoil>>
    {
        Parameters parameters;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter m_b_MSbar;

        UsedParameter m_c;

        UsedParameter m_d_MSbar;

        UsedParameter m_B;

        UsedParameter m_pi;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter f_B;

        UsedParameter f_K;

        UsedParameter lambda_B_p;

        UsedParameter a_1;

        UsedParameter a_2;

        // Mean life times
        UsedParameter tau;

        // Estimation of subleading contributions
        UsedParameter lambda_psd;

        UsedParameter sl_phase_psd;

        // spectator quark charge
        double e_q;

        // spectator quark flavor
        char q;

        bool cp_conjugate;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            parameters(p),
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_c(p["mass::c"], u),
            m_d_MSbar(p["mass::s(2GeV)"], u),
            m_B(p["mass::B_" + o.get("q", "d")], u),
            m_pi(p["mass::pi^" + o.get("c", "+")], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            f_B(p["decay-constant::B_" + o.get("q", "u")], u),
            f_K(p["decay-constant::pi"], u),
            lambda_B_p(p["lambda_B_p"], u),
            a_1(p["B->K::a_1@1GeV"], u),
            a_2(p["B->K::a_2@1GeV"], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            lambda_psd(p["B->Pll::Lambda_pseudo@LargeRecoil"], u),
            sl_phase_psd(p["B->Pll::sl_phase_pseudo@LargeRecoil"], u),
            e_q(-1.0/3.0),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false")))
        {
            form_factors = FormFactorFactory<PToP>::create("B->pi@" + o.get("form-factors", "BCL2008"), p);

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
        }

        complex<double> calT(const double & s) const
        {
            // charges of down- and up-type quarks
            static const double e_d = -1.0 / 3.0;
            static const double e_u = +2.0 / 3.0;

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
            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_ud())) / (model->ckm_tb() * conj(model->ckm_td()));
            if (cp_conjugate)
                lambda_hat_u = std::conj(lambda_hat_u);
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            // Compute the QCDF Integrals
            double invm1_psd = 3.0 * (1.0 + a_1 + a_2); // <ubar^-1>
            QCDFIntegrals::Results qcdf_0 = QCDFIntegrals::dilepton_massless_case(s, m_B, m_pi, mu, 0.0, 0.0, a_1, a_2);
            QCDFIntegrals::Results qcdf_c = QCDFIntegrals::dilepton_charm_case(s, m_c_pole, m_B, m_pi, mu, 0.0, 0.0, a_1, a_2);
            QCDFIntegrals::Results qcdf_b = QCDFIntegrals::dilepton_bottom_case(s, m_b_PS, m_B, m_pi, mu, 0.0, 0.0, a_1, a_2);

            // inverse of the "negative" moment of the B meson LCDA
            // cf. [BFS2001], Eq. (54), p. 15
            double omega_0 = lambda_B_p, lambda_B_p_inv = 1.0 / lambda_B_p;
            complex<double> lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

            /* Y(s) for the up and the top sector */
            // cf. [BFS2001], Eq. (10), p. 4
            complex<double> Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5();
            complex<double> Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64.0 / 3.0 * wc.c6());
            complex<double> Y_top_ = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

            // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
            // then replace b pole mass by the PS mass.
            complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole);
            Y_top += Y_top_b * CharmLoops::h(mu, s, m_b_PS);
            Y_top += Y_top_0 * CharmLoops::h(mu, s);
            Y_top += Y_top_;
            // cf. [BFS2004], Eq. (43), p. 24
            complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

            /* Effective wilson coefficients */
            // cf. [BFS2001], below Eq. (9), p. 4
            complex<double> c7eff = wc.c7() - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6();
            // cf. [BFS2001], below Eq. (26), p. 8
            complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

            /* top sector */
            // cf. [BHP2007], Eq. (B.2) and [BFS2001], Eqs. (14), (15), p. 5, in comparison with \delta_{2,3} = 1
            complex<double> C0_top_psd = 1.0 * (c7eff + wc.c7prime() + m_B / (2.0 * m_b_PS) * Y_top);
            // cf. [BHP2007], Eq. (B.2) and [BFS2004], Eq. (45), p. 24
            // the correct sign in front of C_7^eff is plus, as one can see by
            // comparison with [BF2001], Eq. (63)
            complex<double> C1f_top_psd = 1.0 * (c7eff + wc.c7prime()) * (8.0 * std::log(m_b_PS / mu) + 2.0 * L - 4.0 * (1.0 - mu_f() / m_b_PS));
            // cf. [BHP2007], Eq. (B.2) and [BFS2001], Eqs. (38), p. 9
            complex<double> C1nf_top_psd = -(+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                        + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                        + c8eff * CharmLoops::F89_massless(s, m_b_PS)));

            /* parallel, up sector */
            // cf. [BHP2007], Eq. (B.2) and [BFS2004], comment before Eq. (43), p. 24
            complex<double> C0_up_psd = 1.0 * m_B / (2.0 * m_b_PS) * Y_up;
            // C1f_up_par = 0, cf. second-to-last paragraph in Sec A.1, p. 24
            // cf. [BFS2004], last paragraph in Sec A.1, p. 24
            // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            // Use here FF_massive - FF_massless because FF_massless is defined with an extra '-'
            // compared to [S2004]
            complex<double> C1nf_up_psd = -(+1.0 / QCD::casimir_f) * (
                    (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                    + (m_B / (2.0 * m_b_PS)) * (
                        wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                        + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

            // compute the factorizing contributions
            complex<double> C_psd = C0_top_psd + lambda_hat_u * C0_up_psd
                + a_mu * (C1f_top_psd + C1nf_top_psd + lambda_hat_u * C1nf_up_psd);

            /* parallel, top sector */
            // T0_top_par_p = 0, cf. [BFS2001], Eq. (17), p. 6
            // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
            complex<double> T0_top_psd_m = +e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv;
            // cf. [BHP2007], Eq. (B.2)
            complex<double> T1f_top_psd_p  = -(c7eff + wc.c7prime()) * (4.0 * m_B / energy) * invm1_psd * lambda_B_p_inv;
            // T1f_top_par_m = 0, cf. [BFS2001], Eq. (22), p. 7
            // cf. [BFS2001], Eq. (25), p. 7
            complex<double> T1nf_top_psd_p = -m_B / m_b_PS * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2001], Eq. (26), pp. 7-8
            complex<double> T1nf_top_psd_m = -e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                    + 6.0 * m_B / m_b_PS * (
                        (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                        + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                        + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                        -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

            /* parallel, up sector */
            // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
            // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
            complex<double> T0_up_psd_m = -e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv;
            // cf. [BFS2004], Eq. (50), p. 25
            complex<double> T1nf_up_psd_p = -e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv;
            // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
            complex<double> T1nf_up_psd_m = -e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


            // Compute the nonfactorizing contributions
            complex<double> T_psd = a_mu_f * (T1f_top_psd_p + T1nf_top_psd_p + lambda_hat_u * T1nf_up_psd_p)
                + (T0_top_psd_m + lambda_hat_u * T0_up_psd_m + a_mu_f * (T1nf_top_psd_m + lambda_hat_u * T1nf_up_psd_m));

            // Subleading weak annihilation and hard spectator interaction contributions have only been
            // computed for calT_perp, not for calT_par ~ calT_psd.

            // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
            complex<double> result;
            result = xi_pseudo(s) * C_psd
                + power_of<2>(M_PI) / 3.0 * (f_B * f_K) / m_B  * T_psd;

            return result;
        }

        /* Form factors */
        // cf. [BF2001], Eq. (22)
        double xi_pseudo(const double & s) const
        {
            return form_factors->f_p(s);
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

        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double lam(const double & s) const
        {
            return lambda(m_B() * m_B(), m_pi() * m_pi(), s);
        }

        double energy(const double & s) const
        {
            return (m_B() * m_B() + m_pi() * m_pi() - s) / (2.0 * m_B());
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_A(const WilsonCoefficients<BToS> & wc, const double &) const
        {
            return wc.c10() + wc.c10prime();
        }

        double F_Tkin(const double & s) const
        {
            double result = 2.0 * std::sqrt(lam(s)) * beta_l(s) / (m_B() + m_pi());
            result *= form_factors->f_t(s) / form_factors->f_p(s);
            return result;
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_T(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Tkin(s) * wc.cT();
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_T5(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Tkin(s) * wc.cT5();
        }

        double F_Skin(const double & s) const
        {
            double result = 0.5 * (power_of<2>(m_B()) - power_of<2>(m_pi())) / (m_b_MSbar - m_d_MSbar);
            result *= form_factors->f_0(s) / form_factors->f_p(s);
            return result;
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_S(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Skin(s) * (wc.cS() + wc.cSprime());
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_P(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return F_Skin(s) * (wc.cP() + wc.cPprime()) + m_l() * (wc.c10() + wc.c10prime()) *
                    ((m_B() * m_B() - m_pi() * m_pi()) / s * (form_factors->f_0(s) / form_factors->f_p(s) - 1.0) - 1.0);
        }

        // cf. [BHP2007], Eq. (3.2), p. 4
        std::complex<double> F_V(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            std::complex<double> result = wc.c9() + wc.c9prime();
            result += 2.0 * m_b_PS() / m_B() / xi_pseudo(s) *
                      (calT(s) + lambda_psd / m_B * std::polar(1.0, sl_phase_psd()));
            result += 8.0 * m_l / (m_B() + m_pi()) * form_factors->f_t(s) / form_factors->f_p(s) * wc.cT();
            return result;
        }

        // cf. [BHP2007], Eqs. (4.2), (4.4), (4.5), p. 5
        double N(const double & s) const
        {
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_td()));

            return power_of<2>(g_fermi * alpha_e() * lambda_t) * std::sqrt(lam(s)) * beta_l(s) * xi_pseudo(s) * xi_pseudo(s) /
                    (512.0 * power_of<5>(M_PI) * power_of<3>(m_B()));
        }

        // cf. [BHP2007], Eq. (4.2)
        double a_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::norm(F_S(wc, s)) + std::norm(F_P(wc, s)));
            result += 0.25 * lam(s) * (std::norm(F_A(wc, s)) + std::norm(F_V(wc, s)));
            result += 2.0 * m_l * (m_B() * m_B() - m_pi() * m_pi() + s) * std::real(F_P(wc, s) * std::conj(F_A(wc, s)));
            result += 4.0 * m_l * m_l * m_B() * m_B() * std::norm(F_A(wc, s));

            return N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.3)
        double b_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::real(F_S(wc, s) * std::conj(F_T(wc, s)))
                                 + std::real(F_P(wc, s) * std::conj(F_T5(wc, s))));
            result += m_l * (std::sqrt(lam(s)) * beta_l(s) * std::real(F_S(wc, s) * std::conj(F_V(wc, s)))
                             + (m_B() * m_B() - m_pi() * m_pi() + s) * std::real(F_T5(wc, s) * std::conj(F_A(wc, s))));

            return 2.0 * N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.4)
        double c_l(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double result = s * (power_of<2>(beta_l(s)) * std::norm(F_T(wc, s)) + std::norm(F_T5(wc, s)));
            result -= 0.25 * lam(s) * power_of<2>(beta_l(s)) * (std::norm(F_A(wc, s)) + std::norm(F_V(wc, s)));
            result += 2.0 * m_l * std::sqrt(lam(s)) * beta_l(s) * std::real(F_T(wc, s) * std::conj(F_V(wc, s)));
            return N(s) * result;
        }

        // cf. [BHP2007], Eq. (4.8)
        double unnormalized_decay_width(const double & s) const
        {
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            return 2.0 * (a_l(wc, s) + c_l(wc, s) / 3.0);
        }

        double differential_branching_ratio(const double & s) const
        {
            return unnormalized_decay_width(s) * tau() / hbar();
        }

        // cf. [BHP2007], Eq. (4.9)
        double differential_flat_term_numerator(const double & s) const
        {
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            return 2.0 * (a_l(wc, s) + c_l(wc, s));
        }

        double differential_forward_backward_asymmetry_numerator(const double & s) const
        {
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            return b_l(wc, s);
        }
    };

    BToPiDilepton<LargeRecoil>::BToPiDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPiDilepton<LargeRecoil>>(new Implementation<BToPiDilepton<LargeRecoil>>(parameters, options, *this))
    {
    }

    BToPiDilepton<LargeRecoil>::~BToPiDilepton()
    {
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_A(const double & s) const
    {
        return _imp->F_A(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_V(const double & s) const
    {
        return _imp->F_V(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_S(const double & s) const
    {
        return _imp->F_S(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_P(const double & s) const
    {
        return _imp->F_P(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_T(const double & s) const
    {
        return _imp->F_T(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    std::complex<double>
    BToPiDilepton<LargeRecoil>::F_T5(const double & s) const
    {
        return _imp->F_T5(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    double
    BToPiDilepton<LargeRecoil>::a_l(const double & s) const
    {
        return _imp->a_l(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    double
    BToPiDilepton<LargeRecoil>::b_l(const double & s) const
    {
        return _imp->b_l(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    double
    BToPiDilepton<LargeRecoil>::c_l(const double & s) const
    {
        return _imp->c_l(_imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate), s);
    }

    double
    BToPiDilepton<LargeRecoil>::two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        auto wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        // cf. [BHP2007], Eq. (4.1)
        return _imp->a_l(wc, s) + _imp->b_l(wc, s) * c_theta_l + _imp->c_l(wc, s) * c_theta_l * c_theta_l;
    }

    double
    BToPiDilepton<LargeRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToPiDilepton<LargeRecoil>::differential_flat_term(const double & s) const
    {
        return _imp->differential_flat_term_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToPiDilepton<LargeRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return _imp->differential_forward_backward_asymmetry_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToPiDilepton<LargeRecoil>::differential_ratio_muons_electrons(const double & s) const
    {
        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            br_electrons = BToPiDilepton<LargeRecoil>::differential_branching_ratio(s);
        }

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            br_muons = BToPiDilepton<LargeRecoil>::differential_branching_ratio(s);
        }

        return br_muons / br_electrons;
    }

    // Integrated Observables
    double
    BToPiDilepton<LargeRecoil>::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> f = std::bind(std::mem_fn(&BToPiDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        return integrate(f, 64, s_min, s_max);
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> f = std::bind(&BToPiDilepton<LargeRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        double br = integrate(f, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrate(f, 64, s_min, s_max);

        return (br + br_bar) / 2.0;
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> f = std::bind(&BToPiDilepton<LargeRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        double br = integrate(f, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrate(f, 64, s_min, s_max);

        return (br - br_bar) / (br + br_bar);
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate(num, 64, s_min, s_max);
        denom_integrated += integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    // todo caching of denominator?
    double
    BToPiDilepton<LargeRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LargeRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate(num, 64, s_min, s_max);
        denom_integrated += integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LargeRecoil>::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToPiDilepton<LargeRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

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
}