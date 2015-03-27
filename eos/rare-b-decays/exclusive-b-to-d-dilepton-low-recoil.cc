/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014 Danny van Dyk
 * Copyright (c) 2010, 2011 Christian Wacker
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
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-low-recoil.hh>
#include <eos/rare-b-decays/hard-scattering.hh>
#include <eos/rare-b-decays/long-distance.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/log.hh>
#include <eos/utils/memoise.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/power_of.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/qcd.hh>
#include <eos/utils/save.hh>

#include <cmath>
#include <functional>

namespace eos
{
    //using namespace eos::btovll;
    using std::norm;

    
    /*
     * Decay: B -> K l lbar at Low Recoil
     */
    template <>
    struct Implementation<BToPiDilepton<LowRecoil>>
    {
        Parameters parameters;

        std::shared_ptr<Model> model;

        UsedParameter hbar;

        UsedParameter m_b_MSbar;

        UsedParameter m_d_MSbar;

        UsedParameter m_B;

        UsedParameter m_pi;

        UsedParameter m_l;

        UsedParameter mu;

        UsedParameter alpha_e;

        UsedParameter g_fermi;

        UsedParameter lambda_pseudo;

        UsedParameter sl_phase_pseudo;

        // Mean life time
        UsedParameter tau;

        bool cp_conjugate;

        bool ccbar_resonance;

        std::shared_ptr<FormFactors<PToP>> form_factors;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            parameters(p),
            model(Model::make(o.get("model", "SM"), p, o)),
            hbar(p["hbar"], u),
            m_b_MSbar(p["mass::b(MSbar)"], u),
            m_d_MSbar(p["mass::d(2GeV)"], u),
            m_B(p["mass::B_" + o.get("q", "u")], u),
            m_pi(p["mass::pi^" + o.get("c", "+")], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            mu(p["mu"], u),
            alpha_e(p["QED::alpha_e(m_b)"], u),
            g_fermi(p["G_Fermi"], u),
            lambda_pseudo(p["B->Pll::Lambda_pseudo@LowRecoil"], u),
            sl_phase_pseudo(p["B->Pll::sl_phase_pseudo@LowRecoil"], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            cp_conjugate(destringify<bool>(o.get("cp-conjugate", "false"))),
            ccbar_resonance(destringify<bool>(o.get("ccbar-resonance", "false")))
        {
            form_factors = FormFactorFactory<PToP>::create("B->pi@" + o.get("form-factors", "BCL2008"), p);

            if (! form_factors.get())
                throw InternalError("Form factors not found!");

            std::string spectator_quark = o.get("q", "d");
            if ((spectator_quark != "d") && (spectator_quark != "u"))
                throw InternalError("Unsupported spectator quark");

            u.uses(*form_factors);
            u.uses(*model);
        }

        // We use the PS mass except for kappa
        double m_b_PS() const
        {
            // Actually use m_b_PS at mu_PS = 2.0 GeV
            return model->m_b_ps(2.0);
        }

        // cf. [GP2004], Eq. (56)
        complex<double> c7eff(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            return ShortDistanceLowRecoil::c7eff(s, mu(), model->alpha_s(mu), m_b_PS(), true, wc);
        }

        // cf. [GP2004], Eq. (55), p. 10
        complex<double> c9eff(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_ud())) / (model->ckm_tb() * conj(model->ckm_td()));

            if (cp_conjugate)
            {
                lambda_hat_u = conj(lambda_hat_u);
            }

            return ShortDistanceLowRecoil::c9eff(s, mu(), model->alpha_s(mu), m_b_PS(), model->m_c_msbar(mu), true, ccbar_resonance, lambda_hat_u, wc);
        }

        double kappa() const
        {
            // cf. [BHvD2010], Eq. (3.8), p. 8
            // Use m_b_MSbar(m_b_MSbar) instead m_b_MSbar(mu), as we want kappa up to NLO only.
            return (1.0 - 2.0 * model->alpha_s(mu) / (3.0 * M_PI) * std::log(mu / m_b_MSbar));
        }

        // this is rho_1^+
        double rho_1(const WilsonCoefficients<BToS> & wc, const double & s) const
        {
            double alpha_s = model->alpha_s(mu());

            return std::norm(kappa() * (2.0 * (m_b_MSbar + lambda_pseudo()) * m_B() / s) * (c7eff(wc, s) + wc.c7prime())
                    + 0.5 * alpha_s / m_B * std::polar(lambda_pseudo(), sl_phase_pseudo()) + (c9eff(wc, s) + wc.c9prime()))
                    + std::norm(wc.c10() + wc.c10prime());
        }

        // speed of the lepton
        double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * power_of<2>(m_l()) / s);
        }

        // phase-space function
        double lam(const double & s) const
        {
            return lambda(m_B() * m_B(), m_pi() * m_pi(), s);
        }

        // cf. [BHP2007], Eq. (3.2), p. 3
        std::complex<double> F_A(const WilsonCoefficients<BToS> & wc, const double &) const
        {
            return wc.c10() + wc.c10prime();
        }

        // kinematic part of F_T and F_T5
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

        // kinematic part of F_S and F_P
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
            std::complex<double> result =  c9eff(wc, s) + wc.c9prime();
            result += kappa() * (2.0 * (m_b_MSbar + lambda_pseudo()) * m_B() / s) * (c7eff(wc, s) + wc.c7prime())
                      + 0.5 * model->alpha_s(mu) / m_B * std::polar(lambda_pseudo(), sl_phase_pseudo());
            result += 8.0 * m_l / (m_B() + m_pi()) * form_factors->f_t(s) / form_factors->f_p(s) * wc.cT();
            return result;
        }

        // cf. [BHP2007], Eqs. (4.2), (4.4), (4.5), p. 5
        double N(const double & s) const
        {
            double lambda_t = abs(model->ckm_tb() * conj(model->ckm_td()));

            return power_of<2>(g_fermi * alpha_e() * lambda_t) * std::sqrt(lam(s)) * beta_l(s) * power_of<2>(form_factors->f_p(s)) /
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

        double unnormalized_decay_width(const double & s) const
        {
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(cp_conjugate);

            return 2.0 * (a_l(wc, s) + c_l(wc, s) / 3.0);
        }

        double differential_branching_ratio(const double & s) const
        {
            return unnormalized_decay_width(s) * tau() / hbar();
        }

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

    BToPiDilepton<LowRecoil>::BToPiDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToPiDilepton<LowRecoil>>(new Implementation<BToPiDilepton<LowRecoil>>(parameters, options, *this))
    {
    }

    BToPiDilepton<LowRecoil>::~BToPiDilepton()
    {
    }

    double
    BToPiDilepton<LowRecoil>::real_c9eff(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return real(_imp->c9eff(wc, s));
    }

    double
    BToPiDilepton<LowRecoil>::imag_c9eff(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return imag(_imp->c9eff(wc, s));
    }

    double
    BToPiDilepton<LowRecoil>::real_c7eff(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return real(_imp->c7eff(wc, s));
    }

    double
    BToPiDilepton<LowRecoil>::imag_c7eff(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return imag(_imp->c7eff(wc, s));
    }

    double
    BToPiDilepton<LowRecoil>::a_l(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return _imp->a_l(wc, s);
    }

    double
    BToPiDilepton<LowRecoil>::b_l(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return _imp->b_l(wc, s);
    }

    double
    BToPiDilepton<LowRecoil>::c_l(const double & s) const
    {
        WilsonCoefficients<BToS> wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return _imp->c_l(wc, s);
    }

    // Two Differential Observables
    double
    BToPiDilepton<LowRecoil>::two_differential_decay_width(const double & s, const double & c_theta_l) const
    {
        auto wc = _imp->model->wilson_coefficients_b_to_s(_imp->cp_conjugate);

        return _imp->a_l(wc, s) + _imp->b_l(wc, s) * c_theta_l + _imp->c_l(wc, s) * c_theta_l * c_theta_l;
    }

    // Differential Observables
    double
    BToPiDilepton<LowRecoil>::differential_branching_ratio(const double & s) const
    {
        return _imp->differential_branching_ratio(s);
    }

    double
    BToPiDilepton<LowRecoil>::differential_flat_term(const double & s) const
    {
        return _imp->differential_flat_term_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToPiDilepton<LowRecoil>::differential_forward_backward_asymmetry(const double & s) const
    {
        return _imp->differential_forward_backward_asymmetry_numerator(s) / _imp->unnormalized_decay_width(s);
    }

    double
    BToPiDilepton<LowRecoil>::differential_ratio_muons_electrons(const double & s) const
    {
        double br_electrons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::e"]());
            br_electrons = BToPiDilepton<LowRecoil>::differential_branching_ratio(s);
        }

        double br_muons;
        {
            Save<Parameter, double> save_m_l(_imp->m_l, _imp->parameters["mass::mu"]());
            br_muons = BToPiDilepton<LowRecoil>::differential_branching_ratio(s);
        }

        return br_muons / br_electrons;
    }

    // Integrated Observables
    double
    BToPiDilepton<LowRecoil>::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max);
    }

    double
    BToPiDilepton<LowRecoil>::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToPiDilepton<LowRecoil>::differential_branching_ratio),
                this, std::placeholders::_1);

        return integrate(integrand, 64, s_min, s_max);
    }

    double
    BToPiDilepton<LowRecoil>::integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> integrand = std::bind(&BToPiDilepton<LowRecoil>::differential_branching_ratio,
                this, std::placeholders::_1);

        double br = integrate(integrand, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double br_bar = integrate(integrand, 64, s_min, s_max);

        return (br + br_bar) / 2.0;
    }

    double
    BToPiDilepton<LowRecoil>::integrated_flat_term(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LowRecoil>::integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::differential_flat_term_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate(num, 64, s_min, s_max);
        denom_integrated += integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LowRecoil>::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LowRecoil>::integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->cp_conjugate, false);
        std::function<double (const double &)> num = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::differential_forward_backward_asymmetry_numerator),
                _imp, std::placeholders::_1);
        std::function<double (const double &)> denom = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double num_integrated = integrate(num, 64, s_min, s_max);
        double denom_integrated = integrate(denom, 64, s_min, s_max);

        _imp->cp_conjugate = true;

        num_integrated += integrate(num, 64, s_min, s_max);
        denom_integrated += integrate(denom, 64, s_min, s_max);

        return num_integrated / denom_integrated;
    }

    double
    BToPiDilepton<LowRecoil>::integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const
    {
        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&BToPiDilepton<LowRecoil>::differential_branching_ratio),
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

    double
    BToPiDilepton<LowRecoil>::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->cp_conjugate, false);

        std::function<double (const double &)> integrand = std::bind(std::mem_fn(&Implementation<BToPiDilepton<LowRecoil>>::unnormalized_decay_width),
                _imp, std::placeholders::_1);

        double gamma = integrate(integrand, 64, s_min, s_max);
        _imp->cp_conjugate = true;
        double gamma_bar = integrate(integrand, 64, s_min, s_max);

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }
}