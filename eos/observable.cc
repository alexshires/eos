/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
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

#include <eos/observable.hh>
#include <eos/form-factors/form-factor-adapter.hh>
#include <eos/form-factors/baryonic-impl.hh>
#include <eos/form-factors/mesonic-impl.hh>
#include <eos/b-decays/b-to-pi-l-nu.hh>
#include <eos/b-decays/properties.hh>
#include <eos/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/exclusive-B-to-K-dilepton.hh>
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-low-recoil.hh>
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-large-recoil.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/exclusive-ratio-K-pi.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/utils/concrete_observable.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

#include <map>

namespace eos
{
    ObservableNameError::ObservableNameError(const std::string & name) :
        Exception("Observable name '" + name + "' is malformed")
    {
    }

    ObservableFactory::ObservableFactory()
    {
    }

    ObservableFactory::~ObservableFactory()
    {
    }

    template <typename Decay_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, std::make_tuple()));
    }

    template <typename Decay_, typename Tuple_, typename ... Args_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            double (Decay_::* function)(const Args_ & ...) const,
            const Tuple_ & kinematics_names)
    {
        std::string sname(name);

        return std::make_pair(sname, make_concrete_observable_factory(sname, function, kinematics_names));
    }

    template <typename Transition_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            const char * process,
            double (FormFactors<Transition_>::* function)(const double &) const)
    {
        std::string sname(name), sprocess(process);

        return std::make_pair(sname, new FormFactorAdapterFactory<Transition_>(sname, sprocess, function));
    }

    template <typename Transition_>
    std::pair<std::string, ObservableFactory *> make_observable(const char * name,
            const char * process,
            double (FormFactors<Transition_>::* numerator)(const double &) const,
            double (FormFactors<Transition_>::* denominator)(const double &) const)
    {
        std::string sname(name), sprocess(process);

        return std::make_pair(sname, new FormFactorRatioAdapterFactory<Transition_>(sname, sprocess, numerator, denominator));
    }

    ObservablePtr
    Observable::make(const std::string & _name, const Parameters & parameters, const Kinematics & kinematics, const Options & _options)
    {
        static const std::map<std::string, ObservableFactory *> simple_observables
        {
            /* B Meson Properties */
            make_observable("B::M_B^*-M_B",
                    &BMesonProperties::mass_splitting_j1_j0),

            /* Form Factor for the Exclusive Decays */

            // B -> pi Form Factors
            make_observable("B->pi::f_+(s)", "B->pi",
                    &FormFactors<PToP>::f_p),

            make_observable("B->pi::f_+'(s)", "B->pi",
                    &FormFactors<PToP>::f_p_d1),

            make_observable("B->pi::f_+''(s)", "B->pi",
                    &FormFactors<PToP>::f_p_d2),

            make_observable("B->pi::f_T(s)", "B->pi",
                    &FormFactors<PToP>::f_t),

            make_observable("B->pi::f_0(s)", "B->pi",
                    &FormFactors<PToP>::f_0),

            // B -> pi Form Factors (auxillary variable for determining the
            // LCSR/SVZ threshold parameters)
            make_observable("B->pi::M_B(LCSR)@DKMMO2008",
                    &AnalyticFormFactorBToPiDKMMO2008::MB_lcsr,
                    std::make_tuple("s")),

            make_observable("B->pi::M_B(SVZ)@DKMMO2008",
                    &AnalyticFormFactorBToPiDKMMO2008::MB_svz),

            // B -> K Form Factors
            make_observable("B->K::f_+(s)", "B->K",
                    &FormFactors<PToP>::f_p),

            make_observable("B->K::f_T(s)", "B->K",
                    &FormFactors<PToP>::f_t),

            make_observable("B->K::f_0(s)", "B->K",
                    &FormFactors<PToP>::f_0),

            // B -> K^* Form Factors
            make_observable("B->K^*::V(s)", "B->K^*",
                    &FormFactors<PToV>::v),

            make_observable("B->K^*::A_0(s)", "B->K^*",
                    &FormFactors<PToV>::a_0),

            make_observable("B->K^*::A_1(s)", "B->K^*",
                    &FormFactors<PToV>::a_1),

            make_observable("B->K^*::A_2(s)", "B->K^*",
                    &FormFactors<PToV>::a_2),

            make_observable("B->K^*::A_12(s)", "B->K^*",
                    &FormFactors<PToV>::a_12),

            make_observable("B->K^*::V(s)/A_1(s)", "B->K^*",
                    &FormFactors<PToV>::v, &FormFactors<PToV>::a_1),

            make_observable("B->K^*::A_2(s)/A_1(s)", "B->K^*",
                    &FormFactors<PToV>::a_2, &FormFactors<PToV>::a_1),

            // Lambda_b -> Lambda Form Factors
            make_observable("Lambda_b->Lambda::f_perp^V(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_v),

            make_observable("Lambda_b->Lambda::f_perp^A(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_perp_a),

            make_observable("Lambda_b->Lambda::f_long^V(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_v),

            make_observable("Lambda_b->Lambda::f_long^A(s)", "Lambda_b->Lambda",
                    &FormFactors<OneHalfPlusToOneHalfPlus>::f_long_a),

            /* Exclusive B Decays */

            // B -> pi l nu
            make_observable("B->pilnu::dBR/ds",
                    &BToPiLeptonNeutrino::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->pilnu::BR",
                    &BToPiLeptonNeutrino::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->pilnu::zeta",
                    &BToPiLeptonNeutrino::integrated_zeta,
                    std::make_tuple("s_min", "s_max")),

            /* Exclusive Rare B Decays */

            // B_q -> ll
            make_observable("B_q->ll::BR",
                    &BToDilepton::branching_ratio_time_zero),

            make_observable("B_q->ll::BR@Untagged",
                    &BToDilepton::branching_ratio_untagged_integrated),

            // B -> K^* gamma
            make_observable("B->K^*gamma::BR",
                    &BToKstarGamma::branching_ratio),

            make_observable("B->K^*gamma::BRavg",
                    &BToKstarGamma::branching_ratio_cp_averaged),

            make_observable("B->K^*gamma::A_CP",
                    &BToKstarGamma::cp_asymmetry),

            make_observable("B->K^*gamma::S_K^*gamma",
                    &BToKstarGamma::s_kstar_gamma),

            make_observable("B->K^*gamma::C_K^*gamma",
                    &BToKstarGamma::c_kstar_gamma),

            make_observable("B->K^*gamma::A_I",
                    &BToKstarGamma::isospin_asymmetry),

            // B -> K ll, Large Recoil
            make_observable("B->Kll::d^2Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)")),

            make_observable("B->Kll::dBR/ds@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::A_FB(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::BRavg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_CP@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::Gamma@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_Havg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_flat_term_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FB@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FBavg@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::b_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::b_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::c_l,
                    std::make_tuple("s")),

            // B -> K ll, Low Recoil
            make_observable("B->Kll::d^2Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::two_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)")),

            make_observable("B->Kll::dBR/ds@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->Kll::F_H(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_flat_term,
                    std::make_tuple("s")),

            make_observable("B->Kll::A_FB(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->Kll::BR@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::BRavg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_CP@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::Gamma@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_H@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_flat_term,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::F_Havg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_flat_term_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FB@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::A_FBavg@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->Kll::a_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::a_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::b_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::b_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::c_l@LowRecoil",
                    &BToKDilepton<LowRecoil>::c_l,
                    std::make_tuple("s")),

            make_observable("B->Kll::Re{c9eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::real_c9eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Im{c9eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::imag_c9eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Re{c7eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::real_c7eff,
                    std::make_tuple("s")),

            make_observable("B->Kll::Im{c7eff}@LowRecoil",
                    &BToKDilepton<LowRecoil>::imag_c7eff,
                    std::make_tuple("s")),

            // B -> K^* ll, Large Recoil
            make_observable("B->K^*ll::xi_perp(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::xi_perp,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::xi_para(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::xi_para,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::d^4Gamma@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::four_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)", "cos(theta_k)", "phi")),

            make_observable("B->K^*ll::dBR/ds@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_I(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_isospin_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^2(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^3(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^4(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^5(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^re(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_re,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^im(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transverse_asymmetry_im,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_4(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_p_prime_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_5(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_p_prime_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_6(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_p_prime_6,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_L(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_T(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_transversal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1s(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_1s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1c(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_1c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2s(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_2s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2c(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_2c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3norm(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_3_normalized,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3normavg(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_3_normalized_cp_averaged,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_4(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_5(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6s(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_6s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6c(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_6c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_7(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_7,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_8(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_8,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_9,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9norm(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_9_normalized,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9normavg(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_j_9_normalized_cp_averaged,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_FBavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BR@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BRavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_CP@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_L@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_Lavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_longitudinal_polarisation_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_T@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_Tavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transversal_polarisation_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2avg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_2_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^3@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^4@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^5@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^re@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_re,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^im@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_transverse_asymmetry_im,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_4@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_p_prime_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_5@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_p_prime_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_6@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_p_prime_6,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^1(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_h_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^2(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_h_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^3(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_h_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^4(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_h_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^5(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_h_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^1@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_h_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^2@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_h_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^3@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_h_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^4@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_h_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^5@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_h_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::s_0^A_FB@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::a_fb_zero_crossing),

            make_observable("B->K^*ll::Gamma@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1s@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_1s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1c@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2s@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_2s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2c@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3norm@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3normavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_3_normalized_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_4@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_5@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6s@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6c@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_6c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_7@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_7,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_8@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_8,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_9,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9norm@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9normavg@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_j_9_normalized_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_9@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_a_9,
                    std::make_tuple("s_min", "s_max")),

            // B -> K^* ll, Low Recoil
            make_observable("B->K^*ll::d^4Gamma@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::four_differential_decay_width,
                    std::make_tuple("s", "cos(theta_l)", "cos(theta_k)", "phi")),

            make_observable("B->K^*ll::dBR/ds@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_forward_backward_asymmetry,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^4(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^5(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^re(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_re,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_T^im(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transverse_asymmetry_im,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_4(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_p_prime_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_5(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_p_prime_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::P'_6(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_p_prime_6,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_L(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_longitudinal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::F_T(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_transversal_polarisation,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^1(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_h_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_h_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_h_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^4(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_h_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::H_T^5(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_h_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1s(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_1s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_1c(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_1c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2s(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_2s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_2c(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_2c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3norm(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_3_normalized,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_3normavg(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_3_normalized_cp_averaged,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_4(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_4,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_5(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_5,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6s(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_6s,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_6c(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_6c,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_7(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_7,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_8(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_8,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_9,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9norm(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_9_normalized,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::J_9normavg(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_j_9_normalized_cp_averaged,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::rho_1(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::rho_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::rho_2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::rho_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::A_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_FBavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Abar_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_unnormalized_forward_backward_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_FB@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_forward_backward_asymmetry_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BR@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::BRavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_branching_ratio_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_CP@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_L@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_Lavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_T@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transversal_polarisation,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::F_Tavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transversal_polarisation_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nF_L@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_longitudinal_polarisation_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^2avg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_2_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_3_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nA_T^4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_4_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^5@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^re@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_re,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_T^im@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_transverse_asymmetry_im,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_p_prime_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_5@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_p_prime_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::P'_6@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_p_prime_6,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_1_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_2_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::nH_T^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_3_naive,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::H_T^5@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_h_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Re{Y}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::real_y,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Im{Y}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::imag_y,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Re{C_9^eff}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::real_c9eff,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::Im{C_9^eff}(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::imag_c9eff,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^1(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_1,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^2(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_2,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^3(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_3,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^mix(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_cp_asymmetry_mix,
                    std::make_tuple("s")),

            make_observable("B->K^*ll::a_CP^1@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::a_CP^2@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_2,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::a_CP^3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_asymmetry_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Gamma+Gammabar@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_cp_summed_decay_width,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::Gamma-Gammabar@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_unnormalized_cp_asymmetry_1,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1s@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_1s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_1c@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_1c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2s@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_2s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_2c@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_2c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_3,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3norm@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_3_normalized,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_3normavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_3_normalized_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_4@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_4,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_5@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_5,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6s@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_6s,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_6c@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_6c,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_7@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_7,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_8@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_8,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_9,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9norm@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_9_normalized,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::J_9normavg@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_j_9_normalized_cp_averaged,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::A_9@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_a_9,
                    std::make_tuple("s_min", "s_max")),

            // Lambda_b -> Lambda l^+ l^-, Large Recoil
            make_observable("Lambda_b->Lambdall::dBR/ds@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^l(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^h(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^c(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::F_0(s)@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::BR@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^l@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^h@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^c@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::F_0@LargeRecoil",
                    &LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero,
                    std::make_tuple("s_min", "s_max")),

            // Lambda_b -> Lambda l^+ l^-, Low Recoil
            make_observable("Lambda_b->Lambdall::dBR/ds@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^l(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^h(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::A_FB^c(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::F_0(s)@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::differential_fzero,
                    std::make_tuple("s")),

            make_observable("Lambda_b->Lambdall::BR@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^l@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^h@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::A_FB^c@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined,
                    std::make_tuple("s_min", "s_max")),

            make_observable("Lambda_b->Lambdall::F_0@LowRecoil",
                    &LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero,
                    std::make_tuple("s_min", "s_max")),

            /* Inclusive Decays */

            // B->X_s ll, HLMW2005
            make_observable("B->X_sll::dBR/ds@HLMW2005",
                    &BToXsDilepton<HLMW2005>::differential_branching_ratio,
                    std::make_tuple("s")),

            make_observable("B->X_sll::BR@HLMW2005",
                    &BToXsDilepton<HLMW2005>::integrated_branching_ratio,
                    std::make_tuple("s_min", "s_max")),

            // B->X_s gamma
            make_observable("B->X_sgamma::BR@Minimal",
                    &BToXsGamma<Minimal>::integrated_branching_ratio),

            // B->X_s gamma, NLO implementation
            make_observable("B->X_sgamma::BR(E_min)@NLO",
                    &BToXsGamma<NLO>::integrated_branching_ratio,
                    std::make_tuple("E_min")),

            make_observable("B->X_sgamma::E_1(E_min)@NLO",
                    &BToXsGamma<NLO>::photon_energy_moment_1,
                    std::make_tuple("E_min")),

            make_observable("B->X_sgamma::E_2(E_min)@NLO",
                    &BToXsGamma<NLO>::photon_energy_moment_2,
                    std::make_tuple("E_min")),


            /* pimumu */

            make_observable("B->Pill::BR@LowRecoil",
                    &BToPiDilepton<LowRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min","s_max")),
            make_observable("B->Pill::dBR/ds@LowRecoil",
                    &BToPiDilepton<LowRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),
            make_observable("B->Pill::BR@LargeRecoil",
                    &BToPiDilepton<LargeRecoil>::integrated_branching_ratio,
                    std::make_tuple("s_min","s_max")),
            make_observable("B->Pill::dBR/ds@LargeRecoil",
                    &BToPiDilepton<LargeRecoil>::differential_branching_ratio,
                    std::make_tuple("s")),

            /* lepton unviersality tests */ 
            make_observable("B->Kll::R_K(s)@LargeRecoil",
                    &BToKDilepton<LargeRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),
            make_observable("B->Kll::R_K@LargeRecoil",
                    &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),
            make_observable("B->Kll::R_K(s)@LowRecoil",
                    &BToKDilepton<LowRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),
            make_observable("B->Kll::R_K@LowRecoil",
                    &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),
            make_observable("B->K^*ll::R_K*(s)@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::R_K*@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),
            make_observable("B->K^*ll::R_K*(s)@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::differential_ratio_muons_electrons,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::R_K*@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                    std::make_tuple("s_min", "s_max")),

            make_observable("B->K^*ll::a_perp@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::mag_amp_perp,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::a_para@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::mag_amp_para,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::a_long@LargeRecoil",
                    &BToKstarDilepton<LargeRecoil>::mag_amp_long,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::a_perp@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::mag_amp_perp,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::a_para@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::mag_amp_para,
                    std::make_tuple("s")),
            make_observable("B->K^*ll::a_long@LowRecoil",
                    &BToKstarDilepton<LowRecoil>::mag_amp_long,
                    std::make_tuple("s")),

            //ratio
            make_observable("RPiK::dBRds@LowRecoil",
                    &RatioPiToKDilepton<LowRecoil>::diff_branching_ratio,
                    std::make_tuple("s")),
            make_observable("RPiK::BR@LowRecoil",
                    &RatioPiToKDilepton<LowRecoil>::int_branching_ratio,
                    std::make_tuple("smin","smax")),
            make_observable("RPiK::dBRds@LargeRecoil",
                    &RatioPiToKDilepton<LowRecoil>::diff_branching_ratio,
                    std::make_tuple("s")),
            make_observable("RPiK::BR@LargeRecoil",
                    &RatioPiToKDilepton<LowRecoil>::int_branching_ratio,
                    std::make_tuple("smin","smax")),


        };

        Options options;
        std::string name(_name);

        std::string::size_type pos;
        while (std::string::npos != (pos = name.rfind(',')))
        {
            std::string::size_type sep(name.find('=', pos + 1));
            if (std::string::npos == sep)
                throw ObservableNameError(_name);

            std::string key(name.substr(pos + 1, sep - pos - 1));
            std::string value(name.substr(sep + 1));

            options.set(key, value);
            name.erase(pos);
        }

        auto i(simple_observables.find(name));
        if (simple_observables.end() == i)
            return ObservablePtr();

        return i->second->make(parameters, kinematics, options + _options);
    }
}
