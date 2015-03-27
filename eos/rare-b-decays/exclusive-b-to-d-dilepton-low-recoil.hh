/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014 Danny van Dyk
 * Copyright (c) 2010 Christian Wacker
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_D_DILEPTON_LOW_RECOIL_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_D_DILEPTON_LOW_RECOIL_HH 1

#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /*
     * Decay: B -> K l l at Low Recoil
     */
    template <>
    class BToPiDilepton<LowRecoil> :
        public ParameterUser,
        public PrivateImplementationPattern<BToPiDilepton<LowRecoil>>
    {
        public:
            BToPiDilepton(const Parameters & parameters, const Options & options);
            ~BToPiDilepton();

            // Effective Short-Distance Couplings
            double real_c9eff(const double & s) const;
            double imag_c9eff(const double & s) const;
            double real_c7eff(const double & s) const;
            double imag_c7eff(const double & s) const;

            // Angular Observables
            double a_l(const double & s) const;
            double b_l(const double & s) const;
            double c_l(const double & s) const;

            // Two Differential Observables
            double two_differential_decay_width(const double & s, const double & c_theta_l) const;

            // Differential Observables
            double differential_branching_ratio(const double & s) const;
            double differential_flat_term(const double & s) const;
            double differential_forward_backward_asymmetry(const double & s) const;
            double differential_ratio_muons_electrons(const double & s) const;

            // Integrated Observables
            double integrated_decay_width(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio(const double & s_min, const double & s_max) const;
            double integrated_branching_ratio_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_cp_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_flat_term(const double & s_min, const double & s_max) const;
            double integrated_flat_term_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const;
            double integrated_forward_backward_asymmetry_cp_averaged(const double & s_min, const double & s_max) const;
            double integrated_ratio_muons_electrons(const double & s_min, const double & s_max) const;
    };
}

#endif
