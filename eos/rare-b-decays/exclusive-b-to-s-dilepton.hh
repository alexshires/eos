/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_S_DILEPTON_HH 1

#include <eos/utils/power_of.hh>

#include <array>

namespace eos
{
    namespace btovll
    {
        struct Amplitudes
        {
            complex<double> a_long_right, a_long_left;
            complex<double> a_perp_right, a_perp_left;
            complex<double> a_par_right, a_par_left;
            complex<double> a_timelike;
        };

        struct AngularCoefficients
        {
            double j1s, j1c;
            double j2s, j2c;
            double j3;
            double j4;
            double j5;
            double j6s, j6c;
            double j7;
            double j8;
            double j9;
        };

        inline AngularCoefficients array_to_angular_coefficients(const std::array<double, 12> & arr)
        {
            AngularCoefficients a_c = { arr[0], arr[1], arr[2], arr[3], arr[4],  arr[5],
                arr[6], arr[7], arr[8], arr[9], arr[10], arr[11] };

            return a_c;
        }

        inline double decay_width(const AngularCoefficients & a_c)
        {
            // cf. [BHvD2010], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline std::array<double, 12> angular_coefficients_array(const Amplitudes & amplitudes, const double & s, const double & m_l)
        {
            // cf. [BHvD2010], p. 26, eqs. (A1)-(A11)
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l) / s;
            double beta2 = 1.0 - z;
            double beta = sqrt(beta2);

            // j1s
            result[0] = 3.0 / 4.0 * (
                  (2.0 + beta2) / 4.0 * (norm(amplitudes.a_perp_left) + norm(amplitudes.a_perp_right) + norm(amplitudes.a_par_left) + norm(amplitudes.a_par_right))
                  + z * real(amplitudes.a_perp_left * conj(amplitudes.a_perp_right) + amplitudes.a_par_left * conj(amplitudes.a_par_right)));
            // j1c
            result[1] = 3.0 / 4.0 * (
                  norm(amplitudes.a_long_left) + norm(amplitudes.a_long_right)
                  + z * (norm(amplitudes.a_timelike) + 2.0 * real(amplitudes.a_long_left * conj(amplitudes.a_long_right))));
            // j2s
            result[2] = 3.0 * beta2 / 16.0 * (
                  norm(amplitudes.a_perp_left) + norm(amplitudes.a_perp_right) + norm(amplitudes.a_par_left) + norm(amplitudes.a_par_right));
            // j2c
            result[3] = -3.0 * beta2 / 4.0 * (
                  norm(amplitudes.a_long_left) + norm(amplitudes.a_long_right));
            // j3
            result[4] = 3.0 / 8.0 * beta2 * (
                  norm(amplitudes.a_perp_left) + norm(amplitudes.a_perp_right) - norm(amplitudes.a_par_left) - norm(amplitudes.a_par_right));
            // j4
            result[5] = 3.0 / (4.0 * sqrt(2.0)) * beta2 * real(
                  amplitudes.a_long_left * conj(amplitudes.a_par_left) + amplitudes.a_long_right * conj(amplitudes.a_par_right));
            // j5
            result[6] = 3.0 * sqrt(2.0) / 4.0 * beta * real(
                  amplitudes.a_long_left * conj(amplitudes.a_perp_left) - amplitudes.a_long_right * conj(amplitudes.a_perp_right));
            // j6s
            result[7] = 3.0 / 2.0 * beta * real(
                  amplitudes.a_par_left * conj(amplitudes.a_perp_left) - amplitudes.a_par_right * conj(amplitudes.a_perp_right));
            // j6c
            result[8] = 0.0;
            // j7
            result[9] = 3.0 * sqrt(2.0) / 4.0 * beta * imag(
                  amplitudes.a_long_left * conj(amplitudes.a_par_left) - amplitudes.a_long_right * conj(amplitudes.a_par_right));
            // j8
            result[10] = 3.0 / 4.0 / sqrt(2.0) * beta2 * (imag(
                  amplitudes.a_long_left * conj(amplitudes.a_perp_left)) + imag(amplitudes.a_long_right * conj(amplitudes.a_perp_right)));
            // j9
            result[11] = 3.0 / 4.0 * beta2 * (imag(
                  conj(amplitudes.a_par_left) * amplitudes.a_perp_left) + imag(conj(amplitudes.a_par_right) * amplitudes.a_perp_right));

            return result;
        }
    }



    //base class for B->Xll decays
    //
    //  abstract out all common methods
    //
    class BToPDilepton {
        public :
            //need some virtual functions??
            
            // Angular Observables
            double a_l(const double & s) const;
            double b_l(const double & s) const;
            double c_l(const double & s) const;

            // unnormalised width -= basic funciton for calculating
            //double unnormalized_decay_width(const double & s) const
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

    } ;
}

#endif
