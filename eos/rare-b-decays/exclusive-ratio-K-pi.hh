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

#ifndef EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_K_PI_HH
#define EOS_GUARD_SRC_RARE_B_DECAYS_EXCLUSIVE_B_TO_K_PI_HH 1

#include <eos/rare-b-decays/decays.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-low-recoil.hh>
#include <eos/rare-b-decays/exclusive-B-to-K-dilepton.hh>
#include <eos/rare-b-decays/exclusive-b-to-d-dilepton-large-recoil.hh>

namespace eos
{
    template <>
    class RatioPiToKDilepton<LowRecoil> //:
        public ParameterUser,
        //public PrivateImplementationPattern<RatioPiToKDilepton<LowRecoil>>
    {
        public:
            RatioPiToKDilepton(const Parameters & parameters, const Options & options);
            ~RatioPiToKDilepton();

            //ratio
            double diff_branching_ratio(const double &s) const ;
            double int_branching_ratio(const double &smin, const double &smax) const ;

            //member variables
            BToPiDilepton<LowRecoil> pimumu ;
            BToKDilepton<LowRecoil> kmumu ;
    } ;

    template <>
    class RatioPiToKDilepton<LargeRecoil> //:
        public ParameterUser,
        //public PrivateImplementationPattern<RatioPiToKDilepton<LargeRecoil>>
    {
        public:
            RatioPiToKDilepton(const Parameters & parameters, const Options & options);
            ~RatioPiToKDilepton();

            //ratio
            double diff_branching_ratio(const double &s) const ;
            double int_branching_ratio(const double &smin, const double &smax) const ;

            //member variables
            BToPiDilepton<LargeRecoil> pimumu ;
            BToKDilepton<LargeRecoil> kmumu ;
    } ;

}

#endif /*  */
