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
#include <eos/rare-b-decays/exclusive-ratio-K-pi.hh>
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
     * Ratio of B->pill v.s. B->Kll at Low Recoil
     */

    //implementation struct
    //
    RatioPiToKDilepton<LowRecoil>::RatioPiToKDilepton(const Parameters & parameters, 
                                                      const Options & options) :
        pimumu(parameters, options),
        kmumu(parameters, options)
    {
        ; //
    }


    double RatioPiToKDilepton<LowRecoil>::diff_branching_ratio(const double &s) const
    {
        double pi = pimumu.differential_branching_ratio(s) ;
        double k = kmumu.differential_branching_ratio(s) ;
        return k > 0 ? pi / k : 0 ;
    }

    double RatioPiToKDilepton<LowRecoil>::int_branching_ratio(const double &smin, const double &smax)  const
    {
        double pi = pimumu.integrated_branching_ratio(smin,smax) ;
        double k = kmumu.integrated_branching_ratio(smin,smax) ;
        return k > 0 ? pi / k : 0 ;
    }

    RatioPiToKDilepton<LargeRecoil>::RatioPiToKDilepton(const Parameters & parameters, 
                                                      const Options & options) :
        pimumu(parameters, options),
        kmumu(parameters, options)
    {
        ; //
    }


    double RatioPiToKDilepton<LargeRecoil>::diff_branching_ratio(const double &s) const
    {
        double pi = pimumu.differential_branching_ratio(s) ;
        double k = kmumu.differential_branching_ratio(s) ;
        return k > 0 ? pi / k : 0 ;
    }
    double RatioPiToKDilepton<LargeRecoil>::int_branching_ratio(const double &smin, const double &smax)  const
    {
        double pi = pimumu.integrated_branching_ratio(smin,smax) ;
        double k = kmumu.integrated_branching_ratio(smin,smax) ;
        return k > 0 ? pi / k : 0 ;
    }

}

