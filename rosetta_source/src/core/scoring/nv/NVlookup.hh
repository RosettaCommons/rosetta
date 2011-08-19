// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/NV/NVlookup.hh
/// @brief  Neighbor Vector algorithm lookup table processing class declaration
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nv_NVlookup_hh
#define INCLUDED_core_scoring_nv_NVlookup_hh

//unit headers
#include <utility/vector1.hh>

#include <core/types.hh>

//project headers
// AUTO-REMOVED #include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>
#include <core/types.hh>


//STL header

#include <map>
#include <vector>

//Auto Headers
#include <iostream>



namespace core {
namespace scoring {
namespace nv {

class NVlookup {
public:
	NVlookup(std::string filename);

	core::Real get_potentials(char &name, core::Real &score) const;

private:

	std::map<char, numeric::interpolation::spline::SplineGenerator > lookup_table_;

};


} //NV
} //scoring
} //core

#endif
