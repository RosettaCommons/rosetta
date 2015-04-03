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

#include <core/scoring/nv/NVlookup.fwd.hh>

//project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//STL header
#include <map>
#include <vector>


namespace core {
namespace scoring {
namespace nv {

class NVlookup : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~NVlookup();
	NVlookup(std::string filename);

	core::Real get_potentials(core::chemical::AA const & aa_type, core::Real const & score) const;

private:

	void set_up_spline_from_data(core::chemical::AA const & aa_type, utility::vector1<core::Real> const & bin_centers, utility::vector1<core::Real> const & data);

	//The spline for each amino acid is stored in an array indexed by core::chemical::AA enum
	//numeric::interpolation::spline::InterpolatorOP lookup_table_[core::chemical::num_canonical_aas];

	utility::vector1<numeric::interpolation::spline::InterpolatorOP> lookup_table_;

	//std::map<char, numeric::interpolation::spline::InterpolatorOP > lookup_table_;

};


} //NV
} //scoring
} //core

#endif
