// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author  watkins

#ifndef INCLUDED_protocols_farna_fragments_TorsionSet_HH
#define INCLUDED_protocols_farna_fragments_TorsionSet_HH

#include <protocols/farna/fragments/TorsionSet.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace farna {
namespace fragments {

class TorsionSet {
public:
	TorsionSet & operator =( TorsionSet const & src );
	TorsionSet( core::Size const size, core::Size const position );

	//make this private?
	ObjexxFCL::FArray2D <core::Real>   torsions;  // dimensions: (NUM_RNA_TORSIONS, SRange(0, size) );
	ObjexxFCL::FArray1D <std::string> torsion_source_name;  // dimensions: ( SRange(0, size) );
	ObjexxFCL::FArray1D <char>   secstruct;

	ObjexxFCL::FArray3D <core::Real>   non_main_chain_sugar_coords;
	bool  non_main_chain_sugar_coords_defined;

	inline
	core::Size get_size() const { return size_; }

	inline
	core::Size get_index_in_vall() const { return index_in_vall_; }

private:
	core::Size size_;
	core::Size index_in_vall_;

};

}
}
}

#endif

