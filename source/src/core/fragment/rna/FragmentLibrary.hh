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

#ifndef INCLUDED_core_fragment_rna_FragmentLibrary_HH
#define INCLUDED_core_fragment_rna_FragmentLibrary_HH

#include <core/fragment/rna/FragmentLibrary.fwd.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <core/types.hh>
#include <vector>
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace fragment {
namespace rna {

class FragmentLibrary : public utility::pointer::ReferenceCount  {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~FragmentLibrary();

	core::Real get_fragment_torsion(
		core::Size const num_torsion,
		Size const which_frag,
		core::Size const offset );

	TorsionSet const & get_fragment_torsion_set( core::Size const which_frag ) const;

	void  add_torsion( TorsionSet const & torsion_set );

	void  add_torsion(
		FullAtomRNA_Fragments const & vall,
		core::Size const position,
		core::Size const size
	);

	core::Size get_align_depth() const;

private:
	std::vector< TorsionSet > align_torsions_;

};


} //fragments
} //denovo
} //protocols

#endif

