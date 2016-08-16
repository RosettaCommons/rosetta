// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_minimizer.hh
/// @brief protocols that are specific to RNA_FragmentMover
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_FragmentMover_HH
#define INCLUDED_protocols_rna_RNA_FragmentMover_HH

// Unit headers
#include <protocols/farna/fragments/RNA_Fragments.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <vector>
#include <map>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace protocols::farna::fragments;

namespace protocols {
namespace farna {
namespace movers {

/// @brief The RNA de novo structure modeling protocol
class RNA_FragmentMover: public protocols::moves::Mover {

public:

	/// @brief Construct the protocol object given the RNA fragment library to use.
	RNA_FragmentMover( RNA_Fragments const & all_rna_fragments,
		protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map );

	/// @brief Copy constructor
	RNA_FragmentMover(RNA_FragmentMover const & object_to_copy);

	~RNA_FragmentMover();

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	// virtual protocols::moves::MoverOP fresh_instance() const;

	core::Size
	random_fragment_insertion( core::pose::Pose & pose, Size const & frag_size );

	// is this defunct now? I think so.
	void
	set_frag_size(
		Size const fragment_size
	);

private:

	void
	update_insert_map( core::pose::Pose const & pose );

	RNA_Fragments const & rna_fragments_;
	protocols::toolbox::AtomLevelDomainMapCOP atom_level_domain_map_;

	std::map < Size, Size > insert_map_;
	Size num_insertable_residues_;
	Size insert_map_frag_size_;
	Size frag_size_;

}; // class RNA_FragmentMover


} //movers
} //farna
} //protocols

#endif
