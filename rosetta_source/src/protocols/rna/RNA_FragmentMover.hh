// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_FragmentMover_HH
#define INCLUDED_protocols_rna_RNA_FragmentMover_HH

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rna/RNA_Fragments.fwd.hh>
#include <protocols/rna/AllowInsert.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>
#include <vector>
#include <map>

//Auto Headers
#include <iostream>

namespace protocols {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_FragmentMover: public protocols::moves::Mover {

public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_FragmentMover( RNA_FragmentsOP all_rna_fragments,
										 protocols::rna::AllowInsertOP allow_insert );

	// is this defunct now? I think so.
	RNA_FragmentMover( RNA_FragmentsOP all_rna_fragments,
										 ObjexxFCL::FArray1D<bool> const & allow_insert,
										 core::pose::Pose const & pose );


	~RNA_FragmentMover();

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

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

	protocols::rna::RNA_FragmentsOP rna_fragments_;
	protocols::rna::AllowInsertOP allow_insert_;

	std::map < Size, Size > insert_map_;
	Size num_insertable_residues_;
	Size insert_map_frag_size_;
	Size frag_size_;

}; // class RNA_FragmentMover


} //rna
} // protocols

#endif
