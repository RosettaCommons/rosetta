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
/// @author James Thompson
/// @author Srivatsan Raman


#ifndef INCLUDED_protocols_rna_RNA_FragmentMover_hh
#define INCLUDED_protocols_rna_RNA_FragmentMover_hh

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <protocols/rna/RNA_FragmentsClasses.fwd.hh>
#include <protocols/rna/RNA_FragmentsClasses.hh>

//Oooh.
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
#include <vector>

//Auto Headers
#include <iostream>


namespace protocols {
namespace rna {

/// @brief The RNA de novo structure modeling protocol
class RNA_FragmentMover: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_FragmentMover( RNA_FragmentsOP & all_rna_fragments, FArray1D_bool const & allow_insert );

	~RNA_FragmentMover();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const {
		return new RNA_FragmentMover(*this);
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void set_frag_size( Size const fragment_size );

	void
	insert_fragment(
		core::pose::Pose & pose,
		Size const position,
		rna::TorsionSet const torsion_set
	);

	//Hmmm, this probably reproduces code.
	void
	update_insert_map( core::pose::Pose const & pose );


	void
	random_fragment_insertion( core::pose::Pose & pose, Size const & frag_size );
	// undefined, commentin out to make PyRosetta compile
	// void initialize_allow_insert( core::pose::Pose & pose );

private:

	RNA_FragmentsOP all_rna_fragments_;

	//Hmmm, this probably reproduces code.
	// Move to some kind of RNA Fragment Mover class?
	FArray1D <bool> const allow_insert_;
	std::map < Size, Size > insert_map_;
	Size num_insertable_residues_;
	Size insert_map_frag_size_;

	Size frag_size_;

}; // class RNA_FragmentMover



} //rna
} // protocols

#endif
