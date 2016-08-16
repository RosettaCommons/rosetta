// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTree.hh
/// @brief  class to hold inverse rotamers associated with theozyme / cstfile
/// @author Florian Richter, flosopher@gmail.com, march 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTree_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTree_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTree.fwd.hh>

// package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers

#include <utility/vector1.fwd.hh>


#if defined PYROSETTA || defined WIN32
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.hh>
#include <core/scoring/constraints/Constraint.hh>
#endif


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


///abstract base class for the inverse rotamer tree
///is abstract so that in the future inverse rots could
///be generated in different ways than throug the enz cst machinery
class InvrotTree : public utility::pointer::ReferenceCount {

public:
	typedef core::Size Size;

	InvrotTree();

	virtual ~InvrotTree();

	core::scoring::constraints::ConstraintCOP
	get_constraint_for_target_state( Size target_state ) const;

	/// @brief the main function, generate the constraints
	void
	generate_inverse_rotamer_constraints(
		core::pose::Pose const & pose,
		AllowedSeqposForGeomCstCOP geomcst_seqpos
	);

	/// @brief convenience access function
	/// to the inverse rotamers in the tree
	/// note: the returned invrots can contain
	/// duplications, as every unique definition
	/// of the tree is returned
	utility::vector1< InvrotCollectorCOP >
	collect_all_inverse_rotamers() const;


	/// @brief visualization
	/// this can lead to several files
	/// being written, depending on the
	/// ambiguities in the tree
	void
	dump_invrots_tree_as_multimodel_pdbs( std::string filename_base) const;


	/// @brief important function 1:
	/// populate the tree
	/// should be called before folding a chain
	/// obviously
	/// see TheozymeInvrotTree below for example
	virtual
	void
	generate_targets_and_inverse_rotamers() = 0;


	/// @brief important function 2:
	/// after a chain/pose has been generated,
	/// check whether it's compatible with
	/// the inverse rotamers in the tree
	/// this can change the pose again
	virtual
	bool
	check_pose_tree_compatibility(
		core::pose::Pose & pose ) const = 0;

	Size
	num_target_states() const {
		return invrot_targets_.size(); }

protected:

	void
	clear_target_states();

	void
	add_to_targets( InvrotTargetOP invrot_target );

private:

	// there can be multiple possible targets, i.e. different rotamers
	// for a ligand
	utility::vector1< InvrotTargetOP > invrot_targets_;
	utility::vector1< core::scoring::constraints::ConstraintCOP > invrot_tree_constraints_;

};

/// @brief This class takes an EnzConstraintIO object
/// and generates a 3D model of the theozyme in it,
/// where inverse rotamers are generated for every block
/// in the cstfile. It gets complicated when ambiguous
/// interactions are specified and other residues are
/// interacting (possibly in ambiguous fashion) with
/// these ambiguous residues...
class TheozymeInvrotTree : public InvrotTree  {
public:

	TheozymeInvrotTree( EnzConstraintIOCOP enzcst_io ); //initialize from enz constraint object, all the inverse rotamers are built

	~TheozymeInvrotTree();

public:


	/// @brief
	/// this function generates the 'target' for the inverse rotamers,
	/// i.e. in the case of enzdes, simply the coordinates for the ligand
	/// if in the future somebody wants to use this to build inverse rotamers
	/// against a protein, this function can be reimplemented.
	/// called from the constructor, is this really the best idea?
	virtual
	void
	generate_targets_and_inverse_rotamers();

	bool
	check_pose_tree_compatibility(
		core::pose::Pose & pose ) const;

private:

	EnzConstraintIOCOP enzcst_io_;

	//maybe this also needs access to all the nodes?

}; // class TheozymeInvrotTree


}
}
}

#endif
