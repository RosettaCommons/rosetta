// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_insertion/FusePosesNtoCMover.hh
/// @brief  hh file for FusePosesNtoCMover
/// @author Florian Richter, flosopher@gmail.com, february 2013

#ifndef INCLUDED_devel_domain_insertion_FusePosesNtoCMover_hh
#define INCLUDED_devel_domain_insertion_FusePosesNtoCMover_hh

// Unit headers
#include <devel/domain_insertion/FusePosesNtoCMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// package headers

// Project headers
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


// Utility headers
//#include <utility/pointer/ReferenceCount.hh>

// C++ headers

//#include <utility/vector1.fwd.hh>


namespace devel {
namespace domain_insertion {


/// @brief  mover to stitch two poses together, with the added caveat
///         that these can be multimeric
///         finds the best superposition given residue ranges for the
///         two poses (based on calpha rmsd and simple clash check)
///         then stitches that superposition and calls a relax mover
class FusePosesNtoCMover : public protocols::moves::Mover {

public:

	typedef protocols::moves::Mover parent;

	typedef core::Size Size;
	typedef core::Real Real;

	FusePosesNtoCMover();

	FusePosesNtoCMover( FusePosesNtoCMover const & other );

	~FusePosesNtoCMover() override;

	protocols::moves::MoverOP
	clone() const override;

	
	std::string
	get_name() const override;

	
	void
	apply( core::pose::Pose & pose ) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;


	/// @brief function to set up a fold tree
	/// that insulates movement to the regions
	/// around the fusepoint
	void
	setup_relax_fold_tree(
		core::pose::Pose & pose
	);


	/// @brief function to setup a relax mover
	/// that relaxes the stitched together pose
	/// alternatively this mover could be supplied
	/// through RS, but it probably needs to know
	/// something about how the pose got fused in
	/// the first place, not sure how to best communicate
	/// this in RS
	protocols::moves::MoverOP
	generate_default_relax_mover(
		core::pose::Pose const & pose
	) const;


	/// @brief converts input given through tag
	/// into actual seqpos
	void
	setup_candidate_pdbpos(
		utility::vector1< Size > & pdbpos,
		std::string const tag_input
	) const;

	/// @brief figures out how to superimpose the
	/// two poses given the positions and the superpose_window_ parameter
	core::id::AtomID_Map< core::id::AtomID >
	generate_superposition_map(
		core::pose::Pose const & pose,
		core::pose::Pose const & fuse_pose,
		Size const pose_pdbpos,
		Size const fuse_pose_pdbpos
	) const;

	/// @brief counts the number of hardsphere clashes given
	/// two poses and where they should be fused
	Size
	count_hardsphere_clashes(
		core::pose::Pose const & pose,
		core::pose::Pose const & fuse_pose
	) const;


	/// @brief stitch fuse_pose onto the c-terminus of pose
	void
	fuse_poses(
		core::pose::Pose & pose,
		core::pose::Pose const & fuse_pose
	);


	/// @brief helper function to shave off
	/// residues from the pose that are outside
	/// of the fusion
	core::pose::PoseOP
	truncate_pose_at_fusion_site(
		core::pose::Pose const & pose,
		Size const fusion_pdbpos,
		bool nterm_truncation  //whether to shave off residues upstream or downstream of the fusion pos
	) const;

	/// @brief
	/// analyze the fusion with respect to a few things
	void
	analyze_fused_pose(
		core::pose::Pose & nterm_half,
		core::pose::Pose & cterm_half,
		core::pose::Pose & fusion_raw,
		core::pose::Pose & fusion_rlx,
		std::string const identifier
	) const;


private:

	/// @brief small helper function to navigate
	/// between pdb ids and conformation ids
	void
	setup_chain_mappings(
		core::pose::Pose const & pose );

private:

	//the pose that gets fused to the one coming in through apply
	//convention: this will become the c-terminal end of the output pose
	core::pose::PoseOP fuse_pose_;

	core::scoring::ScoreFunctionOP sfxn_, sfxn_nocb_; //nocb means no chainbreak score term

	//in case the poses to be fused are multimeric
	utility::vector1< char > chains_to_use_; //chains in pdb lettering, from tag

	//which chains (in conformation numbering) are supposed to be fused
	//note pair.first is for fuse_pose
	utility::vector1< std::pair< Size, Size > > chain_mappings_;

	//the (absolute) seqpos where fusions have been done,
	//should have same size as chain_mappings_ vector above
	utility::vector1< Size > fuse_positions_;

	//the regions (in terms of upper and lower seqpos) around the fusion,
	//for each chain / fuseposition
	//used for two things: 1. what's being treated flexible during relaxing
	//2. what's being compared between fused_pose_raw and fused_pose
	utility::vector1< std::pair< Size, Size > > fusion_regions_;

	//which jumps to move during default relax mover
	utility::vector1< Size > jumps_to_move_;

	//the positions in both poses considered for the fusion
	//note: these are pdbpos because depending on how many chains are used,
	//they need to be mapped to several seqpos
	utility::vector1< Size > candidate_pdbpos_pose_, candidate_pdbpos_fuse_pose_;

	//a few parameters regarding superposition of the two
	//poses
	Size superpose_window_;  //what span to do the superposition on
	Size relax_window_; //how much beyond the SS element that the superposition happened is movable during default relax
	Size max_number_allowed_clashes_; //how many hardsphere clashes to allow for the pose to be relaxed
	Real hardsphere_clash_limit_; //what counts as a clash
	Real max_allowed_rms_;  //the maximum allowed superposition rms allowed for the pose to be relaxed

	//raw numbers read in from tag, converted to real pose seqpos
	//when apply is called
	std::string nterm_span_tag_;

	// the mover to be called to relax the raw
	// fused pose
	protocols::moves::MoverOP relax_mover_;
	bool rs_specified_relax_mover_;

	//useful for developing
	bool debugmode_;


};


/// @brief a mover that sets up the fold tree such that
/// movement in a coiled coil gets propagated downstream
/// of the coil without changing the downstream conformation,
/// i.e. the pseudo rigid body dofs of the downstream protein portion
/// change but not the conformation
class SetupCoiledCoilFoldTreeMover : public protocols::moves::Mover {

public:

	SetupCoiledCoilFoldTreeMover();

	SetupCoiledCoilFoldTreeMover( SetupCoiledCoilFoldTreeMover const & other );

	~SetupCoiledCoilFoldTreeMover() override;

	
	std::string
	get_name() const override;

	protocols::moves::MoverOP
	clone() const override;

	void
	apply( core::pose::Pose & pose ) override;

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & ) override;


private:

	utility::vector1< char > coiled_coil_chains_; //chains in pdb lettering, from tag
	core::Size chain2_cutpos_; //the fold tree will have one cut, this is supposed to be pdb numbering
	bool add_chainbreak_variants_; //whether to change the residue types in addition to changing the fold tree

};

/// @helper function to compare two pairs a certain way, necessary to sort above chain mappings
bool
custom_pair_compare(
	std::pair< FusePosesNtoCMover::Size, FusePosesNtoCMover::Size > const & pair1,
	std::pair< FusePosesNtoCMover::Size, FusePosesNtoCMover::Size > const & pair2
)
{
	return (pair1.second < pair2.second );
}


/// @brief helper function to keep code readable
void
set_superpose_residues_in_atom_map(
	core::id::AtomID_Map< core::id::AtomID > & atom_map,
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	core::Size const seqpos1,
	core::Size const seqpos2,
	bool debugout
);


/// @brief
/// %^#$ annoying helper function because the rmsd code wants two different data structures...
std::map< core::id::AtomID, core::id::AtomID >
convert_AtomID_Map_to_std_map(
	core::id::AtomID_Map< core::id::AtomID > const & atom_map );


}
}

#endif
