// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/enzdes/EnzdesRemodelProtocol.hh
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, march 2009


#ifndef INCLUDED_devel_enzdes_EnzdesRemodelProtocol_hh
#define INCLUDED_devel_enzdes_EnzdesRemodelProtocol_hh


#include <protocols/enzdes/EnzdesFlexBBProtocol.hh>

#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>
#include <protocols/forge/remodel/ResidueVicinityRCG.hh>

// utility includes
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>

namespace devel{
namespace enzdes{

class EnzdesRemodelProtocol;
typedef utility::pointer::shared_ptr< EnzdesRemodelProtocol > EnzdesRemodelProtocolOP;
typedef utility::pointer::weak_ptr< EnzdesRemodelProtocol const > EnzdesRemodelProtocolCAP;


class EnzdesRemodelProtocol : public protocols::enzdes::EnzdesFlexBBProtocol
{

public:
	EnzdesRemodelProtocol();

	~EnzdesRemodelProtocol();

	void apply( core::pose::Pose & pose);
	virtual std::string get_name() const;

private:

}; //class EnzdesRemodelProtocol

class EnzdesRemodelMover;
typedef utility::pointer::shared_ptr< EnzdesRemodelMover > EnzdesRemodelMoverOP;
typedef utility::pointer::weak_ptr< EnzdesRemodelMover const > EnzdesRemodelMoverCAP;

class EnzdesRemodelMover : public protocols::moves::Mover
{

public:

	//constructor / destructor
	EnzdesRemodelMover(); //empty constructor for parser stuff

	EnzdesRemodelMover( EnzdesRemodelMover const & other ); //copy constructor

	EnzdesRemodelMover(
		protocols::enzdes::EnzdesFlexBBProtocolOP enz_prot,
		core::pack::task::PackerTaskCOP task,
		protocols::enzdes::EnzdesFlexibleRegionCOP flex_region
	);

	~EnzdesRemodelMover();

	protocols::moves::MoverOP
	clone() const;

	void apply( core::pose::Pose & pose);
	virtual std::string get_name() const;

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & datamap,
		Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		Pose const & pose
	);

	void set_task( core::pack::task::PackerTaskCOP task );

	void set_user_provided_ss( utility::vector1< std::string > const & user_ss );

	protocols::filters::FilterCOP
	setup_packer_neighbor_graph_filter( core::pose::Pose const & pose ) const;

	std::string
	generate_secstruct_string( core::pose::Pose & pose) const;

	/// @brief puts in a random ligand conformer that has a lower score with
	/// the remodel region than the conformer currently in the pose
	void
	apply_random_lowE_ligconf(
		core::pose::Pose & pose ) const;

	/// @brief gets the atom indeces corresponding
	/// to the given atom_names in the given restype_set
	static
	void
	translate_atomnames_to_restype_set_atomids(
		core::pose::Pose const & pose,
		core::chemical::ResidueTypeSetCAP restype_set,
		core::Size const seqpos,
		utility::vector1< std::string > const & atom_names,
		utility::vector1< core::Size > & atom_ids
	);

	/// @brief translates information from enzdes loops file to
	/// rvinfos
	static
	protocols::forge::remodel::ResidueVicinityInfoOP
	translate_res_interactions_to_rvinfos(
		core::pose::Pose const & pose,
		core::Size const targ_pos,
		core::Size const example_loop_seqpos,
		core::chemical::ResidueTypeSetCAP restype_set,
		protocols::toolbox::match_enzdes_util::ResInteractions const & res_ints
	);

  void
  set_seq_mapping( core::id::SequenceMappingOP seq_map_in);

	void
	set_keep_existing_aa_identities( bool setting );

  core::id::SequenceMappingCOP get_seq_mapping() const;

	utility::vector1< core::Size > get_flex_region( ) const;

	void set_target_inverse_rotamers(utility::vector1< std::list < core::conformation::ResidueCOP > > & inv_rot);

	void
	set_max_allowed_score_increase( core::Real sc_increase );

protected:

	/// @brief function to do proper setup if this mover is
	/// called from the parser
	void
	initialize( core::pose::Pose & pose );

	bool
	remodel_pose(
		core::pose::Pose & pose,
		std::string secstruct
	);


	/// @ do fa minimization and some kinematic moves
	/// to get a little bit of fa refinement
	/// remodel constraints will be on in this stage
	bool
	refine_pose(
		core::pose::Pose & pose
	);

	/// @brief function to set up filters according to the initial conformation
	/// @brief of a region to be remodeled
	void
	examine_initial_conformation(
		core::pose::Pose & pose
	);

	void
	setup_cached_observers(
		core::pose::Pose &// pose
	);

	void
	remove_cached_observers(
		core::pose::Pose & pose
	);

	void
	process_length_change(
		core::pose::Pose & pose,
		core::id::SequenceMappingCOP smap
	);

	void
	setup_postdesign_filters(
		core::pose::Pose const & //pose
	);


	/// @brief setup remodel constraint generators
	void
	setup_rcgs(
		core::pose::Pose & pose,
		protocols::forge::components::VarLengthBuild & vlb
	);

	/// @brief setup remodel constraint generators to hit
	/// inverse rotamer backbones
	void
	setup_rcgs_from_inverse_rotamers(
		protocols::forge::components::VarLengthBuild & vlb
	);

	/// @brief generate inverse rotamers according to missing
	/// catalytic interactions
	void
	create_target_inverse_rotamers(
		core::pose::Pose & pose );

	/// @brief very experimental for now
	/// @brief if there are missing interactions in the pose
	/// @brief try to find them with the secondary matcher.
	/// @brief returns true if found, or false otherwise
	bool
	secmatch_after_remodel(
		core::pose::Pose & pose );


private:

	protocols::enzdes::EnzdesFlexBBProtocolOP enz_prot_;
	core::pack::task::PackerTaskCOP orig_task_;
	protocols::enzdes::EnzdesFlexibleRegionCOP flex_region_;
	core::scoring::ScoreFunctionOP vlb_sfxn_;

	std::set< core::Size > remodel_positions_;
	std::set< core::Size > other_design_positions_;
	std::set< core::Size > other_repack_positions_;

	//how many times to remodel this loop
	core::Size remodel_trials_;

	//are we trying to introduce new specific interactions?
	bool remodel_secmatch_;

	//try to make the foldtree at the end of apply
	//look like the input foldtree
	bool reinstate_initial_foldtree_;

	//keep residues not getting added at their
	//original identities instead of making them
	//ala. somewhat arbitrary implementation, still experimental
	//plus secmatch might still change them
	bool keep_existing_aa_identities_;
	utility::vector1< core::chemical::ResidueTypeCAP > init_aa_;

	//which of the regions marked as flexible to remodel
	core::Size region_to_remodel_;

	//filters to be applied after loop closure but before design, to determine whether
	//a newly generated loop is worth designing
	protocols::filters::FilterCollectionOP predesign_filters_;

	//filters to be applied after after design, to determine whether
	//a newly generated loop is better than the original conformation
	protocols::filters::FilterCollectionOP postdesign_filters_;


	//sequence mapping to keep track of the various insertions/deletions
	//that happen in the course of a run
	core::id::SequenceMappingOP start_to_current_smap_;

	//target backbones for remodel secondary matching
	utility::vector1< std::list< core::conformation::ResidueCOP > > target_inverse_rotamers_;

	//should the existing conformation of an inverse rotamer residue be
	//included in the inverse rotamers for that residue?
	bool include_existing_conf_as_invrot_target_;

	//positions of catalytic res that are not being remodeled
	//but are part of an interaction where the other residue
	//is being remodeled
	utility::vector1< core::Size > non_remodel_match_pos_;

	utility::vector1< protocols::forge::remodel::RemodelConstraintGeneratorOP > rcgs_;

	core::Real ss_similarity_probability_;
	core::Real max_allowed_score_increase_; //how much additional score are we allowing the remodeled pose?
	utility::vector1< std::string > user_provided_ss_;

}; //class EnzdesRemodelMover

} //namespace enzdes
} //namespace protocols


#endif // INCLUDED_protocols_enzdes_EnzdesFixBBProtocol_HH
