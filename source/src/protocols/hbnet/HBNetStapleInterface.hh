// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file <protocols/hbnet/HBNetStapleInterface.hh
/// @brief inherits from HBNet; protocol for designing h-bond networks at Prot/Prot interfaces;
/// @author Scott Boyken (sboyken@gmail.com)

#ifndef INCLUDED_protocols_hbnet_HBNetStapleInterface_hh
#define INCLUDED_protocols_hbnet_HBNetStapleInterface_hh

// HBNet INCLUDES
#include <protocols/hbnet/HBNet.hh>
//#include <protocols/hbnet/HBNet.fwd.hh>
#include <protocols/hbnet/HBNetStapleInterface.fwd.hh>

// BASIC INCLUDES
#include <basic/datacache/DataMap.hh>

// UTILTIY INCLUDES
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// CORE INCLUDES
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AtomType.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>

// PROTOCOLS INCLUDES
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/upstream/ProteinUpstreamBuilder.hh>
#include <protocols/match/upstream/ScaffoldBuildPoint.hh>
#include <protocols/match/upstream/UpstreamBuilder.fwd.hh>

namespace protocols {
namespace hbnet {

class HBNetStapleInterface : public HBNet
{

public:
	//constructors
	HBNetStapleInterface();
	HBNetStapleInterface( std::string const name );
	//Constructor from code
	HBNetStapleInterface( core::scoring::ScoreFunctionCOP scorefxn,
		core::Size max_unsat,
		core::Size min_network_size=3,
		core::Real hb_threshold=-0.5,
		core::Size max_network_size=15,
		std::string des_residues="STRKHYWNQDE",
		bool find_native=false,
		bool only_native=false,
		bool keep_existing=false,
		bool extend_existing=false,
		bool only_extend=false
	);

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	//destructor
	~HBNetStapleInterface();

	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;
	void setup_packer_task_and_starting_residues( core::pose::Pose const & pose ) override;
	void prepare_output() override;
	std::string print_additional_info_for_net( hbond_net_struct & i, core::pose::Pose const & pose ) override;
	std::string print_additional_headers() override;
	bool network_meets_initial_criteria( hbond_net_struct const & network ) override;
	// final criteria that reuqires network rotamers placed on pose
	bool network_meets_final_criteria( core::pose::Pose const & pose, hbond_net_struct & network ) override;
	bool state_is_starting_aa_type( core::Size const res, core::Size const rot_id ) override;
	bool pair_meets_starting_criteria( core::Size const res1, core::Size const rot1, core::Size const res2, core::Size const rot2 ) override;
	core::Real scale_twobody_energy( core::Real input_twobody_energy, char res1, char res2 ) override;

	bool same_helix(utility::vector1< std::pair<core::Size,core::Size> > const helix_boundaries, core::Size const r1, core::Size const r2);
	bool interhelical_contact(utility::vector1< std::pair<core::Size,core::Size> > const helix_boundaries, core::Size const r1, core::Size const r2, core::pose::Pose const & pose);
	core::Size get_helix_id( core::Size r1 );

	//void rec_add_staple( std::vector< HBondNetStructOP >::const_iterator netit, HBondNetStructOP new_network, core::Size staple_count );
	void rec_add_staple( std::vector< HBondNetStructOP >::const_iterator netit, std::set< core::Size > net_ids, core::Size staple_count );

	bool network_spans_all_helices( hbond_net_struct const & i );
	core::Size num_helices_w_hbond( hbond_net_struct const & i );
	core::Size num_helices_w_hbond( utility::vector1< HBondResStructCOP > const & residues );
	bool has_pH_His( core::pose::Pose const & pose, hbond_net_struct & i );
	bool residues_are_interface_pairs( core::Size const res1, core::Size const res2 );

	core::Size num_intermolecular_hbonds( hbond_net_struct & i, core::pose::Pose const & pose );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//core::Size symm_num_intermolecular_hbonds( hbond_net_struct & i, core::pose::Pose & pose );

private:
	bool all_helices_;
	bool span_all_helices_;
	bool only_symm_interfaces_;
	bool allow_onebody_networks_;
	bool his_tyr_;
	bool pH_His_;
	bool boundary_his_must_to_hbond_pos_charge_;
	bool only_start_at_interface_pairs_;
	bool use_aa_dependent_weights_;
	core::Size runcount_;
	core::Size min_networks_per_pose_;
	core::Size max_networks_per_pose_;
	core::Size combos_;
	core::Size min_intermolecular_hbonds_;
	core::Size min_helices_contacted_by_network_;
	core::Size min_asp_glu_hbonds_;
	core::Real interf_distance_;
	utility::vector1< core::Size > jump_nums_;
	utility::vector1< std::pair< core::Size, core::Size > > helix_boundaries_; //for helical bundles
	utility::vector1< std::list< Size > > pair_lists_vec_;
};
// end HBNetStapleInterface

} // hbnet namespace
} // protocols namespace

#endif
