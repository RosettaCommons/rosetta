// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song
/// @author Modified by Vikram K. Mulligan (vmullig@uw.edu)

#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/cyclic_peptide/DeclareBondCreator.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/util.hh>

#include <protocols/loops/loop_closure/kinematic_closure/KinematicMover.hh>
#include <protocols/loops/loop_closure/kinematic_closure/KinematicWrapper.hh>

#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <cstdio>

static THREAD_LOCAL basic::Tracer TR( "protocols.cyclic_peptide.DeclareBond" );

namespace protocols {
namespace cyclic_peptide {

DeclareBond::DeclareBond():
	res1_(0),
	atom1_(""),
	res2_(0),
	atom2_(""),
	add_termini_(false),
	run_kic_(false),
	kic_res1_(0),
	kic_res2_(0),
	rebuild_fold_tree_(false)
{}
DeclareBond::~DeclareBond()= default;

void
DeclareBond::set( core::Size const res1,
	std::string const & atom1,
	core::Size const res2,
	std::string const & atom2,
	bool const add_termini,
	bool const run_kic,
	core::Size const kic_res1,
	core::Size const kic_res2,
	bool const rebuild_fold_tree
)
{
	res1_ = res1;
	atom1_ = atom1;
	res2_ = res2;
	atom2_ = atom2;
	add_termini_ = add_termini;
	run_kic_ = run_kic;
	kic_res1_ = kic_res1;
	kic_res2_ = kic_res2;
	rebuild_fold_tree_ = rebuild_fold_tree;
}

void DeclareBond::apply( core::pose::Pose & pose )
{
	using namespace core::chemical;

	for ( core::Size ir=1; ir<=pose.total_residue() ; ++ir ) {
		if ( pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, ir );
		}
		if ( pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
			core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, ir );
		}
	}

	//printf("Stripping termini.\n"); fflush(stdout); //DELETE ME
	if ( atom1_=="N" && (pose.residue(res1_).type().is_alpha_aa() || pose.residue(res1_).type().is_beta_aa()) && pose.residue(res1_).has_variant_type(LOWER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, LOWER_TERMINUS_VARIANT, res1_);
	}
	if ( atom2_=="N" && (pose.residue(res2_).type().is_alpha_aa() || pose.residue(res2_).type().is_beta_aa()) && pose.residue(res2_).has_variant_type(LOWER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, LOWER_TERMINUS_VARIANT, res2_);
	}
	if ( atom1_=="C" && (pose.residue(res1_).type().is_alpha_aa() || pose.residue(res1_).type().is_beta_aa()) && pose.residue(res1_).has_variant_type(UPPER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, UPPER_TERMINUS_VARIANT, res1_);
	}
	if ( atom2_=="C" && (pose.residue(res2_).type().is_alpha_aa() || pose.residue(res2_).type().is_beta_aa()) && pose.residue(res2_).has_variant_type(UPPER_TERMINUS_VARIANT) ) {
		core::pose::remove_variant_type_from_pose_residue(pose, UPPER_TERMINUS_VARIANT, res2_);
	}

	//printf("Declaring bond.\n"); fflush(stdout); //DELETE ME
	pose.conformation().declare_chemical_bond(res1_, atom1_, res2_, atom2_);

	//Rebuild the polymer bond dependent atoms:
	//printf("Rebuilding bond-dependent atoms.\n"); fflush(stdout); //DELETE ME
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(res1_);
	pose.conformation().rebuild_polymer_bond_dependent_atoms_this_residue_only(res2_);

	if ( rebuild_fold_tree_ ) {
		core::pose::Pose const pose_copy(pose); //Make a reference copy of pose (const to prevent accidentally altering it).

		pose.clear();

		for ( Size ires=1; ires<=pose_copy.total_residue(); ++ires ) {
			if ( ires==1 ) {
				pose.append_residue_by_jump(pose_copy.residue(ires),1);
			} else {
				core::Size anchor_rsd(0);
				core::Size anchor_conid(0);
				core::Size icon=1;
				for ( ; icon<=pose_copy.residue_type(ires).n_possible_residue_connections(); ++icon ) {
					if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) != 0 ) {
						if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) < ires ) {
							anchor_rsd = pose_copy.residue(ires).connected_residue_at_resconn(icon);
							anchor_conid = pose_copy.residue(ires).connect_map(icon).connid();
							break;
						}
					}
				}
				if ( anchor_rsd != 0 ) {
					pose.append_residue_by_bond(pose_copy.residue(ires),false, icon, anchor_rsd, anchor_conid);
				} else {
					if ( pose_copy.chain(ires-1)!=pose_copy.chain(ires) ) { //If this is a new chain, connect this by a jump to residue 1, starting a new chain.
						pose.append_residue_by_jump(pose_copy.residue(ires), 1, "", "", true);
					} else {
						pose.append_residue_by_jump(pose_copy.residue(ires), ires-1);
					}
				}
			}
		}

		// add back all the connections
		for ( Size ires=1; ires<=pose_copy.total_residue(); ++ires ) {
			for ( core::Size icon=1; icon<=pose_copy.residue_type(ires).n_possible_residue_connections(); ++icon ) {
				if ( pose_copy.residue(ires).connected_residue_at_resconn(icon) != 0 ) {
					Size anchor_rsd = pose_copy.residue(ires).connected_residue_at_resconn(icon);
					Size anchor_conid = pose_copy.residue(ires).connect_map(icon).connid();

					if ( pose.residue(ires).connected_residue_at_resconn(icon) == 0 ) {
						//if (pose.residue_type(ires).name3() != "CYS") {
						pose.conformation().set_noncanonical_connection(ires, icon, anchor_rsd, anchor_conid);
						TR << "Adding connection Res " << ires << " to Residue " << anchor_rsd << std::endl;
						//}
					}
				}
			}
		}
	}

	//protocols::loops::add_cutpoint_variants( pose );
	for ( core::Size ir=1; ir<=pose.total_residue() ; ++ir ) {
		if ( pose.residue_type(ir).lower_connect_id() != 0 ) {
			if ( pose.residue(ir).connected_residue_at_resconn(pose.residue_type(ir).lower_connect_id()) == 0 ) {
				if ( !pose.residue(ir).has_variant_type(CUTPOINT_LOWER) ) {
					if ( add_termini_ ) {
						core::pose::add_variant_type_to_pose_residue(pose, LOWER_TERMINUS_VARIANT, ir);
					} else {
						core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_LOWER, ir);
					}
				}
			}
		}
		if ( pose.residue_type(ir).upper_connect_id() != 0 ) {
			if ( pose.residue(ir).connected_residue_at_resconn(pose.residue_type(ir).upper_connect_id()) == 0 ) {
				if ( !pose.residue(ir).has_variant_type(CUTPOINT_UPPER) ) {
					if ( add_termini_ ) {
						core::pose::add_variant_type_to_pose_residue(pose, UPPER_TERMINUS_VARIANT, ir);
					} else {
						core::pose::add_variant_type_to_pose_residue(pose, CUTPOINT_UPPER, ir);
					}
				}
			}
		}
	}

	//pose.fold_tree().show(std::cout);

	//Kinematic closure to build the rest of the peptide:
	if ( run_kic_ ) {
		protocols::loops::loop_closure::kinematic_closure::KinematicMoverOP kinmover( new protocols::loops::loop_closure::kinematic_closure::KinematicMover );
		core::scoring::ScoreFunctionOP sfxn;
		sfxn = core::scoring::get_score_function();
		if ( kic_res1_ == 0 ) {
			kic_res1_ = 1;
		}
		if ( kic_res2_ == 0 ) {
			kic_res2_ = pose.total_residue();
		}

		kinmover->set_temperature( 1.0 );
		kinmover->set_vary_bondangles( false );
		kinmover->set_sample_nonpivot_torsions( true );
		kinmover->set_rama_check( true );
		kinmover->set_idealize_loop_first( true );
		kinmover->set_sfxn(sfxn);
		kinmover->set_pivots(kic_res1_, (int)( ( (double)(kic_res1_+kic_res2_) )/2.0 ), kic_res2_);

		pose.update_residue_neighbors();
		kinmover->apply(pose);
	}
}

/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
DeclareBond::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const &
)
{
	res1_ = tag->getOption< Size >( "res1" );
	atom1_ = tag->getOption< std::string >( "atom1" );
	res2_ = tag->getOption< Size >( "res2" );
	atom2_ = tag->getOption< std::string >( "atom2" );
	add_termini_ = tag->getOption< bool >( "add_termini", true );
	rebuild_fold_tree_ = tag->getOption< bool >( "rebuild_fold_tree", false );
	run_kic_ = tag->getOption< bool >( "run_KIC", false);
	kic_res1_ = tag->getOption< Size >( "KIC_res1", 0);
	kic_res2_ = tag->getOption< Size >( "KIC_res2", 0);
}

moves::MoverOP DeclareBond::clone() const { return moves::MoverOP( new DeclareBond( *this ) ); }
moves::MoverOP DeclareBond::fresh_instance() const { return moves::MoverOP( new DeclareBond ); }

protocols::moves::MoverOP
DeclareBondCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeclareBond );
}

std::string
DeclareBondCreator::keyname() const
{
	return DeclareBondCreator::mover_name();
}

std::string
DeclareBondCreator::mover_name()
{
	return "DeclareBond";
}

std::string
DeclareBond::get_name() const {
	return "DeclareBond";
}

} // moves
} // protocols
