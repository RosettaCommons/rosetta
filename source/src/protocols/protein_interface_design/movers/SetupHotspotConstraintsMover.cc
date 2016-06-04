// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/SetupHotspotConstraintsMover.cc
/// @brief
/// @author Jacob Corn (jecorn@u.washington.edu), Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMover.hh>
#include <protocols/protein_interface_design/movers/SetupHotspotConstraintsMoverCreator.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/hotspot_hashing/HotspotStubSet.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <complex>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace protocols::moves;

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.movers.SetupHotspotConstraintsMover" );

std::string
SetupHotspotConstraintsMoverCreator::keyname() const
{
	return SetupHotspotConstraintsMoverCreator::mover_name();
}

protocols::moves::MoverOP
SetupHotspotConstraintsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SetupHotspotConstraintsMover );
}

std::string
SetupHotspotConstraintsMoverCreator::mover_name()
{
	return "SetupHotspotConstraints";
}

SetupHotspotConstraintsMover::SetupHotspotConstraintsMover() :
	protocols::moves::Mover( "SetupHotspotConstraintsMover" ),
	chain_to_design_( 2 ),
	CB_force_constant_( 1.0 ),
	worst_allowed_stub_bonus_( 0.0 ),
	apply_self_energies_( true ),
	bump_cutoff_( 4.0 ),
	apply_ambiguous_constraints_( true ),
	colonyE_( false ),
	stub_energy_fxn_( "backbone_stub_constraint" )
{}

protocols::moves::MoverOP
SetupHotspotConstraintsMover::clone() const{
	return protocols::moves::MoverOP( new SetupHotspotConstraintsMover( *this ) );
}

protocols::moves::MoverOP
SetupHotspotConstraintsMover::fresh_instance() const{
	return protocols::moves::MoverOP( new SetupHotspotConstraintsMover() );
}

SetupHotspotConstraintsMover::SetupHotspotConstraintsMover(
	protocols::hotspot_hashing::HotspotStubSetCOP hotspot_stub_set,
	core::Size const chain_to_design,
	core::Real const & CB_force_constant,
	core::Real const & worst_allowed_stub_bonus,
	bool const apply_self_energies,
	core::Real const & bump_cutoff,
	bool const apply_ambiguous_constraints,
	bool const colonyE,
	std::string const & stub_energy_fxn
) :
	protocols::moves::Mover( "SetupHotspotConstraintMover" ),
	chain_to_design_(chain_to_design),
	CB_force_constant_(CB_force_constant),
	worst_allowed_stub_bonus_(worst_allowed_stub_bonus),
	apply_self_energies_(apply_self_energies),
	bump_cutoff_(bump_cutoff),
	apply_ambiguous_constraints_(apply_ambiguous_constraints),
	colonyE_( colonyE),
	stub_energy_fxn_( stub_energy_fxn)
{
	//  packer_task_ = packer_task->clone();
	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *hotspot_stub_set ) );
}

SetupHotspotConstraintsMover::SetupHotspotConstraintsMover( SetupHotspotConstraintsMover const & init ) :
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( init ),
	chain_to_design_( init.chain_to_design_),
	CB_force_constant_(init.CB_force_constant_),
	worst_allowed_stub_bonus_(init.worst_allowed_stub_bonus_),
	apply_self_energies_(init.apply_self_energies_),
	bump_cutoff_(init.bump_cutoff_),
	apply_ambiguous_constraints_(init.apply_ambiguous_constraints_),
	colonyE_( init.colonyE_ ),
	stub_energy_fxn_( init.stub_energy_fxn_ )
{
	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new protocols::hotspot_hashing::HotspotStubSet( *init.hotspot_stub_set_ ) );
}

void
SetupHotspotConstraintsMover::apply( core::pose::Pose & pose ) {
	if ( colonyE_ ) {
		protocols::hotspot_hashing::HotspotStubSetOP colonyE_set = hotspot_stub_set_->colonyE();
		hotspot_stub_set_ = colonyE_set;
	}
	if ( std::abs(CB_force_constant_) > 1E-9 ) {
		if ( stub_energy_fxn_ == "backbone_stub_constraint" ) {
			hotspot_stub_set_->add_hotspot_constraints_to_pose( pose, chain_to_design_, hotspot_stub_set_,
				CB_force_constant_, worst_allowed_stub_bonus_, apply_self_energies_, bump_cutoff_, apply_ambiguous_constraints_ );
		} else if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
			hotspot_stub_set_->add_hotspot_constraints_to_wholepose( pose, chain_to_design_, hotspot_stub_set_,
				CB_force_constant_, worst_allowed_stub_bonus_, apply_self_energies_, bump_cutoff_, apply_ambiguous_constraints_ );
		}
	} else {
		core::scoring::constraints::ConstraintSetOP empty_constraint_set( new core::scoring::constraints::ConstraintSet );
		pose.constraint_set( empty_constraint_set );
	}
}

std::string
SetupHotspotConstraintsMover::get_name() const {
	return "SetupHotspotConstraintsMover";
}

/// This needs to be parsed before all other movers b/c it changes scorefxns
void
SetupHotspotConstraintsMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	using core::Real;
	chain_to_design_ = tag->getOption<Size>( "redesign_chain", 2 );

	CB_force_constant_ = tag->getOption<Real>( "cb_force", 0.5 );
	worst_allowed_stub_bonus_ = tag->getOption<Real>( "worst_allowed_stub_bonus", 0 );
	apply_self_energies_ = tag->getOption<bool>( "apply_stub_self_energies", 0 );
	bump_cutoff_ = tag->getOption<Real>( "apply_stub_bump_cutoff", 10. );
	apply_ambiguous_constraints_ = tag->getOption<bool>( "pick_best_energy_constraint", 1 );
	core::Real const bb_stub_cst_weight( tag->getOption< core::Real >( "backbone_stub_constraint_weight", 1.0 ) );
	stub_energy_fxn_ = tag->getOption<std::string>( "stubscorefxn", "backbone_stub_constraint" ) ;

	colonyE_ = tag->getOption<bool>( "colonyE", 0 );

	hotspot_stub_set_ = protocols::hotspot_hashing::HotspotStubSetOP( new hotspot_hashing::HotspotStubSet );
	if ( tag->hasOption( "stubfile" ) ) {
		std::string const hotspot_fname( tag->getOption<std::string>( "stubfile", "stubs.pdb" ) );
		hotspot_stub_set_->read_data( hotspot_fname );
	}
	utility::vector1< TagCOP > const branch_tags( tag->getTags() );
	BOOST_FOREACH ( TagCOP const curr_tag, branch_tags ) {
		if ( curr_tag->getName() == "HotspotFiles" ) {
			utility::vector1< TagCOP > const branch_tags2( curr_tag->getTags() );
			BOOST_FOREACH ( TagCOP const curr_tag2, branch_tags2 ) {
				std::string const file_name( curr_tag2->getOption< std::string >( "file_name" ) );
				std::string const nickname( curr_tag2->getOption< std::string >( "nickname" ) );
				core::Size const stub_num( curr_tag2->getOption< core::Size >( "stub_num", 100000 ) );
				hotspot_hashing::HotspotStubSetOP temp_stubset( new hotspot_hashing::HotspotStubSet );
				temp_stubset->read_data( file_name );
				temp_stubset->remove_random_stubs_from_set( temp_stubset->size() - stub_num );
				hotspot_stub_set_->add_stub_set( *temp_stubset );
				TR<<"Read stubset from file "<<file_name<<" and associating it with name "<<nickname<<std::endl;
				TR<<stub_num<<" stubs kept in memory"<<std::endl;
				data.add( "hotspot_library", nickname, temp_stubset );
			}
		} else {
			utility_exit_with_message( curr_tag->getName() + " not recognized by SetupHotspotConstraints, did you mean HotspotFiles?" );
		}
	}

	TR<<"applying " <<  stub_energy_fxn_ <<" constraints to pose with " << " cb_force weight of "<<CB_force_constant_<<", apply ambiguous constraints set to "<<apply_ambiguous_constraints_<< " and colonyE set to " << colonyE_ << std::endl;
	data.add( "constraints" , "hotspot_stubset", hotspot_stub_set_ );

	for ( std::map< std::string, utility::pointer::ReferenceCountOP >::const_iterator it = (data)[ "scorefxns" ].begin(); it!=(data)[ "scorefxns" ].end(); ++it ) {
		using namespace core::scoring;
		ScoreFunctionOP scorefxn( data.get_ptr< ScoreFunction >( "scorefxns", it->first) );
		if ( stub_energy_fxn_ == "backbone_stub_constraint" ) {
			core::Real const weight( scorefxn->get_weight( backbone_stub_constraint ) );
			if ( weight == 0.0 ) {
				scorefxn->set_weight( backbone_stub_constraint, bb_stub_cst_weight );
				TR<<"Setting bacbkone_stub_constraint weight in scorefxn "<<it->first<<" to "<<bb_stub_cst_weight<<std::endl;
			} else {
				TR<<"Skipping resetting of backbone_stub_constraint weight in "<<it->first<<" which is already preset to "<<weight<<std::endl;
			}
		} else if ( stub_energy_fxn_ == "backbone_stub_linear_constraint" ) {
			core::Real const weight( scorefxn->get_weight( backbone_stub_linear_constraint ) );
			if ( weight == 0.0 ) {
				scorefxn->set_weight( backbone_stub_linear_constraint, bb_stub_cst_weight );
				TR<<"Setting backbone_stub_linear_constraint weight in scorefxn "<<it->first<<" to "<<bb_stub_cst_weight<<std::endl;
			} else {
				TR<<"Skipping resetting of backbone_stub_linear_constraint weight in "<<it->first<<" which is already preset to "<<weight<<std::endl;
			}
		} else {
			utility_exit_with_message( "ERROR: unrecognized stub_energy_fxn_. Only support backbone_stub_constraint or backbone_stub_linear_constraint");
		}

	}
	TR.flush();
}

SetupHotspotConstraintsMover::~SetupHotspotConstraintsMover() {}

} //movers
} //protein_interface_design
} //protocols
