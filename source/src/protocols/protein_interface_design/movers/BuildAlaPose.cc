// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/movers/BuildAlaPose.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPoseCreator.hh>

// Package headers
#include <protocols/rosetta_scripts/util.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/types.hh>
#include <core/kinematics/Jump.hh>
#include <protocols/simple_moves/DesignRepackMover.hh>



namespace protocols {
namespace protein_interface_design {
namespace movers {

using namespace core;
using namespace std;
using namespace core::scoring;
using namespace protocols::moves;

static basic::Tracer TR( "protocols.protein_interface_design.movers.BuildAlaPose" );

std::string
BuildAlaPoseCreator::keyname() const
{
	return BuildAlaPoseCreator::mover_name();
}

protocols::moves::MoverOP
BuildAlaPoseCreator::create_mover() const {
	return new BuildAlaPose;
}

std::string
BuildAlaPoseCreator::mover_name()
{
	return "build_Ala_pose";
}

BuildAlaPose::BuildAlaPose() : simple_moves::DesignRepackMover( BuildAlaPoseCreator::mover_name() ) {}
BuildAlaPose::BuildAlaPose(
	bool const partner1,
	bool const partner2,
	core::Real interface_distance_cutoff
) :
	simple_moves::DesignRepackMover( BuildAlaPoseCreator::mover_name() )
{
	repack_partner1_=design_partner1_=partner1;
	repack_partner2_=design_partner2_=partner2;
	interface_distance_cutoff_ = interface_distance_cutoff;
	runtime_assert( interface_distance_cutoff_ >= 0 );
}

BuildAlaPose::~BuildAlaPose() {}

void
BuildAlaPose::apply( pose::Pose & pose )
{
	using namespace core::scoring;
	using namespace core::pack::task::operation;

	allowed_aas_.assign( core::chemical::num_canonical_aas, false );
	allowed_aas_[ core::chemical::aa_from_name(AA_) ] = true;
	//allowed_aas_[ core::chemical::aa_ala ] = true;
	//allowed_aas_[ core::chemical::aa_gly ] = true;
/*	if( repack_partner1_ ^ repack_partner2_ ){
		bool const prevent_chain1( !repack_partner1_ );
		bool const prevent_chain2( !repack_partner2_ );
		core::Size const prevent_chain_begin( pose.conformation().chain_begin( prevent_chain1 ? 1 : 2 ) );
		core::Size const prevent_chain_end( pose.conformation().chain_end( prevent_chain2 ? 2 : 1 ) );
		for( core::Size res=prevent_chain_begin; res<=prevent_chain_end; ++res )
			prevent_repacking_.push_back( res );
	} */
	setup_packer_and_movemap( pose );

	using namespace core::scoring;
	core::scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	pack::pack_rotamers( pose, *scorefxn, task_ );
	(*scorefxn)( pose );
	/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );
}

std::string
BuildAlaPose::get_name() const {
	return BuildAlaPoseCreator::mover_name();
}

void
BuildAlaPose::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &, Movers_map const &, core::pose::Pose const & )
{
	design_partner1_ = tag->getOption<bool>( "partner1", 0 );
	design_partner2_ = tag->getOption<bool>( "partner2", 1 );
	repack_partner1_ = design_partner1_;
	repack_partner2_ = design_partner2_;
	interface_distance_cutoff_ = tag->getOption<core::Real>( "interface_cutoff_distance", 20.0 );
	AA_ = tag->getOption<std::string>( "AA", "ALA" );

	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	TR<<"defined BuildAlaPose mover "<<" for partners "<<( repack_partner1_ ? "1" : "" )<<( repack_partner2_ ? "2": "" )<<" with distance cutoff "<< interface_distance_cutoff_ << " and convert to type " << core::chemical::aa_from_name(AA_) << " NOT WORK for GLY"<< std::endl;
}

protocols::moves::MoverOP
BuildAlaPose::clone() const {
    return( protocols::moves::MoverOP( new BuildAlaPose( *this ) ));
}

} //movers
} //protein_interface_design
} //protocols
