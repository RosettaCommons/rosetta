// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspots.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu)

// Unit headers
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspots.hh>
#include <protocols/protein_interface_design/movers/AddSidechainConstraintsToHotspotsCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <utility/tag/Tag.hh>
#include <protocols/protein_interface_design/movers/PlaceUtils.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <boost/foreach.hpp>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/datacache/DataMap.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/kinematics/FoldTree.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
namespace movers {

using core::pose::Pose;
using namespace protocols::moves;

static thread_local basic::Tracer TR( "protocols.protein_interface_design.movers.AddSidechainConstraintsToHotspots" );
static thread_local basic::Tracer TR_cst( "protocols.protein_interface_design.movers.AddSidechainConstraintsToHotspots_csts" );

std::string
AddSidechainConstraintsToHotspotsCreator::keyname() const
{
	return AddSidechainConstraintsToHotspotsCreator::mover_name();
}

protocols::moves::MoverOP
AddSidechainConstraintsToHotspotsCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddSidechainConstraintsToHotspots );
}

std::string
AddSidechainConstraintsToHotspotsCreator::mover_name()
{
	return "AddSidechainConstraintsToHotspots";
}


AddSidechainConstraintsToHotspots::AddSidechainConstraintsToHotspots() :
	protocols::moves::Mover( AddSidechainConstraintsToHotspotsCreator::mover_name() ),
	chain_( 2 ),
	coord_sdev_( 1.0 )
{ }

AddSidechainConstraintsToHotspots::~AddSidechainConstraintsToHotspots() {}

void
AddSidechainConstraintsToHotspots::apply( Pose & pose )
{
	core::Size const begin( pose.conformation().chain_begin( chain() ) );
	core::Size const end( pose.conformation().chain_end( chain() ) );

	int const num_cutpoints( pose.fold_tree().num_cutpoint() );
	if( num_cutpoints <= 2 && residues().size() == 0 ){
		TR<<"Not enough cutpoints in pose and no residues defined by user. Doing nothing"<<std::endl;
		return;
	}
	for( int i=2; i<=pose.fold_tree().num_cutpoint(); ++i ){
		core::Size const cutpoint = pose.fold_tree().cutpoint( i );
		core::Size const cutpoint_i_1 = pose.fold_tree().cutpoint( i - 1 );
		if( cutpoint - 1 != cutpoint_i_1 ) continue;//only mark residues that are cut on both ends
		if( cutpoint <= end && cutpoint >= begin )
			add_residue( i );
	}
	BOOST_FOREACH( core::Size const residue, residues() )
	{
		using namespace core::scoring::constraints;

		core::scoring::func::HarmonicFuncOP dummy_cst;
		ConstraintCOPs constraint;
		constraint = add_coordinate_constraints( pose, pose.conformation().residue( residue ), chain(), residue, coord_sdev(), dummy_cst );
		BOOST_FOREACH( ConstraintCOP cst, constraint ){
			cst->show_def( TR_cst, pose );
		}
	}
	TR.flush();
}

std::string
AddSidechainConstraintsToHotspots::get_name() const {
	return AddSidechainConstraintsToHotspotsCreator::mover_name();
}

void
AddSidechainConstraintsToHotspots::parse_my_tag( TagCOP const tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		Movers_map const &,
		Pose const & pose )
{
	chain( tag->getOption< core::Size >( "chain", 2 ) );
	coord_sdev( tag->getOption< core::Real >( "coord_sdev", 1.0 ) );
	utility::vector1< core::Size > v1 = core::pose::get_resnum_list( tag, "resnums", pose );
	BOOST_FOREACH( core::Size const r, v1 ){ add_residue( r ); }
}

core::Size
AddSidechainConstraintsToHotspots::chain() const{
	return chain_;
}

void
AddSidechainConstraintsToHotspots::chain( core::Size const c ){
	chain_ = c;
}

core::Real
AddSidechainConstraintsToHotspots::coord_sdev() const{
	return coord_sdev_;
}

void
AddSidechainConstraintsToHotspots::coord_sdev( core::Real const c ){
	coord_sdev_ = c;
}

void
AddSidechainConstraintsToHotspots::add_residue( core::Size const res ){
	residues_.insert( res );
}

std::set< core::Size > const &
AddSidechainConstraintsToHotspots::residues() const{
	return residues_;
}

} //movers
} //protein_interface_design
} //protocols
