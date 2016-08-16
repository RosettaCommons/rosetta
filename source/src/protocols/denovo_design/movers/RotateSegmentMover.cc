// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/RotateSegmentMover.cc
/// @brief Rotates a segment in the pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/RotateSegmentMover.hh>
#include <protocols/denovo_design/movers/RotateSegmentMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh>
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.RotateSegmentMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

RotateSegmentMover::RotateSegmentMover():
	protocols::moves::Mover( "RotateSegmentMover" ),
	range_( -numeric::constants::f::pi, numeric::constants::f::pi ),
	residues_(),
	atoms_(),
	frozen_segments_()
{
}

RotateSegmentMover::~RotateSegmentMover()
{}

void
RotateSegmentMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
}

protocols::moves::MoverOP
RotateSegmentMover::clone() const
{
	return protocols::moves::MoverOP( new RotateSegmentMover( *this ) );
}


moves::MoverOP
RotateSegmentMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RotateSegmentMover );
}

std::string
RotateSegmentMover::get_name() const
{
	return RotateSegmentMover::class_name();
}

std::string
RotateSegmentMover::class_name()
{
	return "RotateSegment";
}

void
RotateSegmentMover::set_rotation_range( core::Real const start, core::Real const stop )
{
	range_ = RotationRange( start, stop );
}

void
RotateSegmentMover::apply( core::pose::Pose & pose )
{
	apply( pose, numeric::random::rg().uniform() );
}

void
RotateSegmentMover::apply( core::pose::Pose & pose, core::Real const random ) const
{
	using components::StructureData;

	// choose a rotation
	core::Real const build_rotation = select_rotation( random );

	StructureData const & sd =
		components::StructureDataFactory::get_instance()->get_from_pose( pose );

	// collect atom ids
	AtomIDs const atom_ids = target_atoms( pose, sd );
	debug_assert( atom_ids.size() == 4 );

	// collect xyz vectors
	utility::vector1< core::Vector > atom_xyzs;
	for ( AtomIDs::const_iterator a=atom_ids.begin(); a!=atom_ids.end(); ++a ) {
		atom_xyzs.push_back( pose.xyz( *a ) );
	}
	debug_assert( atom_xyzs.size() == 4 );

	// start and end of dihedral must be on different chains
	if ( pose.chain( atom_ids.begin()->rsd() ) == pose.chain( atom_ids.rbegin()->rsd() ) ) {
		std::stringstream msg;
		msg << "RotateSegmentMover::apply(): Start (atom=" << *atom_ids.begin() << ", chain="
			<< pose.chain( atom_ids.begin()->rsd() ) << ") and end (atom="
			<< *atom_ids.rbegin() << ", chain=" << pose.chain( atom_ids.rbegin()->rsd() )
			<< ") atoms of the dihedral to be rotated must be on the same chain." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// start and end must be in different movable groups
	if ( sd.movable_group( atom_ids[2].rsd() ) == sd.movable_group( atom_ids[3].rsd() ) ) {
		std::stringstream msg;
		msg << "RotateSegmentMover::apply(): Central atoms (atom=" << atom_ids[2] << ", group="
			<< sd.movable_group( atom_ids[2].rsd() ) << ") and (atom="
			<< atom_ids[3] << ", group=" << sd.movable_group( atom_ids[3].rsd() )
			<< ") of the dihedral to be rotated must be in different groups." << std::endl;
		msg << "StructureData = " << sd << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Vector axis = atom_xyzs[2]-atom_xyzs[3];
	axis.normalize();
	core::Real const current_rotation = calc_dihedral( atom_xyzs[1], atom_xyzs[2], atom_xyzs[3], atom_xyzs[4] );
	core::Real const angle_to_rotate = build_rotation - current_rotation;
	TR << "Rotating current dihedral of " << current_rotation << " by " << angle_to_rotate << " to achieve dihedral of " << build_rotation << std::endl;
	core::Vector const center = (atom_xyzs[2] + atom_xyzs[3]) / 2;

	// prepare FoldTree for rotation
	SegmentNames const fg_roots = boost::assign::list_of
		(sd.segment_name( atom_ids[4].rsd() ))(sd.segment_name( atom_ids[1].rsd() ));
	FoldTreeFromFoldGraphMover make_ft( fg_roots, protocols::loops::Loops() );
	make_ft.apply( pose );

	std::string const parent_seg = sd.segment_name( atom_ids[1].rsd() );
	int const jump_idx = find_jump_rec( pose.fold_tree(), sd.segment( parent_seg ).safe() );
	debug_assert( jump_idx );
	TR << "Grabbed jump " << jump_idx << " of " << pose.fold_tree() << std::endl;

	// get frozen jump indices and make they come off safe residue of segments_[1]
	for ( SegmentNames::const_iterator f=frozen_segments_.begin(); f!=frozen_segments_.end(); ++f ) {
		int const childjump = find_jump_rec( pose.fold_tree(), sd.segment( *f ).safe() );
		debug_assert( childjump );
		if ( jump_idx == childjump ) {
			TR.Warning << "Parent jump and child jump are the same while sliding jump for " << *f
				<< " to have parent " << parent_seg << " --  not doing anything." << std::endl;
			return;
		}
		core::kinematics::Edge const jedge = pose.fold_tree().jump_edge( childjump );
		core::kinematics::FoldTree const new_ft = slide_jump(
			pose.fold_tree(),
			childjump,
			sd.segment(parent_seg).safe(),
			jedge.stop() );

		pose.fold_tree( new_ft );
		//perm.slide_jump( *f, parent_seg );
	}

	core::kinematics::Jump j( pose.jump(jump_idx) );
	j.rotation_by_axis(
		pose.conformation().upstream_jump_stub(jump_idx),
		axis,
		center,
		angle_to_rotate );

	pose.set_jump( jump_idx, j );

	TR << "Pose rotation set in RotatedTomponent to " << build_rotation << std::endl;

	utility::vector1< core::Vector > new_atom_xyzs;
	for ( AtomIDs::const_iterator a=atom_ids.begin(); a!=atom_ids.end(); ++a ) {
		new_atom_xyzs.push_back( pose.xyz( *a ) );
	}
	TR << "Actual rotation is " << calc_dihedral( new_atom_xyzs[1], new_atom_xyzs[2], new_atom_xyzs[3], new_atom_xyzs[4] ) << std::endl;
}

RotateSegmentMover::AtomIDs
RotateSegmentMover::target_atoms(
	core::pose::Pose const & pose,
	components::StructureData const & sd ) const
{
	if ( residues_.size() != 4 ) {
		std::stringstream msg;
		msg << "RotateSegmentMover::target_atoms(): Exactly four residues must be specified to define a proper dihedral."
			<< " The current residue list is: " << residues_ << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( atoms_.size() != 4 ) {
		std::stringstream msg;
		msg << "RotateSegmentMover::target_atoms(): Exactly four atom names must be specified to define a proper dihedral."
			<< " The current residue atom name list is: " << atoms_ << std::endl;
		utility_exit_with_message( msg.str() );
	}

	AtomIDs atom_ids;
	Residues::const_iterator r = residues_.begin();
	Atoms::const_iterator a = atoms_.begin();
	for ( ; ( r != residues_.end() ) && ( a != atoms_.end() ); ++r, ++a ) {
		core::Size const resid = get_resid( sd, *r );
		core::Size const atomid = pose.residue( resid ).type().atom_index( *a );
		atom_ids.push_back( core::id::AtomID( atomid, resid ) );
	}
	debug_assert( atom_ids.size() == 4 );
	return atom_ids;
}

core::Real
RotateSegmentMover::select_rotation( core::Real const random ) const
{
	// RAND_BETWEEN_0_and_1 = ( stop_ - x_ ) / ( stop_ - start_ )
	// therefore,
	// x_ = stop_ - RAND * ( stop_ - start_ )
	return range_.stop - random * ( range_.stop - range_.start );
}

// calculates a dihedral from four coordinates in degrees
core::Real
calc_dihedral(
	core::Vector const & p1,
	core::Vector const & p2,
	core::Vector const & p3,
	core::Vector const & p4 )
{
	core::Vector const b1( (p4-p3).normalize() );
	core::Vector const b2( (p3-p2).normalize() );
	core::Vector const b3( (p2-p1).normalize() );
	core::Vector const n1( b1.cross(b2) );
	core::Vector const n2( b2.cross(b3) );
	core::Vector const m1( n1.cross(b2) );
	core::Real const x( n1.dot(n2) );
	core::Real const y( m1.dot(n2) );
	return atan2(y,x)/(2*numeric::constants::f::pi)*360.0;
}


/////////////// Creator ///////////////

protocols::moves::MoverOP
RotateSegmentMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new RotateSegmentMover );
}

std::string
RotateSegmentMoverCreator::keyname() const
{
	return RotateSegmentMover::class_name();
}

////// RotationRange ////////////
RotationRange::RotationRange( core::Real start_val, core::Real const stop_val ):
	start( start_val ),
	stop( stop_val )
{}

} //protocols
} //denovo_design
} //movers
