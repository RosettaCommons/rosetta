// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ShortBackrubMover.cc
/// @brief implementation of ShortBackrubMover class and functions
/// @author Noah Ollikainen (nollikai@gmail.com)

// Unit headers
#include <protocols/simple_moves/ShortBackrubMover.hh>
#include <protocols/simple_moves/ShortBackrubMoverCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/backrub/BackrubMover.hh>
#include <core/conformation/Residue.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Parser headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/xyz.functions.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.ShortBackrubMover" );
//static numeric::random::RandomGenerator RG(415609);

namespace protocols {
namespace simple_moves {

std::string
ShortBackrubMoverCreator::keyname() const {
	return ShortBackrubMoverCreator::mover_name();
}

protocols::moves::MoverOP
ShortBackrubMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ShortBackrubMover );
}

std::string
ShortBackrubMoverCreator::mover_name() {
	return "ShortBackrubMover";
}

// default constructor
ShortBackrubMover::ShortBackrubMover() : protocols::moves::Mover(),
		backrubmover_(protocols::backrub::BackrubMoverOP( new protocols::backrub::BackrubMover() ))
{
	protocols::moves::Mover::type( "ShortBackrub" );
	backrubmover_->branchopt().read_database();
	backrubmover_->set_min_atoms(3);
	backrubmover_->set_max_atoms(7);
	backrubmover_->set_preserve_detailed_balance(true);
	backrubmover_->set_custom_angle(true);
	resnum_ = 0;
	rotation_std_dev_ = 4.572016;
	randomize_resnum_ = false;
	uniform_backrub_ = false;
}

// constructor that sets the input pose for the BackrubMover
ShortBackrubMover::ShortBackrubMover( core::pose::PoseOP pose ) : protocols::moves::Mover(),
		backrubmover_(protocols::backrub::BackrubMoverOP( new protocols::backrub::BackrubMover() ))
{
	protocols::moves::Mover::type( "ShortBackrub" );
	backrubmover_->branchopt().read_database();
	backrubmover_->set_min_atoms(3);
	backrubmover_->set_max_atoms(7);
	backrubmover_->set_preserve_detailed_balance(true);
	backrubmover_->set_custom_angle(true);
	backrubmover_->set_input_pose(pose);
	resnum_ = 0;
	rotation_std_dev_ = 4.572016;
	randomize_resnum_ = false;
	uniform_backrub_ = false;
}


// copy constructor
ShortBackrubMover::ShortBackrubMover( ShortBackrubMover const & rval ):
	//utility::pointer::ReferenceCount(),
	protocols::moves::Mover( rval ),
	backrubmover_( rval.backrubmover_ ),
	resnum_( rval.resnum_ ),
	rotation_std_dev_( rval.rotation_std_dev_ ),
	randomize_resnum_( rval.randomize_resnum_ ),
	uniform_backrub_( rval.uniform_backrub_ )
{}

// destructor
ShortBackrubMover::~ShortBackrubMover(){}

// clone this object
ShortBackrubMover::MoverOP
ShortBackrubMover::clone() const
{
	return ShortBackrubMover::MoverOP( new protocols::simple_moves::ShortBackrubMover( *this ) );
}

// create this type of object
ShortBackrubMover::MoverOP
ShortBackrubMover::fresh_instance() const
{
	return ShortBackrubMover::MoverOP( new protocols::simple_moves::ShortBackrubMover() );
}

void
ShortBackrubMover::apply( core::pose::Pose & pose )
{
	if (resnum_ == 0)
		randomize_resnum_ = true;

	if (randomize_resnum_) {
		resnum_ = numeric::random::rg().random_range(3, pose.total_residue()-2);
	}

	// only perform move if not adjacent to the end of the chain
	if (core::Size(resnum_-2) > 1 && core::Size(resnum_+2) < pose.total_residue()) {
		utility::vector1<core::Size> pivot_residues;
		utility::vector1<std::string> pivot_atoms;
		backrubmover_->clear_segments();
		pivot_atoms.push_back("CA");
		core::Size start, mid, end;

		// adjust backrub pivots if a proline is at i-1
		if (pose.residue(resnum_-1).name1() == 'P' && pose.residue(resnum_+1).name1() != 'P') {
			pivot_residues.push_back(resnum_-2);
			pivot_residues.push_back(resnum_-1); // proline
			pivot_residues.push_back(resnum_);
			pivot_residues.push_back(resnum_+1);
			start = resnum_ - 2;
			mid = resnum_;
			end = resnum_ + 1;
			backrubmover_->add_mainchain_segments(pivot_residues, pivot_atoms, 3, 34);
		}

		// adjust backrub pivots if a proline is at i+1
		else if (pose.residue(resnum_+1).name1() == 'P' && pose.residue(resnum_-1).name1() != 'P') {
			pivot_residues.push_back(resnum_-1);
			pivot_residues.push_back(resnum_);
			pivot_residues.push_back(resnum_+1); // proline
			pivot_residues.push_back(resnum_+2);
			start = resnum_ - 1;
			mid = resnum_;
			end = resnum_ + 2;
			backrubmover_->add_mainchain_segments(pivot_residues, pivot_atoms, 3, 34);
		}
		else {
			pivot_residues.push_back(resnum_-1);
			pivot_residues.push_back(resnum_);
			pivot_residues.push_back(resnum_+1);
			start = resnum_ - 1;
			mid = resnum_;
			end = resnum_ + 1;
			backrubmover_->add_mainchain_segments(pivot_residues, pivot_atoms, 3, 7);
		}

		core::Vector last_start_o = pose.residue(start).xyz("O");
		core::Vector last_mid_o = pose.residue(mid).xyz("O");

		// if we encountered prolines, there will only be one backrub segment
		// calculate a random rotation angle using a gaussian
		// apply the rotation to the pose
		if (backrubmover_->num_segments() == 1) {
			core::Real rotation_angle = rotation_std_dev_ * numeric::random::rg().gaussian();
			if (uniform_backrub_)
				rotation_angle = 20 - numeric::random::rg().uniform() * 40;
			backrubmover_->set_next_angle(numeric::conversions::radians(rotation_angle));
			backrubmover_->set_next_segment_id(1);
			backrubmover_->apply(pose);
		}

		// if we did not encounter prolines, there will be three backrub segments
		// calculate a random rotation angle using a gaussian
		// apply the rotation to the pose
		// adjust the peptide bonds to minimize the movement of carbonyl oxygens
		else if (backrubmover_->num_segments() == 3) {
			core::Real rotation_angle = rotation_std_dev_ * numeric::random::rg().gaussian();
			if (uniform_backrub_)
				rotation_angle = 20 - numeric::random::rg().uniform() * 40;
			backrubmover_->set_next_angle(numeric::conversions::radians(rotation_angle));
			backrubmover_->set_next_segment_id(2);
			backrubmover_->apply(pose);
			core::Real seg1_dihedral = numeric::dihedral_radians(pose.residue(start).xyz("O"), pose.residue(start).xyz("CA"), pose.residue(mid).xyz("CA"), last_start_o);
			core::Real seg3_dihedral = numeric::dihedral_radians(pose.residue(mid).xyz("O"), pose.residue(mid).xyz("CA"), pose.residue(end).xyz("CA"), last_mid_o);
			backrubmover_->set_next_angle(seg1_dihedral);
			backrubmover_->set_next_segment_id(1);
			backrubmover_->apply(pose);
			backrubmover_->set_next_angle(seg3_dihedral);
			backrubmover_->set_next_segment_id(3);
			backrubmover_->apply(pose);
		}
	}
}

std::string
ShortBackrubMover::get_name() const {
	return "ShortBackrubMover";
}

// setters
void ShortBackrubMover::set_resnum( core::Size resnum ) { resnum_ = resnum; }
void ShortBackrubMover::set_rotation_std_dev( core::Real rotation_std_dev ) { rotation_std_dev_ = rotation_std_dev; }
void ShortBackrubMover::set_randomize_resnum( bool randomize_resnum ) { randomize_resnum_ = randomize_resnum; }
void ShortBackrubMover::set_uniform_backrub( bool uniform_backrub ) { uniform_backrub_ = uniform_backrub; }
void ShortBackrubMover::set_input_pose( core::pose::PoseCOP pose ) {
	backrubmover_->set_input_pose(pose);
}

// getters
core::Size
ShortBackrubMover::get_resnum() const {
	return resnum_;
}
core::Real
ShortBackrubMover::get_rotation_std_dev() const {
	return rotation_std_dev_;
}
bool
ShortBackrubMover::get_randomize_resnum() const {
	return randomize_resnum_;
}
bool
ShortBackrubMover::get_uniform_backrub() const {
	return uniform_backrub_;
}
protocols::backrub::BackrubMoverOP ShortBackrubMover::get_backrubmover() const{
	return backrubmover_;
}

/// @brief parse xml
void
ShortBackrubMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose)
{
	backrubmover_->set_input_pose(core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(pose) ) ));
	backrubmover_->clear_segments();
	backrubmover_->add_mainchain_segments();
	backrubmover_->branchopt().read_database();
}

} // moves
} // protocols
