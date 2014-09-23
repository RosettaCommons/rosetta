// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/simple_moves/RingConformationMover.cc
/// @brief   Method definitions for RingConformationMover.
/// @author  Labonte

// Unit headers
#include <protocols/simple_moves/RingConformationMover.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/RingConformerSet.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>

// Utility headers
//#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <string>
#include <iostream>


// Construct tracers.
static thread_local basic::Tracer TR( "protocols.simple_moves.RingConformationMover" );

// Construct random-number generator.


namespace protocols {
namespace simple_moves {

using namespace core;


// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details  By default, all rings within a given pose will be allowed to move.
RingConformationMover::RingConformationMover(): Mover()
{
	using namespace kinematics;

	// Set default MoveMap.
	MoveMapOP default_movemap( new MoveMap() );
	default_movemap->set_nu(true);

	init(default_movemap);
}

// Copy constructor
RingConformationMover::RingConformationMover(RingConformationMover const & object_to_copy): Mover(object_to_copy)
{
	copy_data(*this, object_to_copy);
}

// Constructor with MoveMap input option
/// @param    <input_movemap>: a MoveMap with desired nu torsions set to true
/// @remarks  Movable cyclic residues will generally be a subset of residues in the MoveMap whose nu
/// torsions are set to true.
RingConformationMover::RingConformationMover(core::kinematics::MoveMapOP input_movemap)
{
	init(input_movemap);
}

// Assignment operator
RingConformationMover &
RingConformationMover::operator=(RingConformationMover const & object_to_copy)
{
	// Abort self-assignment.
	if (this == &object_to_copy) {
		return *this;
	}

	Mover::operator=(object_to_copy);
	copy_data(*this, object_to_copy);
	return *this;
}

// Destructor
RingConformationMover::~RingConformationMover() {}


// Standard Rosetta methods ////////////////////////////////////////////////////
// General methods
void
RingConformationMover::register_options()
{
	// No options to register yet
}

void
RingConformationMover::show(std::ostream & output) const
{
	using namespace std;

	Mover::show(output);  // name, type, tag

	output << "Current MoveMap:" << endl << *movemap_ << endl;
}


// Mover methods
std::string
RingConformationMover::get_name() const
{
	return type();
}

protocols::moves::MoverOP
RingConformationMover::clone() const
{
	return protocols::moves::MoverOP( new RingConformationMover(*this) );
}

protocols::moves::MoverOP
RingConformationMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new RingConformationMover() );
}


/// @details  The mover will create a list of movable residues based on the given MoveMap and select a residue from
/// the list at random.  The torsion angles of a randomly selected ring conformer will be applied to the selected
/// residue.
/// @param    <input_pose>: the structure to be moved
/// @remarks  a work in progess...
void
RingConformationMover::apply(Pose & input_pose)
{
	using namespace std;
	using namespace utility;
	using namespace conformation;
	using namespace chemical;

	show(TR);

	TR << "Getting movable residues...." << endl;

	setup_residue_list(input_pose);

	if (residue_list_.empty()) {
		TR.Warning << "There are no movable cyclic residues available in the given pose." << endl;
		return;
	}

	TR << "Applying " << get_name() << " to pose...." << endl;

	core::uint const i = core::uint(numeric::random::rg().uniform() * residue_list_.size() + 1);
	core::uint const res_num = residue_list_[i];
	Residue const & res = input_pose.residue(res_num);

	TR << "Selected residue " << res_num << ": " << res.name() << endl;

	// TODO: Provide a method for specifying subsets of conformers instead of picking one entirely at random.
	RingConformer const & conformer = res.type().ring_conformer_set()->get_random_conformer();

	TR << "Selected the " << conformer.specific_name << " conformation to apply." << endl;

	TR << "Making move...." << endl;

	input_pose.set_ring_conformation(res_num, conformer);

	TR << "Move complete." << endl;
}


// Accessors/Mutators
kinematics::MoveMapCOP
RingConformationMover::movemap() const
{
	return movemap_;
}

void
RingConformationMover::movemap(kinematics::MoveMapOP new_movemap)
{
	movemap_ = new_movemap;
}


// Private methods /////////////////////////////////////////////////////////////
// Initialize data members from arguments.
void
RingConformationMover::init(core::kinematics::MoveMapOP movemap)
{
	type("RingConformationMover");

	movemap_ = movemap;
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformationMover::copy_data(
		RingConformationMover object_to_copy_to,
		RingConformationMover object_to_copy_from)
{
	object_to_copy_to.movemap_ = object_to_copy_from.movemap_;
	object_to_copy_to.residue_list_ = object_to_copy_from.residue_list_;
}

// Setup list of movable cyclic residues from MoveMap.
void
RingConformationMover::setup_residue_list(core::pose::Pose & pose)
{
	using namespace conformation;

	residue_list_.clear();

	for (core::uint res_num = 1, last_res_num = pose.total_residue(); res_num <= last_res_num; ++res_num) {
		Residue const & residue = pose.residue(res_num);
		if (residue.type().is_cyclic()) {
			if (movemap_->get_nu(res_num) == true) {
				residue_list_.push_back(res_num);
			}
		}
	}
}


// Helper methods //////////////////////////////////////////////////////////////
// Insertion operator (overloaded so that RingConformationMover can be "printed" in PyRosetta).
std::ostream &
operator<<(std::ostream & output, RingConformationMover const & object_to_output)
{
	object_to_output.show(output);
	return output;
}

}  // namespace simple_moves
}  // namespace protocols
