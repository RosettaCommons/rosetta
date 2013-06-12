// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    RingConformationMover.cc
/// @brief   Method definitions for RingConformationMover.
/// @author  labonte

// Unit headers
#include <protocols/simple_moves/carbohydrates/RingConformationMover.hh>
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <string>
#include <iostream>


// Construct tracers.
static basic::Tracer TR("protocols.simple_moves.carbohydrates.RingConformationMover");

using basic::Warning;

// Construct random-number generator.
static numeric::random::RandomGenerator RG(17);  // Does this # mattter?


namespace protocols {
namespace simple_moves {
namespace carbohydrates {

using namespace core;

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
/// @details  By default, all carbohydrate rings within a given pose will be moved.
RingConformationMover::RingConformationMover(): Mover()
{
	using namespace kinematics;

	// Set default MoveMap.
	MoveMapOP default_movemap = new MoveMap();
	default_movemap->set_bb(true);

	init(default_movemap);
}

// Copy constructor
RingConformationMover::RingConformationMover(RingConformationMover const & object_to_copy): Mover(object_to_copy)
{
	copy_data(*this, object_to_copy);
}

// Constructor with MoveMap input option
/// @param    <input_movemap>: a MoveMap with desired backbone torsions set to true
/// @remarks  Movable carbohydrate residues will generally be a subset of residues in the MoveMap whose backbone
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
	return new RingConformationMover(*this);
}

protocols::moves::MoverOP
RingConformationMover::fresh_instance() const
{
	return new RingConformationMover();
}


/// @details  The mover will create a list of movable residues based on the given MoveMap and select a residue from
/// the list at random.  The torsion angles of a randomly selected ring conformer will be applied to the selected
/// residue.
/// @param    <input_pose>: the structure to be moved
/// @remarks  A work in progress...
/// Currently, this algorithm depends on CHI angles being set in the params file to fill in for nu angles, as
/// Rosetta has no mechanism for dealing with ring torsions.  This may cause issues with other protocols using
/// MoveMap, since MoveMap considers part of the ring as BB and part of the ring as CHI, with some overlap.
void
RingConformationMover::apply(Pose & input_pose)
{
	using namespace std;
	using namespace basic::options;
	using namespace utility;
	using namespace conformation;
	using namespace id;
	using namespace pose::carbohydrates;

	if (option[OptionKeys::carbohydrates::lock_rings]) {
		Warning() << "Rings have been locked; no ring conformation moves are being applied to this pose." << endl;
		return;
	}

	show(TR);

	TR << "Getting movable residues...." << endl;

	setup_residue_list(input_pose);

	if (residue_list_.empty()) {
		Warning() << "There are no movable carbohydrate residues available in the given pose." << endl;
		return;
	}

	TR << "Applying " << get_name() << " to pose...." << endl;

	Size i = Size(RG.uniform() * residue_list_.size() + 1);
	Size res_num = residue_list_[i];
	Residue res = input_pose.residue(res_num);

	TR << "Selected residue " << res_num << ": " << res.name() << endl;

	Size ring_size = res.carbohydrate_info()->ring_size();
	ring_conf_def conformer;
	Size j;

	switch (ring_size) {
		case 5:
			//j = Size(RG.uniform() * five_membered_ring_conformers_.size() + 1);
			//conformer = five_membered_ring_conformers_[j];
			Warning() << "Ring flips for 5-membered rings not yet coded!" << endl;
			return;
		case 6:
			j = Size(RG.uniform() * six_membered_ring_conformers_.size() + 1);
			conformer = six_membered_ring_conformers_[j];
			break;
	}

	TR << "Selected the " << conformer.first << " conformation to apply." << endl;

	TR << "Making move...." << endl;

	pair<TorsionType, Size> nu_id;
	for (Size k = 1; k <= ring_size - 2; ++k) {
		nu_id = res.carbohydrate_info()->nu_id(k);
		input_pose.set_torsion(TorsionID(res_num, nu_id.first, nu_id.second), conformer.second[k]);
		align_virtual_atoms_in_carbohydrate_residue(input_pose, res_num);
	}

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
	using namespace std;
	using namespace utility;
	using namespace id;

	type("RingConformationMover");

	movemap_ = movemap;

	// TODO: Create ring_conformer_io.cc file that will include function to read the below data from rosetta_database.

	// TODO: Define furanose ring conformers.

	// TODO: Properly name other conformers.

	// Define pyranose ring conformers.

	// (C-style arrays should normally be avoided, but it makes initializing the vector1 so much easier.)
	// We only need 4 torsion angles to define a 6-membered conformer.
	Angle angles[4];

	// Two chair conformers
	// It's lazy, but for now, I am adding multiple copies of each to compete with entropy of multiple twist-boat states.
	// This will result in a 2:1 chair to twist-boat ratio, which at least matches the CAPRI 27 starting structure.
	for (Size i = 1; i <= 6; ++i) {
		angles[0] = -60.0; angles[1] = 60.0; angles[2] = -60.0; angles[3] = 60.0;
		six_membered_ring_conformers_.push_back(make_pair("1C4", vector1<Angle>(angles, angles + 4)));
		angles[0] = 60.0; angles[1] = -60.0; angles[2] = 60.0; angles[3] = -60.0;
		six_membered_ring_conformers_.push_back(make_pair("4C1", vector1<Angle>(angles, angles + 4)));
	}

	// Six boat conformers
	// We'll assume these are always energy maxima and comment them out for now.  Maybe add option for them later?
	/*angles[0] = 0.0; angles[1] = -60.0; angles[2] = 60.0; angles[3] = 0.0;
	six_membered_ring_conformers_.push_back(make_pair("boat", vector1<Real>(angles, angles + 4)));
	angles[0] = 60.0; angles[1] = 0.0; angles[2] = -60.0; angles[3] = 60.0;
	six_membered_ring_conformers_.push_back(make_pair("B1,4", vector1<Real>(angles, angles + 4)));
	angles[0] = -60.0; angles[1] = 60.0; angles[2] = 0.0; angles[3] = -60.0;
	six_membered_ring_conformers_.push_back(make_pair("boat", vector1<Real>(angles, angles + 4)));

	angles[0] = 0.0; angles[1] = 60.0; angles[2] = -60.0; angles[3] = 0.0;
	six_membered_ring_conformers_.push_back(make_pair("boat", vector1<Real>(angles, angles + 4)));
	angles[0] = -60.0; angles[1] = 0.0; angles[2] = 60.0; angles[3] = -60.0;
	six_membered_ring_conformers_.push_back(make_pair("1,4B", vector1<Real>(angles, angles + 4)));
	angles[0] = 60.0; angles[1] = -60.0; angles[2] = 0.0; angles[3] = 60.0;
	six_membered_ring_conformers_.push_back(make_pair("boat", vector1<Real>(angles, angles + 4)));*/

	// Six twist-boat conformers
	angles[0] = -30.0; angles[1] = -30.0; angles[2] = 60.0; angles[3] = -30.0;
	six_membered_ring_conformers_.push_back(make_pair("4S1", vector1<Real>(angles, angles + 4)));
	angles[0] = -30.0; angles[1] = 60.0; angles[2] = -30.0; angles[3] = -30.0;
	six_membered_ring_conformers_.push_back(make_pair("twist-boat", vector1<Real>(angles, angles + 4)));
	angles[0] = 60.0; angles[1] = -30.0; angles[2] = -30.0; angles[3] = 60.0;
	six_membered_ring_conformers_.push_back(make_pair("twist-boat", vector1<Real>(angles, angles + 4)));

	angles[0] = 30.0; angles[1] = 30.0; angles[2] = -60.0; angles[3] = 30.0;
	six_membered_ring_conformers_.push_back(make_pair("1S4", vector1<Real>(angles, angles + 4)));
	angles[0] = 30.0; angles[1] = -60.0; angles[2] = 30.0; angles[3] = 30.0;
	six_membered_ring_conformers_.push_back(make_pair("twist-boat", vector1<Real>(angles, angles + 4)));
	angles[0] = -60.0; angles[1] = 30.0; angles[2] = 30.0; angles[3] = -60.0;
	six_membered_ring_conformers_.push_back(make_pair("twist-boat", vector1<Real>(angles, angles + 4)));
}

// Copy all data members from <object_to_copy_from> to <object_to_copy_to>.
void
RingConformationMover::copy_data(
		RingConformationMover object_to_copy_to,
		RingConformationMover object_to_copy_from)
{
	object_to_copy_to.movemap_ = object_to_copy_from.movemap_;
	object_to_copy_to.residue_list_ = object_to_copy_from.residue_list_;

	object_to_copy_to.five_membered_ring_conformers_ =
			object_to_copy_from.five_membered_ring_conformers_;
	object_to_copy_to.six_membered_ring_conformers_ =
				object_to_copy_from.six_membered_ring_conformers_;
}

// Setup list of movable carbohydrate residues from MoveMap.
void
RingConformationMover::setup_residue_list(core::pose::Pose & pose)
{
	using namespace conformation;

	residue_list_.clear();

	for (Size res_num = 1, last_res_num = pose.total_residue(); res_num <= last_res_num; ++res_num) {
		Residue const & residue = pose.residue(res_num);
		if (residue.is_carbohydrate()) {
			if (residue.carbohydrate_info()->is_cyclic() && movemap_->get_bb(res_num) == true) {
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

}  // namespace carbohydrates
}  // namespace simple_moves
}  // namespace protocols
