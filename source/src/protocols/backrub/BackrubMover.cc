// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/backrub/BackrubMover.cc
/// @brief implementation of BackrubMover class and functions
/// @author Colin A. Smith (colin.smith@mpibpc.mpg.de)


#include <protocols/backrub/BackrubMover.hh>
#include <protocols/backrub/BackrubMoverCreator.hh>
#include <protocols/backrub/util.hh>

// Protocols Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/select/util.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/DOF_ID_Range.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/types.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/movemap/util.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/IntervalSet.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <iomanip>
#include <set>

// option key includes
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/backrub.OptionKeys.gen.hh>

#include <numeric/trig.functions.hh>
#include <core/id/TorsionID_Range.hh>
#include <core/kinematics/tree/Atom.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/keys/Key3Vector.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <core/select/residue_selector/ResidueSelector.hh> // AUTO IWYU For ResidueSelector


using namespace core;
using namespace core::pose;

static basic::Tracer TR( "protocols.backrub.BackrubMover" );

namespace protocols {
namespace backrub {

BackrubMover::BackrubMover() :
	pivot_atoms_(utility::vector1<std::string>(1, "CA")),
	min_atoms_(3),
	max_atoms_(34),
	max_angle_disp_4_(numeric::conversions::radians(40.)),
	max_angle_disp_7_(numeric::conversions::radians(20.)),
	max_angle_disp_slope_(numeric::conversions::radians(-1./3.)),
	next_segment_id_(0),
	movemap_(/* NULL */),
	preserve_detailed_balance_(false),
	require_mm_bend_(true),
	custom_angle_(false)
{
	init_with_options();
}

BackrubMover::BackrubMover( BackrubMover const & ) = default;

protocols::moves::MoverOP BackrubMover::clone() const { return utility::pointer::make_shared< BackrubMover >( *this ); }
protocols::moves::MoverOP BackrubMover::fresh_instance() const { return utility::pointer::make_shared< BackrubMover >(); }

void
init_backrub_mover_with_options(
	BackrubMover & mover
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::vector1<core::Size> pivot_residues;
	for ( core::Size i = 1; i <= option[ OptionKeys::backrub::pivot_residues ].size(); ++i ) {
		if ( option[ OptionKeys::backrub::pivot_residues ][i] >= 1 ) pivot_residues.push_back(option[ OptionKeys::backrub::pivot_residues ][i]);
	}

	if ( option[in::file::movemap].user() ) {
		core::kinematics::MoveMapOP movemap_ = utility::pointer::make_shared< core::kinematics::MoveMap >();
		movemap_->init_from_file( option[in::file::movemap]() );
		mover.set_movemap(movemap_);
	} else {
		mover.set_pivot_residues(pivot_residues);
	}

	mover.set_pivot_atoms(option[ OptionKeys::backrub::pivot_atoms ]);
	mover.set_min_atoms(option[ OptionKeys::backrub::min_atoms ]);
	mover.set_max_atoms(option[ OptionKeys::backrub::max_atoms ]);
}

void
BackrubMover::init_with_options()
{
	init_backrub_mover_with_options(*this);
}

void
BackrubMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( tag->hasOption("pivot_residues") ) {
		pivot_residue_selector_ = core::pose::get_resnum_selector(tag, "pivot_residues");
		// Clear the command-line input movemap
		if ( option[in::file::movemap].user() ) {
			core::kinematics::MoveMapOP movemap_ = utility::pointer::make_shared< core::kinematics::MoveMap >();
			this->set_movemap(movemap_);
		}
	}

	if ( tag->hasOption("pivot_atoms") ) {
		std::string const pivot_atoms_string( tag->getOption<std::string>("pivot_atoms") );
		this->set_pivot_atoms( utility::string_split( pivot_atoms_string, ',' ) );
	}

	this->set_min_atoms( tag->getOption<core::Size>( "min_atoms", min_atoms_ ) );
	this->set_max_atoms( tag->getOption<core::Size>( "max_atoms", max_atoms_ ) );
	max_angle_disp_4_ = tag->getOption<core::Real>( "max_angle_disp_4", max_angle_disp_4_ );
	max_angle_disp_7_ = tag->getOption<core::Real>( "max_angle_disp_7", max_angle_disp_7_ );
	max_angle_disp_slope_ = tag->getOption<core::Real>( "max_angle_disp_slope", max_angle_disp_slope_ );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", preserve_detailed_balance() ) );
	set_require_mm_bend( tag->getOption<bool>( "require_mm_bend", require_mm_bend() ) );

	if ( tag->hasOption("movemap_factory") ) {
		TR.Debug << "Identified the parser for Movemap Factory" << std::endl;
		parse_movemap_factory(tag, data);
	}

	clear_segments();
	// set_input_pose(PoseCOP( PoseOP( new core::pose::Pose(pose) ) )); // Will be called on first apply() call.
	// add_mainchain_segments(); // Will be called on apply() call

	if ( ! branchopt_.initialized() ) branchopt_.read_database();
}

void
BackrubMover::parse_movemap_factory(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	using namespace core::select::movemap;

	TR << "Defining movemap as per the Movemap Factory Attribute" << std::endl;
	// Including the following option in case a movemap parser is added in the future
	if ( tag->hasOption("chi") || tag->hasOption("bb") ) {
		throw CREATE_EXCEPTION( utility::excn::RosettaScriptsOptionError, "BackrubMover can accept either a MoveMapFactory or bb/chi (as input to a Movemap).");
	} else if ( tag->hasOption("pivot_residues") ) {
		TR.Warning << "*** Pivot residues parser exists. Movemap Factory attribute will be overriden by user-defined pivot residues ***" << std::endl;
	}

	MoveMapFactoryOP mmf = core::select::movemap::parse_movemap_factory(tag, data);
	set_movemap_factory(mmf);
}

void
BackrubMover::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size
)
{
	if ( ! branchopt_.initialized() ) branchopt_.read_database();

	if ( !(num_segments() && get_input_pose() && get_input_pose()->fold_tree() == pose.fold_tree()) ) {

		if ( !(get_input_pose() && get_input_pose()->fold_tree() == pose.fold_tree()) ) {
			set_input_pose(utility::pointer::make_shared< core::pose::Pose >(pose));
		}

		// this code shouldn't get called when the parser is in use
		clear_segments();
		add_mainchain_segments( pose );
	}

	protocols::moves::MonteCarloCOP monte_carlo(metropolis_hastings_mover.monte_carlo());

	if ( monte_carlo->score_function().get_weight(core::scoring::mm_bend) == 0.0 ) {
		TR.Warning << "*** Using BackrubMover without mm_bend Score Term ***" << std::endl;
		if ( require_mm_bend_ ) {
			TR.Warning << "Exit can be prevented by setting require_mm_bend to false" << std::endl;
			utility_exit_with_message("BackrubMover is being used without the mm_bend term.");
		}
	}

	// tell the branch angle optimizer about the score function MMBondAngleResidueTypeParamSet, if any
	if ( monte_carlo->score_function().energy_method_options().bond_angle_residue_type_param_set() ) {
		branchopt_.bond_angle_residue_type_param_set(monte_carlo->score_function().energy_method_options().bond_angle_residue_type_param_set());
	}

	optimize_branch_angles(pose);
}

/// @details
void
BackrubMover::apply(
	Pose & pose
)
{
	// If movemap_factory attribute present, the update the movemap and corresponding pivot_residues
	if ( movemap_factory_ && !(pivot_residue_selector_) ) {
		TR.Debug << "*** Movemap Factory is defined, so setting up movemap ***" << std::endl;
		core::kinematics::MoveMapOP updated_movemap( movemap_factory_->create_movemap_from_pose(pose)->clone() );
		this->set_movemap(updated_movemap);
	}
	// the BranchAngleOptimizer must be initialized by reading a database. Do so now if this has not occurred
	if ( ! branchopt_.initialized() ) branchopt_.read_database();

	if ( !(num_segments() && get_input_pose() && get_input_pose()->fold_tree() == pose.fold_tree()) ) {

		if ( !(get_input_pose() && get_input_pose()->fold_tree() == pose.fold_tree()) ) {
			set_input_pose(utility::pointer::make_shared< core::pose::Pose >(pose));
		}

		clear_segments();
		add_mainchain_segments( pose );
	}

	runtime_assert(segments_.size());

	// pick a random segment
	core::Size segment_id;
	if ( next_segment_id_ ) {
		segment_id = next_segment_id_;
		next_segment_id_ = 0;
	} else {
		segment_id = numeric::random::rg().random_range(1, segments_.size());
	}

	// record data about the move
	last_segment_id_ = segment_id;
	id::AtomID atomid1(segments_[segment_id].start_atomid());
	id::AtomID atomid2(segments_[segment_id].end_atomid());
	if ( atomid2 < atomid1 ) {
		id::AtomID temp(atomid1);
		atomid1 = atomid2;
		atomid2 = temp;
	}
	last_start_atom_name_ = pose.residue_type(atomid1.rsd()).atom_name(atomid1.atomno());
	last_end_atom_name_ = pose.residue_type(atomid2.rsd()).atom_name(atomid2.atomno());

	// the function that calculates the constants would need to be called twice if we didn't save them
	utility::vector0<Real> start_constants;
	utility::vector0<Real> end_constants;

	// pick a random angle
	if ( custom_angle_ ) {
		last_angle_ = next_angle_;
	} else {
		last_angle_ = random_angle(pose, segment_id, start_constants, end_constants);
	}

	TR.Debug << "Applying " << numeric::conversions::degrees(last_angle_) << " degree rotation to segment " << segment_id
		<< std::endl;

	// rotate the segment of the angle is not 0
	if ( last_angle_ != 0 ) {
		rotate_segment(pose, segment_id, last_angle_, start_constants, end_constants);
	}

	update_type();
}


// /// @details
// void
// BackrubMover::test_move(
//  Pose &
// )
// {
//
// }

/// @brief calculate the number of atom tree bonds between the two atoms
/// possibly move this into the Atom class
int tree_distance(
	kinematics::tree::AtomCOP ancestor,
	kinematics::tree::AtomCOP descendent
)
{
	kinematics::tree::Atom const * current = descendent.get();

	for ( int tdist = 0; current; ++tdist ) {
		if ( ancestor.get() == current ) return tdist;
		current = current->parent().get();
	}

	return -1;
}

/// @details
/// The start atom MUST have all three stub atoms defined. However, the end atom
/// may be at the very end of the chain. There must be at least one atom between
/// start and end atoms.
///
/// If the segment is for whatever reason invalid, 0 is returned. Otherwise, the
/// segment is added to the end of the segments vector and the ID to the new
/// segment is returned. If the segment already exists, then nothing is changed
/// and the existing segment ID is returned.
///
/// There has been absolutely no testing of how the code handles jumps! Use this
/// with jumps at your own risk!
Size
BackrubMover::add_segment(
	id::AtomID start_atomid,
	id::AtomID end_atomid,
	Real max_angle_disp // = 0
)
{
	PoseCOP input_pose(get_input_pose());

	runtime_assert( input_pose != nullptr);

	// get references to Atom tree atoms
	kinematics::tree::AtomCOP start_atom(input_pose->atom_tree().atom(start_atomid).get_self_ptr());
	kinematics::tree::AtomCOP end_atom(input_pose->atom_tree().atom(end_atomid).get_self_ptr());

	// calculate initial tree distance
	int tdist(tree_distance(start_atom, end_atom));

	// if no ancestry is detected, check to see if the atoms were reversed
	if ( tdist < 0 ) {
		tdist = tree_distance(end_atom, start_atom);
		if ( tdist >= 2 ) {
			//TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid
			//           << " ordered wrong, reversing" << std::endl;
			kinematics::tree::AtomCOP const temp_atom(start_atom);
			id::AtomID const temp_atomid(start_atomid);
			start_atom = end_atom;
			start_atomid = end_atomid;
			end_atom = temp_atom;
			end_atomid = temp_atomid;
		}
	}

	// if no ancestry is still detected, the segment is invalid
	if ( tdist < 0 ) {
		TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid << " invalid, ignoring"
			<< std::endl;
		return 0;
	}

	// if the tree distance is too short, the segment is invalid
	if ( tdist < 2 ) {
		TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid << " too short, ignoring"
			<< std::endl;
		return 0;
	}

	// check for duplicates
	if ( core::Size segid = segment_id(start_atomid, end_atomid) ) {
		TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid << " already exists, ignoring"
			<< std::endl;
		return segid;
	}

	//TR << "start_atom stub: " << start_atom->stub_atom1_id() << " " << start_atom->stub_atom2_id() << " "
	//   << start_atom->stub_atom3_id() << std::endl;

	if ( !(start_atom->stub_atom1() && start_atom->stub_atom2() && start_atom->stub_atom3()) ) {
		TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid
			<< " has incomplete start stub, ignoring" << std::endl;
		return 0;
	}

	// Find the first two atoms on the path from start to end
	id::AtomID start_atomid1(end_atom->parent()->id());
	id::AtomID start_atomid2(end_atom->id());
	for ( kinematics::tree::Atom const * current_atom = end_atom->parent()->parent().get();
			current_atom != start_atom.get();
			current_atom = current_atom->parent().get() ) {
		start_atomid2 = start_atomid1;
		start_atomid1 = current_atom->id();
	}

	// This catches a rare case that we can't currently handle
	if ( start_atom->stub_atom2()->id() == start_atomid1 ) {
		//TR.Warning << "Backrub segment " << start_atomid << " to " << end_atomid
		//      << " has incompatible stub, ignoring" << std::endl;
		return 0;
	}

	TR.Debug << "Adding Backrub segment " << start_atomid << " to " << end_atomid << " (size " << tdist+1 << ")" << std::endl;

	// use default maximum angular displacement parameters if it was not provided
	if ( max_angle_disp == 0 ) {
		if ( tdist+1 < 7 ) {
			max_angle_disp = max_angle_disp_4_ + (tdist+1 - 4) * max_angle_disp_slope_;
		} else {
			max_angle_disp = max_angle_disp_7_ + (tdist+1 - 7) * max_angle_disp_slope_;
		}
	}

	if ( max_angle_disp <= 0 ) {
		TR.Warning << "Maximum angular displacement is negative! (" << max_angle_disp <<") Using zero instead." << std::endl;
		TR.Warning << "    tdist " << tdist << " max_angle_disp_4 " << max_angle_disp_4_ << " max_angle_disp_7 " << max_angle_disp_7_ << " max_angle_disp_slope " << max_angle_disp_slope_ << std::endl;
		max_angle_disp = 0;
	}

	segments_.push_back(BackrubSegment(start_atomid, start_atomid1, start_atomid2, end_atomid, tdist+1, max_angle_disp));

	core::Size segment_id = segments_.size();

	bond_angle_map_[segments_[segment_id].start_bond_angle_key(*input_pose)] = bond_angle_map_.size() + 1;
	bond_angle_map_[segments_[segment_id].end_bond_angle_key(*input_pose)] = bond_angle_map_.size() + 1;

	return segment_id;
}

core::Size
connected_mainchain_atomids(
	Pose const & pose,
	core::id::AtomID atomid,
	utility::vector1<core::id::AtomID> & atomids
)
{
	// empty the input vector
	atomids.clear();

	// make sure that atomid is a mainchain atom
	core::Size mainchain_index(0);
	{
		chemical::AtomIndices const & mainchain_atoms(pose.residue_type(atomid.rsd()).mainchain_atoms());
		for ( core::Size i = 1; i <= mainchain_atoms.size(); ++i ) {
			if ( mainchain_atoms[i] == atomid.atomno() ) {
				mainchain_index = i;
				break;
			}
		}
	}
	// we weren't given a mainchain atom, return nothing
	if ( ! mainchain_index ) return 0;

	// rewind to the beginning of the chain
	core::Size resnum(atomid.rsd());
	while ( true ) {
		conformation::Residue const & residue(pose.residue(resnum));
		chemical::AtomIndices const & mainchain_atoms(residue.mainchain_atoms());

		if ( !mainchain_atoms.size() ) break;

		// loop over all residue connections in the current residue
		for ( core::Size i = 1; i <= residue.n_possible_residue_connections(); ++i ) {

			// check to see if the connection is to the first mainchain atom
			if ( residue.residue_connect_atom_index(i) == mainchain_atoms.front() && !residue.connection_incomplete(i) ) {

				// check to see if the atom it is connected to is the last mainchain atom of the corresponding residue
				chemical::ResConnID resconid(residue.actual_residue_connection(i));
				core::Size connected_atomno(pose.residue(resconid.resid()).residue_connect_atom_index(resconid.connid()));
				chemical::AtomIndices const & mainchain_atoms2(pose.residue(resconid.resid()).mainchain_atoms());
				if ( ! mainchain_atoms2.empty() && mainchain_atoms2.back() == connected_atomno ) {
					// success, set the new residue number
					resnum = resconid.resid();
					break;
				}
			}
		}

		// get out of the loop if we didn't decrement resnum
		if ( resnum == residue.seqpos() ) break;
	}

	// now iterate over all connected residues and insert their mainchain atomids
	while ( resnum ) {
		conformation::Residue const & residue(pose.residue(resnum));
		chemical::AtomIndices const & mainchain_atoms(residue.mainchain_atoms());

		// set the loop to exit by default
		resnum = 0;

		for ( core::Size i = 1; i <= mainchain_atoms.size(); i++ ) {
			atomids.push_back(id::AtomID(mainchain_atoms[i], residue.seqpos()));
			if ( atomids.back() == atomid ) {
				mainchain_index = atomids.size();
			}
		}

		// loop over all residue connections in the current residue
		for ( core::Size i = 1; i <= residue.n_possible_residue_connections(); ++i ) {

			// check to see if the connection is to the last mainchain atom
			if ( residue.residue_connect_atom_index(i) == mainchain_atoms.back() && !residue.connection_incomplete(i) ) {

				// check to see if the atom it is connected to is the first mainchain atom of the corresponding residue
				chemical::ResConnID resconid(residue.actual_residue_connection(i));
				core::Size connected_atomno(pose.residue(resconid.resid()).residue_connect_atom_index(resconid.connid()));
				chemical::AtomIndices const & mainchain_atoms2(pose.residue(resconid.resid()).mainchain_atoms());
				if ( ! mainchain_atoms2.empty() && mainchain_atoms2.front() == connected_atomno ) {
					// success, set the new residue number
					resnum = resconid.resid();
					break;
				}
			}
		}
	}

	return mainchain_index;
}

core::Size
BackrubMover::add_mainchain_segments(
	utility::vector1<core::id::AtomID> atomids,
	core::Size min_atoms,
	core::Size max_atoms
)
{
	PoseCOP input_pose(get_input_pose());

	runtime_assert(min_atoms >= 3 && max_atoms >= min_atoms);

	std::set<id::AtomID> atomid_set(atomids.begin(), atomids.end());
	core::Size nsegments(0);

	// AMW: cppcheck flags this as possibly inefficient
	while ( !atomid_set.empty() ) {// size()) {

		// get a list of contiguous mainchain atoms connected to the first atomid in the set
		utility::vector1<core::id::AtomID> mainchain_atomids;
		if ( !connected_mainchain_atomids(*input_pose, *atomid_set.begin(), mainchain_atomids) ) {
			// if it wasn't a mainchain atom, delete it from the set
			atomid_set.erase(*atomid_set.begin());
		}

		// iterate over all atoms in the chain
		for ( core::Size i = 1; i <= mainchain_atomids.size(); ++i ) {

			// check to see if the atom is in our set
			if ( atomid_set.count(mainchain_atomids[i]) ) {
				// delete it from the set
				atomid_set.erase(mainchain_atomids[i]);

				// exclude proline N/CA atoms
				if ( input_pose->residue(mainchain_atomids[i].rsd()).aa() == chemical::aa_pro &&
						(input_pose->residue(mainchain_atomids[i].rsd()).atom_name(mainchain_atomids[i].atomno()) == " N  " ||
						input_pose->residue(mainchain_atomids[i].rsd()).atom_name(mainchain_atomids[i].atomno()) == " CA ") ) {
					continue;
				}

				// add any segments for atomids in the set
				for ( core::Size j = i+min_atoms-1; j <= i+max_atoms-1 && j <= mainchain_atomids.size(); ++j ) {
					if ( atomid_set.count(mainchain_atomids[j]) ) {
						// exclude proline N/CA atoms
						if ( input_pose->residue(mainchain_atomids[j].rsd()).aa() == chemical::aa_pro &&
								(input_pose->residue(mainchain_atomids[j].rsd()).atom_name(mainchain_atomids[j].atomno()) == " N  " ||
								input_pose->residue(mainchain_atomids[j].rsd()).atom_name(mainchain_atomids[j].atomno()) == " CA ") ) {
							continue;
						}
						add_segment(mainchain_atomids[i], mainchain_atomids[j]);
					}
				}
			}
		}
	}

	return nsegments;
}

core::Size
BackrubMover::add_mainchain_segments(
	utility::vector1<core::Size> resnums,
	utility::vector1<std::string> atomnames,
	core::Size min_atoms,
	core::Size max_atoms
)
{
	PoseCOP input_pose(get_input_pose());

	runtime_assert(min_atoms >= 3 && max_atoms >= min_atoms);
	core::Size nsegments(0);

	std::sort(resnums.begin(), resnums.end());
	auto new_end(std::unique(resnums.begin(), resnums.end()));
	resnums.erase(new_end, resnums.end());

	TR << "Segment lengths: " << min_atoms << "-" << max_atoms << " atoms" << std::endl;
	TR << "Main chain pivot atoms:";
	for ( core::Size i = 1; i <= pivot_atoms_.size(); ++i ) TR << ' ' << pivot_atoms_[i];
	TR << std::endl;

	// loop over all contiguous segments
	core::Size contiguous_begin = 1;
	while ( contiguous_begin <= resnums.size() ) {

		// find the end of the current segment
		core::Size contiguous_end = contiguous_begin;
		while ( contiguous_end+1 <= resnums.size() && resnums[contiguous_end]+1 == resnums[contiguous_end+1] ) {
			++contiguous_end;
		}

		TR << "Adding backrub segments for residues " << resnums[contiguous_begin] << "-" << resnums[contiguous_end]
			<< std::endl;

		// create a vector of input atomids
		utility::vector1<core::id::AtomID> mainchain_atomids;
		for ( core::Size i = contiguous_begin; i <= contiguous_end; ++i ) {
			core::Size resnum(resnums[i]);
			if ( resnum >= 1 && resnum <= input_pose->size() ) {
				for ( core::Size j = 1; j <= atomnames.size(); j++ ) {
					// KAB - Added if statement to only add atom to vector if it actually exists in the residue
					// This allows for HETATMS without C-alphas (or other pivot atoms)
					//   to be excluded automatically
					if ( input_pose->residue(resnum).type().has( atomnames[j] ) ) {
						mainchain_atomids.push_back(core::id::AtomID(input_pose->residue(resnum).atom_index(atomnames[j]), resnum ));
					}
				}
			}
		}

		// enumerate segments for the current group of contiguous residues
		nsegments += add_mainchain_segments(mainchain_atomids, min_atoms, max_atoms);

		// advance to after the current segment
		contiguous_begin = contiguous_end+1;
	}

	TR << "Total Segments Added: " << num_segments() << std::endl;

	return nsegments;
}

core::Size
BackrubMover::add_mainchain_segments( core::pose::Pose const & pose )
{
	utility::vector1<core::Size> resnums(pivot_residues(pose));
	if ( resnums.size() == 0 ) {
		// use all residues if none defined
		for ( core::Size i = 1; i <= get_input_pose()->size(); ++i ) resnums.push_back(i);
	}

	// add segments to the backrub mover
	return add_mainchain_segments(resnums, pivot_atoms_, min_atoms_, max_atoms_);
}

core::Size
BackrubMover::add_mainchain_segments_from_options( core::pose::Pose const & pose )
{
	init_with_options();

	return add_mainchain_segments( pose );
}

void
BackrubMover::optimize_branch_angles(
	Pose & pose
)
{
	for ( auto & iter : bond_angle_map_ ) {

		BackrubSegment::BondAngleKey const & bond_angle_key(iter.first);
		if ( bond_angle_key.key1().valid() && bond_angle_key.key2().valid() && bond_angle_key.key3().valid() ) {
			TR.Debug << "Optimizing angles for:" << bond_angle_key.key1() << bond_angle_key.key2() << bond_angle_key.key3()
				<< std::endl;
			branchopt_.optimize_angles(pose, bond_angle_key.key1(), bond_angle_key.key2(), bond_angle_key.key3(), preserve_detailed_balance_);
		}
	}
}

utility::vector1<core::Size>
BackrubMover::pivot_residues( core::pose::Pose const & pose ) const
{
	utility::vector1<core::Size> pivot_residues( pivot_residues_ );
	if ( pivot_residue_selector_ ) {
		pivot_residues.append( core::select::get_residues_from_subset( pivot_residue_selector_->apply( pose ) ) );
	}
	return pivot_residues;
}

void
BackrubMover::set_pivot_residues(
	utility::vector1<core::Size> const & pivot_residues
)
{
	pivot_residues_ = pivot_residues;
}

void
BackrubMover::set_pivot_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP pivot_residues
) {
	pivot_residue_selector_ = pivot_residues;
}


void
BackrubMover::set_movemap_factory( core::select::movemap::MoveMapFactoryCOP mmf ){
	movemap_factory_ = mmf;
}

void
BackrubMover::set_movemap(core::kinematics::MoveMapCOP movemap){
	pivot_residues_ = get_pivot_residues_from_movemap(movemap);
}

utility::vector1<std::string> const &
BackrubMover::pivot_atoms() const
{
	return pivot_atoms_;
}

void
BackrubMover::set_pivot_atoms(
	utility::vector1<std::string> const & pivot_atoms
)
{
	pivot_atoms_ = pivot_atoms;
}

core::Size
BackrubMover::min_atoms() const
{
	return min_atoms_;
}

void
BackrubMover::set_min_atoms(
	core::Size min_atoms
)
{
	min_atoms_ = min_atoms;
}

core::Size
BackrubMover::max_atoms() const
{
	return max_atoms_;
}

void
BackrubMover::set_max_atoms(
	core::Size max_atoms
)
{
	max_atoms_ = max_atoms;
}

core::Real
BackrubMover::max_angle_disp_4() const
{
	return max_angle_disp_4_;
}

void
BackrubMover::set_max_angle_disp_4(
	core::Real max_angle_disp_4
)
{
	max_angle_disp_4_ = max_angle_disp_4;
}

core::Real
BackrubMover::max_angle_disp_7() const
{
	return max_angle_disp_7_;
}

void
BackrubMover::set_max_angle_disp_7(
	core::Real max_angle_disp_7
)
{
	max_angle_disp_7_ = max_angle_disp_7;
}

core::Real
BackrubMover::max_angle_disp_slope() const
{
	return max_angle_disp_slope_;
}

void
BackrubMover::set_max_angle_disp_slope(
	core::Real max_angle_disp_slope
)
{
	max_angle_disp_slope_ = max_angle_disp_slope;
}

bool
BackrubMover::custom_angle() const
{
	return custom_angle_;
}

void
BackrubMover::set_custom_angle(
	bool custom_angle
)
{
	custom_angle_ = custom_angle;
}

bool
BackrubMover::preserve_detailed_balance() const
{
	return preserve_detailed_balance_;
}

void
BackrubMover::set_preserve_detailed_balance(
	bool preserve_detailed_balance
)
{
	preserve_detailed_balance_ = preserve_detailed_balance;
}

bool
BackrubMover::require_mm_bend() const
{
	return require_mm_bend_;
}

void
BackrubMover::set_require_mm_bend(
	bool require_mm_bend
)
{
	require_mm_bend_ = require_mm_bend;
}

utility::vector1<core::id::TorsionID_Range>
BackrubMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}

utility::vector1<core::id::DOF_ID_Range>
BackrubMover::dof_id_ranges(
	core::pose::Pose & pose
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	// code for declaring branching atom ranges hasn't been written yet
	runtime_assert(preserve_detailed_balance_);

	std::set<core::id::DOF_ID_Range> range_set;

	for ( core::Size i = 1; i <= segments_.size(); ++i ) {

		BackrubSegment & segment(segments_[i]);

		// get references to Atom tree atoms
		kinematics::tree::AtomCOP start_atom(pose.atom_tree().atom(segment.start_atomid()).get_self_ptr());
		kinematics::tree::AtomCOP start_atom1(pose.atom_tree().atom(segment.start_atomid1()).get_self_ptr());
		kinematics::tree::AtomCOP start_atom2(pose.atom_tree().atom(segment.start_atomid2()).get_self_ptr());
		kinematics::tree::AtomCOP end_atom(pose.atom_tree().atom(segment.end_atomid()).get_self_ptr());
		kinematics::tree::AtomCOP end_atom1(end_atom->get_nonjump_atom(0));
		kinematics::tree::AtomCOP end_atom2(end_atom1 ? end_atom1->get_nonjump_atom(0) : kinematics::tree::AtomCOP( nullptr ) );

		// Only proceed if stub_atom2 != start_atom1
		if ( start_atom->stub_atom2()->id() != segment.start_atomid1() ) {

			// get overall bond angle parameters for the mainchain bond angle
			Real start_Ktheta(0);
			Real start_theta0(0);
			Real start_energy0(0);
			branchopt_.overall_params(pose, start_atom->stub_atom2()->id(), segment.start_atomid(), segment.start_atomid1(),
				start_Ktheta, start_theta0, start_energy0);

			// limit bond angles to those having probabilites greater than 5% at a kT of 0.6
			Real start_angle_dev(std::sqrt(-0.6*std::log(0.05)/start_Ktheta));

			range_set.insert(id::DOF_ID_Range(id::DOF_ID(segment.start_atomid1(), id::THETA),
				pi - start_theta0 - start_angle_dev,
				pi - start_theta0 + start_angle_dev));

			id::DOF_ID const start_dofid1(id::DOF_ID(start_atom1->id(), id::PHI));
			range_set.insert(id::DOF_ID_Range(start_dofid1, -pi, pi));

			id::DOF_ID const start_dofid2(id::DOF_ID(start_atom2->parent()->get_nonjump_atom(0)->id(), id::PHI));
			range_set.insert(id::DOF_ID_Range(start_dofid2, -pi, pi));
		}

		// Only proceed if end_atom1 exists
		if ( end_atom1 ) {

			// get overall bond angle parameters for the mainchain bond angle
			Real end_Ktheta(0);
			Real end_theta0(0);
			Real end_energy0(0);
			branchopt_.overall_params(pose, end_atom->parent()->id(), segment.end_atomid(), end_atom1->id(),
				end_Ktheta, end_theta0, end_energy0);

			// limit bond angles to those having probabilites greater than 5% at a kT of 0.6
			Real end_angle_dev(std::sqrt(-0.6*std::log(0.05)/end_Ktheta));

			range_set.insert(id::DOF_ID_Range(id::DOF_ID(end_atom1->id(), id::THETA),
				pi - end_theta0 - end_angle_dev,
				pi - end_theta0 + end_angle_dev));

			// Set torsion angle 1 by incrementing the oldest sibling PHI DOF
			id::DOF_ID const end_dofid1(id::DOF_ID(end_atom1->id(), id::PHI));
			range_set.insert(id::DOF_ID_Range(end_dofid1, -pi, pi));

			// Only proceed if end_atom2 exists
			if ( end_atom2 ) {

				// Set torsion angle 2 by incrementing the oldest sibling PHI DOF
				id::DOF_ID const end_dofid2(id::DOF_ID(end_atom2->id(), id::PHI));
				range_set.insert(id::DOF_ID_Range(end_dofid2, -pi, pi));
			}
		}
	}

	utility::vector1<core::id::DOF_ID_Range> range_vector(range_set.begin(), range_set.end());

	return range_vector;
}

Real
BackrubMover::random_angle(
	Pose & pose,
	core::Size segment_id,
	utility::vector0<Real> & start_constants,
	utility::vector0<Real> & end_constants
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	BackrubSegment & segment(segments_[segment_id]);

	// get references to Atom tree atoms
	kinematics::tree::AtomCOP start_atom(pose.atom_tree().atom(segment.start_atomid()).get_self_ptr());
	kinematics::tree::AtomCOP start_atom1(pose.atom_tree().atom(segment.start_atomid1()).get_self_ptr());
	kinematics::tree::AtomCOP start_atom2(pose.atom_tree().atom(segment.start_atomid2()).get_self_ptr());
	kinematics::tree::AtomCOP end_atom(pose.atom_tree().atom(segment.end_atomid()).get_self_ptr());
	kinematics::tree::AtomCOP end_atom1(end_atom->get_nonjump_atom(0));
	kinematics::tree::AtomCOP end_atom2(end_atom1 ? end_atom1->get_nonjump_atom(0) : kinematics::tree::AtomCOP(nullptr) );

	//TR << "Start Atom:" << segment.start_atomid() << std::endl;
	//TR << "End Atom:" << segment.end_atomid() << std::endl;

	numeric::IntervalSet<Real> tau_intervals_bond_angle(-pi, pi);

	// Only proceed if stub_atom2 != start_atom1
	if ( start_atom->stub_atom2()->id() != segment.start_atomid1() ) {

		// get overall bond angle parameters for the mainchain bond angle
		Real start_Ktheta(0);
		Real start_theta0(0);
		Real start_energy0(0);
		branchopt_.overall_params(pose, start_atom->stub_atom2()->id(), segment.start_atomid(), segment.start_atomid1(),
			start_Ktheta, start_theta0, start_energy0);

		// limit bond angles to those having probabilites greater than 5% at a kT of 0.6
		Real start_angle_dev(std::sqrt(-0.6*std::log(0.05)/start_Ktheta));

		// get bond angle constraint interval & constants for analytically calculating updated DOF values
		backrub_rotation_constants(start_atom->stub_atom3(), start_atom->stub_atom2(), start_atom, start_atom1, start_atom2,
			end_atom, start_constants, start_theta0 - start_angle_dev, start_theta0 + start_angle_dev,
			&tau_intervals_bond_angle);
		//TR << "Start Bond Angle Intervals: " << tau_intervals_bond_angle << std::endl;
	}

	// Only proceed if end_atom1 exists
	if ( end_atom1 ) {

		// get overall bond angle parameters for the mainchain bond angle
		Real end_Ktheta(0);
		Real end_theta0(0);
		Real end_energy0(0);
		branchopt_.overall_params(pose, end_atom->parent()->id(), segment.end_atomid(), end_atom1->id(),
			end_Ktheta, end_theta0, end_energy0);

		// limit bond angles to those having probabilites greater than 5% at a kT of 0.6
		Real end_angle_dev(std::sqrt(-0.6*std::log(0.05)/end_Ktheta));

		// get bond angle constraint interval & constants for analytically calculating updated DOF values
		numeric::IntervalSet<Real> tau_intervals_bond_angle_end;
		backrub_rotation_constants(end_atom->parent()->parent(), end_atom->parent(), end_atom, end_atom1, end_atom2,
			start_atom, end_constants, end_theta0 - end_angle_dev, end_theta0 + end_angle_dev,
			&tau_intervals_bond_angle_end);
		//TR << "End Bond Angle Intervals: " << tau_intervals_bond_angle_end << std::endl;

		// calculate overall bond angle constraining interval
		tau_intervals_bond_angle = tau_intervals_bond_angle & tau_intervals_bond_angle_end;
		//TR << "Bond Angle Intervals: " << tau_intervals_bond_angle << std::endl;
	}

	// calculate rotation angle interval, the overall constraining interval, and the total length
	Real const angle_disp(segments_[segment_id].max_angle_disp());
	numeric::IntervalSet<Real> tau_intervals_rotation_angle(-angle_disp, angle_disp);
	//TR << "Rotation Angle Intervals: " << tau_intervals_rotation_angle << std::endl;
	numeric::IntervalSet<Real> tau_intervals(tau_intervals_bond_angle & tau_intervals_rotation_angle);
	//TR << "Intervals: " << tau_intervals << std::endl;
	Real const tau_intervals_length(tau_intervals.length());

	// return early if there are no angles to pick from
	if ( tau_intervals_length == 0 ) {
		return 0;
	}

	Real angle(0);
	numeric::IntervalSet<Real> tau_intervals_rotation_angle_p;
	Real const threshold(numeric::random::rg().uniform());
	for ( core::Size i = 0; i < 100000; i++ ) {

		angle = tau_intervals.random_point(numeric::random::rg());

		Real min_angle_p = angle - angle_disp;
		Real max_angle_p = angle + angle_disp;
		min_angle_p = numeric::nearest_angle_radians(min_angle_p, 0.);
		max_angle_p = numeric::nearest_angle_radians(max_angle_p, 0.);
		if ( min_angle_p < max_angle_p ) {
			tau_intervals_rotation_angle_p.set(min_angle_p, max_angle_p);
		} else {
			tau_intervals_rotation_angle_p.set(-pi, max_angle_p, min_angle_p, pi);
		}
		//TR << "Rotation Angle Intervals': " << tau_intervals_rotation_angle_p << std::endl;
		//TR << "Intervals': " << (tau_intervals_bond_angle & tau_intervals_rotation_angle_p) << std::endl;

		Real tau_intervals_length_p = (tau_intervals_bond_angle & tau_intervals_rotation_angle_p).length();

		if ( tau_intervals_length_p == 0 || tau_intervals_length/tau_intervals_length_p >= threshold ) break;
	}

	return angle;
}

/// @details
/// The code currently does not do any optimization of branching atom bond angles.
void
BackrubMover::rotate_segment(
	Pose & pose,
	core::Size segment_id,
	Real angle,
	utility::vector0<Real> & start_constants,
	utility::vector0<Real> & end_constants
)
{
	Real static const pi(numeric::NumericTraits<Real>::pi());

	BackrubSegment & segment(segments_[segment_id]);

	// get references to Atom tree atoms
	kinematics::tree::AtomCOP start_atom(pose.atom_tree().atom(segment.start_atomid()).get_self_ptr());
	kinematics::tree::AtomCOP start_atom1(pose.atom_tree().atom(segment.start_atomid1()).get_self_ptr());
	kinematics::tree::AtomCOP start_atom2(pose.atom_tree().atom(segment.start_atomid2()).get_self_ptr());
	kinematics::tree::AtomCOP end_atom(pose.atom_tree().atom(segment.end_atomid()).get_self_ptr());
	kinematics::tree::AtomCOP end_atom1(end_atom->get_nonjump_atom(0));
	kinematics::tree::AtomCOP end_atom2(end_atom1 ? end_atom1->get_nonjump_atom(0) : kinematics::tree::AtomCOP(nullptr) );

	/*
	PointPosition end_atom_xyz(pose.xyz(end_atom->id()));
	PointPosition end_atom1_xyz;
	if (end_atom1) end_atom1_xyz = pose.xyz(end_atom1->id());
	PointPosition end_atom2_xyz;
	if (end_atom2) end_atom2_xyz = pose.xyz(end_atom2->id());
	*/

	// Get constants for analytically calculating updated DOF values
	// Only recalculate the constants if they haven't already been calculated
	if ( start_constants.size() == 0 ) {
		backrub_rotation_constants(start_atom->stub_atom3(), start_atom->stub_atom2(), start_atom, start_atom1, start_atom2,
			end_atom, start_constants);
	}

	// Get angles before and after the rotation
	Real start_0_bondangle(0), start_0_torsion1(0), start_0_torsion2(0);
	Real start_bondangle(0), start_torsion1(0), start_torsion2(0);
	backrub_rotation_angles(start_constants, 0, start_0_bondangle, start_0_torsion1, start_0_torsion2);
	backrub_rotation_angles(start_constants, angle, start_bondangle, start_torsion1, start_torsion2);

	//TR << "Start: Delta bondangle: " << start_bondangle-start_0_bondangle
	//   << " Delta torsion1: " << start_torsion1 - start_0_torsion1
	//  << " Delta torsion2: " << start_torsion2 - start_0_torsion2 << std::endl;

	// Set bond angles directly
	pose.set_dof(id::DOF_ID(segment.start_atomid1(), id::THETA), pi - start_bondangle);

	// Set torsion angles by incrementing the oldest sibling PHI DOF
	//id::DOF_ID const start_dofid1(id::DOF_ID(start_atom1->parent()->get_nonjump_atom(0)->id(), id::PHI));
	id::DOF_ID const start_dofid1(id::DOF_ID(start_atom1->id(), id::PHI));
	pose.set_dof(start_dofid1, pose.dof(start_dofid1) - start_0_torsion1 + start_torsion1);
	id::DOF_ID const start_dofid2(id::DOF_ID(start_atom2->parent()->get_nonjump_atom(0)->id(), id::PHI));
	//id::DOF_ID const start_dofid2(id::DOF_ID(start_atom2->id(), id::PHI));
	pose.set_dof(start_dofid2, pose.dof(start_dofid2) - start_0_torsion2 + start_torsion2);

	//TR << start_atom1->id() << "\t" << start_atom1->parent()->get_nonjump_atom(0)->id() << std::endl;
	//TR << start_atom2->id() << "\t" << start_atom2->parent()->get_nonjump_atom(0)->id() << std::endl;

	// Only proceed if stub_atom2 != start_atom1
	if ( start_atom->stub_atom2()->id() != segment.start_atomid1() ) {

		// optimize branching atom angles around the start pivot
		if ( !preserve_detailed_balance_ ) {
			branchopt_.optimize_angles(pose, start_atom->stub_atom2()->id(), segment.start_atomid(), segment.start_atomid1());
		}
	}

	// Only proceed if end_atom1 exists
	if ( end_atom1 ) {

		// Get constants for analytically calculating updated DOF values
		// Only recalculate the constants if they haven't already been calculated
		if ( end_constants.size() == 0 ) {
			backrub_rotation_constants(end_atom->parent()->parent(), end_atom->parent(), end_atom, end_atom1, end_atom2,
				start_atom, end_constants);
		}

		Real end_bondangle(0), end_torsion1(0), end_torsion2(0);
		Real end_0_bondangle(0), end_0_torsion1(0), end_0_torsion2(0);
		backrub_rotation_angles(end_constants, 0, end_0_bondangle, end_0_torsion1, end_0_torsion2);
		backrub_rotation_angles(end_constants, angle, end_bondangle, end_torsion1, end_torsion2);

		//TR << "End: Delta bondangle: " << end_bondangle-end_0_bondangle
		//   << " Delta torsion1: " << end_torsion1 - end_0_torsion1
		//  << " Delta torsion2: " << end_torsion2 - end_0_torsion2 << std::endl;

		// Set bond angles directly
		pose.set_dof(id::DOF_ID(end_atom1->id(), id::THETA), pi - end_bondangle);

		// Set torsion angle 1 by incrementing the oldest sibling PHI DOF
		id::DOF_ID const end_dofid1(id::DOF_ID(end_atom1->id(), id::PHI));
		pose.set_dof(end_dofid1, pose.dof(end_dofid1) - end_0_torsion1 + end_torsion1);

		// Only proceed if end_atom2 exists
		if ( end_atom2 ) {

			// Set torsion angle 2 by incrementing the oldest sibling PHI DOF
			id::DOF_ID const end_dofid2(id::DOF_ID(end_atom2->id(), id::PHI));
			pose.set_dof(end_dofid2, pose.dof(end_dofid2) - end_0_torsion2 + end_torsion2);
		}

		// optimize branching atom angles around the end pivot
		if ( !preserve_detailed_balance_ ) {
			branchopt_.optimize_angles(pose, end_atom->parent()->id(), segment.end_atomid(), end_atom1->id());
		}
	}

	/*
	Real end_atom_distsq(pose.xyz(end_atom->id()).distance_squared(end_atom_xyz));
	Real end_atom1_distsq(end_atom1 ? pose.xyz(end_atom1->id()).distance_squared(end_atom1_xyz) : 0);
	Real end_atom2_distsq(end_atom2 ? pose.xyz(end_atom2->id()).distance_squared(end_atom2_xyz) : 0);

	if (end_atom_distsq > 1e-4 || end_atom1_distsq > 1e-4 || end_atom2_distsq > 1e-4) {

	TR.Error << start_atom->atom_id() << "\t" << end_atom->atom_id() << std::endl;
	TR.Error << end_atom_distsq << "\t" << end_atom1_distsq << "\t" << end_atom2_distsq << std::endl;
	//runtime_assert(false);
	}
	*/
}

core::Size
BackrubMover::next_segment_id() const
{
	return next_segment_id_;
}

void
BackrubMover::set_next_segment_id(
	core::Size next_segment_id
)
{
	next_segment_id_ = next_segment_id;
}

core::Size
BackrubMover::last_segment_id() const
{
	return last_segment_id_;
}

std::string
BackrubMover::last_start_atom_name() const
{
	return last_start_atom_name_;
}

std::string
BackrubMover::last_end_atom_name() const
{
	return last_end_atom_name_;
}

core::Real
BackrubMover::next_angle() const
{
	return next_angle_;
}

void
BackrubMover::set_next_angle(
	core::Real next_angle
)
{
	next_angle_ = next_angle;
}

core::Real
BackrubMover::last_angle() const
{
	return last_angle_;
}

/// @details
/// All move types are prefixed with "br". Sections are divided by underscores.
/// The next two sections indicates the names of the atoms at the start and end
/// of the backrub segment. The last section gives the size of the segment.
void
BackrubMover::update_type()
{
	std::stringstream mt;

	std::string start_atom_name(last_start_atom_name_);
	std::string end_atom_name(last_end_atom_name_);
	ObjexxFCL::strip(start_atom_name);
	ObjexxFCL::strip(end_atom_name);

	mt << "br_" << start_atom_name << "_" << end_atom_name << "_"
		<< std::setw(2) << std::setfill('0') << segments_[last_segment_id_].size();

	std::string new_type(mt.str());
	type(new_type);
}

std::string BackrubMover::get_name() const {
	return mover_name();
}

std::string BackrubMover::mover_name() {
	return "Backrub";
}

void BackrubMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	auto ct_gen = complex_type_generator_for_backrub_mover( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Performs purely local moves using rotations around axes defined by two backbone atoms" )
		.write_complex_type_to_schema( xsd );
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
BackrubMover::complex_type_generator_for_backrub_mover( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute(
		"pivot_residues", xs_string,
		"residues for which contiguous stretches can contain segments "
		"(comma separated) can use PDB numbers ([resnum][chain]) or "
		"absolute Rosetta numbers (integer)");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"pivot_atoms", xs_string,
		"main chain atoms usable as pivots (comma separated)",
		"CA");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"min_atoms", xsct_non_negative_integer,
		"minimum backrub segment size (atoms)",
		"3");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_atoms", xsct_non_negative_integer,
		"maximum backrub segment size (atoms)",
		"34");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_angle_disp_4", xsct_real,
		"maximum angular displacement for 4 atom segments (radians)",
		std::to_string(numeric::conversions::radians(40.)));
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_angle_disp_7", xsct_real,
		"maximum angular displacement for 7 atom segments (radians)",
		std::to_string(numeric::conversions::radians(20.)));
	attlist + XMLSchemaAttribute::attribute_w_default(
		"max_angle_disp_slope", xsct_real,
		"maximum angular displacement slope for other atom segments (radians)",
		std::to_string(numeric::conversions::radians(-1./3.)));
	attlist + XMLSchemaAttribute::attribute_w_default(
		"preserve_detailed_balance", xsct_rosetta_bool,
		"if set to true, does not change branching atom angles during apply and "
		"sets ideal branch angles during initialization if used with MetropolisHastings",
		"false");
	attlist + XMLSchemaAttribute::attribute_w_default(
		"require_mm_bend", xsct_rosetta_bool,
		"if true and used with MetropolisHastings, will exit if mm_bend is not in the score function",
		"true");

	core::select::movemap::attributes_for_parse_movemap_factory_default_attr_name( attlist,
		"The name of the pre-defined MoveMapfactory that will be used to alter the"
		" default behaviour of the MoveMap. By default, all backbone, chi, and jump"
		" DOFs are allowed to change. A MoveMapFactory can be used to change which "
		"of those DOFs are actually enabled. The provision of MoveMapFactory can allow"
		"dynamic allocation as the factory can take residues from the residue selector."
		"Be warned that combining a MoveMapFactory with a Movemap can result in "
		"unexpected behaviour. The Movemap provided as a subelement of this element "
		"will be generated, and then the DoF modifications specified in the MoveMap "
		"Factory will be applied afterwards. Note that if residues are defined with "
		"the pivot_residues tag, they will override residues defined by both the movemap"
		" or the movemap factory");

	XMLSchemaSimpleSubelementList subelements;
	rosetta_scripts::append_subelement_for_parse_movemap_factory_legacy( xsd, subelements );

	//XMLSchemaComplexTypeGeneratorOP ct_gen( utility::pointer::make_shared< XMLSchemaComplexTypeGenerator >() );
	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	ct_gen->complex_type_naming_func( & moves::complex_type_name_for_mover )
		.add_attributes( attlist )
		.set_subelements_repeatable( subelements )
		.add_optional_name_attribute();
	return ct_gen;
}

std::string BackrubMoverCreator::keyname() const {
	return BackrubMover::mover_name();
}

protocols::moves::MoverOP
BackrubMoverCreator::create_mover() const {
	return utility::pointer::make_shared< BackrubMover >();
}

void BackrubMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BackrubMover::provide_xml_schema( xsd );
}


/// @details
/// PM1 & PM2 are the parent and grandparent atoms (respectively) of the pivot
/// atom, P. PP1 and PP2 are the child and grandchiled atoms (respectively) of
/// the pivot atom. PM2 and PP2 are optional and may be NULL. REF is the other
/// (reference) pivot atom that defines the rotation axis.
///
/// The first 9 constants returned represent A1-A3, B1-B3, & C1-C3 as described
/// in Betancourt 2005. The last 6 constants allow calculation of the signs of
/// phi and psi. They could be called B4-B6 & C4-C6. For a given tau angle, phi
/// is negative if the following is true:
///
/// B4 < B5 * cos(B6 + tau)
///
/// Similarly, psi is negative if the following is true:
///
/// C4 < C5 * cos(C6 + tau)
void
backrub_rotation_constants(
	core::kinematics::tree::AtomCOP PM2_atom,
	core::kinematics::tree::AtomCOP PM1_atom,
	core::kinematics::tree::AtomCOP P_atom,
	core::kinematics::tree::AtomCOP PP1_atom,
	core::kinematics::tree::AtomCOP PP2_atom,
	core::kinematics::tree::AtomCOP REF_atom,
	utility::vector0<Real> & constants,
	core::Real const alpha_min, // = 0
	core::Real const alpha_max, // = 0
	numeric::IntervalSet<core::Real> *tau_intervals // = NULL
)
{
	using numeric::conversions::radians;
	using numeric::sin_cos_range;
	using numeric::in_sin_cos_range;

	Real static const pi(numeric::NumericTraits<Real>::pi());

	constants.resize(15);
	for ( int i = 0; i < 15; i++ ) {
		constants[i] = 0;
	}

	PointPosition const & N(PM1_atom->xyz());
	PointPosition const & CA(P_atom->xyz());
	PointPosition const & C(PP1_atom->xyz());
	PointPosition const & CAref(REF_atom->xyz());

	PointPosition const V( (CA - N).normalize() );
	PointPosition const U( (C - CA).normalize() );
	PointPosition const R( (CAref - CA).normalize() );

	Real const sigma_v = std::acos(sin_cos_range(-dot(V, R))); // N->CA
	Real const sigma_u = std::acos(sin_cos_range(-dot(U, R))); // CA->C

	Real const alpha = std::acos(sin_cos_range(-dot(V, U))); // N-CA-C

	Real const sin_sigma_v = std::sin(sigma_v);
	Real const cos_sigma_v = std::cos(sigma_v);
	Real const sin_sigma_u = std::sin(sigma_u);
	Real const cos_sigma_u = std::cos(sigma_u);
	Real const cos_alpha = std::cos(alpha);

	Real const tau_u = ((dot(U, cross(V, R)) < 0) ? -1 : 1) *
		std::acos(sin_cos_range((cos_alpha + cos_sigma_v * cos_sigma_u) /
		(sin_sigma_v * sin_sigma_u)));

	Real const a1 = constants[0] = sin_sigma_v * sin_sigma_u * std::cos(tau_u);
	Real const b1 = constants[1] = -sin_sigma_v * sin_sigma_u * std::sin(tau_u);
	Real const c1 = constants[2] = -cos_sigma_v * cos_sigma_u;

	if ( PM2_atom ) {

		PointPosition const & Cm1(PM2_atom->xyz());

		PointPosition const W( (N - Cm1).normalize() );

		Real const sigma_w = std::acos(sin_cos_range(-dot(W, R))); // C(-1)->N

		Real const gamma = std::acos(sin_cos_range(-dot(W, V))); // C(-1)-N-CA

		Real const sin_sigma_w = std::sin(sigma_w);
		Real const cos_sigma_w = std::cos(sigma_w);
		Real const sin_gamma = std::sin(gamma);
		Real const cos_gamma = std::cos(gamma);

		Real const tau_w = ((dot(U, cross(W, R)) < 0) ? -1 : 1) *
			std::acos(sin_cos_range((-dot(W, U) + cos_sigma_w * cos_sigma_u) /
			(sin_sigma_w * sin_sigma_u)));

		constants[3] = (sin_sigma_w * sin_sigma_u * std::cos(tau_w) + a1 * cos_gamma) / sin_gamma;
		constants[4] = (-sin_sigma_w * sin_sigma_u * std::sin(tau_w) + b1 * cos_gamma) / sin_gamma;
		constants[5] = (-cos_sigma_u * cos_sigma_w + c1 * cos_gamma) / sin_gamma;

		PointPosition const Wn( (cross(W, V)).normalize() ); // V-W plane normal

		Real const sigma_wn = std::acos(sin_cos_range(-dot(Wn, R)));

		Real const sin_sigma_wn = std::sin(sigma_wn);
		Real const cos_sigma_wn = std::cos(sigma_wn);

		Real const tau_wn = ((dot(U, cross(Wn, R)) < 0) ? -1 : 1) *
			std::acos(sin_cos_range((-dot(Wn, U) + cos_sigma_wn * cos_sigma_u) /
			(sin_sigma_wn * sin_sigma_u)));

		constants[9] = cos_sigma_wn * cos_sigma_u;
		constants[10] = sin_sigma_wn * sin_sigma_u;
		constants[11] = tau_wn;
	}

	if ( PP2_atom ) {

		PointPosition const & Np1(PP2_atom->xyz());

		PointPosition const W1( (Np1 - C).normalize() );

		Real const sigma_w1 = std::acos(sin_cos_range(-dot(W1, R))); // C->N(+1)

		Real const beta = std::acos(sin_cos_range(-dot(U, W1))); // CA-C-N(+1)

		Real const sin_sigma_w1 = std::sin(sigma_w1);
		Real const cos_sigma_w1 = std::cos(sigma_w1);
		Real const sin_beta = std::sin(beta);
		Real const cos_beta = std::cos(beta);

		Real const tau_w1 = ((dot(W1, cross(V, R)) < 0) ? -1 : 1) *
			std::acos(sin_cos_range((-dot(V, W1) + cos_sigma_v * cos_sigma_w1) /
			(sin_sigma_v * sin_sigma_w1)));

		constants[6] = (sin_sigma_v * sin_sigma_w1 * std::cos(tau_w1) + a1 * cos_beta) / sin_beta;
		constants[7] = (-sin_sigma_v * sin_sigma_w1 * std::sin(tau_w1) + b1 * cos_beta) / sin_beta;
		constants[8] = (-cos_sigma_v * cos_sigma_w1 + c1 * cos_beta) / sin_beta;

		PointPosition const W1n( (cross(U, W1)).normalize() ); // U-W1 plane normal

		Real const sigma_w1n = std::acos(sin_cos_range(-dot(W1n, R)));

		Real const sin_sigma_w1n = std::sin(sigma_w1n);
		Real const cos_sigma_w1n = std::cos(sigma_w1n);

		Real const tau_w1n = ((dot(W1n, cross(V, R)) < 0) ? -1 : 1) *
			std::acos(sin_cos_range((-dot(V, W1n) + cos_sigma_v * cos_sigma_w1n) /
			(sin_sigma_v * sin_sigma_w1n)));

		constants[12] = cos_sigma_w1n * cos_sigma_v;
		constants[13] = sin_sigma_w1n * sin_sigma_v;
		constants[14] = tau_w1n;
	}

	if ( tau_intervals != nullptr && alpha_min >= 0 && alpha_max > alpha_min ) {

		if ( alpha_min == 0 && alpha_max == pi ) {

			// easiest base case: unconstrained
			tau_intervals->set(-pi, pi);

		} else {

			numeric::IntervalSet<Real> min_intervals;
			numeric::IntervalSet<Real> max_intervals;

			Real const min_value = (cos_sigma_v * cos_sigma_u + std::cos(alpha_min)) / (sin_sigma_v * sin_sigma_u);

			if ( in_sin_cos_range(min_value, 0.00000001) ) {

				Real const acos_min = std::acos(sin_cos_range(min_value));
				Real tau1 = -acos_min - tau_u;
				Real tau2 = acos_min - tau_u;

				tau1 = numeric::nearest_angle_radians(tau1, 0.);
				tau2 = numeric::nearest_angle_radians(tau2, 0.);

				if ( tau1 > tau2 ) {
					Real const temp = tau1;
					tau1 = tau2;
					tau2 = temp;
				}

				Real const alpha_mid = std::acos(sin_cos_range(sin_sigma_v * sin_sigma_u *
					std::cos(tau_u + (tau1 + tau2)/2) - cos_sigma_v * cos_sigma_u));

				if ( alpha_mid > alpha_min ) {
					min_intervals.set(tau1, tau2);
				} else {
					min_intervals.set(-pi, tau1, tau2, pi);
				}
			} else if ( alpha > alpha_min ) {
				min_intervals.set(-pi, pi);
			} else {
				//TR.Error << "No tau angles meet minimum alpha bond angle" << std::endl;
			}

			/*
			if (min_intervals.endpoints().size() && !min_intervals.is_inside(0)) {
			std::cout << "Min not inside " << tau_u << " "
			<< (in_sin_cos_range(min_value) ? std::acos(sin_cos_range(min_value)) : -1) << " "
			<< min_value << std::endl;
			std::cout << min_intervals << std::endl;
			}
			*/

			Real const max_value = (cos_sigma_v * cos_sigma_u + std::cos(alpha_max)) / (sin_sigma_v * sin_sigma_u);

			if ( in_sin_cos_range(max_value, 0.00000001) ) {

				Real const acos_max = std::acos(sin_cos_range(max_value));
				Real tau1 = -acos_max - tau_u;
				Real tau2 = acos_max - tau_u;

				tau1 = numeric::nearest_angle_radians(tau1, 0.);
				tau2 = numeric::nearest_angle_radians(tau2, 0.);

				if ( tau1 > tau2 ) {
					Real const temp = tau1;
					tau1 = tau2;
					tau2 = temp;
				}

				Real const alpha_mid = std::acos(sin_cos_range(sin_sigma_v * sin_sigma_u *
					std::cos(tau_u + (tau1 + tau2)/2) - cos_sigma_v * cos_sigma_u));

				if ( alpha_mid < alpha_max ) {
					max_intervals.set(tau1, tau2);
				} else {
					max_intervals.set(-pi, tau1, tau2, pi);
				}
			} else if ( alpha < alpha_max ) {
				max_intervals.set(-pi, pi);
			} else {
				//TR.Error << "No tau angles meet maximum alpha bond angle" << std::endl;
			}

			/*
			if (max_intervals.endpoints().size() && !max_intervals.is_inside(0)) {
			std::cout << "Max not inside " << tau_u << " "
			<< (in_sin_cos_range(max_value) ? std::acos(sin_cos_range(max_value)) : -1) << " "
			<< max_value << std::endl;
			std::cout << max_intervals << std::endl;
			}
			*/

			*tau_intervals = min_intervals & max_intervals;
		}
	}
}

/// @details
/// tau is the angular displacement
void
backrub_rotation_angles(
	utility::vector0<Real> const & constants,
	Real const tau,
	Real & bondange,
	Real & torsion1,
	Real & torsion2
)
{
	using numeric::sin_cos_range;

	Real const sin_tau = std::sin(tau);
	Real const cos_tau = std::cos(tau);

	bondange = std::acos(sin_cos_range(constants[0]*cos_tau + constants[1]*sin_tau + constants[2]));

	Real const sin_bondange = std::sin(bondange);

	torsion1 = acos(sin_cos_range((constants[3]*cos_tau + constants[4]*sin_tau + constants[5]) / sin_bondange));
	if ( constants[9] < constants[10]*std::cos(constants[11] + tau) ) torsion1 = -torsion1;

	torsion2 = acos(sin_cos_range((constants[6]*cos_tau + constants[7]*sin_tau + constants[8]) / sin_bondange));
	if ( constants[12] < constants[13]*std::cos(constants[14] + tau) ) torsion2 = -torsion2;
}

} // moves
} // protocols
