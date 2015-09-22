// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.cc
/// @brief  Helper class for generalized closure of arbitrary segments that could go through side-chains (e.g. disulfides).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.hh>
#include <protocols/generalized_kinematic_closure/GeneralizedKIC.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <core/id/AtomID.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>


#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/BBGaussianMover.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace generalized_kinematic_closure {
namespace perturber {

static THREAD_LOCAL basic::Tracer TR( "protocols.generalized_kinematic_closure.perturber.GeneralizedKICperturber" );

/// @brief Creator for GeneralizedKICperturber.
GeneralizedKICperturber::GeneralizedKICperturber():
	utility::pointer::ReferenceCount(),
	bbgmover_(),
	bin_transition_calculator_(),
	effect_(no_effect),
	inputvalues_real_(),
	residues_(),
	atoms_(),
	iterations_(1),
	must_switch_bins_(false),
	bin_("")
	//TODO -- make sure above data are copied properly when duplicating this mover.
{}

/// @brief Copy constructor for GeneralizedKICperturber.
///
GeneralizedKICperturber::GeneralizedKICperturber( GeneralizedKICperturber const &src ):
	utility::pointer::ReferenceCount(),
	bbgmover_( ), //Cloned later
	bin_transition_calculator_( ), //Cloned later
	effect_(src.effect_),
	inputvalues_real_(src.inputvalues_real_),
	residues_(src.residues_),
	atoms_(src.atoms_),
	iterations_(src.iterations_),
	must_switch_bins_(src.must_switch_bins_),
	bin_(src.bin_)
{
	if ( src.bbgmover_ ) bbgmover_ = utility::pointer::dynamic_pointer_cast< protocols::simple_moves::BBGaussianMover >(src.bbgmover_->clone());
	if ( src.bin_transition_calculator_ ) bin_transition_calculator_ = utility::pointer::dynamic_pointer_cast< core::scoring::bin_transitions::BinTransitionCalculator >(src.bin_transition_calculator_->clone());
}

/// @brief Destructor for GeneralizedKICperturber mover.
GeneralizedKICperturber::~GeneralizedKICperturber() {}

/// @brief Clone function for GeneralizedKICperturber:
/// @details Returns an owning pointer to a copy of this perturber.
GeneralizedKICperturberOP GeneralizedKICperturber::clone() const
{
	return GeneralizedKICperturberOP( new GeneralizedKICperturber(*this) );
}

/// @brief Returns the name of this class ("GeneralizedKICperturber").
std::string GeneralizedKICperturber::get_name() const{
	return "GeneralizedKICperturber";
}


/// @brief Returns the enum type for the effect of a pertuber based on a perturber name.
///        Returns unknown_effect if can't find a match for the name.
perturber_effect GeneralizedKICperturber::get_perturber_effect_from_name( std::string const &name ) const
{
	for ( core::Size i=1; i<end_of_effect_list; ++i ) {
		if ( get_perturber_effect_name( i ) == name ) return static_cast<perturber_effect>(i);
	}
	return unknown_effect;
}


/// @brief Returns the name of a perturber given the enum type.
///        Returns "unknown_effect" if no such effect exists.
std::string GeneralizedKICperturber::get_perturber_effect_name( core::Size &effect ) const
{
	std::string returnstring = "";

	switch(effect) {
	case no_effect :
		returnstring = "no_effect";
		break;
	case set_dihedral :
		returnstring = "set_dihedral";
		break;
	case set_bondangle :
		returnstring = "set_bondangle";
		break;
	case set_bondlength :
		returnstring = "set_bondlength";
		break;
	case set_backbone_bin :
		returnstring = "set_backbone_bin";
		break;
	case randomize_dihedral :
		returnstring = "randomize_dihedral";
		break;
	case randomize_alpha_backbone_by_rama :
		returnstring = "randomize_alpha_backbone_by_rama";
		break;
	case randomize_backbone_by_bins :
		returnstring = "randomize_backbone_by_bins";
		break;
	case perturb_dihedral :
		returnstring = "perturb_dihedral";
		break;
	case perturb_dihedral_bbg :
		returnstring = "perturb_dihedral_bbg";
		break;
	case perturb_backbone_by_bins :
		returnstring = "perturb_backbone_by_bins";
		break;
	case sample_cis_peptide_bond :
		returnstring = "sample_cis_peptide_bond";
		break;
	default :
		returnstring = "unknown_effect";
		break;
	}

	return returnstring;
}

/// @brief Sets the effect of this perturber.
void GeneralizedKICperturber::set_perturber_effect( perturber_effect const &effect )
{
	runtime_assert_string_msg(effect > 0 && effect < end_of_effect_list, "The perturber effect type is not recognized.");
	effect_ = effect;

	return;
}


/// @brief Sets the effect of this perturber using the perturber effect name.
///        Exits with an error message if the name is unknown.
void GeneralizedKICperturber::set_perturber_effect( std::string const &effectname )
{
	core::Size effect = get_perturber_effect_from_name(effectname);
	runtime_assert_string_msg(effect < end_of_effect_list, "Perturber effect type " + effectname + "not recognized.  Error in GeneralizedKICperturber::set_perturber_effect.");
	set_perturber_effect( (perturber_effect)effect );
	return;
}

/// @brief Initializes the BinTransitionCalculator object and loads a bin_params file.
///
void GeneralizedKICperturber::load_bin_params( std::string const &bin_params_file )
{
	using namespace core::scoring::bin_transitions;

	//Create the object, if it doesn't exist.
	if ( !bin_transition_calculator_ ) {
		if ( TR.visible() ) TR << "Creating BinTransitionCalculator." << std::endl;
		bin_transition_calculator_=BinTransitionCalculatorOP( new BinTransitionCalculator );
	}

	if ( TR.visible() ) TR << "Loading bin_params file " << bin_params_file << "." << std::endl;
	bin_transition_calculator_->load_bin_params(bin_params_file);

	return;
}

/// @brief Applies the perturbation to the vectors of desired torsions, desired angles, and desired bond lengths.
///
/// @details
///
/// @param[in] original_pose - The original pose.
/// @param[in] loop_pose - A pose consisting of just the loop to be perturbed, plus one residue on each side establishing the frame.
/// @param[in] residue_map - Mapping of (loop residue, original pose residue).
/// @param[in] tail_residue_map - Mapping of (tail residue in loop pose, original pose tail residue).
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in,out] torsions - Desired torsions for each atom; can be set or altered by the apply() function.
/// @param[in,out] bondangles - Desired bond angles for each atom; can be set or altered by the apply() function.
/// @param[in,out] bondlengths - Desired bond lengths for each atom; can be set or altered by the apply() function.
void GeneralizedKICperturber::apply(
	core::pose::Pose const &original_pose, //The original pose
	core::pose::Pose const &loop_pose, //A pose consisting of just the loop to be perturbed, plus one residue on each side establishing the frame
	utility::vector1< std::pair< core::Size, core::Size > > const &residue_map, //mapping of (loop residue, original pose residue)
	utility::vector1< std::pair< core::Size, core::Size > > const &tail_residue_map, //mapping of (tail residue in loop pose, original pose tail residue)
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1< core::Real > &torsions, //desired torsions for each atom (input/output)
	utility::vector1< core::Real > &bondangles, //desired bond angle for each atom (input/output)
	utility::vector1< core::Real > &bondlengths //desired bond length for each atom (input/output)
) const {

	//utility::vector1 < core::Size > residues_loopindexed; //The list of residues, with indices based on loop_pose.  If necessary, will be generated.
	utility::vector1 < utility::vector1 < core::id::AtomID > > AtomIDs_loopindexed; //The list of lists of AtomIDs, with indices based on loop_pose.  If necessary, will be generated.

	switch(effect_) {
	case set_dihedral :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_set_dihedral(AtomIDs_loopindexed, atomlist, inputvalues_real_, torsions, 0);
		break;
	case set_bondangle :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_set_bondangle(AtomIDs_loopindexed, atomlist, inputvalues_real_, bondangles);
		break;
	case set_bondlength :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_set_bondlength(AtomIDs_loopindexed, atomlist, inputvalues_real_, bondlengths);
		break;
	case set_backbone_bin :
		apply_set_backbone_bin( original_pose, loop_pose, residues_, atomlist, residue_map, tail_residue_map, torsions );
		break;
	case randomize_dihedral :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_set_dihedral(AtomIDs_loopindexed, atomlist, inputvalues_real_, torsions, 1); //We recycle the apply_set_dihedral() function to avoid code duplication
		break;
	case randomize_alpha_backbone_by_rama :
		apply_randomize_alpha_backbone_by_rama(original_pose, loop_pose, residues_, atomlist, residue_map, tail_residue_map, torsions);
		break;
	case randomize_backbone_by_bins :
		apply_randomize_backbone_by_bins( original_pose, loop_pose, residues_, atomlist, residue_map, tail_residue_map, torsions );
		break;
	case perturb_dihedral :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_set_dihedral(AtomIDs_loopindexed, atomlist, inputvalues_real_, torsions, 2); //We recycle the apply_set_dihedral() function to avoid code duplication
		break;
	case perturb_dihedral_bbg :
		reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
		apply_perturb_dihedral_bbg(original_pose, loop_pose, residues_, atomlist, residue_map, torsions);
		break;
	case perturb_backbone_by_bins :
		apply_perturb_backbone_by_bins( original_pose, loop_pose, residues_, atomlist, residue_map, tail_residue_map, torsions );
		break;
	case sample_cis_peptide_bond :
		apply_sample_cis_peptide_bond(loop_pose, atomlist, residues_, residue_map, tail_residue_map, torsions);
		break;
	default :
		break;
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Given an index in the original pose and a mapping from loop to pose,
/// return the index in the loop.
core::Size GeneralizedKICperturber::get_loop_index (
	core::Size const original_pose_index,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
) const {
	for ( core::Size i=1, imax=residue_map.size(); i<=imax; ++i ) {
		if ( residue_map[i].second == original_pose_index ) return residue_map[i].first;
	}

	utility_exit_with_message("Residue does not exist in loop.  Exiting from GeneralizedKICperturber::get_loop_index with error status.");

	return 0;
}

/// @brief Given a list of lists of atoms (as std::pair <residue_index, atom_name>)
/// where residue indices are based on the original pose, and a mapping of
/// original pose residue ID values to loop residue ID values, generate a new list
/// of lists of AtomIDs, where the residue indices are based on the loop pose.
void GeneralizedKICperturber::reindex_AtomIDs (
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //input
	utility::vector1 < utility::vector1 < core::id::AtomID > > &AtomIDs_reindexed, //output
	core::pose::Pose const &original_pose //input -- for reference
) const {
	using namespace core::id;
	AtomIDs_reindexed.clear();
	for ( core::Size i=1, imax=atoms_.size(); i<=imax; ++i ) {
		utility::vector1 < AtomID > newAtomIDset;
		for ( core::Size j=1, jmax=atoms_[i].size(); j<=jmax; ++j ) { //Looping through atoms in the set
			runtime_assert_string_msg( atoms_[i][j].rsd() <= original_pose.n_residue(), "GeneralizedKICperturber::reindex_AtomIDs can't find the specified residue." );
			runtime_assert_string_msg( original_pose.residue( atoms_[i][j].rsd() ).has( atoms_[i][j].atom() ),
				"GeneralizedKICperturber::reindex_AtomIDs can't find atom " + atoms_[i][j].atom() + " in the specified residue." );
			newAtomIDset.push_back(
				AtomID(
				original_pose.residue( atoms_[i][j].rsd() ).atom_index( atoms_[i][j].atom() ),
				get_loop_index(atoms_[i][j].rsd(), residue_map)
				)
			);
		} //Looping through AtomIDs in the set
		AtomIDs_reindexed.push_back(newAtomIDset);
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTIONS FOR SPECIFIC EFFECTS                              //
////////////////////////////////////////////////////////////////////////////////

/// @brief Applies a set_dihedral perturbation to a list of torsions.
/// @details  Can also be used to randomize dihedral values.
/// @param[in] dihedrallist - List of sets of atoms defining dihedrals, indexed based on the loop_pose.
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in] inputvalues_real -- Vector of input values (one value or one for each dihedral to be set).
/// @param[in,out] torsions - Desired torsions for each atom; set by this function.
/// @param[in] effect - Should the specified torsions be set (0), randomized (1), or perturbed (2)?
void GeneralizedKICperturber::apply_set_dihedral (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &dihedrallist, //List of sets of atoms defining dihedrals, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < core::Real > const &inputvalues_real,
	utility::vector1< core::Real > &torsions, //desired torsions for each atom (input/output)
	core::Size const effect //0=set, 1=randomize, 2=perturb
) const {
	//TR << "Applying set_dihedral perturbation effect." << std::endl;

	runtime_assert_string_msg( dihedrallist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_dihedral, since no atoms were provided as input." );
	if ( effect==0 ) runtime_assert_string_msg( inputvalues_real.size() > 0, "Could not set dihedral value with GeneralizedKICperturber::apply_set_dihedral, since no value for the dihedral angle was provided as input." );
	if ( effect==2 ) runtime_assert_string_msg( inputvalues_real.size() > 0, "Could not perturb dihedrals with GeneralizedKICperturber::apply_set_dihedral, since no value for the dihedral angle perturbation magnitude provided as input." );


	bool separate_values = false; //Have separate dihedral values been provided for each dihedral in the list, or are we setting everything to one value?
	if ( effect==0 || effect==2 ) {
		std::string str1=" set";
		std::string str2=" to";
		std::string str3="Setting ";
		std::string str4=" to ";
		if ( effect==2 ) {
			str1=" perturbed"; str2=" by"; str3="Perturbing "; str4=" by ";
		}
		if ( inputvalues_real.size()>=dihedrallist.size() ) {
			separate_values = true;
			if ( inputvalues_real.size()>dihedrallist.size() ) {
				if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values for set_dihedral pertruber effect exceeds the number of torsions to be" << str1 << ".  Using only the first " << dihedrallist.size() << " values." << std::endl;
			}
		} else if ( inputvalues_real.size()!=1 ) {
			separate_values = false;
			if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values does not match the number of dihedral angles specified.  All angles will be" << str1 << str2 << " the first value." << std::endl;
		} else { //inputvalues_real.size()==1 and dihedrallist.size() > 1
			if ( TR.Debug.visible() ) TR.Debug << str3 << "all specified dihedral angles" << str4 << inputvalues_real[1] << "." << std::endl;
		}
	}

	for ( core::Size i=1, imax=dihedrallist.size(); i<=imax; ++i ) { //Loop through each dihedral angle given as input
		runtime_assert_string_msg( dihedrallist[i].size()==4 || dihedrallist[i].size()==2, "Error in GeneralizedKICperturber::apply_set_dihedral().  Either two or four AtomIDs must be provided to define a dihedral angle." );

		//If four atoms are specified, that uniquely defines the dihedral angle.
		//If only two are specified, we assume that the upstream and downstream atoms of those two are the other two atoms needed to define the dihedral.
		utility::vector1 < core::id::AtomID > dihed_atoms;
		if ( dihedrallist[i].size()==4 ) dihed_atoms=dihedrallist[i];
		else if ( dihedrallist[i].size()==2 ) {
			for ( core::Size j=2, jmax=atomlist.size()-2; j<=jmax; ++j ) { //Find the upstream and downstream atoms:
				if ( dihedrallist[i][1]==atomlist[j].first && dihedrallist[i][2]==atomlist[j+1].first ) {
					dihed_atoms.push_back( atomlist[j-1].first );
					dihed_atoms.push_back( atomlist[j].first );
					dihed_atoms.push_back( atomlist[j+1].first );
					dihed_atoms.push_back( atomlist[j+2].first );
					break;
				} else if ( dihedrallist[i][2]==atomlist[j].first && dihedrallist[i][1]==atomlist[j+1].first ) {
					dihed_atoms.push_back( atomlist[j+2].first );
					dihed_atoms.push_back( atomlist[j+1].first );
					dihed_atoms.push_back( atomlist[j].first );
					dihed_atoms.push_back( atomlist[j-1].first );
					break;
				}
			}
		}

		runtime_assert_string_msg(dihed_atoms.size()==4, "Error in GeneralizedKICperturber::apply_set_dihedral().  A dihedral angle was specified that was not found in the atoms making up the loop to be closed.");

		//Find this dihedral angle in the atom list:
		core::Size torsion_index = 0;
		for ( core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j ) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if ( dihed_atoms[2]==atomlist[j].first ) {
				if ( j<jmax && dihed_atoms[3]==atomlist[j+1].first && dihed_atoms[4]==atomlist[j+2].first && dihed_atoms[1]==atomlist[j-1].first ) { //dihedral going forward
					torsion_index=j;
				} else if ( j>3 && dihed_atoms[1]==atomlist[j+1].first && dihed_atoms[3]==atomlist[j-1].first && dihed_atoms[4]==atomlist[j-2].first ) { //dihedral going backward
					torsion_index=j-1;
				}
			}
			if ( torsion_index > 0 ) break;
		}
		runtime_assert_string_msg(torsion_index > 0, "Error in GeneralizedKICperturber::apply_set_dihedral.  The dihedral angle specified was not found in the chain of atoms to be closed.");
		if ( effect==2 ) {
			torsions[torsion_index] += numeric::random::rg().gaussian() * (separate_values ? inputvalues_real[i] : inputvalues_real[1]); //Add a randomly chosen value from a gaussian distribution of specified bredth.
		} else if ( effect==1 ) { //randomizing torsions
			torsions[torsion_index] = numeric::random::rg().uniform()*360.0; //Set the desired torsion to the user-specified value.
		} else if ( effect==0 ) { //setting torsions
			torsions[torsion_index] = (separate_values ? inputvalues_real[i] : inputvalues_real[1]); //Set the desired torsion to the user-specified value.
		}
	}

	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	if ( TR.visible() ) TR.flush();

	return;
}

void GeneralizedKICperturber::init_bbgmover(
	core::pose::Pose const &loop_pose,
	utility::vector1< std::pair< core::Size, core::Size > > const &residue_map
) {
	if ( !bbgmover_ ) {
		bbgmover_ = simple_moves::BBGaussianMoverOP( new simple_moves::BBGaussianMover() );

		Size looplength = loop_pose.n_residue();
		//reset movemap, default is empty (no move)
		//the first(0) and the last(N-1) res are anchors
		//1 and N-2 res are pivots
		//all controled by residue list
		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );

		for ( Size i=1, imax=residues_.size(); i<=imax; i++ ) {
			for ( Size j=1, jmax=residue_map.size(); j<=jmax; j++ ) {
				//std::cout << residue_map[j].first << "<-->" << residue_map[j].second << std::endl;
				if ( residue_map[j].second == residues_[i] ) {
					mm->set_bb(residue_map[j].first, true);
					//std::cout << residue_map[j].first << " set to be true!" << std::endl;
				}
			}
		}

		bbgmover_->init_kic_loop( looplength, mm );
	}
}

/// @brief Applies a perturb_dihedral_bbg perturbation to a list of torsions.
/// @details  Backbone Gaussian Perturbation
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.  Note that
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_perturb_dihedral_bbg(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	using namespace protocols::generalized_kinematic_closure;

	if ( TR.visible() ) { TR << "Applying apply_perturb_dihedral_bbg perturbation effect." << std::endl;  TR.flush(); } //DELETE ME

	runtime_assert_string_msg( residues.size() > 0 , "Residues must be specified for the apply_perturb_dihedral_bbg generalized kinematic closure perturber." );

	core::Size nres = original_pose.n_residue();
	//Make a copy of the loop_pose
	core::pose::Pose loop_pose_copy = loop_pose;

	//call bbg to perturb the loop_pose_copy, then copy the torsions back
	//do we really need copying this? --in case foldtree propagate to downstream
	if ( bbgmover_ ) {
		//loop_pose_copy.dump_pdb("before.pdb");
		bbgmover_->apply(loop_pose_copy);
		bbgmover_->last_proposal_density_ratio();
		//loop_pose_copy.dump_pdb("after.pdb");
		//runtime_assert_string_msg(false, "exit");
	} else {
		runtime_assert_string_msg( residues.size() > 0 , "BBGmover is not initialized correctly!");
	}

	for ( core::Size ir=1, irmax=residues.size(); ir<=irmax; ++ir ) {
		runtime_assert_string_msg( residues[ir] <= nres, "Unable to apply apply_perturb_dihedral_bbg perturbation.  Residue list includes residues that are not in the original pose." );

		//Check for alpha amino acids:
		if ( !original_pose.residue(residues[ir]).type().is_alpha_aa() ) {
			if ( TR.Warning.visible() ) {
				TR.Warning << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_perturb_dihedral_bbg, but this residue is not an alpha-amino acid.  Skipping." << std::endl;
				TR.Warning.flush();
			}
			continue;
		}

		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		//TR << "Current loop index is " << loopindex << std::endl; TR.flush(); //DELETE ME

		//Finding and setting phi and psi:
		utility::vector1 < std::string > at1list;
		utility::vector1 < std::string > at2list;
		utility::vector1 < std::string > angname;
		at1list.push_back("N"); at1list.push_back("CA");
		at2list.push_back("CA"); at2list.push_back("C");
		angname.push_back("phi"); angname.push_back("psi");
		for ( core::Size j=1; j<=2; ++j ) {
			for ( core::Size ia=4, iamax=atomlist.size()-3; ia<=iamax; ++ia ) {
				if ( atomlist[ia].first.rsd()!=loopindex ) continue;
				if ( atomlist[ia].first.atomno() == loop_pose_copy.residue(loopindex).atom_index(at1list[j]) ) {
					core::Real angleval=0.0;
					if ( ia<iamax && atomlist[ia+1].first.rsd()==loopindex && atomlist[ia+1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going forward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia+2].first),
							angleval
						);
						torsions[ia]=angleval;
					} else if ( ia>4 && atomlist[ia-1].first.rsd()==loopindex && atomlist[ia-1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going backward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia-2].first),
							angleval
						);
						torsions[ia-1]=angleval;
					} else {
						if ( TR.Warning.visible() ) {
							TR.Warning << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_perturb_dihedral_bbg, but its " << angname[j] << " dihedral angle was not part of the chain of atoms to be closed." << std::endl;  TR.Warning.flush();
						}
					}
					break;
				}
			}
		}

	}

	return;
}

/// @brief Applies a set_bondangle perturbation to a list of bond angles.
///
/// @details
///
/// @param[in] bondanglelist - List of sets of atoms defining bond angles, indexed based on the loop_pose.
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in] inputvalues_real - List of real-valued input values (one for each bondangle to be set OR one single one).
/// @param[in,out] bondangles - Desired bond angles for each atom; set by this function.
void GeneralizedKICperturber::apply_set_bondangle (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &bondanglelist, //List of sets of atoms defining bond angles, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < core::Real > const &inputvalues_real,
	utility::vector1< core::Real > &bondangles //desired bond angles for each atom (input/output)
) const {
	//TR << "Applying set_bondangle perturbation effect." << std::endl;

	runtime_assert_string_msg( bondanglelist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondangle, since no atoms were provided as input." );
	runtime_assert_string_msg( inputvalues_real.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondangle, since no value for the bond angle was provided as input." );

	bool separate_values = false; //Have separate bond angle values been provided for each bond angle in the list, or are we setting everything to one value?
	if ( inputvalues_real.size()>=bondanglelist.size() ) {
		separate_values = true;
		if ( inputvalues_real.size()>bondanglelist.size() ) {
			if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values for set_bondangle pertruber effect exceeds the number of bond angles to be set.  Using only the first " << bondanglelist.size() << " values." << std::endl;
		}
	} else if ( inputvalues_real.size()!=1 ) {
		separate_values = false;
		if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values does not match the number of bond angles specified.  All angles will be set to the first value." << std::endl;
	} else { //inputvalues_real.size()==1 and bondanglelist.size() > 1
		if ( TR.Debug.visible() ) TR.Debug << "Setting all specified bond angles to " << inputvalues_real[1] << "." << std::endl;
	}

	for ( core::Size i=1, imax=bondanglelist.size(); i<=imax; ++i ) { //Loop through each bond angle given as input
		runtime_assert_string_msg( bondanglelist[i].size()==3 || bondanglelist[i].size()==1, "Error in GeneralizedKICperturber::apply_set_bondangle.  Either one or three AtomIDs must be provided to define a bond angle." );

		utility::vector1 < core::id::AtomID > ang_atoms;
		if ( bondanglelist[i].size()==3 ) ang_atoms = bondanglelist[i]; //If three atoms are specified for this bond angle, use them.
		else if ( bondanglelist[i].size()==1 ) { //If only one atom is specified for this bond angle, we assume that the upstream and downstream atoms in the atom list are the other two atoms defining the angle.
			for ( core::Size j=2, jmax=atomlist.size()-1; j<=jmax; ++j ) { //Find the upstream and downstream atoms:
				if ( bondanglelist[i][1]==atomlist[j].first ) {
					ang_atoms.push_back( atomlist[j-1].first );
					ang_atoms.push_back( atomlist[j].first );
					ang_atoms.push_back( atomlist[j+1].first );
					break;
				}
			}
		}

		//Find this bond angle in the atom list:
		core::Size bondangle_index = 0;
		for ( core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j ) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if ( ang_atoms[2]==atomlist[j].first ) {
				if ( ang_atoms[1]==atomlist[j-1].first && ang_atoms[3]==atomlist[j+1].first ) { //bond angle going forward
					bondangle_index=j;
				} else if ( ang_atoms[1]==atomlist[j+1].first && ang_atoms[3]==atomlist[j-1].first ) { //bond angle going backward
					bondangle_index=j;
				}
			}
			if ( bondangle_index > 0 ) break;
		}
		runtime_assert_string_msg(bondangle_index > 0, "Error in GeneralizedKICperturber::apply_set_bondangle.  The bond angle specified was not found in the chain of atoms to be closed.");
		bondangles[bondangle_index] = (separate_values ? inputvalues_real[i] : inputvalues_real[1]); //Set the desired bond angle to the user-specified value.
	}

	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	if ( TR.visible() ) TR.flush();

	return;
}

/// @brief Applies a set_bondlength perturbation to a list of bond lengths.
///
/// @details
///
/// @param[in] bondlengthlist - List of sets of atoms defining bond lengths, indexed based on the loop_pose.
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in] inputvalues_real - List of input values (a single one or one for each bondlength to be set).
/// @param[in,out] bondlengths - Desired bond lengths for each atom; set by this function.
void GeneralizedKICperturber::apply_set_bondlength (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &bondlengthlist, //List of sets of atoms defining bond lengths, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < core::Real > const &inputvalues_real,
	utility::vector1< core::Real > &bondlengths //desired bond lengths for each atom (input/output)
) const {
	//TR << "Applying set_bondlength perturbation effect." << std::endl;

	runtime_assert_string_msg( bondlengthlist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondlength, since no atoms were provided as input." );
	runtime_assert_string_msg( inputvalues_real.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondlength, since no value for the bond length was provided as input." );

	bool separate_values = false; //Have separate bond length values been provided for each bond length in the list, or are we setting everything to one value?
	if ( inputvalues_real.size()>=bondlengthlist.size() ) {
		separate_values = true;
		if ( inputvalues_real.size()>bondlengthlist.size() ) {
			if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values for set_bondlength pertruber effect exceeds the number of bond lengths to be set.  Using only the first " << bondlengthlist.size() << " values." << std::endl;
		}
	} else if ( inputvalues_real.size()!=1 ) {
		separate_values = false;
		if ( TR.Warning.visible() ) TR.Warning << "Warning! Number of input values does not match the number of bond lengths specified.  All lengths will be set to the first value." << std::endl;
	} else { //inputvalues_real.size()==1 and bondlengthlist.size() > 1
		if ( TR.Debug.visible() ) TR.Debug << "Setting all specified bond lengths to " << inputvalues_real[1] << "." << std::endl;
	}

	for ( core::Size i=1, imax=bondlengthlist.size(); i<=imax; ++i ) { //Loop through each bond length given as input
		runtime_assert_string_msg( bondlengthlist[i].size()==2, "Error in GeneralizedKICperturber::apply_set_bondlength.  Two AtomIDs must be provided to define a bond length." );

		//Find this bond length in the atom list:
		core::Size bondlength_index = 0;
		for ( core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j ) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if ( bondlengthlist[i][1]==atomlist[j].first ) {
				if ( j<jmax && bondlengthlist[i][2]==atomlist[j+1].first ) { //bond length going forward
					bondlength_index=j;
				} else if ( j>3 && bondlengthlist[i][2]==atomlist[j-1].first ) { //bond length going backward
					bondlength_index=j-1;
				}
			}
			if ( bondlength_index > 0 ) break;
		}
		runtime_assert_string_msg(bondlength_index > 0, "Error in GeneralizedKICperturber::apply_set_bondlength.  The bond length specified was not found in the chain of atoms to be closed.");
		bondlengths[bondlength_index] = (separate_values ? inputvalues_real[i] : inputvalues_real[1]); //Set the desired bond length to the user-specified value.
	}

	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	if ( TR.visible() ) TR.flush();

	return;
} //set_bondlength

/// @brief Applies a set_backbone_bin perturbation to the list of torsions.
/// @details  Sets the mainchain torsion bin for each of a list of residues, then picks random mainchain torsion
/// values from within the bin.
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in] tail_residue_map - A vector of pairs of (loop pose index of tail residue, original pose index of tail residue).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_set_backbone_bin (
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1 < std::pair < core::Size, core::Size > > const &,//tail_residue_map, //Mapping of (tail residue in loop_pose, tail residue in original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	//Initial checks:
	if ( !bin_transition_calculator_ ) {
		utility_exit_with_message( "In GeneralizedKICperturber::apply_set_backbone_bin(): The BinTransitionCalculator object has not been created!  This should be impossible -- consult a developer or a mortician." );
	}
	if ( !bin_transition_calculator_->bin_params_loaded() ) {
		utility_exit_with_message( "In GeneralizedKICperturber::apply_set_backbone_bin(): The BinTransitionCalculator object has not loaded a bin params file!" );
	}
	if ( bin_=="" ) utility_exit_with_message( "No bin has been specified!  Failing in GeneralizedKICperturber::apply_set_backbone_bin()." );

	core::Size const rescount( residues.size() );
	if ( rescount<1 ) utility_exit_with_message( "In GeneralizedKICperturber::apply_set_backbone_bin(): no residues have been specified!  Failing gracelessly." );
	core::Size const nres( original_pose.n_residue() );

	utility::vector1 <core::Size> loop_indices; //The indices in the loop of the defined residues.
	loop_indices.resize( rescount, 0 );
	for ( core::Size ir=1, irmax=rescount; ir<=irmax; ++ir ) { //Loop through all residues in the residues list
		//Confirm that this residue exists in the pose.
		if ( !(residues[ir]>0 && residues[ir]<=nres) ) utility_exit_with_message("In GeneralizedKICperturber::apply_set_backbone_bin(): At least one of the residues specified does not exist in the pose!" );
		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		loop_indices[ir]=loopindex;
		if ( ir>1 ) {
			for ( core::Size jr=1, jrmax=ir-1; jr<=jrmax; ++jr ) {
				if ( loop_indices[jr]==loopindex ) utility_exit_with_message( "In GeneralizedKICperturber::apply_set_backbone_bin(): The same residue cannot be specified twice!" ) ;
			}
		}
	}

	utility::vector1 < utility::vector1 < core::Real > > mainchain_torsions; //Vector of mainchain torsion values that will be populated by the BinTransitionCalculator object.
	bin_transition_calculator_->random_mainchain_torsions_from_bin( bin_, loop_pose.conformation(), loop_indices, mainchain_torsions ); //Get the BinTransitionCalculator to draw a random sequence of torsions from the specified bin.

	//Generate list of atoms defining the mainchain torsions to set:
	utility::vector1 < utility::vector1 <core::id::AtomID> > dihedral_list;
	utility::vector1 < core::Real > inputvalues_real;
	for ( core::Size i=1, imax=loop_indices.size(); i<=imax; ++i ) {
		for ( core::Size j=1, jmax=mainchain_torsions[i].size(); j<=jmax; ++j ) {
			utility::vector1 < core::id::AtomID > cur_atomid_list;
			cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(loop_indices[i]).mainchain_atoms()[j], loop_indices[i] ) );
			if ( j<jmax ) {
				cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(loop_indices[i]).mainchain_atoms()[j+1], loop_indices[i] ) );
			} else { //For the last torsion, we have to look up what the last mainchain atom is bonded to.
				core::conformation::Residue const &cur_res( loop_pose.conformation().residue(loop_indices[i]) );
				core::Size const other_res( cur_res.connect_map( cur_res.type().upper_connect_id() ).resid() ); //Index of the other residue that this one is connected to.
				core::Size const other_conn_id( cur_res.connect_map( cur_res.type().upper_connect_id() ).connid() ); //Index of the connection on the other residue.
				cur_atomid_list.push_back(
					core::id::AtomID (
					loop_pose.conformation().residue(other_res).residue_connect_atom_index( other_conn_id  ), //Other atom that this one is connected to at its upper connection
					other_res //Other residue that this one is connected to at its upper connection.
					)
				);
			}
			dihedral_list.push_back( cur_atomid_list );
			inputvalues_real.push_back( mainchain_torsions[i][j] );
		}
	}

	apply_set_dihedral ( dihedral_list, atomlist, inputvalues_real, torsions, 0); //Recycle this function.

	return;
}

/// @brief Applies a randomize_alpha_backbone_by_rama perturbation to the list of torsions.
/// @details This checks whether each residue is an alpha-amino acid.
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.  Note that
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in] tail_residue_map - A vector of pairs of (loop pose index of tail residue, original pose index of tail residue).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1 < std::pair < core::Size, core::Size > > const &/*tail_residue_map*/, //Mapping of (tail residue in loop_pose, tail residue in original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	using namespace protocols::generalized_kinematic_closure;

	//TR << "Applying randomize_alpha_backbone_by_rama perturbation effect." << std::endl;  TR.flush(); //DELETE ME

	runtime_assert_string_msg( residues.size() > 0 , "Residues must be specified for the randomize_alpha_backbone_by_rama generalized kinematic closure perturber." );

	core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

	core::Size nres = original_pose.n_residue();

	//Make a copy of the loop_pose
	core::pose::Pose loop_pose_copy = loop_pose; //TODO -- switch this to just a conformation.

	for ( core::Size ir=1, irmax=residues.size(); ir<=irmax; ++ir ) {
		runtime_assert_string_msg( residues[ir] <= nres, "Unable to apply randomize_alpha_backbone_by_rama perturbation.  Residue list includes residues that are not in the original pose." );

		//Check for alpha amino acids:
		if ( !original_pose.residue(residues[ir]).type().is_alpha_aa() ) {
			if ( TR.Warning.visible() ) {
				TR.Warning << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama, but this residue is not an alpha-amino acid.  Skipping." << std::endl;
				TR.Warning.flush();
			}
			continue;
		}

		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		//TR << "Current loop index is " << loopindex << std::endl; TR.flush(); //DELETE ME

		//Randomize phi and psi for this residue:
		core::Real rama_phi=0;
		core::Real rama_psi=0;
		rama.random_phipsi_from_rama(loop_pose_copy.aa(loopindex), rama_phi, rama_psi);
		// This check necessery because the aa is passed
		if ( loop_pose_copy.residue( loopindex ).has_property( "D_AA" ) ) {
			rama_phi *= -1;
			rama_psi *= -1;
		}
		general_set_phi(loop_pose_copy, loopindex, rama_phi);
		general_set_psi(loop_pose_copy, loopindex, rama_psi);

		//Finding and setting phi and psi:
		utility::vector1 < std::string > at1list;
		utility::vector1 < std::string > at2list;
		utility::vector1 < std::string > angname;
		at1list.push_back("N"); at1list.push_back("CA");
		at2list.push_back("CA"); at2list.push_back("C");
		angname.push_back("phi"); angname.push_back("psi");
		for ( core::Size j=1; j<=2; ++j ) {
			for ( core::Size ia=4, iamax=atomlist.size()-3; ia<=iamax; ++ia ) {
				if ( atomlist[ia].first.rsd()!=loopindex ) continue;
				if ( atomlist[ia].first.atomno() == loop_pose_copy.residue(loopindex).atom_index(at1list[j]) ) {
					core::Real angleval=0.0;
					if ( ia<iamax && atomlist[ia+1].first.rsd()==loopindex && atomlist[ia+1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going forward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia+2].first),
							angleval
						);
						torsions[ia]=angleval;
					} else if ( ia>4 && atomlist[ia-1].first.rsd()==loopindex && atomlist[ia-1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going backward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia-2].first),
							angleval
						);
						torsions[ia-1]=angleval;
					} else {
						if ( TR.Debug.visible() ) TR.Debug << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama, but its " << angname[j] << " dihedral angle was not part of the chain of atoms to be closed." << std::endl;
					}
					break;
				}
			}
		}

	}

	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();
	if ( TR.visible() ) TR.flush();

	return;
} //apply_randomize_alpha_backbone_by_rama

/// @brief Applies a randomize_backbone_by_bins perturbation to the list of torsions.
/// @details  This randomly assigns torsion bins based on transition probabilities of i->i+1 residues, then picks mainchain torsion angles
/// within each bin randomly or in a biased manner (if Ramachandran information is available, for example, in the case of alpha-amino acids).
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in] tail_residue_map - A vector of pairs of (loop pose index of tail residue, original pose index of tail residue).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_randomize_backbone_by_bins(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1 < std::pair < core::Size, core::Size > > const &/*tail_residue_map*/, //Mapping of (tail residue in loop_pose, tail residue in original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {

	runtime_assert_string_msg( bin_transition_calculator_, "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): The BinTransitionCalculator object has not been created!  This should be impossible -- consult a developer or a mortician." );
	runtime_assert_string_msg( bin_transition_calculator_->bin_params_loaded(), "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): The BinTransitionCalculator object has not loaded a bin params file!" );

	core::Size const rescount( residues.size() ); //Number of residues that have been defined for this perturber
	core::Size const nres( original_pose.n_residue() ); //Number of residues in the original pose

	runtime_assert_string_msg(rescount>0, "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): At least one residue must be defined!");

	//Confirm that the residues are contiguous, and conventionally joined:
	utility::vector1 <core::Size> loop_indices; //The indices in the loop of the defined residues.
	loop_indices.resize( rescount, 0 );
	for ( core::Size ir=1, irmax=rescount; ir<=irmax; ++ir ) { //Loop through all residues in the residues list
		//Confirm that this residue exists in the pose.
		runtime_assert_string_msg( residues[ir]>0 && residues[ir]<=nres, "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): At least one of the residues specified does not exist in the pose!" );

		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		loop_indices[ir]=loopindex;
		if ( ir>1 ) { for ( core::Size jr=1, jrmax=ir-1; jr<=jrmax; ++jr ) { runtime_assert_string_msg( loop_indices[jr]!=loopindex, "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): The same residue cannot be specified twice!" ); } }
		//TR << "Current loop index is " << loopindex << std::endl; TR.flush(); //DELETE ME
	}
	std::string const errmsg( "In GeneralizedKICperturber::apply_randomize_backbone_by_bins(): The selected residues are not regularly connected through mainchain polymer bonds!" );
	for ( core::Size ir=1, irmax=rescount; ir<=irmax; ++ir ) {
		if ( ir==1 ) {
			if ( irmax>1 ) runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir+1] ), errmsg + "  Problem is with start residue." );
		} else if ( ir==irmax ) {
			if ( irmax>1 ) runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir-1] ), errmsg + "  Problem is with end residue." );
		} else {
			runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir-1] ) && loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir+1] ) , errmsg + "  Problem is with a middle residue." );
		}
	}

	utility::vector1 < utility::vector1 < core::Real > > mainchain_torsions; //Vector of mainchain torsion values that will be populated by the BinTransitionCalculator object.
	bin_transition_calculator_->random_mainchain_torsions_from_bins( loop_pose.conformation(), loop_indices, mainchain_torsions ); //Get the BinTransitionCalculator to generate a random sequence of bins, and then draw a random sequence of torsions from it.

	//Generate list of atoms defining the mainchain torsions to set:
	utility::vector1 < utility::vector1 <core::id::AtomID> > dihedral_list;
	utility::vector1 < core::Real > inputvalues_real;
	for ( core::Size i=1, imax=loop_indices.size(); i<=imax; ++i ) {
		for ( core::Size j=1, jmax=mainchain_torsions[i].size(); j<=jmax; ++j ) {
			utility::vector1 < core::id::AtomID > cur_atomid_list;
			cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(loop_indices[i]).mainchain_atoms()[j], loop_indices[i] ) );
			if ( j<jmax ) {
				cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(loop_indices[i]).mainchain_atoms()[j+1], loop_indices[i] ) );
			} else { //For the last torsion, we have to look up what the last mainchain atom is bonded to.
				core::conformation::Residue const &cur_res( loop_pose.conformation().residue(loop_indices[i]) );
				core::Size const other_res( cur_res.connect_map( cur_res.type().upper_connect_id() ).resid() ); //Index of the other residue that this one is connected to.
				core::Size const other_conn_id( cur_res.connect_map( cur_res.type().upper_connect_id() ).connid() ); //Index of the connection on the other residue.
				cur_atomid_list.push_back(
					core::id::AtomID (
					loop_pose.conformation().residue(other_res).residue_connect_atom_index( other_conn_id  ), //Other atom that this one is connected to at its upper connection
					other_res //Other residue that this one is connected to at its upper connection.
					)
				);
			}
			dihedral_list.push_back( cur_atomid_list );
			inputvalues_real.push_back( mainchain_torsions[i][j] );
		}
	}

	apply_set_dihedral ( dihedral_list, atomlist, inputvalues_real, torsions, 0); //Recycle this function.

	return;
}

/// @brief Applies a perturb_backbone_by_bins perturbation to the list of torsions.
/// @details
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in] tail_residue_map - A vector of pairs of (loop pose index of tail residue, original pose index of tail residue).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_perturb_backbone_by_bins(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1 < std::pair < core::Size, core::Size > > const &/*tail_residue_map*/, //Mapping of (tail residue in loop_pose, tail residue in original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {

	core::pose::Pose loop_pose_copy(loop_pose);

	core::Size const n_iterations( iterations() ); //Number of iterations for this mover
	runtime_assert_string_msg( n_iterations>0,
		"In GeneralizedKICperturber::apply_perturb_backbon_by_bins():  The number of iterations for this mover must be greater than zero.");

	runtime_assert_string_msg( bin_transition_calculator_, "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): The BinTransitionCalculator object has not been created!  This should be impossible -- consult a developer or a mortician." );
	runtime_assert_string_msg( bin_transition_calculator_->bin_params_loaded(), "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): The BinTransitionCalculator object has not loaded a bin params file!" );

	core::Size const rescount( residues.size() ); //Number of residues that have been defined for this perturber
	core::Size const nres( original_pose.n_residue() ); //Number of residues in the original pose

	runtime_assert_string_msg(rescount>0, "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): At least one residue must be defined!");

	//Confirm that the residues are contiguous, and conventionally joined:
	utility::vector1 <core::Size> loop_indices; //The indices in the loop of the defined residues.
	loop_indices.resize( rescount, 0 );
	for ( core::Size ir=1, irmax=rescount; ir<=irmax; ++ir ) { //Loop through all residues in the residues list
		//Confirm that this residue exists in the pose.
		runtime_assert_string_msg( residues[ir]>0 && residues[ir]<=nres, "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): At least one of the residues specified does not exist in the pose!" );

		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		loop_indices[ir]=loopindex;
		if ( ir>1 ) { for ( core::Size jr=1, jrmax=ir-1; jr<=jrmax; ++jr ) { runtime_assert_string_msg( loop_indices[jr]!=loopindex, "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): The same residue cannot be specified twice!" ); } }
		//TR << "Current loop index is " << loopindex << std::endl; TR.flush(); //DELETE ME
	}
	std::string const errmsg( "In GeneralizedKICperturber::apply_perturb_backbone_by_bins(): The selected residues are not regularly connected through mainchain polymer bonds!" );
	for ( core::Size ir=1, irmax=rescount; ir<=irmax; ++ir ) {
		if ( ir==1 ) {
			if ( irmax>1 ) runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir+1] ), errmsg + "  Problem is with start residue." );
		} else if ( ir==irmax ) {
			if ( irmax>1 ) runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir-1] ), errmsg + "  Problem is with end residue." );
		} else {
			runtime_assert_string_msg( loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir-1] ) && loop_pose.residue(loop_indices[ir]).is_polymer_bonded( loop_indices[ir+1] ) , errmsg + "  Problem is with a middle residue." );
		}
	}

	for ( core::Size irepeat=1; irepeat<=n_iterations; ++irepeat ) { //Loop through the iterations
		utility::vector1 < core::Real > mainchain_torsions; //Vector of mainchain torsion values for a single residue that will be populated by the BinTransitionCalculator object.
		core::Size res_index( loop_indices[ numeric::random::rg().random_range( 1, loop_indices.size() ) ] ); //Index (in the loop) of the residue that will be set my the mover this iteration.  This is randomly drawn from the range [startres, endres] inclusive.

		//Set the torsions for the loop pose copy to the current torsion values:
		protocols::generalized_kinematic_closure::set_loop_pose ( loop_pose_copy, atomlist, torsions ); //Stripped-down version of the torsion-setting function.  Should be reasonably fast.

		//Actually generate the mainchain torsion values:
		bin_transition_calculator_->random_mainchain_torsions_using_adjacent_bins( loop_pose_copy.conformation(), res_index, must_switch_bins(), mainchain_torsions) ; //pose.conformation() and res_index are const inputs; mainchain_torsions is the output

		//Generate list of atoms defining the mainchain torsions to set:
		utility::vector1 < utility::vector1 <core::id::AtomID> > dihedral_list;
		for ( core::Size j=1, jmax=mainchain_torsions.size(); j<=jmax; ++j ) {
			utility::vector1 < core::id::AtomID > cur_atomid_list;
			cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(res_index).mainchain_atoms()[j], res_index ) );
			if ( j<jmax ) {
				cur_atomid_list.push_back( core::id::AtomID( loop_pose.conformation().residue(res_index).mainchain_atoms()[j+1], res_index ) );
			} else { //For the last torsion, we have to look up what the last mainchain atom is bonded to.
				core::conformation::Residue const &cur_res( loop_pose.conformation().residue(res_index) );
				core::Size const other_res( cur_res.connect_map( cur_res.type().upper_connect_id() ).resid() ); //Index of the other residue that this one is connected to.
				core::Size const other_conn_id( cur_res.connect_map( cur_res.type().upper_connect_id() ).connid() ); //Index of the connection on the other residue.
				cur_atomid_list.push_back(
					core::id::AtomID (
					loop_pose.conformation().residue(other_res).residue_connect_atom_index( other_conn_id  ), //Other atom that this one is connected to at its upper connection
					other_res //Other residue that this one is connected to at its upper connection.
					)
				);
			}
			dihedral_list.push_back( cur_atomid_list );
		}

		apply_set_dihedral ( dihedral_list, atomlist, mainchain_torsions, torsions, 0); //Recycle this function.
	} //Iterate until the number of iterations is reached

	return;
}

/// @brief Applies a sample_cis_peptide_bond perturbation to the list of torsions.
/// @details This checks whether each residue specified is an alpha- or beta-amino acid.  If it is, it samples the cis version of the omega angle (if omega is in the chain of atoms).
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residues - A vector of the indices of residues affected by this perturber.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in] tail_residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_sample_cis_peptide_bond(
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1 < std::pair < core::Size, core::Size > > const &/*tail_residue_map*/, //Mapping of (tail residue in loop_pose, tail residue in original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	using namespace protocols::generalized_kinematic_closure;

	core::Size const nres=residues.size();
	runtime_assert_string_msg( nres>0, "The sample_cis_peptide_bond perturber requires at least one residue to be specified.  Error in GeneralizdKICperturber::apply_sample_cis_peptide_bond." );

	//Cis probability:
	core::Real cis_prob = 0.1;
	if ( inputvalues_real_.size()>0 ) {
		cis_prob=inputvalues_real_[1];
		if ( inputvalues_real_.size()>1 ) {
			if ( TR.Warning.visible() ) TR.Warning << "Warning!  Multiple input values were passed to a sample_cis_peptide_bond perturber.  The first value is being used as the probability of a cis peptide bond." << std::endl;
		}
	}

	for ( core::Size ir=1; ir<=nres; ++ir ) { //Loop through all specified residues.
		core::Size const curres = get_loop_index(residues[ir], residue_map);
		if ( loop_pose.residue(curres).type().is_alpha_aa() || loop_pose.residue(curres).type().is_beta_aa() ) { //If this is either an alpha- or a beta-amino acid.
			core::Size omegaindex=0;
			for ( core::Size ia=4, iamax=atomlist.size()-3; ia<=iamax; ++ia ) { //Loop through the atom list and find the appropriate omega value
				if ( atomlist[ia].first.rsd()!=curres ) continue; //First find an atom with the current residue number.
				//TR << "Found residue." << std::endl; TR.flush(); //DELETE ME
				if ( loop_pose.residue(atomlist[ia].first.rsd()).atom_index("C") == atomlist[ia].first.atomno() ) { //If we've found the carbonyl carbon
					//TR << "Found carbonyl carbon." << std::endl; TR.flush(); //DELETE ME
					if ( loop_pose.residue(atomlist[ia+1].first.rsd()).atom_index("N") == atomlist[ia+1].first.atomno() ) { //Omega going forward
						omegaindex=ia;
						break;
					} else if ( loop_pose.residue(atomlist[ia-1].first.rsd()).atom_index("N") == atomlist[ia-1].first.atomno() ) { //Omega going backward
						omegaindex=ia-1;
						break;
					}
				}
			}

			if ( omegaindex!=0 ) { //If omega was found.
				if ( numeric::random::rg().uniform() < cis_prob ) { //Die-roll to decide whether to set this to a cis peptide bond
					torsions[omegaindex]=0.0; //Set to a cis peptide bond.
				}
			} else { //Else if omega wasn't found:
				if ( TR.Warning.visible() ) TR.Warning<<"Warning!  No omega angle was found for residue " << residues[ir] << " in the chain of atoms to close by GeneralizedKIC." << std::endl;
			}

		} else { //If this is neither an alpha- nor a beta-amino acid.
			if ( TR.Warning.visible() ) TR.Warning << "Warning!  Residue " << residues[ir] << " was passed to a sample_cis_peptide_bond perturber, but this residue is neither an alpha- nor a beta-amino acid.  Skipping." << std::endl;
		}
	}

	if ( TR.visible() ) TR.flush();
	if ( TR.Warning.visible() ) TR.Warning.flush();
	if ( TR.Debug.visible() ) TR.Debug.flush();

	return;
}

} //namespace perturber
} //namespace generalized_kinematic_closure
} //namespace protocols
