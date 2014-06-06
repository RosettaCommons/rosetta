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

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace generalized_kinematic_closure {
namespace perturber {

static numeric::random::RandomGenerator RG(40031);  // <- Magic number, do not change it!
static basic::Tracer TR("protocols.generalized_kinematic_closure.perturber.GeneralizedKICperturber");

///@brief Creator for GeneralizedKICperturber.
GeneralizedKICperturber::GeneralizedKICperturber():
		//utility::pointer::ReferenceCount(),
		effect_(no_effect)
		//TODO -- make sure above data are copied properly when duplicating this mover.
{}

///@brief Destructor for GeneralizedKICperturber mover.
GeneralizedKICperturber::~GeneralizedKICperturber() {}

///@brief Returns the name of this class ("GeneralizedKICperturber").
std::string GeneralizedKICperturber::get_name() const{
	return "GeneralizedKICperturber";
}

///
/// @brief Returns the enum type for the effect of a pertuber based on a perturber name.
///        Returns unknown_effect if can't find a match for the name.
perturber_effect GeneralizedKICperturber::get_perturber_effect_from_name( std::string const &name ) const
{
	for(core::Size i=1; i<end_of_effect_list; ++i) {
		if(get_perturber_effect_name( i ) == name) return (perturber_effect)i;
	}
	return unknown_effect;
}

///
/// @brief Returns the name of a perturber given the enum type.
///        Returns "unknown_effect" if no such effect exists.
std::string GeneralizedKICperturber::get_perturber_effect_name( core::Size &effect ) const
{
	std::string returnstring = "";

	switch(effect) {
		case no_effect:
			returnstring = "no_effect";
			break;
		case set_dihedral:
			returnstring = "set_dihedral";
			break;
		case set_bondangle:
			returnstring = "set_bondangle";
			break;
		case set_bondlength:
			returnstring = "set_bondlength";
			break;
		case randomize_dihedral:
			returnstring = "randomize_dihedral";
			break;
		case randomize_alpha_backbone_by_rama:
			returnstring = "randomize_alpha_backbone_by_rama";
			break;
		case perturb_dihedral:
			returnstring = "perturb_dihedral";
			break;
		case sample_cis_peptide_bond:
			returnstring = "sample_cis_peptide_bond";
			break;
		default:
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

///
/// @brief Sets the effect of this perturber using the perturber effect name.
///        Exits with an error message if the name is unknown.
void GeneralizedKICperturber::set_perturber_effect( std::string const &effectname )
{
	core::Size effect = get_perturber_effect_from_name(effectname);
	runtime_assert_string_msg(effect < end_of_effect_list, "Perturber effect type " + effectname + "not recognized.  Error in GeneralizedKICperturber::set_perturber_effect.");
	set_perturber_effect( (perturber_effect)effect );
	return;
}

/// @brief Applies the perturbation to the vectors of desired torsions, desired angles, and desired bond lengths.
///
/// @detailed
///
/// @param[in] original_pose - The original pose.
/// @param[in] loop_pose - A pose consisting of just the loop to be perturbed, plus one residue on each side establishing the frame.
/// @param[in] residue_map - Mapping of (loop residue, original pose residue).
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in,out] torsions - Desired torsions for each atom; can be set or altered by the apply() function.
/// @param[in,out] bondangles - Desired bond angles for each atom; can be set or altered by the apply() function.
/// @param[in,out] bondlengths - Desired bond lengths for each atom; can be set or altered by the apply() function.
void GeneralizedKICperturber::apply (
	core::pose::Pose const &original_pose, //The original pose
	core::pose::Pose const &loop_pose, //A pose consisting of just the loop to be perturbed, plus one residue on each side establishing the frame
	utility::vector1< std::pair< core::Size, core::Size > > const &residue_map, //mapping of (loop residue, original pose residue)
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1< core::Real > &torsions, //desired torsions for each atom (input/output)
	utility::vector1< core::Real > &bondangles, //desired bond angle for each atom (input/output)
	utility::vector1< core::Real > &bondlengths //desired bond length for each atom (input/output)
) const {

	//utility::vector1 < core::Size > residues_loopindexed; //The list of residues, with indices based on loop_pose.  If necessary, will be generated.
	utility::vector1 < utility::vector1 < core::id::AtomID > > AtomIDs_loopindexed; //The list of lists of AtomIDs, with indices based on loop_pose.  If necessary, will be generated.

	switch(effect_) {
		case set_dihedral:
			reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
			apply_set_dihedral(AtomIDs_loopindexed, atomlist, torsions, 0);
			break;
		case set_bondangle:
			reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
			apply_set_bondangle(AtomIDs_loopindexed, atomlist, bondangles);
			break;
		case set_bondlength:
			reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
			apply_set_bondlength(AtomIDs_loopindexed, atomlist, bondlengths);
			break;
		case randomize_dihedral:
			reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
			apply_set_dihedral(AtomIDs_loopindexed, atomlist, torsions, 1); //We recycle the apply_set_dihedral() function to avoid code duplication
			break;
		case randomize_alpha_backbone_by_rama:
			apply_randomize_alpha_backbone_by_rama(original_pose, loop_pose, residues_, atomlist, residue_map, torsions);
			break;
		case perturb_dihedral:
			reindex_AtomIDs(residue_map, AtomIDs_loopindexed, original_pose);
			apply_set_dihedral(AtomIDs_loopindexed, atomlist, torsions, 2); //We recycle the apply_set_dihedral() function to avoid code duplication
			break;
		case sample_cis_peptide_bond:
			apply_sample_cis_peptide_bond(loop_pose, atomlist, residues_, residue_map, torsions);
			break;
		default:
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
	for(core::Size i=1, imax=residue_map.size(); i<=imax; ++i) {
		if(residue_map[i].second == original_pose_index) return residue_map[i].first;
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
	for(core::Size i=1, imax=atoms_.size(); i<=imax; ++i) {
		utility::vector1 < AtomID > newAtomIDset;
		for(core::Size j=1, jmax=atoms_[i].size(); j<=jmax; ++j) { //Looping through atoms in the set
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
/// @param[in,out] torsions - Desired torsions for each atom; set by this function.
/// @param[in] effect - Should the specified torsions be set (0), randomized (1), or perturbed (2)?
void GeneralizedKICperturber::apply_set_dihedral (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &dihedrallist, //List of sets of atoms defining dihedrals, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1< core::Real > &torsions, //desired torsions for each atom (input/output)
	core::Size const effect //0=set, 1=randomize, 2=perturb
) const {
	//TR << "Applying set_dihedral perturbation effect." << std::endl;

	runtime_assert_string_msg( dihedrallist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_dihedral, since no atoms were provided as input." );
	if(effect==0) runtime_assert_string_msg( inputvalues_real_.size() > 0, "Could not set dihedral value with GeneralizedKICperturber::apply_set_dihedral, since no value for the dihedral angle was provided as input." );
	if(effect==2) runtime_assert_string_msg( inputvalues_real_.size() > 0, "Could not perturb dihedrals with GeneralizedKICperturber::apply_set_dihedral, since no value for the dihedral angle perturbation magnitude provided as input." );


	bool separate_values = false; //Have separate dihedral values been provided for each dihedral in the list, or are we setting everything to one value?
	if(effect==0 || effect==2) {
		std::string str1=" set";
		std::string str2=" to";
		std::string str3="Setting ";
		std::string str4=" to ";
		if(effect==2) {
			str1=" perturbed"; str2=" by"; str3="Perturbing "; str4=" by ";
		}
		if(inputvalues_real_.size()>=dihedrallist.size()) {
			separate_values = true;
			if(inputvalues_real_.size()>dihedrallist.size()) TR.Warning << "Warning! Number of input values for set_dihedral pertruber effect exceeds the number of torsions to be" << str1 << ".  Using only the first " << dihedrallist.size() << " values." << std::endl;
		} else if (inputvalues_real_.size()!=1) {
			separate_values = false;
			TR.Warning << "Warning! Number of input values does not match the number of dihedral angles specified.  All angles will be" << str1 << str2 << " the first value." << std::endl;
		} else { //inputvalues_real_.size()==1 and dihedrallist.size() > 1
			TR << str3 << "all specified dihedral angles" << str4 << inputvalues_real_[1] << "." << std::endl;
		}
	}

	for( core::Size i=1, imax=dihedrallist.size(); i<=imax; ++i ) { //Loop through each dihedral angle given as input
		runtime_assert_string_msg( dihedrallist[i].size()==4 || dihedrallist[i].size()==2, "Error in GeneralizedKICperturber::apply_set_dihedral.  Either two or four AtomIDs must be provided to define a dihedral angle." );

		//If four atoms are specified, that uniquely defines the dihedral angle.
		//If only two are specified, we assume that the upstream and downstream atoms of those two are the other two atoms needed to define the dihedral.
		utility::vector1 < core::id::AtomID > dihed_atoms;
		if(dihedrallist[i].size()==4) dihed_atoms=dihedrallist[i];
		else if (dihedrallist[i].size()==2) {
			for(core::Size j=2, jmax=atomlist.size()-2; j<=jmax; ++j) { //Find the upstream and downstream atoms:
				if(dihedrallist[i][1]==atomlist[j].first && dihedrallist[i][2]==atomlist[j+1].first) {
					dihed_atoms.push_back( atomlist[j-1].first );
					dihed_atoms.push_back( atomlist[j].first );
					dihed_atoms.push_back( atomlist[j+1].first );
					dihed_atoms.push_back( atomlist[j+2].first );
					break;
				} else if(dihedrallist[i][2]==atomlist[j].first && dihedrallist[i][1]==atomlist[j+1].first) {
					dihed_atoms.push_back( atomlist[j+2].first );
					dihed_atoms.push_back( atomlist[j+1].first );
					dihed_atoms.push_back( atomlist[j].first );
					dihed_atoms.push_back( atomlist[j-1].first );
					break;
				}
			}
		}

		//Find this dihedral angle in the atom list:
		core::Size torsion_index = 0;
		for(core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if(dihed_atoms[2]==atomlist[j].first) {
				if( j<jmax && dihed_atoms[3]==atomlist[j+1].first && dihed_atoms[4]==atomlist[j+2].first && dihed_atoms[1]==atomlist[j-1].first ) { //dihedral going forward
					torsion_index=j;
				} else if ( j>3 && dihed_atoms[1]==atomlist[j+1].first && dihed_atoms[3]==atomlist[j-1].first && dihed_atoms[4]==atomlist[j-2].first ) { //dihedral going backward
					torsion_index=j-1;
				}
			}
			if(torsion_index > 0) break;
		}
		runtime_assert_string_msg(torsion_index > 0, "Error in GeneralizedKICperturber::apply_set_dihedral.  The dihedral angle specified was not found in the chain of atoms to be closed.");
		if (effect==2) {
			torsions[torsion_index] += RG.gaussian() * (separate_values ? inputvalues_real_[i] : inputvalues_real_[1]); //Add a randomly chosen value from a gaussian distribution of specified bredth.
		} else if(effect==1) { //randomizing torsions
			torsions[torsion_index] = RG.uniform()*360.0; //Set the desired torsion to the user-specified value.
		} else if (effect==0) { //setting torsions
			torsions[torsion_index] = (separate_values ? inputvalues_real_[i] : inputvalues_real_[1]); //Set the desired torsion to the user-specified value.
		}
	}
	
	return;
}

/// @brief Applies a set_bondangle perturbation to a list of bond angles.
///
/// @detailed
///
/// @param[in] bondanglelist - List of sets of atoms defining bond angles, indexed based on the loop_pose.
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in,out] bondangles - Desired bond angles for each atom; set by this function.
void GeneralizedKICperturber::apply_set_bondangle (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &bondanglelist, //List of sets of atoms defining bond angles, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1< core::Real > &bondangles //desired bond angles for each atom (input/output)
) const {
	//TR << "Applying set_bondangle perturbation effect." << std::endl;

	runtime_assert_string_msg( bondanglelist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondangle, since no atoms were provided as input." );
	runtime_assert_string_msg( inputvalues_real_.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondangle, since no value for the bond angle was provided as input." );

	bool separate_values = false; //Have separate bond angle values been provided for each bond angle in the list, or are we setting everything to one value?
	if(inputvalues_real_.size()>=bondanglelist.size()) {
		separate_values = true;
		if(inputvalues_real_.size()>bondanglelist.size()) TR.Warning << "Warning! Number of input values for set_bondangle pertruber effect exceeds the number of bond angles to be set.  Using only the first " << bondanglelist.size() << " values." << std::endl;
	} else if (inputvalues_real_.size()!=1) {
		separate_values = false;
		TR.Warning << "Warning! Number of input values does not match the number of bond angles specified.  All angles will be set to the first value." << std::endl;
	} else { //inputvalues_real_.size()==1 and bondanglelist.size() > 1
		TR << "Setting all specified bond angles to " << inputvalues_real_[1] << "." << std::endl;
	}

	for( core::Size i=1, imax=bondanglelist.size(); i<=imax; ++i ) { //Loop through each bond angle given as input
		runtime_assert_string_msg( bondanglelist[i].size()==3 || bondanglelist[i].size()==1, "Error in GeneralizedKICperturber::apply_set_bondangle.  Either one or three AtomIDs must be provided to define a bond angle." );

		utility::vector1 < core::id::AtomID > ang_atoms;
		if(bondanglelist[i].size()==3) ang_atoms = bondanglelist[i]; //If three atoms are specified for this bond angle, use them.
		else if (bondanglelist[i].size()==1) { //If only one atom is specified for this bond angle, we assume that the upstream and downstream atoms in the atom list are the other two atoms defining the angle.
			for(core::Size j=2, jmax=atomlist.size()-1; j<=jmax; ++j) { //Find the upstream and downstream atoms:
				if(bondanglelist[i][1]==atomlist[j].first) {
					ang_atoms.push_back( atomlist[j-1].first );
					ang_atoms.push_back( atomlist[j].first );
					ang_atoms.push_back( atomlist[j+1].first );
					break;
				}
			}
		}

		//Find this bond angle in the atom list:
		core::Size bondangle_index = 0;
		for(core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if(ang_atoms[2]==atomlist[j].first) {
				if( ang_atoms[1]==atomlist[j-1].first && ang_atoms[3]==atomlist[j+1].first ) { //bond angle going forward
					bondangle_index=j;
				} else if ( ang_atoms[1]==atomlist[j+1].first && ang_atoms[3]==atomlist[j-1].first ) { //bond angle going backward
					bondangle_index=j;
				}
			}
			if(bondangle_index > 0) break;
		}
		runtime_assert_string_msg(bondangle_index > 0, "Error in GeneralizedKICperturber::apply_set_bondangle.  The bond angle specified was not found in the chain of atoms to be closed.");
		bondangles[bondangle_index] = (separate_values ? inputvalues_real_[i] : inputvalues_real_[1]); //Set the desired bond angle to the user-specified value.
	}

	return;
}

/// @brief Applies a set_bondlength perturbation to a list of bond lengths.
///
/// @detailed
///
/// @param[in] bondlengthlist - List of sets of atoms defining bond lengths, indexed based on the loop_pose.
/// @param[in] atomlist - List of atoms (residue indices are based on the loop_pose).
/// @param[in,out] bondlengths - Desired bond lengths for each atom; set by this function.
void GeneralizedKICperturber::apply_set_bondlength (
	utility::vector1 < utility::vector1 < core::id::AtomID > > const &bondlengthlist, //List of sets of atoms defining bond lengths, indexed based on the loop_pose.
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1< core::Real > &bondlengths //desired bond lengths for each atom (input/output)
) const {
	//TR << "Applying set_bondlength perturbation effect." << std::endl;

	runtime_assert_string_msg( bondlengthlist.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondlength, since no atoms were provided as input." );
	runtime_assert_string_msg( inputvalues_real_.size() > 0, "Could not apply GeneralizedKICperturber::apply_set_bondlength, since no value for the bond length was provided as input." );

	bool separate_values = false; //Have separate bond length values been provided for each bond length in the list, or are we setting everything to one value?
	if(inputvalues_real_.size()>=bondlengthlist.size()) {
		separate_values = true;
		if(inputvalues_real_.size()>bondlengthlist.size()) TR.Warning << "Warning! Number of input values for set_bondlength pertruber effect exceeds the number of bond lengths to be set.  Using only the first " << bondlengthlist.size() << " values." << std::endl;
	} else if (inputvalues_real_.size()!=1) {
		separate_values = false;
		TR.Warning << "Warning! Number of input values does not match the number of bond lengths specified.  All lengths will be set to the first value." << std::endl;
	} else { //inputvalues_real_.size()==1 and bondlengthlist.size() > 1
		TR << "Setting all specified bond lengths to " << inputvalues_real_[1] << "." << std::endl;
	}

	for( core::Size i=1, imax=bondlengthlist.size(); i<=imax; ++i ) { //Loop through each bond length given as input
		runtime_assert_string_msg( bondlengthlist[i].size()==2, "Error in GeneralizedKICperturber::apply_set_bondlength.  Two AtomIDs must be provided to define a bond length." );

		//Find this bond length in the atom list:
		core::Size bondlength_index = 0;
		for(core::Size j=3, jmax=atomlist.size()-2; j<=jmax; ++j) { //Loop through all atoms in the chain of atoms to be closed (excluding the anchors).
			if(bondlengthlist[i][1]==atomlist[j].first) {
				if(j<jmax && bondlengthlist[i][2]==atomlist[j+1].first ) { //bond length going forward
					bondlength_index=j;
				} else if (j>3 && bondlengthlist[i][2]==atomlist[j-1].first ) { //bond length going backward
					bondlength_index=j-1;
				}
			}
			if(bondlength_index > 0) break;
		}
		runtime_assert_string_msg(bondlength_index > 0, "Error in GeneralizedKICperturber::apply_set_bondlength.  The bond length specified was not found in the chain of atoms to be closed.");
		bondlengths[bondlength_index] = (separate_values ? inputvalues_real_[i] : inputvalues_real_[1]); //Set the desired bond length to the user-specified value.
	}

	return;
}

/// @brief Applies a randomize_alpha_backbone_by_rama perturbation to the list of torsions.
/// @details This checks whether each residue is an alpha-amino acid.
/// @param[in] original_pose - The input pose.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] residues - A vector of the indices of residues affected by this perturber.  Note that 
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
		utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	using namespace protocols::generalized_kinematic_closure;

	//TR << "Applying randomize_alpha_backbone_by_rama perturbation effect." << std::endl;  TR.flush(); //DELETE ME

	runtime_assert_string_msg( residues.size() > 0 , "Residues must be specified for the randomize_alpha_backbone_by_rama generalized kinematic closure perturber." );

	core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran();

	core::Size nres = original_pose.n_residue();

	//Make a copy of the loop_pose
	core::pose::Pose loop_pose_copy = loop_pose;

	for(core::Size ir=1, irmax=residues.size(); ir<=irmax; ++ir) {
		runtime_assert_string_msg( residues[ir] <= nres, "Unable to apply randomize_alpha_backbone_by_rama perturbation.  Residue list includes residues that are not in the original pose." );

		//Check for alpha amino acids:
		if(!original_pose.residue(residues[ir]).type().is_alpha_aa())	{
			TR.Warning << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama, but this residue is not an alpha-amino acid.  Skipping." << std::endl;
			TR.Warning.flush();
			continue;
		}

		//Get this residue's index in the loop pose:
		core::Size const loopindex = get_loop_index(residues[ir], residue_map);
		//TR << "Current loop index is " << loopindex << std::endl; TR.flush(); //DELETE ME

		//Randomize phi and psi for this residue:
		core::Real rama_phi=0;
		core::Real rama_psi=0;
		rama.random_phipsi_from_rama(loop_pose_copy.aa(loopindex), rama_phi, rama_psi);
		general_set_phi(loop_pose_copy, loopindex, rama_phi);
		general_set_psi(loop_pose_copy, loopindex, rama_psi);

		//Finding and setting phi and psi:
		utility::vector1 < std::string > at1list;
		utility::vector1 < std::string > at2list;
		utility::vector1 < std::string > angname;
		at1list.push_back("N"); at1list.push_back("CA");
		at2list.push_back("CA"); at2list.push_back("C");
		angname.push_back("phi"); angname.push_back("psi");
		for(core::Size j=1; j<=2; ++j) {
			for(core::Size ia=4, iamax=atomlist.size()-3; ia<=iamax; ++ia) {
				if(atomlist[ia].first.rsd()!=loopindex) continue;
				if(atomlist[ia].first.atomno() == loop_pose_copy.residue(loopindex).atom_index(at1list[j])) {
					core::Real angleval=0.0;
					if(ia<iamax && atomlist[ia+1].first.rsd()==loopindex && atomlist[ia+1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going forward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia+2].first),
							angleval
						);
						torsions[ia]=angleval;
					} else if (ia>4 && atomlist[ia-1].first.rsd()==loopindex && atomlist[ia-1].first.atomno()==loop_pose_copy.residue(loopindex).atom_index(at2list[j]) ) { //Torsion going backward
						numeric::dihedral_degrees (
							loop_pose_copy.xyz(atomlist[ia+1].first),
							loop_pose_copy.xyz(atomlist[ia].first),
							loop_pose_copy.xyz(atomlist[ia-1].first),
							loop_pose_copy.xyz(atomlist[ia-2].first),
							angleval
						);
						torsions[ia-1]=angleval;
					} else {
						TR.Warning << "Warning! Residue " << residues[ir] << " was passed to GeneralizedKICperturber::apply_randomize_alpha_backbone_by_rama, but its " << angname[j] << " dihedral angle was not part of the chain of atoms to be closed." << std::endl;  TR.Warning.flush();
					}
					break;
				}
			}
		}

	}

	return;
}

/// @brief Applies a sample_cis_peptide_bond perturbation to the list of torsions.
/// @details This checks whether each residue specified is an alpha- or beta-amino acid.  If it is, it samples the cis version of the omega angle
/// (if omega is in the chain of atoms).  Note that this sets the omega value to 0 with some probability specified in the value passed to this
/// perturber.  It assumes that the omega values have all been initialized to 180.
/// @param[in] loop_pose - A pose that is just the loop to be closed (possibly with other things hanging off of it).
/// @param[in] atomlist - A vector of pairs of AtomID, xyz coordinate.  Residue indices are based on the loop pose, NOT the original pose.
/// @param[in] residues - A vector of the indices of residues affected by this perturber.  Note that 
/// @param[in] residue_map - A vector of pairs of (loop pose index, original pose index).
/// @param[in,out] torsions - A vector of desired torsions, some of which are randomized by this function.
void GeneralizedKICperturber::apply_sample_cis_peptide_bond(
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair < core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 <core::Size> const &residues,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map, //Mapping of (loop_pose, original_pose).
	utility::vector1< core::Real > &torsions //desired torsions for each atom (input/output)
) const {
	using namespace protocols::generalized_kinematic_closure;

	core::Size const nres=residues.size();
	runtime_assert_string_msg( nres>0, "The sample_cis_peptide_bond perturber requires at least one residue to be specified.  Error in GeneralizdKICperturber::apply_sample_cis_peptide_bond." );

	//Cis probability:
	core::Real cis_prob = 0.1;
	if(inputvalues_real_.size()>0) {
		cis_prob=inputvalues_real_[1];
		if(inputvalues_real_.size()>1) {
			TR.Warning << "Warning!  Multiple input values were passed to a sample_cis_peptide_bond perturber.  The first value is being used as the probability of a cis peptide bond." << std::endl;
		}
	}

	for(core::Size ir=1; ir<=nres; ++ir) { //Loop through all specified residues.
		core::Size const curres = get_loop_index(residues[ir], residue_map);
		if(loop_pose.residue(curres).type().is_alpha_aa() || loop_pose.residue(curres).type().is_beta_aa() ) { //If this is either an alpha- or a beta-amino acid.
			core::Size omegaindex=0;
			for(core::Size ia=4, iamax=atomlist.size()-3; ia<=iamax; ++ia) { //Loop through the atom list and find the appropriate omega value
				if(atomlist[ia].first.rsd()!=curres) continue; //First find an atom with the current residue number.
				//TR << "Found residue." << std::endl; TR.flush(); //DELETE ME
				if(loop_pose.residue(atomlist[ia].first.rsd()).atom_index("C") == atomlist[ia].first.atomno()) { //If we've found the carbonyl carbon
					//TR << "Found carbonyl carbon." << std::endl; TR.flush(); //DELETE ME
					if(loop_pose.residue(atomlist[ia+1].first.rsd()).atom_index("N") == atomlist[ia+1].first.atomno()) { //Omega going forward
						omegaindex=ia;
						break;
					} else if( loop_pose.residue(atomlist[ia-1].first.rsd()).atom_index("N") == atomlist[ia-1].first.atomno() ) { //Omega going backward
						omegaindex=ia-1;
						break;
					}
				}
			}

			if(omegaindex!=0) { //If omega was found.
				if(RG.uniform() < cis_prob) { //Die-roll to decide whether to set this to a cis peptide bond
					torsions[omegaindex]=0.0; //Set to a cis peptide bond.
				}
			} else { //Else if omega wasn't found:
				TR.Warning<<"Warning!  No omega angle was found for residue " << residues[ir] << " in the chain of atoms to close by GeneralizedKIC." << std::endl;
			}

		} else { //If this is neither an alpha- nor a beta-amino acid.
			TR.Warning << "Warning!  Residue " << residues[ir] << " was passed to a sample_cis_peptide_bond perturber, but this residue is neither an alpha- nor a beta-amino acid.  Skipping." << std::endl;
		}
	}

	TR.flush();
	TR.Warning.flush();

	return;
}

} //namespace perturber
} //namespace generalized_kinematic_closure
} //namespace protocols
