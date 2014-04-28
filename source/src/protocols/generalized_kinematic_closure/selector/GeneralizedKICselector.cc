// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.cc
/// @brief  Helper class for GeneralizedKIC defining how solutions are chosen that pass filters.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

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
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace generalized_kinematic_closure {
namespace selector {

static basic::Tracer TR("protocols.generalized_kinematic_closure.selector.GeneralizedKICselector");
static numeric::random::RandomGenerator RG(9980002);  // <- Magic number, do not change it!

///@brief Creator for GeneralizedKICselector.
GeneralizedKICselector::GeneralizedKICselector():
		selectortype_(no_selector),
		selector_sfxn_(NULL)
		//utility::pointer::ReferenceCount(),
		//TODO -- make sure above data are copied properly when duplicating this mover.
{}

///@brief Destructor for GeneralizedKICselector mover.
GeneralizedKICselector::~GeneralizedKICselector() {}

///@brief Returns the name of this class ("GeneralizedKICselector").
std::string GeneralizedKICselector::get_name() const{
	return "GeneralizedKICselector";
}

///
/// @brief Given a selector type, return its name.  Returns "unknown_selector" if not recognized.
std::string GeneralizedKICselector::get_selector_type_name( core::Size const selector_type ) const {
	std::string returnstring = "";
	switch(selector_type) {
		case no_selector:
			returnstring = "no_selector";
			break;
		case random_selector:
			returnstring = "random_selector";
			break;
		case lowest_energy_selector:
			returnstring = "lowest_energy_selector";
			break;
		default:
			returnstring = "unknown_selector";
			break;
	}
	return returnstring;
}

///
/// @brief Given the name of a selector type, return the selector type enum.  Returns unknown_selector if not recognized.
selector_type GeneralizedKICselector::get_selector_type_by_name( std::string const &selectorname ) const {
	for(core::Size i=1, imax=end_of_selector_list; i<imax; ++i) {
		if(get_selector_type_name(i)==selectorname) return (selector_type)i;
	}
	return unknown_selector;
}

///
/// @brief Sets the selector type for this selector.
void GeneralizedKICselector::set_selector_type( selector_type const &stype) {
	runtime_assert_string_msg(stype > 0 && stype < end_of_selector_list, "Selector type not recognized.  Error in GeneralizedKICselector::set_selector_type().");
	selectortype_ = stype;
	return;
}

///
/// @brief Sets the selector type for this selector by name.
void GeneralizedKICselector::set_selector_type( std::string const &stypename) {
	selector_type stype = get_selector_type_by_name(stypename);
	runtime_assert_string_msg( stype < end_of_selector_list, "Selector type " + stypename + " not recognized.  Error in GeneralizedKICselector::set_selector_type()." );
	selectortype_ = stype;
	return;
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////


/// @brief Applies a selector type to choose a solution and set a loop pose.
/// @details
/// @param[in,out] pose -- The loop to be closed.  This function puts it into its new, closed conformation.
/// @param[in] original_pose -- The original pose.  Can be used for reference by selectors.
/// @param[in] residue_map -- Mapping of (loop residue, original pose residue).
/// @param[in] atomlist -- The list of (AtomID, original XYZ coordinates of atoms) representing the chain that was closed.
/// @param[in] torsions -- Matrix of [closure attempt #][solution #][torsion #] with torsion values for each torsion angle in the chain.  A selector will pick one solution.
/// @param[in] bondangles -- Matrix of [closure attempt #][solution #][angle #] with bond angle values for each bond angle in the chain.  A selector will pick one solution.
/// @param[in] bondlengths -- Matrix of [closure attempt #][solution #][bondlength #] with bond length for each bond in the chain.  A selector will pick one solution.
/// @param[in] nsol_for_attempt -- List of the number of solutions for each attempt.
/// @param[in] total_solutions -- Total number of solutions found.
void GeneralizedKICselector::apply (
	core::pose::Pose &pose,
	core::pose::Pose const &original_pose, //The original pose
	utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map, //mapping of (loop residue, original pose residue)
	utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, //torsions for each atom 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, //bond angle for each atom
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, //bond length for each atom
	utility::vector1 <core::Size> const &nsol_for_attempt,
	core::Size const total_solutions
) const {

	TR << "Choosing GeneralizedKIC solution." << std::endl;
	if(total_solutions < 1) {
		TR.Warning << "Warning!  No solutions passed to GeneralizedKICselector::apply.  The loop pose could not be updated!  No solution chosen!" << std::endl;
	}

	//Indices that specify where the solution will be found:
	core::Size chosen_attempt_number=0;
	core::Size chosen_solution = 0;

	switch( selectortype_ ) {
	case random_selector:
		apply_random_selector( nsol_for_attempt, total_solutions, chosen_attempt_number, chosen_solution );
		break;
	case lowest_energy_selector:
		apply_lowest_energy_selector(
			nsol_for_attempt,
			total_solutions,
			chosen_attempt_number,
			chosen_solution,
			selector_sfxn_,
			residue_map,
			atomlist,
			torsions,
			bondangles,
			bondlengths,
			pose,
			original_pose
		);
		break;
	default:
		TR.Warning << "Warning!  No selector specified for GeneralizedKICselector::apply.  The loop pose could not be updated!  No solution chosen!" << std::endl;
		return;
		break;
	}

	if(chosen_attempt_number!=0 && chosen_solution!=0) {
		set_loop_pose ( pose, atomlist, torsions[chosen_attempt_number][chosen_solution], bondangles[chosen_attempt_number][chosen_solution], bondlengths[chosen_attempt_number][chosen_solution]);
	} else {
		utility_exit_with_message("Internal program error.  For some reason, GeneralizedKICselector::apply() did not do its job.  Something has gone wrong that shouldn't have.  Contact a developer or an exorcist.");
	}

	return;
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE APPLY FUNCTIONS FOR EACH FILTER                           //
////////////////////////////////////////////////////////////////////////////////

/// @brief Applies a random_selector selector.
/// @details This picks a solution randomly from the solutions that passed filters.
void GeneralizedKICselector::apply_random_selector(
	utility::vector1<core::Size> const &nsol_for_attempt,
	core::Size const total_solutions,
	core::Size &chosen_attempt_number,
	core::Size &chosen_solution
) const {
	core::Size solutionnumber = static_cast<core::Size>(RG.random_range(1, total_solutions)); //Pick a random solution

	core::Size accumulator=0;
	for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) {
		accumulator+=nsol_for_attempt[i];
		if(accumulator>=solutionnumber)	{ //If we've passed the random number we picked, the solution is in this interval
			accumulator -= nsol_for_attempt[i]; //Back up and look for it.
			chosen_attempt_number = i;
			for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) {
				++accumulator;
				if(accumulator==solutionnumber) {
					chosen_solution=j;
					break;
				}
			}
			break;
		}
	}
	
	assert(chosen_attempt_number!=0 && chosen_solution!=0); //We should have picked an attempt number and a solution within that attempt at this point.

	return;
}

/// @brief Applies a lowest_energy_selector selector.
/// @details This picks the lowest-energy solution, as scored with sfxn.  It's a good idea to use a modified
/// scorefunction for this (something that just has the backbone conformation and H-bonding terms, for
/// example, since side-chains will not be repacked by default prior to invoking this selector).
void GeneralizedKICselector::apply_lowest_energy_selector(
	utility::vector1<core::Size> const &nsol_for_attempt,
	core::Size const /*total_solutions*/,
	core::Size &chosen_attempt_number,
	core::Size &chosen_solution,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, 
	core::pose::Pose const &ref_loop_pose,
	core::pose::Pose const &ref_pose
) const {
	using namespace protocols::generalized_kinematic_closure;

	core::scoring::ScoreFunctionOP my_sfxn=sfxn;
	if(!my_sfxn) my_sfxn=core::scoring::getScoreFunction(); //Get the default scorefunction if one has not been supplied.
	
	core::Real lowest_energy = 0.0;
	core::Size lowest_energy_attempt = 0;
	core::Size lowest_energy_solution = 0;

	//Copies of the loop pose and the full pose:
	core::pose::Pose fullpose = ref_pose;
	core::pose::Pose looppose = ref_loop_pose;

	for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) { //Loop through all attempts
		if(nsol_for_attempt[i]==0) continue;
		for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) { //Loop through all solutions for this attempt
			set_loop_pose( looppose, atomlist, torsions[i][j], bondangles[i][j], bondlengths[i][j]);
			copy_loop_pose_to_original( fullpose, looppose, residue_map);
			(*my_sfxn)(fullpose);
			TR << "Scoring solution " << j << " from closure attempt " << i << ".  E = " << fullpose.energies().total_energy() << std::endl;
			if(lowest_energy_attempt==0 || fullpose.energies().total_energy() < lowest_energy) {
				lowest_energy=fullpose.energies().total_energy();
				lowest_energy_attempt=i;
				lowest_energy_solution=j;
			}
		}
	}

	TR << "Lowest energy found = " << lowest_energy << std::endl;

	chosen_attempt_number = lowest_energy_attempt;
	chosen_solution = lowest_energy_solution;

	TR.flush();

	return;
}

} //namespace selector
} //namespace generalized_kinematic_closure
} //namespace protocols
