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
#include <numeric/constants.hh>
#include <numeric/random/random.hh>
#include <numeric/model_quality/rms.hh>
#include <ObjexxFCL/FArray2D.hh>

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

static thread_local basic::Tracer TR( "protocols.generalized_kinematic_closure.selector.GeneralizedKICselector" );

///@brief Creator for GeneralizedKICselector.
GeneralizedKICselector::GeneralizedKICselector():
		selectortype_(no_selector),
		selector_sfxn_(/* NULL */),
		boltzmann_kbt_(1.0)
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
		case boltzmann_energy_selector:
			returnstring = "boltzmann_energy_selector";
			break;
		case lowest_rmsd_selector:
			returnstring = "lowest_rmsd_selector";
			break;
		case lowest_delta_torsion_selector:
			returnstring = "lowest_delta_torsion_selector";
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
/// @param[in] tail_residue_map -- Mapping of (tail residue index in pose, tail residue index in original_pose).
/// @param[in] atomlist -- The list of (AtomID, original XYZ coordinates of atoms) representing the chain that was closed.
/// @param[in] torsions -- Matrix of [closure attempt #][solution #][torsion #] with torsion values for each torsion angle in the chain.  A selector will pick one solution.
/// @param[in] bondangles -- Matrix of [closure attempt #][solution #][angle #] with bond angle values for each bond angle in the chain.  A selector will pick one solution.
/// @param[in] bondlengths -- Matrix of [closure attempt #][solution #][bondlength #] with bond length for each bond in the chain.  A selector will pick one solution.
/// @param[in] nsol_for_attempt -- List of the number of solutions for each attempt.
/// @param[in] total_solutions -- Total number of solutions found.
/// @param[in] pre_selectoin_mover -- Pointer to a mover applied to each solution before applying the selector.
/// @param[in] preselection_mover_exists -- Boolean that determines whether a mover has been specified.
void GeneralizedKICselector::apply (
	core::pose::Pose &pose,
	core::pose::Pose const &original_pose, //The original pose
	utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map, //mapping of (loop residue, original pose residue)
	utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map, //mapping of (tail residue index in pose, tail residue index in original_pose)
	utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist, //list of atoms (residue indices are based on the loop_pose)
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, //torsions for each atom 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, //bond angle for each atom
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, //bond length for each atom
	utility::vector1 <core::Size> const &nsol_for_attempt,
	core::Size const total_solutions,
	protocols::moves::MoverOP pre_selection_mover,
	bool const preselection_mover_exists
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
		apply_lowest_energy_selector( nsol_for_attempt, total_solutions, chosen_attempt_number,
			chosen_solution, selector_sfxn_, residue_map, tail_residue_map, atomlist, torsions,	bondangles,
			bondlengths, pose, original_pose, boltzmann_kbt_, pre_selection_mover, preselection_mover_exists, false
		);
		break;
	case boltzmann_energy_selector:
		apply_lowest_energy_selector( nsol_for_attempt, total_solutions, chosen_attempt_number,
			chosen_solution, selector_sfxn_, residue_map, tail_residue_map, atomlist, torsions,	bondangles,
			bondlengths, pose, original_pose, boltzmann_kbt_, pre_selection_mover, preselection_mover_exists, true
		); //Recycle this function from lowest_energy_selector to avoid code duplication
		break;
	case lowest_rmsd_selector:
		apply_lowest_rmsd_selector( nsol_for_attempt, total_solutions, chosen_attempt_number,
			chosen_solution, residue_map, atomlist, torsions,	bondangles, bondlengths, pose, preselection_mover_exists
		);
		break;
	case lowest_delta_torsion_selector:
		apply_lowest_delta_torsion_selector( nsol_for_attempt, total_solutions, chosen_attempt_number,
			chosen_solution, atomlist, torsions, pose, preselection_mover_exists 
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
	core::Size solutionnumber = static_cast<core::Size>(numeric::random::rg().random_range(1, total_solutions)); //Pick a random solution

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
/// example, since side-chains will not be repacked by default prior to invoking this selector).  If the
/// use_boltzmann parameter is set to true, the function randomly selects a solution weighted by exp(-E/kbt).
void GeneralizedKICselector::apply_lowest_energy_selector(
	utility::vector1<core::Size> const &nsol_for_attempt,
	core::Size const /*total_solutions*/,
	core::Size &chosen_attempt_number,
	core::Size &chosen_solution,
	core::scoring::ScoreFunctionOP sfxn,
	utility::vector1 <std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 <std::pair <core::Size, core::Size> > const &tail_residue_map,
	utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, 
	utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, 
	core::pose::Pose const &ref_loop_pose,
	core::pose::Pose const &ref_pose,
	core::Real const &boltzmann_kbt,
	protocols::moves::MoverOP pre_selection_mover,
	bool const preselection_mover_exists, 
	bool const use_boltzmann
) const {
	using namespace protocols::generalized_kinematic_closure;

	core::scoring::ScoreFunctionOP my_sfxn=sfxn;
	if(!my_sfxn) my_sfxn=core::scoring::get_score_function(); //Get the default scorefunction if one has not been supplied.
	
	//Variables for finding the lowest energy:
	core::Real lowest_energy = 0.0;
	core::Size lowest_energy_attempt = 0;
	core::Size lowest_energy_solution = 0;

	//Variables for randomly selecting weighted by exp(-E/kbt):
	utility::vector1 < core::Real > boltzmann_factors;
	core::Real boltzmann_factor_accumulator = 0.0;

	//Copies of the loop pose and the full pose:
	core::pose::Pose fullpose = ref_pose;
	core::pose::Pose looppose = ref_loop_pose;

	for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) { //Loop through all attempts
		if(nsol_for_attempt[i]==0) continue;
		for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) { //Loop through all solutions for this attempt
			fullpose = ref_pose;
			looppose = ref_loop_pose;
			set_loop_pose( looppose, atomlist, torsions[i][j], bondangles[i][j], bondlengths[i][j]);
			copy_loop_pose_to_original( fullpose, looppose, residue_map, tail_residue_map);
			if(preselection_mover_exists) {
				TR.Debug << "Applying preselection mover to solution " << j << " from closure attempt " << i << "." << std::endl;
				pre_selection_mover->apply(fullpose);
			}
			(*my_sfxn)(fullpose);
			if(!use_boltzmann) { //If we're just finding the lowest-energy solution:
				TR.Debug << "Scoring solution " << j << " from closure attempt " << i << ".  E = " << fullpose.energies().total_energy() << std::endl;
				if(lowest_energy_attempt==0 || fullpose.energies().total_energy() < lowest_energy) {
					lowest_energy=fullpose.energies().total_energy();
					lowest_energy_attempt=i;
					lowest_energy_solution=j;
				}
			} else { //If we're randomly picking weighted by Boltzmann factors:
				boltzmann_factors.push_back( exp(-fullpose.energies().total_energy() / boltzmann_kbt) );
				boltzmann_factor_accumulator += boltzmann_factors[boltzmann_factors.size()]; //For normalization
				TR.Debug << "Scoring solution " << j << " from closure attempt " << i << ".  E = " << fullpose.energies().total_energy() << "  exp(-E/kbt) = " << boltzmann_factors[boltzmann_factors.size()] << std::endl;
			}
		}
	}

	if(!use_boltzmann) { //If we're picking the lowest-energy solution:
		TR.Debug << "Lowest energy found = " << lowest_energy << std::endl;
		chosen_attempt_number = lowest_energy_attempt;
		chosen_solution = lowest_energy_solution;
	} else { //If we're choosing randomly
		bool breaknow=false;
		while(!breaknow) {
			//Pick a random number from 0 to 1:
			core::Real const randnum = numeric::random::rg().uniform();
			core::Size counter=0;
			core::Real accumulator = 0.0;
			//Loop through all solutions, accumulate the normalized Boltzmann factors, and break when we reach the bin
			//that contains the random number chosen above.  (Since the bin widths correspond to the normalized Boltzmann
			//factors, the probability of picking a given bin is proportional to the Boltzmann factors).
			for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) {
				for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) {
					++counter;
					accumulator+=boltzmann_factors[counter]/boltzmann_factor_accumulator;
					if(accumulator > randnum) { //We've found the appropriate bin
						chosen_attempt_number=i;
						chosen_solution=j;
						breaknow=true;
						break;
					}
				}
				if(breaknow) break;
			}
			if(!breaknow) TR.Debug << "Boltzmann energy selector did not converge.  Trying again.  Final accumulator= " << accumulator << " randnum=" << randnum << std::endl; //DELETE
		}
	}

	TR.flush();
	TR.Debug.flush();
	TR.Warning.flush();

	return;
}

/// @brief Applies a lowest_rmsd_selector selector.
/// @details This picks the solution with the lowest RMSD from the starting pose.
void GeneralizedKICselector::apply_lowest_rmsd_selector( 
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const /*total_solutions*/,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution,
		utility::vector1 <std::pair <core::Size, core::Size> > const &/*residue_map*/,
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondangles, 
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &bondlengths, 
		core::pose::Pose const &ref_loop_pose,
		bool const preselection_mover_exists
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace ObjexxFCL;
	using namespace numeric::model_quality;

	if(preselection_mover_exists) TR.Warning << "Warning!  A preselection mover was specified, but this is incompatible with the lowest_rmsd_selector." << std::endl; 

	//Copies of the loop pose and the full pose:
	core::pose::Pose looppose = ref_loop_pose;

	//Vars for finding lowest RMSD:
	core::Real rmsd_current = 0.0;
	core::Real rmsd_lowest = 0.0;
	core::Size attempt_lowest = 0;
	core::Size solution_lowest = 0;

	//FArray2D of chain atoms in the reference pose (for RMSD calculation):
	FArray2D<core::Real> ref_chain;
	ref_chain.redimension( 3, atomlist.size()-6 );
	for(core::Size i=4, imax=atomlist.size()-3; i<=imax; ++i) {
		ref_chain(1,i-3)=atomlist[i].second[0];
		ref_chain(2,i-3)=atomlist[i].second[1];
		ref_chain(3,i-3)=atomlist[i].second[2];
	}

	//FArray2D of chain atoms for the current pose (for RMSD calculation):
	FArray2D<core::Real> cur_chain;
	cur_chain.redimension( 3, atomlist.size()-6 );
	for(core::Size i=4, imax=atomlist.size()-3; i<=imax; ++i) { //Initialize to zero, just to be safe.
		cur_chain(1,i-3)=0.0;
		cur_chain(2,i-3)=0.0;
		cur_chain(3,i-3)=0.0;
	}

	for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) { //Loop through all attempts
		for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) { //Loop through all solutions to this attempt
			//Set the loop pose to the current solution:
			set_loop_pose( looppose, atomlist, torsions[i][j], bondangles[i][j], bondlengths[i][j] );

			//Populate the cur_chain FArray2D
			for(core::Size ii=4, iimax=atomlist.size()-3; ii<=iimax; ++ii) {
				core::Size const rsd = atomlist[ii].first.rsd();
				core::Size const atomno = atomlist[ii].first.atomno();
				cur_chain(1,ii-3)=looppose.residue(rsd).xyz(atomno)[0] ;
				cur_chain(2,ii-3)=looppose.residue(rsd).xyz(atomno)[1] ;
				cur_chain(3,ii-3)=looppose.residue(rsd).xyz(atomno)[2] ;
			}

			//Calculate RMSD here for all atoms in the chain of atoms to be closed:
			rmsd_current = rms_wrapper(atomlist.size()-6, ref_chain, cur_chain);
			TR.Debug << "Attempt " << i << " solution " << j << " rmsd=" << rmsd_current << std::endl;

			if(attempt_lowest==0 || rmsd_current < rmsd_lowest) {
				rmsd_lowest=rmsd_current;
				attempt_lowest=i;
				solution_lowest=j;
			}
		}
	}

	chosen_attempt_number=attempt_lowest;
	chosen_solution=solution_lowest;

	TR.flush();
	TR.Debug.flush();
	TR.Warning.flush();

	return;
}

/// @brief Applies a lowest_delta_torsion_selector selector.
/// @details This picks the solution with the lowest delta torsion from the starting pose.
void GeneralizedKICselector::apply_lowest_delta_torsion_selector( 
		utility::vector1<core::Size> const &nsol_for_attempt,
		core::Size const /*total_solutions*/,
		core::Size &chosen_attempt_number,
		core::Size &chosen_solution,
		utility::vector1 <std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
		utility::vector1 <utility::vector1 <utility::vector1<core::Real> > > const &torsions, 
		core::pose::Pose const &ref_pose,
		bool const preselection_mover_exists
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace ObjexxFCL;
	using namespace numeric::model_quality;

	if(preselection_mover_exists) TR.Warning << "Warning!  A preselection mover was specified, but this is incompatible with the lowest_delta_torsion_selector." << std::endl; 

	//Vars for finding lowest RMSD:
	core::Real dtor_current = 0.0;
	core::Real dtor_lowest = 0.0;
	core::Size attempt_lowest = 0;
	core::Size solution_lowest = 0;

	for(core::Size i=1, imax=nsol_for_attempt.size(); i<=imax; ++i) { //Loop through all attempts
		for(core::Size j=1; j<=nsol_for_attempt[i]; ++j) { //Loop through all solutions to this attempt
			//Calculate delta torsion
			dtor_current = 0;
			//actually only 6 tors change, but we don't know the pivot at this point, do we?
			//be careful! don't perturb pivot atoms
			for(core::Size k=2, kmax=(torsions[i][j].size()-2); k<=kmax; ++k){ //loop torsion angles
				core::Real tor = ref_pose.conformation().torsion_angle(atomlist[k-1].first, atomlist[k].first, atomlist[k+1].first, atomlist[k+2].first);
				core::Real dtor = tor - numeric::conversions::radians(torsions[i][j][k]);
				if (dtor>numeric::constants::f::pi) dtor -= numeric::constants::f::pi_2;
				if (dtor<-numeric::constants::f::pi) dtor += numeric::constants::f::pi_2;
				//std::cout << k << "->" << numeric::conversions::degrees(dtor) << " " << torsions[i][j][k] << std::endl;
				dtor_current += dtor*dtor;
			}
			TR << "Attempt " << i << " solution " << j << " dtor=" << dtor_current << std::endl;

			if(attempt_lowest==0 || dtor_current < dtor_lowest) {
				dtor_lowest=dtor_current;
				attempt_lowest=i;
				solution_lowest=j;
			}
		}
	}

	chosen_attempt_number=attempt_lowest;
	chosen_solution=solution_lowest;

	TR.flush();

	return;
}

} //namespace selector
} //namespace generalized_kinematic_closure
} //namespace protocols
