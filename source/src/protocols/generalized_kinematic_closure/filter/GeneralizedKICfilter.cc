// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.cc
/// @brief  Helper class defining success filters for generalized closure of arbitrary segments that could go through side-chains (e.g. disulfides).
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// BOINC includes -- keep these first:
#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC

// Unit Headers
#include <protocols/generalized_kinematic_closure/filter/GeneralizedKICfilter.hh>
#include <protocols/generalized_kinematic_closure/util.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/ScoringManager.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/id/AtomID.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>
#include <numeric/kinematic_closure/kinematic_closure_helpers.hh>
#include <numeric/conversions.hh>

#include <core/scoring/bin_transitions/BinTransitionCalculator.hh>
#include <core/scoring/bin_transitions/BinTransitionData.hh>

#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace generalized_kinematic_closure {
namespace filter {

static THREAD_LOCAL basic::Tracer TR( "protocols.generalized_kinematic_closure.filter.GeneralizedKICfilter" );

/// @brief Constructor for GeneralizedKICfilter.
///
GeneralizedKICfilter::GeneralizedKICfilter():
	filtertype_(no_filter),
	filter_params_real_(),
	filter_params_size_(),
	filter_params_bool_(),
	filter_params_string_(),
	bin_transition_calculator_(),
	bin_(""),
	resnum_(0),
	rama_threshold_(0.3),
	attach_boinc_ghost_observer_(false)
	//TODO -- make sure above data are copied properly when duplicating this mover.
{}

/// @brief Copy constructor for GeneralizedKICfilter.
///
GeneralizedKICfilter::GeneralizedKICfilter( GeneralizedKICfilter const &src ):
	utility::pointer::ReferenceCount(),
	filtertype_(src.filtertype_),
	filter_params_real_( src.filter_params_real_ ),
	filter_params_size_( src.filter_params_size_ ),
	filter_params_bool_( src.filter_params_bool_ ),
	filter_params_string_( src.filter_params_string_ ),
	bin_transition_calculator_( ), //CLONE this, below
	bin_( src.bin_ ),
	resnum_( src.resnum_ ),
	rama_threshold_( src.rama_threshold_ ),
	attach_boinc_ghost_observer_(src.attach_boinc_ghost_observer_)
	//TODO -- make sure above data are copied properly when duplicating this mover.
{
	if ( src.bin_transition_calculator_ ) bin_transition_calculator_ = utility::pointer::dynamic_pointer_cast< core::scoring::bin_transitions::BinTransitionCalculator >(src.bin_transition_calculator_->clone());
}

/// @brief Destructor for GeneralizedKICfilter mover.
///
GeneralizedKICfilter::~GeneralizedKICfilter() {}

/// @brief Clone operator to create a pointer to a fresh GeneralizedKICfilter object that copies this one.
///
GeneralizedKICfilterOP GeneralizedKICfilter::clone() const {
	return GeneralizedKICfilterOP( new GeneralizedKICfilter( *this ) );
}

/// @brief Returns the name of this class ("GeneralizedKICfilter").
std::string GeneralizedKICfilter::get_name() const{
	return "GeneralizedKICfilter";
}


/// @brief Given a filter type, return its name.  Returns "unknown_filter" if not recognized.
std::string GeneralizedKICfilter::get_filter_type_name( core::Size const filter_type ) const {
	std::string returnstring = "";
	switch(filter_type) {
	case no_filter :
		returnstring = "no_filter";
		break;
	case loop_bump_check :
		returnstring = "loop_bump_check";
		break;
	case atom_pair_distance :
		returnstring = "atom_pair_distance";
		break;
	case backbone_bin :
		returnstring = "backbone_bin";
		break;
	case alpha_aa_rama_check :
		returnstring = "alpha_aa_rama_check";
		break;
	case rama_prepro_check :
		returnstring = "rama_prepro_check";
		break;
	default :
		returnstring = "unknown_filter";
		break;
	}
	return returnstring;
}


/// @brief Given the name of a filter type, return the filter type enum.  Returns unknown_filter if not recognized.
filter_type GeneralizedKICfilter::get_filter_type_by_name( std::string const &filtername ) const {
	for ( core::Size i=1, imax=end_of_filter_list; i<imax; ++i ) {
		if ( get_filter_type_name(i)==filtername ) return (filter_type)i;
	}
	return unknown_filter;
}


/// @brief Sets the filter type for this filter.
void GeneralizedKICfilter::set_filter_type( filter_type const &ftype) {
	runtime_assert_string_msg(ftype > 0 && ftype < end_of_filter_list, "Filter type not recognized.  Error in GeneralizedKICfilter::set_filter_type().");
	filtertype_ = ftype;
	filter_params_real_.clear();
	filter_params_size_.clear();
	filter_params_bool_.clear();
	return;
}


/// @brief Sets the filter type for this filter by name.
void GeneralizedKICfilter::set_filter_type( std::string const &ftypename) {
	filter_type ftype = get_filter_type_by_name(ftypename);
	runtime_assert_string_msg( ftype < end_of_filter_list, "Filter type " + ftypename + " not recognized.  Error in GeneralizedKICfilter::set_filter_type()." );
	filtertype_ = ftype;
	filter_params_real_.clear();
	filter_params_size_.clear();
	filter_params_bool_.clear();
	return;
}


/// @brief Gets the filter type name for THIS filter.
std::string GeneralizedKICfilter::get_this_filter_type_name () const {
	return get_filter_type_name( filtertype_ );
}


/// @brief Add a real-valued filter parameter.
void GeneralizedKICfilter::add_filter_param( std::string const &param_name, core::Real const &value )
{
	filter_params_real_.push_back( std::pair< std::string, core::Real >( param_name, value ) );
	return;
}


/// @brief Add a integer-valued filter parameter.
void GeneralizedKICfilter::add_filter_param( std::string const &param_name, core::Size const value )
{
	filter_params_size_.push_back( std::pair< std::string, core::Size >( param_name, value ) );
	return;
}


/// @brief Add a Boolean-valued filter parameter.
void GeneralizedKICfilter::add_filter_param( std::string const &param_name, bool const value )
{
	filter_params_bool_.push_back( std::pair< std::string, bool >( param_name, value ) );
	return;
}


/// @brief Add a string-valued filter parameter.
void GeneralizedKICfilter::add_filter_param( std::string const &param_name, std::string const &value )
{
	filter_params_string_.push_back( std::pair< std::string, std::string >( param_name, value ) );
	return;
}

/// @brief Get a real-valued filter parameter.
/// @details Returns false if the parameter couldn't be found.
bool GeneralizedKICfilter::get_filter_param( std::string const &param_name, core::Real &outvalue ) const
{
	outvalue = 0.0; // Dummy value to suppress uninitialize warning (-Werror=maybe-uninitialized)

	core::Size const listsize = filter_params_real_.size();
	if ( !listsize ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( filter_params_real_[i].first==param_name ) {
			outvalue=filter_params_real_[i].second;
			return true;
		}
	}

	return false;
}

/// @brief Get a integer-valued filter parameter.
/// @details Returns false if the parameter couldn't be found.
bool GeneralizedKICfilter::get_filter_param( std::string const &param_name, core::Size &outvalue ) const
{
	outvalue = 0; // Dummy value to suppress uninitialize warning (-Werror=maybe-uninitialized)

	core::Size const listsize = filter_params_size_.size();
	if ( !listsize ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( filter_params_size_[i].first==param_name ) {
			outvalue=filter_params_size_[i].second;
			return true;
		}
	}

	return false;
}

/// @brief Get a Boolean-valued filter parameter.
/// @details Returns false if the parameter couldn't be found.
bool GeneralizedKICfilter::get_filter_param( std::string const &param_name, bool &outvalue ) const
{
	outvalue = false; // Dummy value to suppress uninitialize warning (-Werror=maybe-uninitialized)

	core::Size const listsize = filter_params_bool_.size();
	if ( !listsize ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( filter_params_bool_[i].first==param_name ) {
			outvalue=filter_params_bool_[i].second;
			return true;
		}
	}

	return false;
}

/// @brief Get a string-valued filter parameter.
/// @details Returns false if the parameter couldn't be found.
bool GeneralizedKICfilter::get_filter_param( std::string const &param_name, std::string &outvalue ) const
{
	core::Size const listsize = filter_params_string_.size();
	if ( !listsize ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( filter_params_string_[i].first==param_name ) {
			outvalue=filter_params_string_[i].second;
			return true;
		}
	}

	return false;
}

/// @brief Initializes the BinTransitionCalculator object and loads a bin_params file.
///
void GeneralizedKICfilter::load_bin_params(
	std::string const &bin_params_file
) {
	using namespace core::scoring::bin_transitions;

	//Create the object, if it doesn't exist.
	if ( !bin_transition_calculator_ ) {
		if ( TR.visible() ) TR << "Creating BinTransitionCalculator." << std::endl;
		bin_transition_calculator_=BinTransitionCalculatorOP( new BinTransitionCalculator );
	}

	if ( TR.visible() ) TR << "Loading bin_params file " << bin_params_file << "." << std::endl;
	bin_transition_calculator_->load_bin_params(bin_params_file);

	if ( TR.visible() ) TR.flush();

	return;
}

/// @brief Apply this filter to ONE of the kinematic closure solutions produced by the bridgeObjects function,
/// and return pass or fail.
/// @details
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool GeneralizedKICfilter::apply (
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 < core::Real > const &torsions,
	utility::vector1 < core::Real > const &bondangles,
	utility::vector1 < core::Real > const &bondlengths
) const {

	utility::vector1 < utility::vector1 < core::id::AtomID > > AtomIDs_loopindexed; //The list of lists of AtomIDs, with indices based on loop_pose.  If necessary, will be generated.

	switch(filtertype_) {
	case loop_bump_check :
		return apply_loop_bump_check( original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles, bondlengths);
	case atom_pair_distance :
		return apply_atom_pair_distance( original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles, bondlengths);
	case backbone_bin :
		return apply_backbone_bin( original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles, bondlengths);
	case alpha_aa_rama_check :
		return apply_alpha_aa_rama_check( original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles, bondlengths );
	case rama_prepro_check :
		return apply_rama_prepro_check( original_pose, loop_pose, residue_map, tail_residue_map, atomlist, torsions, bondangles, bondlengths );
	default :
		break;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Given an index in the original pose and a mapping from loop to pose,
/// return the index in the loop.
core::Size GeneralizedKICfilter::get_loop_index (
	core::Size const original_pose_index,
	utility::vector1 < std::pair < core::Size, core::Size > > const &residue_map
) const {
	for ( core::Size i=1, imax=residue_map.size(); i<=imax; ++i ) {
		if ( residue_map[i].second == original_pose_index ) return residue_map[i].first;
	}

	utility_exit_with_message("Residue does not exist in loop.  Exiting from GeneralizedKICperturber::get_loop_index with error status.");

	return 0;
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE APPLY FUNCTIONS FOR EACH FILTER                           //
////////////////////////////////////////////////////////////////////////////////

/// @brief Applies the loop_bump_check filter, which checks for clashes between the atoms in the chain
/// to be closed and the rest of the structure (or for clashes within these atoms).
/// @details Returns "true" for pass and "false" for fail.  Does NOT check for clashes with tail residues.
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool GeneralizedKICfilter::apply_loop_bump_check(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 < core::Real > const &torsions,
	utility::vector1 < core::Real > const &bondangles,
	utility::vector1 < core::Real > const &bondlengths
) const {
	using namespace protocols::generalized_kinematic_closure;
	using namespace core::id;

	core::Real const multiplier = 0.6; //TEMPORARY -- multiplier to reduce the stringency of the bump check.

	core::Size const nres=original_pose.size(); //Total number of residues in the original pose.

	core::pose::Pose pose(loop_pose); //Make a copy of the loop pose
	set_loop_pose (pose, atomlist, torsions, bondangles, bondlengths);

	//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
	if ( attach_boinc_ghost_observer() ) {
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( pose );
		//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
		//std::cerr.flush();
	}
#endif


	//Make a copy of the list of atoms in the chain to be closed:
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > atomlist_prime = atomlist;
	//Remove the first three and last three entries in atomlist_prime.  These are atoms outside of the loop used to establish frame:
	for ( core::Size i=1; i<=3; ++i ) {
		atomlist_prime.erase(atomlist_prime.end()-1);
		atomlist_prime.erase(atomlist_prime.begin());
	}
	//Loop through atomlist_prime and find alpha- or beta-amino acids with carbonyl carbons in the list.  Add the carbonyl oxygens.
	//Also find alpha- or beta-amino acids with CA carbons and one of (N, CM/C, CB).  Add the missing one of (N, CM/C, CB) to the list.
	for ( core::Size ia=1, iamax=atomlist_prime.size(); ia<=iamax; ++ia ) {
		core::Size const ia_res = atomlist_prime[ia].first.rsd();
		core::Size const ia_atomno = atomlist_prime[ia].first.atomno();
		if ( pose.residue(ia_res).type().is_alpha_aa() || pose.residue(ia_res).type().is_beta_aa() ) { //if alpha or beta amino acid
			if ( pose.residue(ia_res).atom_name(ia_atomno)=="C" && pose.residue(ia_res).has("O") ) {
				atomlist_prime.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> > (AtomID( pose.residue(ia_res).atom_index("O"), ia_res), loop_pose.residue(ia_res).xyz("O")  ) );
			} else if ( pose.residue(ia_res).atom_name(ia_atomno)=="CA" ) { //Adding whatever's missing from the (N, CM/C, CB) atoms surrounding a CA.
				utility::vector1 <bool> atsfound;
				atsfound.resize(3,false);
				utility::vector1 <std::string> ats;
				ats.resize(3, "");
				ats[1]="N"; ats[2]=(pose.residue(ia_res).type().is_beta_aa()?"CM":"C"); ats[3]="C";
				for ( core::Size i=1; i<=3; ++i ) {
					if ( ia>1 && atomlist_prime[ia-1].first.rsd()==ia_res && pose.residue(ia_res).atom_name(atomlist_prime[ia-1].first.atomno()) == ats[i] ) atsfound[i]=true;
					if ( ia<iamax && atomlist_prime[ia+1].first.rsd()==ia_res && pose.residue(ia_res).atom_name(atomlist_prime[ia+1].first.atomno()) == ats[i] ) atsfound[i]=true;
				}
				core::Size numfound = 0;
				for ( core::Size i=1; i<=3; ++i ) { if ( atsfound[i] ) ++numfound; } //Confirm that two of the three were found
				if ( numfound==2 ) {
					for ( core::Size i=1; i<=3; ++i ) {
						if ( !atsfound[i] && pose.residue(ia_res).has(ats[i]) ) {
							atomlist_prime.push_back( std::pair<AtomID, numeric::xyzVector<core::Real> > (AtomID( pose.residue(ia_res).atom_index(ats[i]), ia_res ), loop_pose.residue(ia_res).xyz(ats[i])) ); //Add the missing one of the three to the atom list
							//TR.Debug << "Adding the " << ats[i] << " atom for residue number " << ia_res << " in the chain of residues to be closed to the bump check." << std::endl; //DELETE ME
							break;
						}
					}
				}
			} //else if this is a CA atom
		} //if alpha or beta amino acid
		//TODO -- figure out what heavyatoms to add in the general case (not an alpha- or beta-amino acid).
	}


	//Loop through atoms in atomlist_prime list:
	for ( core::Size ia=1, iamax=atomlist_prime.size(); ia<=iamax; ++ia ) {
		core::Size const ia_res=atomlist_prime[ia].first.rsd();
		core::Size const ia_atomno=atomlist_prime[ia].first.atomno();
		core::Real const ia_radius = pose.residue(ia_res).type().atom_type(ia_atomno).lj_radius();

		//First, check internal clashes:
		for ( core::Size ja=4; ja<ia-1; ++ja ) {
			core::Size const ja_res=atomlist_prime[ja].first.rsd();
			if ( ja_res==ia_res || ja_res+1==ia_res || ja_res==ia_res+1 ) continue; //Don't bother checking for intra-residue clashes, or clashes between adjacent residues.
			core::Size const ja_atomno=atomlist_prime[ja].first.atomno();
			core::Real const ja_radius = pose.residue(ja_res).type().atom_type(ja_atomno).lj_radius();
			if ( pose.residue(ja_res).xyz(ja_atomno).distance_squared( pose.residue(ia_res).xyz(ia_atomno) )  < pow((ja_radius+ia_radius)*multiplier, 2) ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "GeneralizedKICfilter::apply_loop_bump_check filter failed due to internal clash within loop that was closed." << std::endl;
				}
				return false;
			}
		}

		//Next, check external clashes:
		for ( core::Size jr=1; jr<=nres; ++jr ) { //Loop through all residues of the main chain.
			if ( original_pose_residue_is_in_residue_map(jr, residue_map) ) continue; //If this is a loop residue, we're already considering clashes to it, so skip loop residues.
			if ( original_pose_residue_is_in_residue_map(jr, tail_residue_map) ) continue; //We won't consider clashes with tail residues.

			//If this residue is connected to the residue containing the current loop atom, skip it.
			core::Size const ia_original_res = get_original_pose_rsd(ia_res, residue_map);
			bool is_connected_to_loop_rsd = false;
			for ( core::Size i=1, imax=original_pose.residue(jr).n_possible_residue_connections(); i<=imax; ++i ) {
				if ( original_pose.residue(jr).connected_residue_at_resconn(i)==ia_original_res ) {
					is_connected_to_loop_rsd=true;
					break;
				}
			}
			if ( is_connected_to_loop_rsd ) continue;

			utility::vector1<core::Size> atoms_to_consider = original_pose.residue(jr).mainchain_atoms();
			if ( (original_pose.residue(jr).type().is_alpha_aa() || original_pose.residue(jr).type().is_beta_aa()) && original_pose.residue(jr).has("CB") ) {
				atoms_to_consider.push_back( original_pose.residue(jr).atom_index("CB") ); //Also consider clashes with the beta carbon, if present.
			}
			for ( core::Size ja=1, jamax=atoms_to_consider.size(); ja<=jamax; ++ja ) { //Loop through the mainchain atoms.
				core::Real const ja_radius = original_pose.residue(jr).type().atom_type( atoms_to_consider[ja] ).lj_radius();
				if ( original_pose.residue(jr).xyz( atoms_to_consider[ja] ).distance_squared( pose.residue(ia_res).xyz(ia_atomno) ) < pow((ja_radius+ia_radius)*multiplier, 2) ) {
					if ( TR.Debug.visible() ) {
						TR.Debug << "GeneralizedKICfilter::apply_loop_bump_check filter failed due to external clash between loop that was closed and the rest of the structure." << std::endl;
					}
					return false;
				}
			}
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "GeneralizedKICfilter::apply_loop_bump_check() filter passed." << std::endl;
		TR.Debug.flush();
	}
	return true;
}

/// @brief Applies the atom_pair_distance filter, checking that the distance between two atoms is less than
/// a given threshold (or greater than a given threshold if the user so specifies with the "greater_than"
/// option).
/// @details Returns "true" for pass and "false" for fail.  The user can set the following options:
/// "distance" (real-valued, mandatory)
/// "atom1" (string-valued, mandatory)
/// "atom2" (string-valued, mandatory)
/// "res1" (integer-valued, mandatory, based on original pose numbering)
/// "res2" (integer-valued, mandatory, based on original pose numbering)
/// "greater_than" (boolean, optional, false by default)
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool GeneralizedKICfilter::apply_atom_pair_distance(
	core::pose::Pose const & original_pose,
	core::pose::Pose const & loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const & residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const & tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const & atomlist,
	utility::vector1 < core::Real > const & torsions,
	utility::vector1 < core::Real > const & bondangles,
	utility::vector1 < core::Real > const & bondlengths
) const {

	//First, check that all necessary params have been set and store necessary information in local vars:
	std::string at1, at2;
	core::Size res1, res2;
	core::Real dist_cutoff, dist_cutoff_sq;
	bool greaterthan = false;
	runtime_assert_string_msg(get_filter_param("atom1", at1), "In GeneralizedKICfilter::apply_atom_pair_distance(): the atom_pair_distance filter cannot be used without specifying atom1.");
	runtime_assert_string_msg(get_filter_param("atom2", at2), "In GeneralizedKICfilter::apply_atom_pair_distance(): the atom_pair_distance filter cannot be used without specifying atom2.");
	runtime_assert_string_msg(get_filter_param("res1", res1), "In GeneralizedKICfilter::apply_atom_pair_distance(): the atom_pair_distance filter cannot be used without specifying res1.");
	runtime_assert_string_msg(get_filter_param("res2", res2), "In GeneralizedKICfilter::apply_atom_pair_distance(): the atom_pair_distance filter cannot be used without specifying res2.");
	runtime_assert_string_msg(get_filter_param("distance", dist_cutoff), "In GeneralizedKICfilter::apply_atom_pair_distance(): the atom_pair_distance filter cannot be used without specifying a distance cutoff (\"distance\" parameter).");
	get_filter_param("greater_than", greaterthan);
	dist_cutoff_sq = dist_cutoff*dist_cutoff;

	runtime_assert_string_msg(original_pose_residue_is_in_residue_map(res1, residue_map) || original_pose_residue_is_in_residue_map(res2, residue_map),
		"In GeneralizedKICfilter::apply_atom_pair_distance(): at least one of the residues passed to the atom_pair_distance filter must be in the loop to be closed." );

	runtime_assert_string_msg(original_pose.residue(res1).has(at1),
		"In GeneralizedKICfilter::apply_atom_pair_distance(): the residue specified with \"res1\" does not contain the atom specified with \"atom1\"." );
	runtime_assert_string_msg(original_pose.residue(res2).has(at2),
		"In GeneralizedKICfilter::apply_atom_pair_distance(): the residue specified with \"res2\" does not contain the atom specified with \"atom2\"." );


	core::pose::Pose looppose(loop_pose); //Make a copy of the loop pose
	core::pose::Pose fullpose(original_pose); //Make a copy of the full pose
	set_loop_pose (looppose, atomlist, torsions, bondangles, bondlengths); //Set the loop conformation using the torsions, bondangles, and bondlengths vectors.
	copy_loop_pose_to_original( fullpose, looppose, residue_map, tail_residue_map); //Copy the loop conformation to the full pose.

	//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
	if ( attach_boinc_ghost_observer() ) {
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( fullpose );
		protocols::boinc::Boinc::update_graphics_current_ghost( fullpose );
		//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
		//std::cerr.flush();
	}
#endif


	if ( fullpose.residue(res1).xyz(at1).distance_squared( fullpose.residue(res2).xyz(at2) ) > dist_cutoff_sq ) {
		if ( !greaterthan ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "GeneralizedKICfilter::apply_atom_pair_distance filter() failed." << std::endl;
				TR.Debug.flush();
			}
			return false;
		}
	} else {
		if ( greaterthan ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "GeneralizedKICfilter::apply_atom_pair_distance filter() failed." << std::endl;
				TR.Debug.flush();
			}
			return false;
		}
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "GeneralizedKICfilter::apply_atom_pair_distance filter() passed." << std::endl;
		TR.Debug.flush();
	}
	return true;
} //apply_atom_pair_distance

/// @brief Applies the backbone_bin filter, checking that a given residue lies within a defined
/// mainchain torsion bin and failing if it does not.
/// @details Returns "true" for pass and "false" for fail.  The user needs to have set a bin
/// transition probabilities file, a bin, and a residue.
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool GeneralizedKICfilter::apply_backbone_bin(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &,//tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 < core::Real > const &torsions,
	utility::vector1 < core::Real > const &bondangles,
	utility::vector1 < core::Real > const &bondlengths
) const {
	//Initial checks:
	if ( resnum() < 1 || resnum() > original_pose.size() ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_backbone_bin(): Could not apply filter.  The residue was not specified, or is out of range." );
	}
	if ( !original_pose_residue_is_in_residue_map(resnum(), residue_map) ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_backbone_bin(): The residue must be part of the loop being closed." );
	}
	if ( !bin_transition_calculator_ ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_backbone_bin(): No bin transition calculator object was set up." );
	}
	if ( bin_=="" ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_backbone_bin(): No bin was specified." );
	}

	core::Size const curres( get_loop_index( resnum(), residue_map ) );

	core::pose::Pose pose(loop_pose); //Make a copy of the loop pose
	set_loop_pose (pose, atomlist, torsions, bondangles, bondlengths); //Set the loop pose to the current solution

	//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
	if ( attach_boinc_ghost_observer() ) {
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( pose );
		//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
		//std::cerr.flush();
	}
#endif

	bool const inbin (bin_transition_calculator_->is_in_bin( pose.residue(curres), bin_ ));
	if ( TR.visible() ) {
		if ( inbin ) TR << "The backbone_bin filter reports that residue " << resnum() << " is in bin " << bin_ << ".  Passing." << std::endl;
		else  {
			TR << "The backbone_bin filter reports that residue " << resnum() << " is not in bin " << bin_ << ".  Failing and rejecting solution." << std::endl;
			for ( core::Size j=1, jmax=pose.residue(curres).mainchain_torsions().size(); j<=jmax; ++j ) {
				TR << "\tMainchain torsion " << j << ":\t" << pose.residue(curres).mainchain_torsion(j) << std::endl;
			}
		}
		TR.flush();
	}
	return inbin;
} //apply_backbone_bin

/// @brief Calculates Ramachandran energy for an alpha-amino acid based on its phi/psi values.
/// @details Returns "true" for pass (below threshold) and "false" for fail.
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool
GeneralizedKICfilter::apply_alpha_aa_rama_check (
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &,//tail_residue_map,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 < core::Real > const &torsions,
	utility::vector1 < core::Real > const &bondangles,
	utility::vector1 < core::Real > const &bondlengths
) const {
	//Initial checks:
	if ( resnum() < 1 || resnum() > original_pose.size() ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_alpha_aa_rama_check(): Could not apply filter.  The residue was not specified, or is out of range." );
	}
	if ( !original_pose_residue_is_in_residue_map(resnum(), residue_map) ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_alpha_aa_rama_check(): The residue must be part of the loop being closed." );
	}
	if ( !original_pose.residue(resnum()).type().is_alpha_aa() ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_alpha_aa_rama_check(): Cannot apply to a non-alpha amino acid!" );
	}

	core::Size const curres( get_loop_index( resnum(), residue_map ) );

	core::pose::Pose pose(loop_pose); //Make a copy of the loop pose
	set_loop_pose (pose, atomlist, torsions, bondangles, bondlengths); //Set the loop pose to the current solution

	//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
	if ( attach_boinc_ghost_observer() ) {
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( pose );
		//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
		//std::cerr.flush();
	}
#endif

	core::scoring::Ramachandran const & rama = core::scoring::ScoringManager::get_instance()->get_Ramachandran(); //Get the Rama scoring function
	core::Real rama_out(0.0), drama_dphi(0.0), drama_dpsi(0.0);
	rama.eval_rama_score_residue_nonstandard_connection(pose, pose.residue(curres), rama_out, drama_dphi, drama_dpsi);

	rama_out = rama_out * 0.25;  //Multiplied by talaris2014 scorefunction weights.  TODO: let the user set the scorefunction?

	bool const rama_passed( rama_out < rama_cutoff_energy() );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Res" << original_pose.residue(resnum()).name3() << resnum() << " Rama_score=" << rama_out << " Rama_cutoff=" << rama_cutoff_energy();
		if ( rama_passed ) {
			TR.Debug << " PASSED" << std::endl;
		} else {
			TR.Debug << " FAILED" << std::endl;
		}
	}
	return rama_passed;
} //apply_alpha_aa_rama_check

/// @brief Calculates RamaPrePro energy for a residue based on its mainchain torsion values.
/// @details Returns "true" for pass (below threshold) and "false" for fail.
/// @param[in] original_pose -- The full, initial pose.
/// @param[in] loop_pose -- A pose consisting of just the loop to be closed.
/// @param[in] residue_map -- The mapping of (residue index in loop_pose, residue index in original_pose).
/// @param[in] tail_residue_map -- The mapping of (tail residue index in loop_pose, tail residue index in original_pose).
/// @param[in] atomlist -- A list of atoms making the chain that was closed by bridgeObjects, with residue indices corresponding to loop_pose.
/// @param[in] torsions -- A vector of dihedral angles that the bridgeObjects function spat out.
/// @param[in] bondangles -- A vector of bond angles that the bridgeObjects function spat out.
/// @param[in] bondlengths -- A vector of bond lengths that the bridgeObjects function spat out.
bool
GeneralizedKICfilter::apply_rama_prepro_check(
	core::pose::Pose const &original_pose,
	core::pose::Pose const &loop_pose,
	utility::vector1 < std::pair <core::Size, core::Size> > const &residue_map,
	utility::vector1 < std::pair <core::Size, core::Size> > const &/*tail_residue_map*/,
	utility::vector1 < std::pair <core::id::AtomID, numeric::xyzVector<core::Real> > > const &atomlist,
	utility::vector1 < core::Real > const &torsions,
	utility::vector1 < core::Real > const &bondangles,
	utility::vector1 < core::Real > const &bondlengths
) const {
	//Initial checks:
	if ( resnum() < 1 || resnum() > original_pose.size() ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_rama_prepro_check(): Could not apply filter.  The residue was not specified, or is out of range." );
	}
	if ( !original_pose_residue_is_in_residue_map(resnum(), residue_map) ) {
		utility_exit_with_message( "In GeneralizedKICfilter::apply_rama_prepro_check(): The residue must be part of the loop being closed." );
	}

	core::Size const curres( get_loop_index( resnum(), residue_map ) );

	core::pose::Pose pose(loop_pose); //Make a copy of the loop pose
	set_loop_pose (pose, atomlist, torsions, bondangles, bondlengths); //Set the loop pose to the current solution

	//If this is the BOINC graphics build, and we're using the ghost pose observer, attach the observer now:
#ifdef BOINC_GRAPHICS
	if ( attach_boinc_ghost_observer() ) {
		protocols::boinc::Boinc::attach_graphics_current_pose_ghost_observer( pose );
		protocols::boinc::Boinc::update_graphics_current_ghost( pose );
		//std::cerr << "GenKIC attached a BOINC ghost observer." << std::endl;
		//std::cerr.flush();
	}
#endif

	core::scoring::RamaPrePro const & rama = core::scoring::ScoringManager::get_instance()->get_RamaPrePro(); //Get the Rama scoring function

	core::Real rama_out(0.0);
	utility::vector1 < core::Real > gradient; //Dummy var needed for next function call.
	core::conformation::Residue const &this_residue( pose.residue(curres) );
	utility::vector1 < core::Real > mainchain_torsions( this_residue.mainchain_torsions().size() - 1 );
	for ( core::Size i=1, imax=mainchain_torsions.size(); i<=imax; ++i ) mainchain_torsions[i] = this_residue.mainchain_torsions()[i];
	core::Size const &that_residue_index( this_residue.residue_connection_partner( this_residue.upper_connect().index() ) );

	rama.eval_rpp_rama_score(
		pose.conformation(),
		this_residue.type().get_self_ptr(),
		pose.residue_type(that_residue_index).get_self_ptr(),
		mainchain_torsions,
		rama_out,
		gradient,
		false /*Don't return gradient*/
	);

	rama_out = rama_out * 0.45;  //Multiplied by beta_nov15 scorefunction weights.  TODO: let the user set the scorefunction?

	bool const rama_passed( rama_out < rama_cutoff_energy() );

	if ( TR.Debug.visible() ) {
		TR.Debug << "Res" << original_pose.residue(resnum()).name3() << resnum() << " RamaPrePro_score=" << rama_out << " RamaPrePro_cutoff=" << rama_cutoff_energy();
		if ( rama_passed ) {
			TR.Debug << " PASSED" << std::endl;
		} else {
			TR.Debug << " FAILED" << std::endl;
		}
	}
	return rama_passed;
} //apply_rama_prepro_check

void
GeneralizedKICfilter::define_valid_filter_name_enumeration( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaRestriction genkic_filter_name;
	genkic_filter_name.name( "genkic_filter_name" );
	genkic_filter_name.base_type( xs_string );
	genkic_filter_name.add_restriction( xsr_enumeration, "no_filter" );
	genkic_filter_name.add_restriction( xsr_enumeration, "loop_bump_check" );
	genkic_filter_name.add_restriction( xsr_enumeration, "atom_pair_distance" );
	genkic_filter_name.add_restriction( xsr_enumeration, "backbone_bin" );
	genkic_filter_name.add_restriction( xsr_enumeration, "alpha_aa_rama_check" );
	genkic_filter_name.add_restriction( xsr_enumeration, "rama_prepro_check" );
	xsd.add_top_level_element( genkic_filter_name );

}

} //namespace filter
} //namespace generalized_kinematic_closure
} //namespace protocols
