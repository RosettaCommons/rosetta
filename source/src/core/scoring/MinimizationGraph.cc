// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/MinimizationGraph.cc
/// @brief  Minimization graph class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/MinimizationGraph.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>

#include <core/kinematics/MinimizerMapBase.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>

// Numeric headers

// Boost Headers
#include <utility/graph/unordered_object_pool.hpp>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>

#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh> // AUTO IWYU For LongRangeTwoBodyEnergy
#include <boost/pool/object_pool.hpp> // AUTO IWYU For unordered_object_pool::construct


namespace core {
namespace scoring {

using namespace utility::graph;

///////// Minimization Node Class /////////////

MinimizationNode::MinimizationNode( Graph * owner, Size index ) :
	parent( owner, index )
{}

MinimizationNode::~MinimizationNode() = default;


void MinimizationNode::print() const
{
	std::cout << "MinimizationNode::print() deferring to parent::print()" << std::endl;
	parent::print();
}

/// @brief copy mmember data from source node
///
/// invoked by copy ctor and operator= methods from Graph base class
void MinimizationNode::copy_from( parent const * source )
{
	auto const & mn_source = static_cast< MinimizationNode const & > ( * source );

	res_min_data_ = mn_source.res_min_data_; // deep copy

	onebody_enmeths_ = mn_source.onebody_enmeths_;
	active_1benmeths_ = mn_source.active_1benmeths_;
	active_1benmeths_std_ = mn_source.active_1benmeths_std_;
	active_1benmeths_ext_ = mn_source.active_1benmeths_ext_;
	dof_deriv_1benmeths_ = mn_source.dof_deriv_1benmeths_;
	sfs_dm_req_1benmeths_ = mn_source.sfs_dm_req_1benmeths_;
	sfd_req_1benmeths_ = mn_source.sfd_req_1benmeths_;
	twobody_enmeths_ = mn_source.twobody_enmeths_;
	active_intrares2benmeths_ = mn_source.active_intrares2benmeths_;
	active_intrares2benmeths_std_ = mn_source.active_intrares2benmeths_std_;
	active_intrares2benmeths_ext_ = mn_source.active_intrares2benmeths_ext_;
	dof_deriv_2benmeths_ = mn_source.dof_deriv_2benmeths_;
	sfs_dm_req_2benmeths_ = mn_source.sfs_dm_req_2benmeths_;
	sfd_req_2benmeths_ = mn_source.sfd_req_2benmeths_;
	sfs_drs_req_enmeths_ = mn_source.sfs_drs_req_enmeths_;
	weight_ = mn_source.weight_;
	dweight_ = mn_source.dweight_;
}

Size MinimizationNode::count_static_memory() const
{
	return sizeof ( MinimizationNode );
}

Size MinimizationNode::count_dynamic_memory() const
{
	Size tot = 0;
	tot += parent::count_dynamic_memory(); //recurse to parent
	/// add in the cost of storing the minimization data here
	return tot;
}

bool MinimizationNode::add_onebody_enmeth( OneBodyEnergyCOP enmeth, Residue const & rsd, Pose const & pose, int domain_map_color )
{
	onebody_enmeths_.push_back( enmeth );
	return classify_onebody_enmeth( enmeth, rsd, pose, domain_map_color );
}

/// @details Note that the minimization graph is used both to evaluate the scores for
/// both the "modern" EnergyMethods and the "old" energy methods (those which still
/// define eval_atom_derivative instead of the pairwise-decomposable derivative
/// evaluation methods) so when these older energy methods are added to the minimization
/// graph, they should not simply be rejected.
bool MinimizationNode::add_twobody_enmeth(
	TwoBodyEnergyCOP enmeth,
	Residue const & rsd,
	Pose const & pose,
	EnergyMap const & weights,
	int domain_map_color
) {
	//if ( enmeth->minimize_in_whole_structure_context( pose ) ) return false;

	twobody_enmeths_.push_back( enmeth );
	return classify_twobody_enmeth( enmeth, rsd, pose, weights, domain_map_color );
}

/// @details basically the logic of classify_onebody_enmeth without
/// the decision of whether or not residue should be active
void
MinimizationNode::activate_dof_deriv_one_body_method(
	OneBodyEnergyCOP enmeth,
	pose::Pose const & pose
)
{
	debug_assert( enmeth->defines_dof_derivatives( pose ) );

	if ( std::find( dof_deriv_1benmeths_.begin(), dof_deriv_1benmeths_.end(), enmeth ) !=
			dof_deriv_1benmeths_.end() ) {
		// Do not add the method twice; if it is already in the dof_deriv_1benmeths_
		// list, then it is already in the other lists
		return;
	}

	_activate_1benmeth( enmeth );
	add_dof_deriv_1benmeth( enmeth );
	_check_sfs_and_sfd_for_1benmeth( enmeth, pose );
}

void MinimizationNode::activate_dof_deriv_two_body_method(
	TwoBodyEnergyCOP enmeth,
	pose::Pose const & pose
)
{
	debug_assert( enmeth->defines_intrares_dof_derivatives( pose ) );

	// Do not add the energy method twice. If it is already in the
	// dof_deriv_2beneths_ list, it is already in the other lists
	if ( std::find( dof_deriv_2benmeths_.begin(), dof_deriv_2benmeths_.end(), enmeth ) !=
			dof_deriv_2benmeths_.end() ) {
		return;
	}

	_activate_2benmeth( enmeth );
	add_dof_deriv_2benmeth( enmeth );
	_check_sfs_and_sfd_for_2benmeth( enmeth, pose );
}


void
MinimizationNode::setup_for_minimizing(
	Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	basic::datacache::BasicDataCache & res_data_cache
) {
	for ( auto iter = active_1benmeths_begin(),
			iter_end = active_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing_for_residue( rsd, pose, sfxn, min_map, res_data_cache, res_min_data_ );
	}
	for ( TwoBodyEnergiesIterator iter = twobody_enmeths_.begin(),
			iter_end = twobody_enmeths_.end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing_for_residue( rsd, pose, sfxn, min_map, res_data_cache, res_min_data_ );
	}
}

void MinimizationNode::setup_for_scoring(
	Residue const & rsd,
	basic::datacache::BasicDataCache & residue_data_cache,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) {
	/// 1a 1body energy methods
	for ( auto iter = sfs_dm_req_1benmeths_begin(),
			iter_end = sfs_dm_req_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue( rsd, pose, sfxn, res_min_data_ );
	}

	/// 1b 2body intraresidue contributions
	for ( auto iter = sfs_dm_req_2benmeths_begin(),
			iter_end = sfs_dm_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue( rsd, pose, sfxn, res_min_data_ );
	}

	// 1c EnergyMethods that require sfs-drs calls
	for ( auto iter = sfs_drs_req_enmeths_begin(),
			iter_end = sfs_drs_req_enmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue( rsd, pose, sfxn, residue_data_cache );
	}
}

void MinimizationNode::setup_for_derivatives(
	Residue const & rsd,
	basic::datacache::BasicDataCache & residue_data_cache,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) {
	/// 1a 1body energy methods
	for ( auto iter = sfd_req_1benmeths_begin(),
			iter_end = sfd_req_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives_for_residue( rsd, pose, sfxn, res_min_data_, residue_data_cache );
	}
	/// 1b 2body intraresidue contributions
	for ( auto iter = sfd_req_2benmeths_begin(),
			iter_end = sfd_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives_for_residue( rsd, pose, sfxn, res_min_data_, residue_data_cache );
	}
}

void MinimizationNode::update_active_enmeths_for_residue(
	Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap const & weights,
	int domain_map_color
) {
	active_1benmeths_.clear();
	active_1benmeths_std_.clear();
	active_1benmeths_ext_.clear();
	dof_deriv_1benmeths_.clear();
	sfs_dm_req_1benmeths_.clear();
	sfd_req_1benmeths_.clear();

	active_intrares2benmeths_.clear();
	active_intrares2benmeths_std_.clear();
	active_intrares2benmeths_ext_.clear();
	dof_deriv_2benmeths_.clear();
	sfs_dm_req_2benmeths_.clear();
	sfd_req_2benmeths_.clear();

	sfs_drs_req_enmeths_.clear();

	for ( OneBodyEnergiesIterator iter = onebody_enmeths_.begin(),
			iter_end = onebody_enmeths_.end(); iter != iter_end; ++iter ) {
		classify_onebody_enmeth( *iter, rsd, pose, domain_map_color );
	}
	for ( TwoBodyEnergiesIterator iter = twobody_enmeths_.begin(),
			iter_end = twobody_enmeths_.end(); iter != iter_end; ++iter ) {
		classify_twobody_enmeth( *iter, rsd, pose, weights, domain_map_color );
	}

	for ( auto iter = get_mingraph_owner()->whole_pose_context_enmeths_begin(),
			iter_end = get_mingraph_owner()->whole_pose_context_enmeths_end();
			iter != iter_end; ++iter ) {
		if ( (*iter)->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose ) ) {
			add_sfs_drs_enmeth( *iter );
		}
	}
}

MinimizationGraph const * MinimizationNode::get_mingraph_owner() const
{
	return static_cast< MinimizationGraph const * > ( get_owner() );
}

MinimizationGraph * MinimizationNode::get_mingraph_owner()
{
	return static_cast< MinimizationGraph * > ( get_owner() );
}

void MinimizationNode::add_active_1benmeth_std( OneBodyEnergyCOP enmeth ) { active_1benmeths_.push_back( enmeth ); active_1benmeths_std_.push_back( enmeth );}
void MinimizationNode::add_active_1benmeth_ext( OneBodyEnergyCOP enmeth ) { active_1benmeths_.push_back( enmeth ); active_1benmeths_ext_.push_back( enmeth ); }
void MinimizationNode::add_dof_deriv_1benmeth( OneBodyEnergyCOP enmeth ) { dof_deriv_1benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfs_dm_1benmeth( OneBodyEnergyCOP enmeth ) { sfs_dm_req_1benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfd_1benmeth( OneBodyEnergyCOP enmeth ) { sfd_req_1benmeths_.push_back( enmeth ); }

void MinimizationNode::add_active_2benmeth_std( TwoBodyEnergyCOP enmeth ) { active_intrares2benmeths_.push_back( enmeth ); active_intrares2benmeths_std_.push_back( enmeth );}
void MinimizationNode::add_active_2benmeth_ext( TwoBodyEnergyCOP enmeth ) { active_intrares2benmeths_.push_back( enmeth ); active_intrares2benmeths_ext_.push_back( enmeth ); }
void MinimizationNode::add_dof_deriv_2benmeth( TwoBodyEnergyCOP enmeth ) { dof_deriv_2benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfs_dm_2benmeth( TwoBodyEnergyCOP enmeth ) { sfs_dm_req_2benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfd_2benmeth( TwoBodyEnergyCOP enmeth ) { sfd_req_2benmeths_.push_back( enmeth ); }

void MinimizationNode::add_sfs_drs_enmeth( EnergyMethodCOP enmeth ) { sfs_drs_req_enmeths_.push_back( enmeth ); }

bool
MinimizationNode::classify_onebody_enmeth( OneBodyEnergyCOP enmeth, Residue const & rsd, Pose const & pose, int domain_map_color )
{
	if ( domain_map_color == 0 || enmeth->method_type() == methods::cd_1b ) {
		if ( enmeth->defines_score_for_residue( rsd ) ) {
			_activate_1benmeth( enmeth );

			if ( enmeth->defines_dof_derivatives( pose ) ) {
				add_dof_deriv_1benmeth( enmeth );
			}
			_check_sfs_and_sfd_for_1benmeth( enmeth, pose );
			return true;
		} else {
			return false;
		}
	} else {
		return false; // need to think about this
	}
}

bool MinimizationNode::classify_twobody_enmeth(
	TwoBodyEnergyCOP enmeth,
	Residue const & rsd,
	Pose const & pose,
	EnergyMap const & weights,
	int domain_map_color
)
{
	bool added( false );
	// should we let energy methods setup for scoring for a residue even if they
	// don't define an intrares energy for that residue?  Yes.  These methods
	// may need to initialize per-residue data that is used during two-body scoring,
	// such as water coordinates for LKBall term.
	_check_sfs_and_sfd_for_2benmeth( enmeth, pose );

	// Domain map check only prevents intra-residue sfxn/deriv evaluations; but not
	// SFS or SFD registration.
	if ( domain_map_color == 0 ) {
		if ( enmeth->defines_intrares_energy( weights ) && enmeth->defines_intrares_energy_for_residue( rsd ) ) {
			_activate_2benmeth( enmeth );
			added = true;
		} else {
			added = false;
		}
		if ( enmeth->defines_intrares_dof_derivatives( pose ) ) {
			add_dof_deriv_2benmeth( enmeth );
		}
	}
	return added;
}

void
MinimizationNode::_activate_1benmeth( OneBodyEnergyCOP enmeth )
{
	if ( enmeth->use_extended_residue_energy_interface() ) {
		add_active_1benmeth_ext( enmeth );
	} else {
		add_active_1benmeth_std( enmeth );
	}
}


void
MinimizationNode::_check_sfs_and_sfd_for_1benmeth(
	OneBodyEnergyCOP enmeth,
	pose::Pose const & pose
)
{
	if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose ) ) {
		add_sfs_dm_1benmeth( enmeth );
	} else if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose ) ) {
		add_sfs_drs_enmeth( enmeth );
	}
	if ( enmeth->requires_a_setup_for_derivatives_for_residue_opportunity( pose ) ) {
		add_sfd_1benmeth( enmeth );
	}
}

void
MinimizationNode::_activate_2benmeth( TwoBodyEnergyCOP enmeth ) {
	if ( enmeth->use_extended_intrares_energy_interface() ) {
		add_active_2benmeth_ext( enmeth );
	} else {
		add_active_2benmeth_std( enmeth );
	}
}

void
MinimizationNode::_check_sfs_and_sfd_for_2benmeth(
	TwoBodyEnergyCOP enmeth,
	pose::Pose const & pose
)
{
	if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity_during_minimization( pose ) ) {
		add_sfs_dm_2benmeth( enmeth );
	} else if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose ) ) {
		add_sfs_drs_enmeth( enmeth );
	}
	if ( enmeth->requires_a_setup_for_derivatives_for_residue_opportunity( pose ) ) {
		add_sfd_2benmeth( enmeth );
	}
}


MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_begin() const {
	return active_1benmeths_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_end() const {
	return active_1benmeths_.end();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_std_begin() const {
	return active_1benmeths_std_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_std_end() const {
	return active_1benmeths_std_.end();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_ext_begin() const {
	return active_1benmeths_ext_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::active_1benmeths_ext_end() const {
	return active_1benmeths_ext_.end();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::dof_deriv_1benmeths_begin() const {
	return dof_deriv_1benmeths_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::dof_deriv_1benmeths_end() const {
	return dof_deriv_1benmeths_.end();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::sfs_dm_req_1benmeths_begin() const {
	return sfs_dm_req_1benmeths_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::sfs_dm_req_1benmeths_end() const {
	return sfs_dm_req_1benmeths_.end();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::sfd_req_1benmeths_begin() const {
	return sfd_req_1benmeths_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::sfd_req_1benmeths_end() const {
	return sfd_req_1benmeths_.end();
}


MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_begin() const {
	return active_intrares2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_end() const {
	return active_intrares2benmeths_.end();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_std_begin() const {
	return active_intrares2benmeths_std_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_std_end() const {
	return active_intrares2benmeths_std_.end();
}


MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_ext_begin() const {
	return active_intrares2benmeths_ext_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::active_intrares2benmeths_ext_end() const {
	return active_intrares2benmeths_ext_.end();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::dof_deriv_2benmeths_begin() const {
	return dof_deriv_2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::dof_deriv_2benmeths_end() const {
	return dof_deriv_2benmeths_.end();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfs_dm_req_2benmeths_begin() const {
	return sfs_dm_req_2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfs_dm_req_2benmeths_end() const {
	return sfs_dm_req_2benmeths_.end();
}


MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfd_req_2benmeths_begin() const {
	return sfd_req_2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfd_req_2benmeths_end() const {
	return sfd_req_2benmeths_.end();
}

MinimizationNode::EnergyMethodsIterator
MinimizationNode::sfs_drs_req_enmeths_begin() const {
	return sfs_drs_req_enmeths_.begin();
}

MinimizationNode::EnergyMethodsIterator
MinimizationNode::sfs_drs_req_enmeths_end() const {
	return sfs_drs_req_enmeths_.end();
}

///////// Minimization Edge Class /////////////

/// @brief Minimization edge ctor
MinimizationEdge::MinimizationEdge(
	MinimizationGraph * owner,
	Size n1,
	Size n2
) :
	parent( owner, n1, n2 ),
	weight_( 1.0 ),
	dweight_( 1.0 )
{}

MinimizationEdge::MinimizationEdge(
	MinimizationGraph * owner,
	MinimizationEdge const & example_edge
)
:
	parent( owner, example_edge.get_first_node_ind(), example_edge.get_second_node_ind() ),
	weight_( 1.0 ),
	dweight_( 1.0 )
{
	MinimizationEdge::copy_from( & example_edge ); // Virtual dispatch doesn't work fully in constructor
}

/// @brief virtual dstor; The MinimizationEdge must free the array pool element it
/// holds before it disappears.
MinimizationEdge::~MinimizationEdge() = default;

/// @brief copies data from MinimizationEdge const * source;
///
/// called from the copy ctor and operator= methods defined in the Graph base class
void MinimizationEdge::copy_from( parent const * source )
{
	//MinimizationEdge const * ee = static_cast< MinimizationEdge const * > ( source );
	// down_cast is *supposed* to assert the dynamic cast in debug builds; doesn't work for some reason
	auto const & minedge = static_cast< MinimizationEdge const & > ( * source );

	res_pair_min_data_ = minedge.res_pair_min_data_;
	twobody_enmeths_ = minedge.twobody_enmeths_;
	active_2benmeths_ = minedge.active_2benmeths_;
	active_2benmeths_std_ = minedge.active_2benmeths_std_;
	active_2benmeths_ext_ = minedge.active_2benmeths_ext_;
	sfs_req_2benmeths_ = minedge.sfs_req_2benmeths_;
	sfd_req_2benmeths_ = minedge.sfd_req_2benmeths_;
	weight_ = minedge.weight_;
	dweight_ = minedge.dweight_;
}


/// @brief virtual call to determine the static size of an Edge object
/// dynamic memory use is counted through the recursive count_dynamic_memory()
/// calling path
Size MinimizationEdge::count_static_memory() const
{
	return sizeof ( MinimizationEdge );
}

/// @brief virtual call to determine the amount of dynamic memory allocated by
/// an edge; this function must recurse to the parent class to determine how
/// much memory the parent class is responsible for.  Do not account for
/// the size of the ArrayPool array here; instead, that is accounted for in
/// the MinimizationGraph::count_dynamic_memory method.
Size MinimizationEdge::count_dynamic_memory() const
{
	Size tot = 0;
	tot += parent::count_dynamic_memory(); //recurse to parent
	/// add in cost of storing the residue-pair minimization data here
	return tot;
}

bool MinimizationEdge::add_twobody_enmeth(
	TwoBodyEnergyCOP enmeth,
	Residue const & rsd1,
	Residue const & rsd2,
	pose::Pose const & pose,
	bool res_moving_wrt_eachother
)
{
	twobody_enmeths_.push_back( enmeth );
	return classify_twobody_enmeth( enmeth, rsd1, rsd2, pose, res_moving_wrt_eachother );

}

void
MinimizationEdge::activate_dof_deriv_two_body_method(
	TwoBodyEnergyCOP enmeth,
	pose::Pose const & pose
)
{
	debug_assert( enmeth->defines_intrares_dof_derivatives( pose ) );
	if ( std::find( active_2benmeths_.begin(), active_2benmeths_.end(), enmeth ) !=
			active_2benmeths_.end() ) {
		return;
	}

	_activate_2benmeth( enmeth );
	_check_sfs_and_sfd_for_2benmeth( enmeth, pose );
}


void MinimizationEdge::setup_for_minimizing(
	Residue const & rsd1,
	Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
)
{
	for ( TwoBodyEnergiesIterator iter = active_2benmeths_.begin(),
			iter_end = active_2benmeths_.end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing_for_residue_pair(
			rsd1, rsd2, pose, sfxn, min_map,
			get_minimization_node( 0 )->res_min_data(),
			get_minimization_node( 1 )->res_min_data(),
			res_pair_min_data_ );
	}
}

void MinimizationEdge::setup_for_scoring(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const & pose,
	ScoreFunction const & sfxn
)
{
	for ( auto iter = sfs_req_2benmeths_begin(),
			iter_end = sfs_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue_pair(
			rsd1, rsd2,
			get_minimization_node( 0 )->res_min_data(),
			get_minimization_node( 1 )->res_min_data(),
			pose, sfxn, res_pair_min_data_ );
	}
}

void MinimizationEdge::setup_for_derivatives(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const & pose,
	ScoreFunction const & sfxn
)
{
	for ( auto iter = sfd_req_2benmeths_begin(),
			iter_end = sfd_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives_for_residue_pair(
			rsd1, rsd2,
			get_minimization_node( 0 )->res_min_data(),
			get_minimization_node( 1 )->res_min_data(),
			pose, sfxn, res_pair_min_data_ );
	}

}

void MinimizationEdge::reinitialize_active_energy_methods(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const & pose,
	bool res_moving_wrt_eachother
)
{
	active_2benmeths_.clear();
	active_2benmeths_std_.clear();
	active_2benmeths_ext_.clear();
	sfs_req_2benmeths_.clear();
	sfd_req_2benmeths_.clear();

	for ( TwoBodyEnergiesIterator iter = twobody_enmeths_.begin(),
			iter_end = twobody_enmeths_.end(); iter != iter_end; ++iter ) {
		classify_twobody_enmeth( *iter, rsd1, rsd2, pose, res_moving_wrt_eachother );
	}
}

void MinimizationEdge::add_active_enmeth_std( TwoBodyEnergyCOP enmeth ) { active_2benmeths_.push_back( enmeth ); active_2benmeths_std_.push_back( enmeth );}
void MinimizationEdge::add_active_enmeth_ext( TwoBodyEnergyCOP enmeth ) { active_2benmeths_.push_back( enmeth ); active_2benmeths_ext_.push_back( enmeth ); }
void MinimizationEdge::add_sfs_enmeth( TwoBodyEnergyCOP enmeth ) { sfs_req_2benmeths_.push_back( enmeth ); }
void MinimizationEdge::add_sfd_enmeth( TwoBodyEnergyCOP enmeth ) { sfd_req_2benmeths_.push_back( enmeth ); }

bool MinimizationEdge::classify_twobody_enmeth(
	TwoBodyEnergyCOP enmeth,
	Residue const & rsd1,
	Residue const & rsd2,
	pose::Pose const & pose,
	bool res_moving_wrt_eachother
)
{
	if ( ! enmeth->defines_score_for_residue_pair( rsd1, rsd2, res_moving_wrt_eachother ) ) return false;

	_activate_2benmeth( enmeth );
	_check_sfs_and_sfd_for_2benmeth( enmeth, pose );
	return true;
}

void
MinimizationEdge::_activate_2benmeth( TwoBodyEnergyCOP enmeth ) {
	if ( enmeth->use_extended_residue_pair_energy_interface() ) {
		add_active_enmeth_ext( enmeth );
	} else {
		add_active_enmeth_std( enmeth );
	}
}

void
MinimizationEdge::_check_sfs_and_sfd_for_2benmeth( TwoBodyEnergyCOP enmeth, Pose const & pose )
{
	if ( enmeth->requires_a_setup_for_scoring_for_residue_pair_opportunity( pose ) ) {
		add_sfs_enmeth( enmeth );
	}
	if ( enmeth->requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose ) ) {
		add_sfd_enmeth( enmeth );
	}
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_begin() const {
	return active_2benmeths_.begin();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_end() const {
	return active_2benmeths_.end();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_std_begin() const {
	return active_2benmeths_std_.begin();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_std_end() const {
	return active_2benmeths_std_.end();
}


MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_ext_begin() const {
	return active_2benmeths_ext_.begin();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::active_2benmeths_ext_end() const {
	return active_2benmeths_ext_.end();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::sfs_req_2benmeths_begin() const {
	return sfs_req_2benmeths_.begin();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::sfs_req_2benmeths_end() const {
	return sfs_req_2benmeths_.end();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::sfd_req_2benmeths_begin() const {
	return sfd_req_2benmeths_.begin();
}

MinimizationEdge::TwoBodyEnergiesIterator
MinimizationEdge::sfd_req_2benmeths_end() const {
	return sfd_req_2benmeths_.end();
}

///////// Minimization Graph Class /////////////

MinimizationEdge * MinimizationGraph::find_minimization_edge( Size n1, Size n2)
{
	Edge * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< MinimizationEdge * > ( edge );
		//return utility::down_cast< MinimizationEdge * > ( edge );
	} else {
		return nullptr;
	}
}

MinimizationEdge const * MinimizationGraph::find_minimization_edge( Size n1, Size n2) const
{
	Edge const * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< MinimizationEdge const * > ( edge );
		//return utility::down_cast< MinimizationEdge const * > ( edge );
	} else {
		return nullptr;
	}
}


MinimizationGraph::MinimizationGraph()
:
	parent(),
	minimization_edge_pool_( new boost::unordered_object_pool< MinimizationEdge > ( 256 ) )
{}

/// @details This does not call the base class parent( Size ) constructor since
/// that produces calls to the polymorphic function create_new_node() and polymorphism
/// does not work during constructor intialization.
MinimizationGraph::MinimizationGraph( Size num_nodes )
:
	parent(),
	minimization_edge_pool_( new boost::unordered_object_pool< MinimizationEdge > ( 256 ) )
{
	set_num_nodes( num_nodes );
}

/// @details Notice that this does not call the parent( src ) copy constructor.
/// This is because the copy constructor relies on polymorphic functions which
/// are unavailable during the Graph constructor.  Instead, this function waits
/// until parent construction is complete, and relies on the assigmnent operator.
MinimizationGraph::MinimizationGraph( // NOLINT(bugprone-copy-constructor-init) -- deliberately not calling parent copy constructor
	MinimizationGraph const & src
):
	parent( ),
	minimization_edge_pool_( new boost::unordered_object_pool< MinimizationEdge > ( 256 ) )
{
	parent::operator = ( src );
}

MinimizationGraph::~MinimizationGraph() {
	delete_everything();
	delete minimization_edge_pool_; minimization_edge_pool_ = nullptr;
}


/// @brief assignment operator -- performs a deep copy
MinimizationGraph &
MinimizationGraph::operator = ( MinimizationGraph const & rhs )
{
	parent::operator = ( rhs );
	return *this;
}

void MinimizationGraph::delete_edge( utility::graph::Edge * edge )
{
	debug_assert( dynamic_cast< MinimizationEdge* > (edge) );
	minimization_edge_pool_->destroy( static_cast< MinimizationEdge* > (edge) );
}

Size MinimizationGraph::count_static_memory() const
{
	return sizeof ( MinimizationGraph );
}

Size MinimizationGraph::count_dynamic_memory() const
{
	Size tot(0);
	tot += parent::count_dynamic_memory(); //recurse to parent
	return tot;
}

void MinimizationGraph::add_whole_pose_context_enmeth( EnergyMethodCOP enmeth, core::pose::Pose const & pose )
{
	whole_pose_context_enmeths_.push_back( enmeth );
	if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity_during_regular_scoring( pose ) ) {
		for ( core::Size ii = 1; ii <= num_nodes(); ++ii ) {
			get_minimization_node( ii )->add_sfs_drs_enmeth( enmeth );
		}
	}
}

MinimizationGraph::EnergiesIterator
MinimizationGraph::whole_pose_context_enmeths_begin() const
{
	return whole_pose_context_enmeths_.begin();
}

MinimizationGraph::EnergiesIterator
MinimizationGraph::whole_pose_context_enmeths_end() const
{
	return whole_pose_context_enmeths_.end();
}

void MinimizationGraph::set_fixed_energies( EnergyMap const & vals ) { fixed_energies_ = vals; }
EnergyMap const & MinimizationGraph::fixed_energies() const { return fixed_energies_; }


/// @brief Factory method for node creation
Node*
MinimizationGraph::create_new_node( Size index )
{
	return new MinimizationNode( this, index );
}

/// @brief Factory method for edge creation
Edge*
MinimizationGraph::create_new_edge( Size index1, Size index2 )
{
	return minimization_edge_pool_->construct( this, index1, index2 );
}

/// @brief Factory copy-constructor method for edge creation
Edge*
MinimizationGraph::create_new_edge( Edge const * example_edge )
{
	return minimization_edge_pool_->construct(
		this,
		static_cast< MinimizationEdge const & > (*example_edge)
	);
}


///// NON MEMBER FUNCTIONS////

void
eval_atom_derivatives_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap const & res_weights,
	utility::vector1< DerivVectorPair > & atom_derivs
) {
	for ( auto
			iter = min_node.active_1benmeths_begin(),
			iter_end = min_node.active_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->eval_residue_derivatives(
			rsd, min_node.res_min_data(), pose, res_weights, atom_derivs );
	}
	/// 1b 2body intraresidue contributions
	for ( auto
			iter = min_node.active_intrares2benmeths_begin(),
			iter_end = min_node.active_intrares2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_derivatives(
			rsd, min_node.res_min_data(), pose, res_weights, atom_derivs );
	}

}

void
eval_res_onebody_energies_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap // accumulated into
) {
	for ( auto
			iter = min_node.active_1benmeths_std_begin(),
			iter_end = min_node.active_1benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, emap );
	}
	for ( auto
			iter = min_node.active_1benmeths_ext_begin(),
			iter_end = min_node.active_1benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy_ext( rsd, min_node.res_min_data(), pose, emap );
	}
	/// 1b 2body intraresidue contributions
	for ( auto
			iter = min_node.active_intrares2benmeths_std_begin(),
			iter_end = min_node.active_intrares2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, sfxn, emap );
	}
	for ( auto
			iter = min_node.active_intrares2benmeths_ext_begin(),
			iter_end = min_node.active_intrares2benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy_ext( rsd, min_node.res_min_data(), pose, sfxn, emap );
	}

}

void
eval_atom_derivatives_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	ResSingleMinimizationData const & res1_min_data,
	ResSingleMinimizationData const & res2_min_data,
	pose::Pose const & pose,
	EnergyMap const & respair_weights,
	utility::vector1< DerivVectorPair > & r1atom_derivs,
	utility::vector1< DerivVectorPair > & r2atom_derivs
) {
	for ( auto
			iter = min_edge.active_2benmeths_begin(),
			iter_end = min_edge.active_2benmeths_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_residue_pair_derivatives(
			res1, res2, res1_min_data, res2_min_data, min_edge.res_pair_min_data(),
			pose, respair_weights, r1atom_derivs, r2atom_derivs );
	}
}

void
eval_weighted_atom_derivatives_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	ResSingleMinimizationData const & res1_min_data,
	ResSingleMinimizationData const & res2_min_data,
	pose::Pose const & pose,
	EnergyMap const & respair_weights,
	utility::vector1< DerivVectorPair > & r1atom_derivs,
	utility::vector1< DerivVectorPair > & r2atom_derivs
) {
	//fpd rather then change eval_residue_pair_derivatives interface
	//    we will just modify the energymap
	EnergyMap respair_weight_new = respair_weights;
	respair_weight_new *= min_edge.dweight();

	for ( auto
			iter = min_edge.active_2benmeths_begin(),
			iter_end = min_edge.active_2benmeths_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_residue_pair_derivatives(
			res1, res2, res1_min_data, res2_min_data, min_edge.res_pair_min_data(),
			pose, respair_weight_new, r1atom_derivs, r2atom_derivs );
	}
}

void
eval_res_pair_energy_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
) {
	for ( auto
			iter = min_edge.active_2benmeths_std_begin(),
			iter_end = min_edge.active_2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy(
			res1, res2, pose, sfxn, emap );
	}
	for ( auto
			iter = min_edge.active_2benmeths_ext_begin(),
			iter_end = min_edge.active_2benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy_ext(
			res1, res2, min_edge.res_pair_min_data(), pose, sfxn, emap );
	}
}

Real
eval_dof_deriv_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
) {
	Real deriv( 0 );
	/// 1. eval 1 body derivatives
	/// 1a 1body energy methods
	for ( auto
			iter = min_node.dof_deriv_1benmeths_begin(),
			iter_end = min_node.dof_deriv_1benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_residue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}

	for ( auto
			iter = min_node.dof_deriv_2benmeths_begin(),
			iter_end = min_node.dof_deriv_2benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_intraresidue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}
	return deriv;
}

void
eval_weighted_res_onebody_energies_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap, // accumulated into
	EnergyMap & scratch_emap
) {
	for ( auto
			iter = min_node.active_1benmeths_std_begin(),
			iter_end = min_node.active_1benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( auto
			iter = min_node.active_1benmeths_ext_begin(),
			iter_end = min_node.active_1benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy_ext( rsd, min_node.res_min_data(), pose, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	/// 1b 2body intraresidue contributions
	for ( auto
			iter = min_node.active_intrares2benmeths_std_begin(),
			iter_end = min_node.active_intrares2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( auto
			iter = min_node.active_intrares2benmeths_ext_begin(),
			iter_end = min_node.active_intrares2benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy_ext( rsd, min_node.res_min_data(), pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
}

void
eval_weighted_res_pair_energy_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap,
	EnergyMap & scratch_emap // should be zeros coming in, left zeroed at the end;
) {
	for ( auto
			iter = min_edge.active_2benmeths_std_begin(),
			iter_end = min_edge.active_2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy(
			res1, res2, pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_edge.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( auto
			iter = min_edge.active_2benmeths_ext_begin(),
			iter_end = min_edge.active_2benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy_ext(
			res1, res2, min_edge.res_pair_min_data(), pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_edge.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
}

Real
eval_weighted_dof_deriv_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
) {
	Real deriv( 0 );
	/// 1. eval 1 body derivatives
	/// 1a 1body energy methods
	for ( auto
			iter = min_node.dof_deriv_1benmeths_begin(),
			iter_end = min_node.dof_deriv_1benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_residue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}

	for ( auto
			iter = min_node.dof_deriv_2benmeths_begin(),
			iter_end = min_node.dof_deriv_2benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_intraresidue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}
	return deriv * min_node.weight();
}

void
create_and_store_atom_tree_minimization_graph_asym(
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose
);

void
create_and_store_atom_tree_minimization_graph_symm(
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose
);

void
create_and_store_atom_tree_minimization_graph(
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose
)
{
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		create_and_store_atom_tree_minimization_graph_symm( sfxn, min_map, pose );
	} else {
		create_and_store_atom_tree_minimization_graph_asym( sfxn, min_map, pose );
	}
}

void
create_and_store_atom_tree_minimization_graph_asym(
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose
)
{
	/// 1. Initialize the nodes in the energy graph with one-body
	/// and the intra-residue two-body energy methods

	/// 2. Initialize the short-ranged edges in the energy graph with
	/// the short-ranged two-body energy methods.

	/// 3. Initialize any additional edges with the long-range
	/// two-body energies

	/// 4. Drop any edges from the minimization graph that produce no
	/// energies, and call setup-for-minimization on all edges that remain

	/// 5. Setup whole-structure energies and energy methods that opt-out
	/// of the minimization-graph control over derivative evaluation.

	MinimizationGraphOP g( new MinimizationGraph( pose.size() ) );
	std::list< methods::EnergyMethodCOP > eval_derivs_with_pose_enmeths;
	for ( auto const & all_method : sfxn.all_methods() ) {
		if ( all_method->defines_high_order_terms( pose ) || all_method->minimize_in_whole_structure_context( pose ) ) {
			eval_derivs_with_pose_enmeths.push_back( all_method );
		}
	}

	EnergyMap fixed_energies; // portions of the score function that will not change over the course of minimization.

	kinematics::DomainMap const & domain_map = min_map.domain_map();
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		sfxn.setup_for_minimizing_for_node( * g->get_minimization_node( ii ), pose.residue( ii ), pose.residue_data( ii ),
			min_map, pose, true, fixed_energies );
	}

	g->copy_connectivity( pose.energies().energy_graph() );
	// 2. Now initialize the edges of the minimization graph using the edges of the EnergyGraph;
	// The energy graph should be up-to-date before this occurs
	for ( utility::graph::Graph::EdgeListIter
			edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(),
			ee_edge_iter = pose.energies().energy_graph().edge_list_begin();
			edge_iter != edge_iter_end; ++edge_iter, ++ee_edge_iter ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();
		debug_assert( node1 == (*ee_edge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( node2 == (*ee_edge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering

		auto & minedge( static_cast< MinimizationEdge & > (**edge_iter) );
		// domain map check here?
		bool const res_moving_wrt_eachother(
			domain_map( node1 ) == 0 ||
			domain_map( node2 ) == 0 ||
			domain_map( node1 ) != domain_map( node2 ) );

		sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
			pose.residue( node1 ), pose.residue( node2 ),
			minedge, min_map, pose, res_moving_wrt_eachother, true,
			static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies );
	}

	/// 3. Long range energies need time to get included into the graph, which may require the addition of new edges
	/// 3a: CILR2B energies
	///    i.   Iterate across all ci2b long range energy methods
	///    ii.  Iterate across all residue pairs indicated in the long-range energy containers
	///    iii. If two residues have the same non-zero domain-map coloring, then accumulate their interaction energy and move on
	///    iv.  otherwise, find the corresponding minimization-graph edge for this residue pair
	///    v.   adding a new edge if necessary,
	///    vi.  and prepare the minimization data for this edge

	for ( auto
			iter = sfxn.long_range_energies_begin(),
			iter_end = sfxn.long_range_energies_end();
			iter != iter_end; ++iter ) {

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();
				bool const res_moving_wrt_eachother(
					domain_map( ii ) == 0 ||
					domain_map( jj ) == 0 ||
					domain_map( ii ) != domain_map( jj ) );
				sfxn.setup_for_lr2benmeth_minimization_for_respair(
					pose.residue( ii ), pose.residue( jj ), *iter, *g, min_map, pose,
					res_moving_wrt_eachother, true, rni, fixed_energies );
			}
		}
	}

	/// 4. Call setup_for_minimizing on each edge that has active twobody energies, and drop
	/// all other edges.
	for ( utility::graph::Graph::EdgeListIter edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		auto & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		utility::graph::Graph::EdgeListIter edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, sfxn, min_map );
		} else {
			/// The edge will not contribute anything to scoring during minimization,
			/// so delete it from the graph, so we don't have to pay the expense of traversing
			/// through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 5.  Whole structure energies and energies that are opting out of the MinimizationGraph
	/// routines get a chance to setup for minimizing (using the entire pose as context) and
	for ( std::list< methods::EnergyMethodCOP >::const_iterator
			iter     = eval_derivs_with_pose_enmeths.begin(),
			iter_end = eval_derivs_with_pose_enmeths.end();
			iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing( pose, sfxn, min_map );
		g->add_whole_pose_context_enmeth( *iter, pose );
	}


	//std::cout << "Fixed energies: ";
	//fixed_energies.show_if_nonzero_weight( std::cout, weights_ );
	//std::cout << std::endl;

	g->set_fixed_energies( fixed_energies );
	pose.energies().set_minimization_graph( g );
}

void
create_and_store_atom_tree_minimization_graph_symm(
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map,
	pose::Pose & pose
)
{
	/// 1. Initialize the nodes of the minimization graph
	/// 2. Initialize the edges with the short-ranged two-body energies
	/// 3. Initialize the edges with the long-ranged two-body energies
	/// 4. Run setup-for-minimization on the edges of the mingraph; dropping edges with no active 2b enmeths.
	/// 5. Let whole-structure energies initialize themselves

	using namespace core::conformation::symmetry;
	using namespace core::scoring::symmetry;

	auto & symm_conf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
	SymmetryInfo const & symm_info( * symm_conf.Symmetry_Info() );
	auto & symm_energies( dynamic_cast< SymmetricEnergies & > ( pose.energies()) );

	MinimizationGraphOP g( new MinimizationGraph( pose.size() ) );
	MinimizationGraphOP dg( new MinimizationGraph( pose.size() ) ); // derivative graph

	std::list< methods::EnergyMethodCOP > eval_derivs_with_pose_enmeths;
	for ( auto const & iter : sfxn.all_methods() ) {
		if ( iter->defines_high_order_terms( pose ) || iter->minimize_in_whole_structure_context( pose ) ) {
			eval_derivs_with_pose_enmeths.push_back( iter );
		}
	}

	EnergyMap fixed_energies; // portions of the score function that will not change over the course of minimization.

	/// Accumulate the portions of the scoring function for those residues with non-zero domain-map values
	/// but also let all energy methods register with each node, since some may require a setup-for-scoring opportunity
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {

		bool accumulate_fixed_energies( symm_info.fa_is_independent(ii) &&
			ii > symm_info.num_total_residues_without_pseudo() );

		sfxn.setup_for_minimizing_for_node( * g->get_minimization_node( ii ), pose.residue( ii ),
			pose.residue_data( ii ), min_map, pose, accumulate_fixed_energies, fixed_energies );
		g->get_minimization_node( ii )->weight( symm_info.score_multiply_factor() );
		sfxn.setup_for_minimizing_for_node( * dg->get_minimization_node( ii ), pose.residue( ii ),
			pose.residue_data( ii ), min_map, pose, false, fixed_energies ); // only accumulate once
	}
	g->copy_connectivity(  pose.energies().energy_graph() );
	dg->copy_connectivity( pose.energies().energy_graph() );

	kinematics::DomainMap const & domain_map( min_map.domain_map() );
	for ( utility::graph::Graph::EdgeListIter
			edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(),
			dedge_iter = dg->edge_list_begin(),
			//dedge_iter_end = dg->edge_list_end(),
			ee_edge_iter = pose.energies().energy_graph().edge_list_begin();
			edge_iter != edge_iter_end; ++edge_iter, ++dedge_iter, ++ee_edge_iter ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();
		debug_assert( node1 == (*ee_edge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( node2 == (*ee_edge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( node1 == (*dedge_iter)->get_first_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( node2 == (*dedge_iter)->get_second_node_ind() ); // copy_connectivity should preserve the energy-graph edge ordering
		debug_assert( symm_info.bb_follows( node1 ) == 0 || symm_info.bb_follows( node2 ) == 0 );

		// domain map check here?
		bool const res_moving_wrt_eachother(
			domain_map( node1 ) == 0 ||
			domain_map( node2 ) == 0 ||
			domain_map( node1 ) != domain_map( node2 ) );

		Real edge_weight = symm_info.score_multiply( node1, node2 );
		Real edge_dweight = symm_info.deriv_multiply( node1, node2 );

		// fd, 7/17, removing this separate logic for old_sym_min and new_sym_min
		//  the reason is that the derivative weights serve 2 purposes:
		//    1) handing scaling of jumps in lattice systems
		//    2) handing computation of derivatives on edges where the weight on scoring is 0
		//  the new_sym_min logic (commented out) was meant to address the first case but broke the second case
		//  instead, we will address the first case in SymmetryInfo and restore old_sym_min behavior here

		//bool const new_sym_min( !basic::options::option[ basic::options::OptionKeys::optimization::old_sym_min ]() );
		// if ( new_sym_min ) { // new way
		// } else { // classic way
		if ( edge_weight != 0.0 ) {
			auto & minedge( static_cast< MinimizationEdge & > (**edge_iter) );
			sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
				pose.residue( node1 ), pose.residue( node2 ),
				minedge, min_map, pose, res_moving_wrt_eachother, true,
				static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies, edge_weight );
			minedge.weight( edge_weight );
			minedge.dweight( edge_dweight );
		} else {
			auto & minedge( static_cast< MinimizationEdge & > (**dedge_iter) );
			sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
				pose.residue( node1 ), pose.residue( node2 ),
				minedge, min_map, pose, res_moving_wrt_eachother, false,
				static_cast< EnergyEdge const * > (*ee_edge_iter), fixed_energies ); // edge weight of 1
			minedge.dweight( edge_dweight );
		}
	}

	/// 3. Long range energies need time to get included into the graph, which may require the addition of new edges
	/// 3a: CILR2B energies
	///    i.   Iterate across all ci2b long range energy methods
	///    ii.  Iterate across all residue pairs indicated in the long-range energy containers
	///    iii. If two residues have the same non-zero domain-map coloring, then accumulate their interaction energy and move on
	///    iv.  otherwise, find the corresponding minimization-graph edge for this residue pair
	///    v.   adding a new edge if necessary,
	///    vi.  and prepare the minimization data for this edge

	for ( auto
			iter = sfxn.long_range_energies_begin(),
			iter_end = sfxn.long_range_energies_end();
			iter != iter_end; ++iter ) {
		// NO! add these terms to the minimization graph for scoring, even if not for derivative evaluation
		///if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const jj = rni->upper_neighbor_id();
				bool const res_moving_wrt_eachother(
					domain_map( ii ) == 0 ||
					domain_map( jj ) == 0 ||
					domain_map( ii ) != domain_map( jj ) );

				Real edge_weight = symm_info.score_multiply( ii, jj );

				// fd, 7/17, removing this separate logic for old and new
				//  see comment above

				//bool const new_sym_min( !basic::options::option[ basic::options::OptionKeys::optimization::old_sym_min ]() );
				// if ( new_sym_min ) { // new way
				// } else { // classic way

				if ( edge_weight != 0.0 ) {
					// adjust/add the edge to the scoring graph
					sfxn.setup_for_lr2benmeth_minimization_for_respair(
						pose.residue( ii ), pose.residue( jj ), *iter, *g, min_map, pose,
						res_moving_wrt_eachother, true, rni, fixed_energies, edge_weight );
				} else {
					// adjust/add this edge to the derivative graph
					sfxn.setup_for_lr2benmeth_minimization_for_respair(
						pose.residue( ii ), pose.residue( jj ), *iter, *dg, min_map, pose,
						res_moving_wrt_eachother, false, rni, fixed_energies ); // no edge weight
				}
			}
		}
	}

	/// 4a. drop unused edges from the scoring graph; call setup for minimizing on those remaining
	for ( utility::graph::Graph::EdgeListIter edge_iter = g->edge_list_begin(),
			edge_iter_end = g->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		auto & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		utility::graph::Graph::EdgeListIter edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			//std::cout << " active scoring graph edge: " << node1 << " " << node2 << std::endl;
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, sfxn, min_map );
		} else {
			/// The edge will not contribute anything to scoring during minimization,
			/// so delete it from the graph, so we don't have to pay the expense of traversing
			/// through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 4b. drop unused edges from the derivatives graph; call setup for minimizing on those remaining
	for ( utility::graph::Graph::EdgeListIter edge_iter = dg->edge_list_begin(),
			edge_iter_end = dg->edge_list_end(); edge_iter != edge_iter_end; /* no increment */ ) {
		Size const node1 = (*edge_iter)->get_first_node_ind();
		Size const node2 = (*edge_iter)->get_second_node_ind();

		auto & minedge( static_cast< MinimizationEdge & > (**edge_iter) );

		utility::graph::Graph::EdgeListIter edge_iter_next( edge_iter );
		++edge_iter_next;

		if ( minedge.any_active_enmeths() ) {
			//std::cout << " active deriv graph edge: " << node1 << " " << node2 << std::endl;
			minedge.setup_for_minimizing( pose.residue(node1), pose.residue(node2), pose, sfxn, min_map );
		} else {
			/// The edge will not contribute anything to derivative evaluation during minimization;
			/// it either represents an interaction that is not changing as the result of minimization,
			/// or the interaction is handled by the scoring graph. Delete it from the graph so
			/// we don't have to pay the expense of traversing through it.
			g->delete_edge( *edge_iter );
		}
		edge_iter = edge_iter_next;
	}

	/// 5.  Whole structure energies and energies that are opting out of the MinimizationGraph
	/// routines get a chance to setup for minimizing (using the entire pose as context) and
	for ( std::list< methods::EnergyMethodCOP >::const_iterator
			iter     = eval_derivs_with_pose_enmeths.begin(),
			iter_end = eval_derivs_with_pose_enmeths.end();
			iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing( pose, sfxn, min_map );
		g->add_whole_pose_context_enmeth( *iter, pose );
	}

	//std::cout << "Fixed energies: ";
	//fixed_energies.show_if_nonzero_weight( std::cout, weights() );
	//std::cout << std::endl;

	g->set_fixed_energies( fixed_energies );
	symm_energies.set_minimization_graph( g );
	symm_energies.set_derivative_graph( dg );
}


} //namespace scoring
} //namespace core
