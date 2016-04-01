// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/MinimizationGraph.cc
/// @brief  Minimization graph class implementation
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit Headers
#include <core/scoring/MinimizationGraph.hh>

// Package Headers
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>

// Numeric headers

// Boost Headers
#include <core/graph/unordered_object_pool.hpp>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>
#include <boost/pool/pool.hpp>


namespace core {
namespace scoring {

using namespace graph;

///////// Minimization Node Class /////////////

MinimizationNode::MinimizationNode( Graph * owner, Size index ) :
	parent( owner, index )
{}

MinimizationNode::~MinimizationNode() {}


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
	MinimizationNode const & mn_source = static_cast< MinimizationNode const & > ( * source );

	res_min_data_ = mn_source.res_min_data_; // deep copy

	onebody_enmeths_ = mn_source.onebody_enmeths_;
	active_1benmeths_ = mn_source.active_1benmeths_;
	active_1benmeths_std_ = mn_source.active_1benmeths_std_;
	active_1benmeths_ext_ = mn_source.active_1benmeths_ext_;
	dof_deriv_1benmeths_ = mn_source.dof_deriv_1benmeths_;
	sfs_req_1benmeths_ = mn_source.sfs_req_1benmeths_;
	sfd_req_1benmeths_ = mn_source.sfd_req_1benmeths_;
	twobody_enmeths_ = mn_source.twobody_enmeths_;
	active_intrares2benmeths_ = mn_source.active_intrares2benmeths_;
	active_intrares2benmeths_std_ = mn_source.active_intrares2benmeths_std_;
	active_intrares2benmeths_ext_ = mn_source.active_intrares2benmeths_ext_;
	dof_deriv_2benmeths_ = mn_source.dof_deriv_2benmeths_;
	sfs_req_2benmeths_ = mn_source.sfs_req_2benmeths_;
	sfd_req_2benmeths_ = mn_source.sfd_req_2benmeths_;
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

void
MinimizationNode::setup_for_minimizing(
	Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	kinematics::MinimizerMapBase const & min_map
) {
	for ( OneBodyEnergiesIterator iter = active_1benmeths_begin(),
			iter_end = active_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing_for_residue( rsd, pose, sfxn, min_map, res_min_data_ );
	}
	for ( TwoBodyEnergiesIterator iter = twobody_enmeths_.begin(),
			iter_end = twobody_enmeths_.end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_minimizing_for_residue( rsd, pose, sfxn, min_map, res_min_data_ );
	}
}

void MinimizationNode::setup_for_scoring(
	Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) {
	/// 1a 1body energy methods
	for ( OneBodyEnergiesIterator iter = sfs_req_1benmeths_begin(),
			iter_end = sfs_req_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue( rsd, pose, sfxn, res_min_data_ );
	}
	/// 1b 2body intraresidue contributions
	for ( TwoBodyEnergiesIterator iter = sfs_req_2benmeths_begin(),
			iter_end = sfs_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_scoring_for_residue( rsd, pose, sfxn, res_min_data_ );
	}
}

void MinimizationNode::setup_for_derivatives(
	Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn
) {
	/// 1a 1body energy methods
	for ( OneBodyEnergiesIterator iter = sfd_req_1benmeths_begin(),
			iter_end = sfd_req_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives_for_residue( rsd, pose, sfxn, res_min_data_ );
	}
	/// 1b 2body intraresidue contributions
	for ( TwoBodyEnergiesIterator iter = sfd_req_2benmeths_begin(),
			iter_end = sfd_req_2benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->setup_for_derivatives_for_residue( rsd, pose, sfxn, res_min_data_ );
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
	sfs_req_1benmeths_.clear();
	sfd_req_1benmeths_.clear();

	active_intrares2benmeths_.clear();
	active_intrares2benmeths_std_.clear();
	active_intrares2benmeths_ext_.clear();
	dof_deriv_2benmeths_.clear();
	sfs_req_2benmeths_.clear();
	sfd_req_2benmeths_.clear();

	for ( OneBodyEnergiesIterator iter = onebody_enmeths_.begin(),
			iter_end = onebody_enmeths_.end(); iter != iter_end; ++iter ) {
		classify_onebody_enmeth( *iter, rsd, pose, domain_map_color );
	}
	for ( TwoBodyEnergiesIterator iter = twobody_enmeths_.begin(),
			iter_end = twobody_enmeths_.end(); iter != iter_end; ++iter ) {
		classify_twobody_enmeth( *iter, rsd, pose, weights, domain_map_color );
	}
}

void MinimizationNode::add_active_1benmeth_std( OneBodyEnergyCOP enmeth ) { active_1benmeths_.push_back( enmeth ); active_1benmeths_std_.push_back( enmeth );}
void MinimizationNode::add_active_1benmeth_ext( OneBodyEnergyCOP enmeth ) { active_1benmeths_.push_back( enmeth ); active_1benmeths_ext_.push_back( enmeth ); }
void MinimizationNode::add_dof_deriv_1benmeth( OneBodyEnergyCOP enmeth ) { dof_deriv_1benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfs_1benmeth( OneBodyEnergyCOP enmeth ) { sfs_req_1benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfd_1benmeth( OneBodyEnergyCOP enmeth ) { sfd_req_1benmeths_.push_back( enmeth ); }

void MinimizationNode::add_active_2benmeth_std( TwoBodyEnergyCOP enmeth ) { active_intrares2benmeths_.push_back( enmeth ); active_intrares2benmeths_std_.push_back( enmeth );}
void MinimizationNode::add_active_2benmeth_ext( TwoBodyEnergyCOP enmeth ) { active_intrares2benmeths_.push_back( enmeth ); active_intrares2benmeths_ext_.push_back( enmeth ); }
void MinimizationNode::add_dof_deriv_2benmeth( TwoBodyEnergyCOP enmeth ) { dof_deriv_2benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfs_2benmeth( TwoBodyEnergyCOP enmeth ) { sfs_req_2benmeths_.push_back( enmeth ); }
void MinimizationNode::add_sfd_2benmeth( TwoBodyEnergyCOP enmeth ) { sfd_req_2benmeths_.push_back( enmeth ); }

bool
MinimizationNode::classify_onebody_enmeth( OneBodyEnergyCOP enmeth, Residue const & rsd, Pose const & pose, int domain_map_color )
{
	if ( domain_map_color == 0 || enmeth->method_type() == methods::cd_1b ) {
		if ( enmeth->defines_score_for_residue( rsd ) ) {
			if ( enmeth->use_extended_residue_energy_interface() ) {
				add_active_1benmeth_ext( enmeth );
			} else {
				add_active_1benmeth_std( enmeth );
			}
			if ( enmeth->defines_dof_derivatives( pose ) ) {
				add_dof_deriv_1benmeth( enmeth );
			}
			if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity( pose ) ) {
				add_sfs_1benmeth( enmeth );
			}
			if ( enmeth->requires_a_setup_for_derivatives_for_residue_opportunity( pose ) ) {
				add_sfd_1benmeth( enmeth );
			}
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
	if ( enmeth->requires_a_setup_for_scoring_for_residue_opportunity( pose ) ) {
		/// should we let energy methods setup for scoring for a residue even if they
		/// don't define an intrares energy for that residue?  Yes.  These
		add_sfs_2benmeth( enmeth );
	}
	if ( enmeth->requires_a_setup_for_derivatives_for_residue_opportunity( pose ) ) {
		/// should we let energy methods setup for scoring for a residue even if they
		/// don't define an intrares energy for that residue?  Yes.  These
		add_sfd_2benmeth( enmeth );
	}

	/// Domain map check only prevents intra-residue sfxn/deriv evaluations; but not
	/// SFS or SFD registration.
	if ( domain_map_color == 0 ) {
		if ( enmeth->defines_intrares_energy( weights ) && enmeth->defines_intrares_energy_for_residue( rsd ) ) {
			if ( enmeth->use_extended_intrares_energy_interface() ) {
				add_active_2benmeth_ext( enmeth );
			} else {
				add_active_2benmeth_std( enmeth );
			}
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
MinimizationNode::sfs_req_1benmeths_begin() const {
	return sfs_req_1benmeths_.begin();
}

MinimizationNode::OneBodyEnergiesIterator
MinimizationNode::sfs_req_1benmeths_end() const {
	return sfs_req_1benmeths_.end();
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
MinimizationNode::sfs_req_2benmeths_begin() const {
	return sfs_req_2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfs_req_2benmeths_end() const {
	return sfs_req_2benmeths_.end();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfd_req_2benmeths_begin() const {
	return sfd_req_2benmeths_.begin();
}

MinimizationNode::TwoBodyEnergiesIterator
MinimizationNode::sfd_req_2benmeths_end() const {
	return sfd_req_2benmeths_.end();
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
	copy_from( & example_edge );
}

/// @brief virtual dstor; The MinimizationEdge must free the array pool element it
/// holds before it disappears.
MinimizationEdge::~MinimizationEdge()
{
}

/// @brief copies data from MinimizationEdge const * source;
///
/// called from the copy ctor and operator= methods defined in the Graph base class
void MinimizationEdge::copy_from( parent const * source )
{
	//MinimizationEdge const * ee = static_cast< MinimizationEdge const * > ( source );
	// down_cast is *supposed* to assert the dynamic cast in debug builds; doesn't work for some reason
	MinimizationEdge const & minedge = static_cast< MinimizationEdge const & > ( * source );

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
	for ( TwoBodyEnergiesIterator iter = sfs_req_2benmeths_begin(),
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
	for ( TwoBodyEnergiesIterator iter = sfd_req_2benmeths_begin(),
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

	if ( enmeth->use_extended_residue_pair_energy_interface() ) {
		add_active_enmeth_ext( enmeth );
	} else {
		add_active_enmeth_std( enmeth );
	}

	if ( enmeth->requires_a_setup_for_scoring_for_residue_pair_opportunity( pose ) ) {
		add_sfs_enmeth( enmeth );
	}
	if ( enmeth->requires_a_setup_for_derivatives_for_residue_pair_opportunity( pose ) ) {
		add_sfd_enmeth( enmeth );
	}
	return true;
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
		return 0;
	}
}

MinimizationEdge const * MinimizationGraph::find_minimization_edge( Size n1, Size n2) const
{
	Edge const * edge( find_edge( n1, n2 ) );
	if ( edge ) {
		return static_cast< MinimizationEdge const * > ( edge );
		//return utility::down_cast< MinimizationEdge const * > ( edge );
	} else {
		return 0;
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
MinimizationGraph::MinimizationGraph( MinimizationGraph const & src )
:
	parent( ),
	minimization_edge_pool_( new boost::unordered_object_pool< MinimizationEdge > ( 256 ) )
{
	parent::operator = ( src );
}

MinimizationGraph::~MinimizationGraph() {
	delete_everything();
	delete minimization_edge_pool_; minimization_edge_pool_ = 0;
}


/// @brief assignment operator -- performs a deep copy
MinimizationGraph &
MinimizationGraph::operator = ( MinimizationGraph const & rhs )
{
	parent::operator = ( rhs );
	return *this;
}

void MinimizationGraph::delete_edge( graph::Edge * edge )
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

void MinimizationGraph::add_whole_pose_context_enmeth( EnergyMethodCOP enmeth )
{
	whole_pose_context_enmeths_.push_back( enmeth );
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
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.active_1benmeths_begin(),
			iter_end = min_node.active_1benmeths_end(); iter != iter_end; ++iter ) {
		(*iter)->eval_residue_derivatives(
			rsd, min_node.res_min_data(), pose, res_weights, atom_derivs );
	}
	/// 1b 2body intraresidue contributions
	for ( MinimizationNode::TwoBodyEnergiesIterator
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
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.active_1benmeths_std_begin(),
			iter_end = min_node.active_1benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, emap );
	}
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.active_1benmeths_ext_begin(),
			iter_end = min_node.active_1benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy_ext( rsd, min_node.res_min_data(), pose, emap );
	}
	/// 1b 2body intraresidue contributions
	for ( MinimizationNode::TwoBodyEnergiesIterator
			iter = min_node.active_intrares2benmeths_std_begin(),
			iter_end = min_node.active_intrares2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, sfxn, emap );
	}
	for ( MinimizationNode::TwoBodyEnergiesIterator
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
	for ( MinimizationEdge::TwoBodyEnergiesIterator
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

	for ( MinimizationEdge::TwoBodyEnergiesIterator
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
	for ( MinimizationEdge::TwoBodyEnergiesIterator
			iter = min_edge.active_2benmeths_std_begin(),
			iter_end = min_edge.active_2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy(
			res1, res2, pose, sfxn, emap );
	}
	for ( MinimizationEdge::TwoBodyEnergiesIterator
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
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.dof_deriv_1benmeths_begin(),
			iter_end = min_node.dof_deriv_1benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_residue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}

	for ( MinimizationNode::TwoBodyEnergiesIterator
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
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.active_1benmeths_std_begin(),
			iter_end = min_node.active_1benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy( rsd, pose, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.active_1benmeths_ext_begin(),
			iter_end = min_node.active_1benmeths_ext_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_energy_ext( rsd, min_node.res_min_data(), pose, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	/// 1b 2body intraresidue contributions
	for ( MinimizationNode::TwoBodyEnergiesIterator
			iter = min_node.active_intrares2benmeths_std_begin(),
			iter_end = min_node.active_intrares2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->eval_intrares_energy( rsd, pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_node.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( MinimizationNode::TwoBodyEnergiesIterator
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
	for ( MinimizationEdge::TwoBodyEnergiesIterator
			iter = min_edge.active_2benmeths_std_begin(),
			iter_end = min_edge.active_2benmeths_std_end();
			iter != iter_end; ++iter ) {
		(*iter)->residue_pair_energy(
			res1, res2, pose, sfxn, scratch_emap );
		emap.accumulate( scratch_emap, (*iter)->score_types(), min_edge.weight() );
		scratch_emap.zero( (*iter)->score_types() );
	}
	for ( MinimizationEdge::TwoBodyEnergiesIterator
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
	for ( MinimizationNode::OneBodyEnergiesIterator
			iter = min_node.dof_deriv_1benmeths_begin(),
			iter_end = min_node.dof_deriv_1benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_residue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}

	for ( MinimizationNode::TwoBodyEnergiesIterator
			iter = min_node.dof_deriv_2benmeths_begin(),
			iter_end = min_node.dof_deriv_2benmeths_end(); iter != iter_end; ++iter ) {
		deriv += (*iter)->eval_intraresidue_dof_derivative(
			rsd, min_node.res_min_data(), dof_id, torsion_id, pose, sfxn, weights );
	}
	return deriv * min_node.weight();
}


} //namespace scoring
} //namespace core
