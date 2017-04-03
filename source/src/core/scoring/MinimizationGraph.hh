// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/MinimizationGraph.hh
/// @brief  Minimization graph class declaration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_scoring_MinimizationGraph_hh
#define INCLUDED_core_scoring_MinimizationGraph_hh

// Unit Headers
#include <core/scoring/MinimizationGraph.fwd.hh>

// Package Headers
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/OneBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>

#ifdef WIN32
#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#endif

// Project headers
#include <utility/graph/Graph.hh>
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Basic headers
#include <basic/datacache/BasicDataCache.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <utility/vector1.hh>


namespace core {
namespace scoring {

/// Class MinimizationNode holds the ResSingleMinimizationData information for
/// a single residue in a Pose which is being minimized. The data held in this
/// node will be used in both scoring the residue one-body energies and evaluating
/// atom derivatives during minimization.
class MinimizationNode : public utility::graph::Node
{
public:
	typedef utility::graph::Node               parent;
	typedef methods::EnergyMethodCOP  EnergyMethodCOP;
	typedef methods::OneBodyEnergyCOP OneBodyEnergyCOP;
	typedef methods::TwoBodyEnergyCOP TwoBodyEnergyCOP;

	typedef utility::vector1< OneBodyEnergyCOP > OneBodyEnergies;
	typedef OneBodyEnergies::const_iterator      OneBodyEnergiesIterator;

	typedef utility::vector1< TwoBodyEnergyCOP > TwoBodyEnergies;
	typedef TwoBodyEnergies::const_iterator      TwoBodyEnergiesIterator;

	typedef utility::vector1< EnergyMethodCOP > EnergyMethods;
	typedef EnergyMethods::const_iterator       EnergyMethodsIterator;

	typedef conformation::Residue                Residue;
	typedef pose::Pose                           Pose;


public:
	MinimizationNode( utility::graph::Graph * owner, Size index );
	~MinimizationNode() override;
	void copy_from( parent const * source ) override;

	void print() const override;
	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;

	ResSingleMinimizationData const & res_min_data() const { return res_min_data_; }
	ResSingleMinimizationData & res_min_data() { return res_min_data_; }

	bool add_onebody_enmeth( OneBodyEnergyCOP enmeth, Residue const & rsd, Pose const & pose, int domain_map_color );
	bool add_twobody_enmeth(
		TwoBodyEnergyCOP enmeth,
		Residue const & rsd,
		Pose const & pose,
		EnergyMap const & weights,
		int domain_map_color
	);

	void setup_for_minimizing(
		Residue const & rsd,
		Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map,
		basic::datacache::BasicDataCache & res_data_cache
	);
	void setup_for_scoring( Residue const & rsd, basic::datacache::BasicDataCache & residue_data_cache, Pose const & pose, ScoreFunction const & sfxn );
	void setup_for_derivatives( Residue const & rsd, basic::datacache::BasicDataCache & residue_data_cache, pose::Pose const & pose, ScoreFunction const & sfxn );
	void update_active_enmeths_for_residue(
		Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap const & weights,
		int domain_map_color
	);

	Real weight() const { return weight_; }
	void weight( Real setting ) { weight_ = setting; }
	Real dweight() const { return dweight_; }
	void dweight( Real setting ) { dweight_ = setting; }

protected:
	MinimizationGraph const * get_mingraph_owner() const;
	MinimizationGraph * get_mingraph_owner();

private:
	void add_active_1benmeth_std( OneBodyEnergyCOP enmeth );
	void add_active_1benmeth_ext( OneBodyEnergyCOP enmeth );
	void add_dof_deriv_1benmeth( OneBodyEnergyCOP enmeth );
	void add_sfs_dm_1benmeth( OneBodyEnergyCOP enmeth );
	void add_sfd_1benmeth( OneBodyEnergyCOP enmeth );

	void add_active_2benmeth_std( TwoBodyEnergyCOP enmeth );
	void add_active_2benmeth_ext( TwoBodyEnergyCOP enmeth );
	void add_dof_deriv_2benmeth( TwoBodyEnergyCOP enmeth );
	void add_sfs_dm_2benmeth( TwoBodyEnergyCOP enmeth );
	void add_sfd_2benmeth( TwoBodyEnergyCOP enmeth );

public:
	/// @brief This method is not meant for general use; it's only to be called
	/// by the MinimizationNode and the MinimizationGraph.
	void add_sfs_drs_enmeth( EnergyMethodCOP enmeth );

private:
	bool classify_onebody_enmeth( OneBodyEnergyCOP enmeth, Residue const & rsd, Pose const & pose,int domain_map_color );
	bool classify_twobody_enmeth(
		TwoBodyEnergyCOP enmeth,
		Residue const & rsd,
		Pose const & pose,
		EnergyMap const & weights,
		int domain_map_color
	);

public:
	OneBodyEnergiesIterator active_1benmeths_begin() const;
	OneBodyEnergiesIterator active_1benmeths_end() const;
	OneBodyEnergiesIterator active_1benmeths_std_begin() const;
	OneBodyEnergiesIterator active_1benmeths_std_end() const;
	OneBodyEnergiesIterator active_1benmeths_ext_begin() const;
	OneBodyEnergiesIterator active_1benmeths_ext_end() const;
	OneBodyEnergiesIterator dof_deriv_1benmeths_begin() const;
	OneBodyEnergiesIterator dof_deriv_1benmeths_end() const;
	OneBodyEnergiesIterator sfs_dm_req_1benmeths_begin() const;
	OneBodyEnergiesIterator sfs_dm_req_1benmeths_end() const;
	OneBodyEnergiesIterator sfd_req_1benmeths_begin() const;
	OneBodyEnergiesIterator sfd_req_1benmeths_end() const;

	TwoBodyEnergiesIterator active_intrares2benmeths_begin() const;
	TwoBodyEnergiesIterator active_intrares2benmeths_end() const;
	TwoBodyEnergiesIterator active_intrares2benmeths_std_begin() const;
	TwoBodyEnergiesIterator active_intrares2benmeths_std_end() const;
	TwoBodyEnergiesIterator active_intrares2benmeths_ext_begin() const;
	TwoBodyEnergiesIterator active_intrares2benmeths_ext_end() const;
	TwoBodyEnergiesIterator dof_deriv_2benmeths_begin() const;
	TwoBodyEnergiesIterator dof_deriv_2benmeths_end() const;
	TwoBodyEnergiesIterator sfs_dm_req_2benmeths_begin() const;
	TwoBodyEnergiesIterator sfs_dm_req_2benmeths_end() const;
	TwoBodyEnergiesIterator sfd_req_2benmeths_begin() const;
	TwoBodyEnergiesIterator sfd_req_2benmeths_end() const;

	EnergyMethodsIterator sfs_drs_req_enmeths_begin() const;
	EnergyMethodsIterator sfs_drs_req_enmeths_end() const;

private:

	ResSingleMinimizationData res_min_data_;

	// The list of all acive and inactive 1body energy methods; inactive methods are not used during score evaluation,
	// but they might become active if the residue-type changes at this position.
	OneBodyEnergies onebody_enmeths_;
	// The list of all active 1body energy methods
	OneBodyEnergies active_1benmeths_;
	// one-body energy methods that evaluate a residue energy for this
	// residue using the standard residue_energy() interface
	OneBodyEnergies active_1benmeths_std_;
	// one-body energy methods that evaluate a residue energy for this
	// residue using the residue_energy_ext() interface
	OneBodyEnergies active_1benmeths_ext_;
	// one-body energy methods that define DOF derivatives
	OneBodyEnergies dof_deriv_1benmeths_;
	// one-body energy methods that require a setup-for-scoring-during-minimization (sfs_dm) opportunity
	OneBodyEnergies sfs_dm_req_1benmeths_;
	// one-body energy methods that require a setup-for-derivatives (sfd) opportunity
	OneBodyEnergies sfd_req_1benmeths_;


	// The list of all active and inactive 2body energy methods; inactive methods are not used during score evaluation,
	// but they might become active if the residue-type changes at this position.
	TwoBodyEnergies twobody_enmeths_;

	// The list of all active 2body energy methods that define an intraresidue energy
	TwoBodyEnergies active_intrares2benmeths_;
	// two-body energy methods that define an intra-residue energy for this
	// residue using the eval_intrares_energy() interface
	TwoBodyEnergies active_intrares2benmeths_std_;
	// two-body energy methods that define an intra-residue energy for this
	// residue using the eval_intrares_energy_ext() interface
	TwoBodyEnergies active_intrares2benmeths_ext_;
	// two-body energy methods that define DOF derivatives
	TwoBodyEnergies dof_deriv_2benmeths_;
	// two-body energy methods that require a setup-for-scoring-during-minimization (sfs_dm) opportunity
	TwoBodyEnergies sfs_dm_req_2benmeths_;
	// two-body energy methods that require a setup-for-derivatives (sfs) opportunity
	TwoBodyEnergies sfd_req_2benmeths_;

	// two-body energy methods that require a setup-for-scoring-during-regular-scoring (sfs_drs) opportunity
	EnergyMethods sfs_drs_req_enmeths_;

	Real weight_,dweight_;
};

/// Class MinimizationEdge holds ResPairMinimizationData for a certain pair
/// of interacting residues; this data might be a neighborlist for this residue
/// pair, for example.  The data held in this edge will be used in both scoring
/// the residue-pair energies and evaluating atom derivatives during minimization.
class MinimizationEdge : public utility::graph::Edge
{
public:
	typedef utility::graph::Edge parent;
	typedef methods::TwoBodyEnergyCOP            TwoBodyEnergyCOP;
	typedef utility::vector1< TwoBodyEnergyCOP > TwoBodyEnergies;
	typedef TwoBodyEnergies::const_iterator      TwoBodyEnergiesIterator;
	typedef conformation::Residue                Residue;
	typedef pose::Pose                           Pose;

public:
	MinimizationEdge( MinimizationGraph * owner, Size n1, Size n2 );
	MinimizationEdge( MinimizationGraph * owner, MinimizationEdge const & example_edge );
	~MinimizationEdge() override;

	/// @brief Copy the data held on the example edge, source.
	/// The source edge must be castable to class MinimizationEdge.
	void copy_from( parent const * source ) override;

	ResPairMinimizationData const & res_pair_min_data() const { return res_pair_min_data_; }
	ResPairMinimizationData & res_pair_min_data() { return res_pair_min_data_; }

	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;

	/// @brief Include a particular energy method as part of this edge.  It may not show up
	/// in the active energy methods should this energy method not define an energy for the
	/// residues.
	bool add_twobody_enmeth( TwoBodyEnergyCOP enmeth, Residue const & rsd1, Residue const & rsd2, Pose const & pose, bool residues_mwrt_eachother );

	/// @brief It may be possible to determine that an edge does not need to belong to the
	/// minimization graph if there are no active two-body energy methods; this is a convenience
	/// function that answers quickly if active_2benmths_.begin() == active_2benmeths_.end().
	bool any_active_enmeths() const {
		return ! active_2benmeths_.empty();
	}

	void setup_for_minimizing(
		Residue const & rsd1,
		Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		kinematics::MinimizerMapBase const & min_map
	);

	/// @brief Initialize the active energy methods for score function evaluation
	void setup_for_scoring( Residue const & rsd1, Residue const & rsd2, Pose const & pose, ScoreFunction const & sfxn );

	/// @brief Initialize the active energy methods for derivative evaluation
	void setup_for_derivatives( Residue const & rsd1, Residue const & rsd2, Pose const & pose, ScoreFunction const & sfxn );

	/// @brief Setup the active and inactive energy methods
	void reinitialize_active_energy_methods(
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		bool res_moving_wrt_eachother
	);

private:
	void add_active_enmeth_std( TwoBodyEnergyCOP enmeth );
	void add_active_enmeth_ext( TwoBodyEnergyCOP enmeth );
	void add_sfs_enmeth( TwoBodyEnergyCOP enmeth );
	void add_sfd_enmeth( TwoBodyEnergyCOP enmeth );

	bool classify_twobody_enmeth(
		TwoBodyEnergyCOP enmeth,
		Residue const & rsd1,
		Residue const & rsd2,
		Pose const & pose,
		bool res_moving_wrt_eachother
	);

public:
	TwoBodyEnergiesIterator active_2benmeths_begin() const;
	TwoBodyEnergiesIterator active_2benmeths_end() const;
	TwoBodyEnergiesIterator active_2benmeths_std_begin() const;
	TwoBodyEnergiesIterator active_2benmeths_std_end() const;
	TwoBodyEnergiesIterator active_2benmeths_ext_begin() const;
	TwoBodyEnergiesIterator active_2benmeths_ext_end() const;
	TwoBodyEnergiesIterator sfs_req_2benmeths_begin() const;
	TwoBodyEnergiesIterator sfs_req_2benmeths_end() const;
	TwoBodyEnergiesIterator sfd_req_2benmeths_begin() const;
	TwoBodyEnergiesIterator sfd_req_2benmeths_end() const;


	/// @brief The minimization graph will allow the storage of edge weights, should that
	/// prove useful for any application (e.g. symmetric minimization)
	Real weight() const { return weight_; }
	/// @brief Set the weight for an edge
	void weight( Real setting ) { weight_ = setting; }

	// weights to be used during minization
	Real dweight() const { return dweight_; }
	void dweight( Real setting ) { dweight_ = setting; }

protected:
	/// Downcasts

	inline
	MinimizationGraph const *
	get_minimization_owner() const;

	inline
	MinimizationGraph *
	get_minimization_owner();

	inline
	MinimizationNode *
	get_minimization_node( Size index ) {
		return static_cast< MinimizationNode * > ( get_node( index ) );
	}

	inline
	MinimizationNode const *
	get_minimization_node( Size index ) const {
		return static_cast< MinimizationNode const * > ( get_node( index ) );
	}

private:

	ResPairMinimizationData res_pair_min_data_;

	// The list of all active and inactive two-body energy methods for this position
	TwoBodyEnergies twobody_enmeths_;
	// The list of all two-body energy methods that define a residue pair energy
	// for this residue pair
	TwoBodyEnergies active_2benmeths_;
	// two-body energy methods that define a residue pair energy for this
	// residue using the residue_pair_energy() interface
	TwoBodyEnergies active_2benmeths_std_;
	// two-body energy methods that define an intra-residue energy for this
	// residue using the residue_pair_energy_ext() interface
	TwoBodyEnergies active_2benmeths_ext_;
	// two-body energy methods that require a setup-for-scoring opportunity for
	// this residue pair
	TwoBodyEnergies sfs_req_2benmeths_;
	// two-body energy methods that require a setup-for-derivatives opportunity for
	// this residue pair
	TwoBodyEnergies sfd_req_2benmeths_;

	Real weight_, dweight_;

};

/// @brief Class to hold all the minimization-specific data that's required
/// to efficiently evaluate the score function and its derivatives on a structure
/// of fixed sequence and chemical identity.
class MinimizationGraph : public utility::graph::Graph
{
public:
	typedef utility::graph::Graph                 parent;
	typedef methods::EnergyMethodCOP     EnergyMethodCOP;
	typedef std::list< EnergyMethodCOP > Energies;
	typedef Energies::const_iterator     EnergiesIterator;

public:
	~MinimizationGraph() override;

	MinimizationGraph( Size num_nodes );
	MinimizationGraph();
	MinimizationGraph( MinimizationGraph const & src );

	MinimizationGraph & operator = ( MinimizationGraph const & rhs );

	inline
	MinimizationNode const *
	get_minimization_node( Size index ) const
	{
		return static_cast< MinimizationNode const * > ( get_node( index ));
	}

	inline
	MinimizationNode *
	get_minimization_node( Size index )
	{
		return static_cast< MinimizationNode * > ( get_node( index ));
	}


	MinimizationEdge * find_minimization_edge( Size n1, Size n2);
	MinimizationEdge const * find_minimization_edge( Size n1, Size n2) const;

	void delete_edge( utility::graph::Edge * edge ) override;

	void add_whole_pose_context_enmeth( EnergyMethodCOP enmeth, core::pose::Pose const & pose );
	EnergiesIterator whole_pose_context_enmeths_begin() const;
	EnergiesIterator whole_pose_context_enmeths_end() const;

	void set_fixed_energies( EnergyMap const & );
	EnergyMap const & fixed_energies() const;

protected:
	Size count_static_memory() const override;
	Size count_dynamic_memory() const override;

	utility::graph::Node * create_new_node( Size index ) override;
	utility::graph::Edge * create_new_edge( Size index1, Size index2 ) override;
	utility::graph::Edge * create_new_edge( utility::graph::Edge const * example_edge ) override;

private:
	Energies whole_pose_context_enmeths_;
	EnergyMap fixed_energies_;
	boost::unordered_object_pool< MinimizationEdge > * minimization_edge_pool_;

};


MinimizationGraph const *
MinimizationEdge::get_minimization_owner() const {
	return static_cast< MinimizationGraph const * > (get_owner());
}


MinimizationGraph *
MinimizationEdge::get_minimization_owner() {
	return static_cast< MinimizationGraph * > (get_owner());
}

//// Non member functions for score and derivative evaluation when using a minimization graph

/*void
reinitialize_minedge_for_respair(
MinimizationEdge & min_edge,
conformation::Residue const & rsd1,
conformation::Residue const & rsd2,
ResSingleMinimizationData const & res1_ressingle_min_data,
ResSingleMinimizationData const & res2_ressingle_min_data,
pose::Pose & p,
ScoreFunction const & sfxn,
kinematics::MinimizerMapBase const & minmap
);*/

/*void
setup_for_scoring_for_minnode(
MinimizationNode & min_node,
conformation::Residue const & rsd,
pose::Pose const & pose
);*/

/*void
setup_for_scoring_for_minedge(
MinimizationEdge & min_edge,
conformation::Residue const & rsd1,
conformation::Residue const & rsd2,
pose::Pose const & pose,
ResSingleMinimizationData const & res1_ressingle_min_data,
ResSingleMinimizationData const & res2_ressingle_min_data
);*/

/*void
setup_for_derivatives_for_minnode(
MinimizationNode & min_node,
conformation::Residue const & rsd,
pose::Pose const & pose
);*/

/*void
setup_for_derivatives_for_minedge(
MinimizationEdge & min_edge,
conformation::Residue const & rsd1,
conformation::Residue const & rsd2,
pose::Pose const & pose,
ResSingleMinimizationData const & res1_min_data,
ResSingleMinimizationData const & res2_min_data
);*/

/// @brief Evaluate the derivatives for all atoms on the input residue
/// for the terms that apply to this residue (which are stored on the input
/// minimization node).
void
eval_atom_derivatives_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	EnergyMap const & res_weights,
	utility::vector1< DerivVectorPair > & atom_derivs
);

/// @brief Deprecated
/*void
eval_atom_derivative_for_minnode(
MinimizationNode const & min_node,
Size atom_index,
conformation::Residue const & rsd,
pose::Pose const & pose,
kinematics::DomainMap const & domain_map,
ScoreFunction const & sfxn,
EnergyMap const & res_weights,
Vector & F1, // accumulated into
Vector & F2  // accumulated into
);*/

void
eval_res_onebody_energies_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap // accumulated into
);

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
);

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
);

/// @brief Deprecated
/*void
eval_atom_deriv_for_minedge(
MinimizationEdge const & min_edge,
Size atom_index,
conformation::Residue const & res1,
conformation::Residue const & res2,
ResSingleMinimizationData const & res1_min_data,
ResSingleMinimizationData const & res2_min_data,
pose::Pose const & pose,
kinematics::DomainMap const & domain_map,
ScoreFunction const & sfxn,
EnergyMap const & respair_weights,
Vector & F1, // accumulated into
Vector & F2  // accumulated into
);*/

void
eval_res_pair_energy_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap
);

Real
eval_dof_deriv_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
);

/*void
eval_weighted_atom_derivative_for_minnode(
MinimizationNode const & min_node,
Size atom_index,
conformation::Residue const & rsd,
pose::Pose const & pose,
kinematics::DomainMap const & domain_map,
ScoreFunction const & sfxn,
EnergyMap const & res_weights,
Vector & F1, // accumulated into
Vector & F2  // accumulated into
);*/

void
eval_weighted_res_onebody_energies_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap, // accumulated into
	EnergyMap & scratch_emap // should be zeros coming in, left zeroed at the end;
);


/*void
eval_weighted_atom_deriv_for_minedge(
MinimizationEdge const & min_edge,
Size atom_index,
conformation::Residue const & res1,
conformation::Residue const & res2,
ResSingleMinimizationData const & res1_min_data,
ResSingleMinimizationData const & res2_min_data,
pose::Pose const & pose,
kinematics::DomainMap const & domain_map,
ScoreFunction const & sfxn,
EnergyMap const & respair_weights,
Vector & F1, // accumulated into
Vector & F2  // accumulated into
);*/

void
eval_weighted_res_pair_energy_for_minedge(
	MinimizationEdge const & min_edge,
	conformation::Residue const & res1,
	conformation::Residue const & res2,
	pose::Pose const & pose,
	ScoreFunction const & sfxn,
	EnergyMap & emap,
	EnergyMap & scratch_emap // should be zeros coming in, left zeroed at the end;
);

Real
eval_weighted_dof_deriv_for_minnode(
	MinimizationNode const & min_node,
	conformation::Residue const & rsd,
	pose::Pose const & pose,
	id::DOF_ID const & dof_id,
	id::TorsionID const & torsion_id,
	ScoreFunction const & sfxn,
	EnergyMap const & weights
);


} //namespace scoring
} //namespace core

#endif

