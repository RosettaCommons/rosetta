// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.hh
/// @brief  Score function class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_ScoreFunction_hh
#define INCLUDED_core_scoring_ScoreFunction_hh

// Unit header
#include <core/scoring/ScoreFunction.fwd.hh>

// C/C++ headers
#include <map>
#include <string>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/ContextDependentOneBodyEnergy.fwd.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.fwd.hh>

// Basic headers
#include <basic/datacache/BasicDataCache.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.fwd.hh>


#ifdef WIN32 //VC++ needs full class declaration
#include <core/scoring/methods/EnergyMethod.hh> // WIN32 INCLUDE
#include <core/scoring/methods/EnergyMethodOptions.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/TwoBodyEnergy.hh> // WIN32 INCLUDE
#include <core/scoring/methods/WholeStructureEnergy.hh> // WIN32 INCLUDE
#endif

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

/// @brief This object defines a ScoreFunction, it contains methods for
/// calculating the various scoring components (called ScoreType's) used in
/// Rosetta. It also contains weights that are applied to each of those
/// components. Only scoring components with non-zero weights are calculated.
class ScoreFunction : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ScoreFunction >
{
	/////////////////////////////////////////////////////////////////////////////
	// typedefs
	/////////////////////////////////////////////////////////////////////////////

public:
	typedef utility::vector1< methods::EnergyMethodCOP > AllMethods;

	typedef utility::vector1< methods::ContextIndependentTwoBodyEnergyOP >   CI_2B_Methods;
	typedef utility::vector1< methods::ContextDependentTwoBodyEnergyOP >     CD_2B_Methods;
	typedef utility::vector1< methods::ContextIndependentOneBodyEnergyOP >   CI_1B_Methods;
	typedef utility::vector1< methods::ContextDependentOneBodyEnergyOP >     CD_1B_Methods;
	typedef utility::vector1< methods::ContextIndependentLRTwoBodyEnergyOP > CI_LR_2B_Methods;
	typedef utility::vector1< methods::ContextDependentLRTwoBodyEnergyOP >   CD_LR_2B_Methods;
	typedef utility::vector1< methods::LongRangeTwoBodyEnergyCOP >           LR_2B_Methods;
	typedef utility::vector1< methods::TwoBodyEnergyOP >                     TWO_B_Methods;
	typedef utility::vector1< methods::WholeStructureEnergyOP >              WS_Methods;

	typedef utility::vector1< methods::WholeStructureEnergyOP >::const_iterator WS_MethodIterator;
	typedef utility::vector1< methods::LongRangeTwoBodyEnergyCOP >::const_iterator LR_2B_MethodIterator;
	typedef utility::vector1< methods::TwoBodyEnergyOP >::const_iterator TWO_B_MethodIterator;
	typedef utility::vector1< methods::EnergyMethodCOP >::const_iterator AllMethodsIterator;
	typedef utility::vector1< methods::ContextIndependentLRTwoBodyEnergyOP >::const_iterator CI_LR_2B_MethodIterator;
	typedef utility::vector1< methods::ContextDependentLRTwoBodyEnergyOP >::const_iterator   CD_LR_2B_MethodIterator;

public:

	/// @brief default ctor that initializes its EnergyMethodOptions object from the global option collection
	ScoreFunction();

	/// @brief ctor that initializes its EnergyMethodOptions object from a (possibly local) option collection
	ScoreFunction( utility::options::OptionCollection const & options );

	~ScoreFunction() override;

	static
	void
	list_options_read( utility::options::OptionKeyList & option_list );

	/// self pointers
	inline ScoreFunctionCOP get_self_ptr() const { return shared_from_this(); }
	inline ScoreFunctionOP get_self_ptr() { return shared_from_this(); }
	// inline ScoreFunctionCAP get_self_weak_ptr() const { return ScoreFunctionCAP( shared_from_this() ); }
	// inline ScoreFunctionAP get_self_weak_ptr() { return ScoreFunctionAP( shared_from_this() ); }

private:
	/// @brief The ScoreFunction copy constructor is explicitly private
	/// as using it to make a copy is just too attractive, but discards subclass information.
	/// Use ScoreFunction::clone() instead.
	ScoreFunction( ScoreFunction const & );

	/// @brief The ScoreFunction assignment operator is explicitly private
	/// as using it discards subclass information.
	/// Rework your algorithm to use ScoreFunctionOP's instead.
	ScoreFunction &
	operator=( ScoreFunction const & );

public:

	/// @brief NOT FOR GENERAL USE
	/// Copy the information about src into the current score function.
	/// There are deep interactions with subclasses,
	/// (the subclass information doesn't necessarily get copied)
	/// so this is primarily for advanced scorefunction manipulation.
	/// Normal usage should just use clone() and replace the OP.
	virtual void
	assign( ScoreFunction const & src);

	/// @brief Create a copy of the scorefunction
	/// Virtual to keep subclass information.
	virtual ScoreFunctionOP clone() const;

	/// @brief If you *want* to discard subclass information, the following function is availible
	ScoreFunctionOP clone_as_base_class() const;

	/// @brief Resets the ScoreFunction to default values, reading from the global options collection
	void reset();

	/// @brief Resets the ScoreFunction to default values, reading from a (possibly local) options collection
	void reset( utility::options::OptionCollection const & options );

	/// @brief Randomly perturbs non-zero score function weights
	void perturb_weights();

	/// @brief Serializes the non-zero score function term weights
	/// Format: { term : weight, ... }
	std::string serialize_weights() const;

	/// @brief Initializes this ScoreFunction from the given  <filename>
	void
	add_weights_from_file( std::string const & filename );

	/// @brief Initializes this ScoreFunction from the given  <filename>
	/// no lookup in database directory
	void
	_add_weights_from_file( std::string const & filename, bool patch=false );

	void
	_add_weights_from_stream( std::istream & data, bool patch=false, std::string const & filename="");

	/// @brief Resets everything before reading the  <filename>
	void
	initialize_from_file( std::string const & filename );

	/// @brief Applies a patch from the given  <filename>
	void
	apply_patch_from_file( std::string const & filename );

	/// @brief Given a  <filename> (represented by a std::string), set the
	/// e_table for this ScoreFunction.
	void
	set_etable( std::string const & etable_name );


	void
	set_method_weights( ScoreType const & t, utility::vector1< Real > const & wts );


	/// @brief Returns the EnergyMethodOptions object contained in this
	/// ScoreFunction (const access)
	methods::EnergyMethodOptions const &
	energy_method_options() const
	{
		return *energy_method_options_;
	}


	// NO LONGER IN USE! Users trying to take advantage of this saw some non-intuitive behavior.
	// @brief Returns the EnergyMethodOptions object contained in this ScoreFunction.
	// methods::EnergyMethodOptions &
	// energy_method_options()
	// {
	//  return *energy_method_options_;
	// }

	/// @brief Sets the EnergyMethodOptions object contained in this ScoreFunction.
	/// with appropriate update of all the energy methods.
	void
	set_energy_method_options(
		methods::EnergyMethodOptions const &
		energy_method_options_in
	);

	void
	reset_energy_methods();

	/////////////////////////////////////////////////////////////////////////////
	// score
	/////////////////////////////////////////////////////////////////////////////
	/// @brief Scores the given  <pose>  using this ScoreFunction. Alters the
	/// Energies object within  <pose>, but does not alter this ScoreFunction
	virtual Real
	operator ()( pose::Pose & pose ) const;

	/// @brief Scores the given  <pose>  using this ScoreFunction. Alters the
	/// Energies object within  <pose>, but does not alter this ScoreFunction
	/// @note: Synonym for () operator. Makes code look a little nicer. Doesn't
	/// do anything but call () operator.
	virtual Real
	score( pose::Pose & pose ) const;

	/// @brief Score the structure and store the component energies in the EnergyGraph
	/// without requiring a second evaluation of the short-ranged two body energies.
	/// Note: pose copy operations do not copy the EnergyGraph, so cloning or copying a pose
	/// that has had it's components scored will not copy over the component energies.
	//virtual Real
	//score_components( pose::Pose & pose ) const;

	/// @brief return an object to describe abstractly the methods contained in this
	/// ScoreFunction so that class Energies can ensure that the ScoreFunction is
	/// properly evaluated (ie, no obsolete cashed data is used )
	ScoreFunctionInfoOP
	info() const;

	/// Returns the largest atomic interaction cutoff required by the
	/// EnergyMethods
	Distance
	max_atomic_interaction_cutoff() const;

	/// find which context graphs the energy methods require
	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;


	/// @brief  Determine if two residues could have non-zero interaction energy after possibly changing rotamers/chi
	/// @brief  angles; uses any_lr_residue_pair_energy, combines with max_atomic_interaction_cutoff.
	/// @note  Includes long-range energies like constraints, gb...

	/// @brief Returns true if the  <pose>  positions  <pos1>  and  <pos2>
	/// are neighbors
	bool
	are_they_neighbors( pose::Pose const & pose, Size const pos1, Size const pos2 ) const;

	/////////////////////////////////////////////////////////////////////////////
	// access and/or set the weights
	/////////////////////////////////////////////////////////////////////////////


	/// @brief Returns the weight of the ScoreType  <t>
	///
	/// example(s):
	///     scorefxn[fa_sol]
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.set_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	Real
	operator []( ScoreType const & t ) const
	{
		return weights_[ t ];
	}

	/// @brief Returns the score of the ScoreType  <t>
	Real score_by_scoretype(
		pose::Pose & pose,
		ScoreType const t,
		bool const weighted = true
	) const;

	/// @brief Returns an EnergyMap of the current set of weights
	///
	/// example(s):
	///     scorefxn.weights()
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.set_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	EnergyMap const &
	weights() const
	{
		return weights_;
	}

	/// @brief Returns true if the ScoreType  <t>  has a weight of zero,
	///
	/// example(s):
	///     scorefxn.has_zero_weight(fa_sol)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.has_nonzero_weight
	///     ScoreFunction.set_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	bool
	has_zero_weight( ScoreType const & t ) const
	{
		return ( weights_[t] == Real( 0.0 ) );
	}

	/// @brief Returns true if the ScoreType  <t>  has a non-zero weight
	///
	/// example(s):
	///     scorefxn.has_nonzero_weight(fa_sol)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.has_zero_weight
	///     ScoreFunction.set_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	bool
	has_nonzero_weight( ScoreType const & t ) const
	{
		return ( weights_[t] != Real( 0.0 ) );
	}

	/// @brief Returns a list of the ScoreTypes which are non-zero with
	/// their current weights
	///
	/// example(s):
	///     scorefxn.get_nonzero_weighted_scoretypes()
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.set_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	ScoreTypes
	get_nonzero_weighted_scoretypes() const {
		ScoreTypes scoretypes;
		for ( int i=1; i<= n_score_types; ++i ) {
			if ( has_nonzero_weight( ScoreType(i) ) ) {
				scoretypes.push_back( ScoreType(i) );
			}
		}
		return scoretypes;
	}

	/// @brief Sets the weight for ScoreType  <t>  to  <setting>
	///
	/// example(s):
	///     scorefxn.set_weight(fa_sol,.5)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.weights
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	void
	set_weight( ScoreType const & t, Real const & setting );

	/// @brief Sets the weight for ScoreType  <t>  to  <setting> if weight is originally zero
	///
	/// example(s):
	///     scorefxn.set_weight_if_zero(fa_sol,.5)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.set_weight
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	void
	set_weight_if_zero( ScoreType const & t, Real const & setting);


	/// @brief Increments the weight for ScoreType  <t>  by  <setting>
	///
	/// example(s):
	///     scorefxn.add_to_weight(fa_sol,.5)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.get_weight
	///     ScoreFunction.set_weight
	///     ScoreType
	///     name_from_score_type
	///     score_type_from_name
	void
	add_to_weight( ScoreType const & t, Real const & increment );

	/// @brief Returns the weight for ScoreType  <t>
	///
	/// examples(s):
	///     scorefxn.get_weight(fa_sol)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.set_weight
	///     ScoreType
	///     create_score_function
	///     name_from_score_type
	///     score_type_from_name
	Real
	get_weight( ScoreType const & t ) const;


	/// @brief Adds a scoring method that is not necessarily included in
	/// the core library
	void
	add_extra_method(
		ScoreType const & new_type,
		Real const new_weight,
		methods::EnergyMethod const & new_method
	);

	/// @brief Adds a scoring method that is not necessarily included in
	/// the core library
	void
	add_extra_method(
		std::map< ScoreType, Real > const & new_weights,
		methods::EnergyMethod const & new_method
	);

	/// @brief Initializes a MinimizationGraph and caches it in
	/// Energies object of  <pose>
	/// @note: for use during minimization
	virtual
	void
	setup_for_minimizing(
		pose::Pose & pose,
		kinematics::MinimizerMapBase const & min_map
	) const;

	/// @brief Initialize a single node of a MinimizationGraph with the one-body and two-body
	/// energy methods that are held within this ScoreFunction object.
	virtual
	void
	setup_for_minimizing_for_node(
		MinimizationNode & min_node,
		conformation::Residue const & rsd,
		basic::datacache::BasicDataCache & res_data_cache,
		kinematics::MinimizerMapBase const & min_map,
		pose::Pose & pose, // context
		bool accumulate_fixed_energies,
		EnergyMap & fixed_energies
	) const;

	void
	reinitialize_minnode_for_residue(
		MinimizationNode & min_node,
		conformation::Residue const & rsd,
		basic::datacache::BasicDataCache & res_data_cache,
		kinematics::MinimizerMapBase const & min_map,
		pose::Pose & pose // context
	) const;

	/// @brief Initialize a single MinimizationEdge for a particular part of residues, storing
	/// sr2b energy method pointers on the edge for those sr2b energy methods in this ScoreFunction
	void
	setup_for_minimizing_sr2b_enmeths_for_minedge(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		MinimizationEdge & min_edge,
		kinematics::MinimizerMapBase const & min_map,
		pose::Pose & pose,
		bool const res_moving_wrt_eachother,
		bool accumulate_fixed_energies,
		EnergyEdge const * energy_edge,
		EnergyMap & fixed_energies,
		Real const edge_weight = 1.0  //fpd this is unused so I am not adding 'edge_dweight' parameter
	) const;

	/// @brief Initialize an edge in the MinimizationGraph with a particular long-range two body
	void
	setup_for_lr2benmeth_minimization_for_respair(
		conformation::Residue const & res1,
		conformation::Residue const & res2,
		methods::LongRangeTwoBodyEnergyCOP lr2benergy,
		MinimizationGraph & g,
		kinematics::MinimizerMapBase const & min_map,
		pose::Pose & p,
		bool const res_moving_wrt_eachother,
		bool accumulate_fixed_energies,
		ResidueNeighborConstIteratorOP rni,
		EnergyMap & fixed_energies,
		Real const edge_weight = 1.0,
		Real const edge_dweight = 1.0
	) const;

	/// @brief For external scorers: Let the energy methods prepare for
	/// evaluating their scores on a particular structure
	/// @note: invoked during scoring.
	void
	setup_for_scoring(
		pose::Pose & pose
	) const;


	/// @brief Lets the scoring functions cache anything they need to calculate
	/// energies in a packing step (rotamer_trials or pack_rotamers)
	/// @note: the Etable caches tries for each of the residues, the
	/// hydrogen bond function caches backbone/backbone hydrogen bonds
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const & residues_repacking, utility::vector1< bool > const & residues_designing ) const;

	/// @brief Lets the scoring functions cache anything they need to rapidly
	/// calculate rotamer pair energies used in packing (like a trie, e.g.)
	void
	prepare_rotamers_for_packing( pose::Pose const & pose, conformation::RotamerSetBase & set ) const;

	/// @brief If inside packing, the pose changes conformation, inform the
	/// scoring functions that any data they have cached in the Energies object
	/// is out of date. In particular, this is to update the trie(s) during
	/// rotamer trials.
	void
	update_residue_for_packing( pose::Pose & pose, Size resid ) const;


	//void
	//setup_for_scoring( pose::Pose & pose ) const;

	/// @brief
	virtual
	void
	setup_for_derivatives( pose::Pose & pose ) const;

	/// @brief
	virtual
	void
	finalize_after_derivatives( pose::Pose & pose ) const;


	/// @brief Compute the score for subset of residues
	core::Real
	get_sub_score(
		pose::Pose const & pose,
		utility::vector1< bool > const & residue_mask ) const;

	/// @brief Compute the score for subset of residues
	void
	get_sub_score(
		pose::Pose const & pose,
		utility::vector1< bool > const & residue_mask,
		EnergyMap & emap) const;

	/// @brief Compute the score for subset of residues
	core::Real
	get_sub_score(
		pose::Pose & pose,
		utility::vector1< bool > const & residue_mask ) const;

	/// @brief Compute the score for subset of residues
	void
	get_sub_score(
		pose::Pose & pose,
		utility::vector1< bool > const & residue_mask,
		EnergyMap & emap) const;


	/// @brief Compute the score for subset of residues
	core::Real
	get_sub_score_exclude_res(
		pose::Pose const & pose,
		utility::vector1< core::Size > const & exclude_list ) const;

	/// @brief Compute the score for subset of residues
	void
	get_sub_score_exclude_res(
		pose::Pose const & pose,
		utility::vector1< core::Size > const & exclude_list,
		EnergyMap & emap) const;

	/// @brief Compute the score for subset of residues
	core::Real
	get_sub_score_exclude_res(
		pose::Pose & pose,
		utility::vector1< core::Size > const & exclude_list ) const;

	/// @brief Compute the score for subset of residues
	void
	get_sub_score_exclude_res(
		pose::Pose & pose,
		utility::vector1< core::Size > const & exclude_list,
		EnergyMap & emap) const;


	/// @brief
	virtual
	void
	eval_twobody_neighbor_energies( pose::Pose & pose ) const;

	/// @brief
	virtual
	void
	eval_long_range_twobody_energies( pose::Pose & pose ) const;

	AllMethodsIterator
	all_energies_begin() const;

	AllMethodsIterator
	all_energies_end() const;

	LR_2B_MethodIterator
	long_range_energies_begin() const;

	LR_2B_MethodIterator
	long_range_energies_end() const;

	TWO_B_MethodIterator
	ci_2b_intrares_begin() const;

	TWO_B_MethodIterator
	ci_2b_intrares_end() const;

	TWO_B_MethodIterator
	cd_2b_intrares_begin() const;

	TWO_B_MethodIterator
	cd_2b_intrares_end() const;

	CI_2B_Methods::const_iterator
	ci_2b_begin() const;

	CI_2B_Methods::const_iterator
	ci_2b_end() const;

	CD_2B_Methods::const_iterator
	cd_2b_begin() const;

	CD_2B_Methods::const_iterator
	cd_2b_end() const;

	CI_LR_2B_MethodIterator
	ci_lr_2b_methods_begin() const;

	CI_LR_2B_MethodIterator
	ci_lr_2b_methods_end() const;

	CD_LR_2B_MethodIterator
	cd_lr_2b_methods_begin() const;

	CD_LR_2B_MethodIterator
	cd_lr_2b_methods_end() const;

	WS_MethodIterator
	ws_methods_begin() const;

	WS_MethodIterator
	ws_methods_end() const;

	/// @brief check order of methods
	bool
	check_methods_in_right_order( ScoreType const & score_type_in_first_method,
		ScoreType const & score_type_in_second_method ) const;

	/// @brief
	virtual
	void
	eval_onebody_energies( pose::Pose & pose ) const;

	/// This is now handled lazily by the Energies object itself
	/// void
	/// accumulate_residue_total_energies( pose::Pose & pose ) const;

	/// @brief
	void
	show( std::ostream & out ) const;

	/// @brief Merges in the weights of another score function
	///
	/// example(s):
	///     scorefxn.merge(scorefxn2)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.weights
	///     Energies
	///     create_score_function
	void
	merge( const ScoreFunction & scorefxn_to_be_merged );

	/// @brief Scores  <pose>  and shows the raw and weighted scores for each
	/// non-zero ScoreType
	///
	/// example(s):
	///     scorefxn.show(pose)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.weights
	///     Energies
	///     create_score_function
	void
	show( std::ostream & out,  pose::Pose & pose ) const;

	/// @brief similar output as show( ostream, pose ) but without the pose
	void
	show_pretty( std::ostream & out ) const;


	/// @brief Scores  <pose>  and shows the raw and weighted scores for each
	/// non-zero ScoreType
	/// @note: this function is mostly for convenience in PyRosetta
	///
	/// example(s):
	///     scorefxn.show(pose)
	/// See also:
	///     ScoreFunction
	///     ScoreFunction.weights
	///     Energies
	///     create_score_function
	void
	show(pose::Pose & pose ) const;


	void
	show_line_headers( std::ostream & out ) const;


	void
	show_line( std::ostream & out,  pose::Pose const & pose ) const;

	void
	show_additional( std::ostream & out,  pose::Pose & pose, bool verbose=false ) const;

	/// @brief Accumulates the unweighted one body energies of all context
	/// independent one body energies for  <pose>  Residue  <rsd>  into
	/// EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector

	void name( std::string const & weights_tag ) { name_ = weights_tag; }
	std::string get_name() const { return name_; }

	void
	eval_ci_1b(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted one body energies of all context
	/// dependent one body energies for  <pose>  Residue  <rsd>  into
	/// EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	virtual
	void
	eval_cd_1b(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted context independent two body
	/// interaction energies of  <pose>  between Residue  <rsd1>  and Residue
	/// <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_ci_2b(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted context independent two body
	/// interaction energies of  <pose>  between the backbones of Residue
	/// <rsd1>  and Residue  <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_ci_2b_bb_bb(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted short range context independent two
	/// body interaction energies of  <pose>  between the backbone of Residue
	/// <rsd1>  and the sidechain of Residue  <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_ci_2b_bb_sc(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted short range context dependent two body
	/// interaction energies of  <pose>  between the sidechains of Residue  <rsd1>
	/// and Residue  <rsd2>  into Energymap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_ci_2b_sc_sc(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	/// @brief Accumulate the unweighted short range context dependent two body
	/// interaction energies of  <pose>  between Residue  <rsd1>  and Residue
	/// <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	virtual
	void
	eval_cd_2b(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	/// @brief Accumulates the unweighted short ranged context dependent two body
	/// interaction energies of  <pose>  between the backbones of Residue
	/// <rsd1>  and  <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_cd_2b_bb_bb(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted short ranged context dependent two body
	/// interaction energies of  <pose>  between the backbone of Residue  <rsd1>
	/// and the sidechain of Residue  <rsd2>  into EnergyMap <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_cd_2b_bb_sc(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted short ranged context dependent two body
	/// interaction energies of  <pose>  between the sidechains of Residue
	/// <rsd1>  and Residue  <rsd2>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_cd_2b_sc_sc(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates for rsd the unweighted intra-residue one body energies
	/// for all context dependent and context independent two body terms that
	/// define intra-residue energies
	/// @note: EnergyMap is an EMapVector
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted intra-residue one body energies for all
	/// context independent two body terms that define intra-residue energies of
	/// <pose>  Residue  <rsd>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	void
	eval_ci_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Accumulates the unweighted intra-residue one body energies for all
	/// context dependent two body terms that define intra-residue energies of
	/// <pose>  Residue  <rsd>  into EnergyMap  <emap>
	/// @note: EnergyMap is an EMapVector
	virtual
	void
	eval_cd_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	void
	bump_check_full(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	/// @brief Scores the sidechain from  <pose>  Residue  <rsd1>  against the
	/// backbone of Residue  <rsd2>
	void
	bump_check_backbone(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	void
	evaluate_rotamer_intrares_energies(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		utility::vector1< core::PackerEnergy > & energies
	) const;

	void
	evaluate_rotamer_intrares_energy_maps(
		conformation::RotamerSetBase const & set,
		pose::Pose const & pose,
		utility::vector1< EnergyMap > & emaps
	) const;


	void
	evaluate_rotamer_pair_energies(
		conformation::RotamerSetBase const & set1,
		conformation::RotamerSetBase const & set2,
		pose::Pose const & pose,
		ObjexxFCL::FArray2D< core::PackerEnergy > & energy_table
	) const;


	void
	evaluate_rotamer_background_energies(
		conformation::RotamerSetBase const & set1,
		conformation::Residue const & residue2,
		pose::Pose const & pose,
		utility::vector1< core::PackerEnergy > & energy_vector
	) const;

	bool
	any_lr_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;


	virtual
	void
	eval_npd_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		Vector & F1,
		Vector & F2
	) const;


	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose
	) const;


	ScoreTypes const &
	score_types_by_method_type( methods::EnergyMethodType const & t ) const
	{
		return score_types_by_method_type_[ t ];
	}

	/// @brief convenience access to all ci 1b score types
	ScoreTypes const &
	ci_1b_types() const
	{
		return score_types_by_method_type_[ methods::ci_1b ];
	}

	/// @brief convenience access to all cd 1b score types
	ScoreTypes const &
	cd_1b_types() const
	{
		return score_types_by_method_type_[ methods::cd_1b ];
	}

	/// @brief convenience access to all ci 2b score types
	ScoreTypes const &
	ci_2b_types() const
	{
		return score_types_by_method_type_[ methods::ci_2b ];
	}

	/// @brief convenience access to all cd 2b score types
	ScoreTypes const &
	cd_2b_types() const
	{
		return score_types_by_method_type_[ methods::cd_2b ];
	}


	ScoreTypes const &
	ci_lr_2b_types() const
	{
		return score_types_by_method_type_[ methods::ci_lr_2b ];
	}

	ScoreTypes const &
	cd_lr_2b_types() const
	{
		return score_types_by_method_type_[ methods::cd_lr_2b ];
	}

	ScoreTypes const &
	whole_structure_types() const
	{
		return score_types_by_method_type_[ methods::ws ];
	}

	bool
	ready_for_nonideal_scoring() const;


public:

	AllMethods const & all_methods() const;

	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////

private:


	void
	add_method( methods::EnergyMethodOP method );


	void
	remove_method( methods::EnergyMethodOP method );

	/// private
	bool
	check_methods() const;

	/// @brief sets up the derived method vectors from the primary vectors
	void
	initialize_methods_arrays();

	/// @brief check if any of the 2body methods define intra-residue energies, given
	/// our current weight set
	void
	update_intrares_energy_status();

protected:

	bool
	any_intrares_energies() const
	{
		return any_intrares_energies_;
	}

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
private:

	/// the scorefxn name
	std::string name_;

	/// the scorefxn weights
	EnergyMap weights_;

	/// @brief  Options that control choice/configuration of EnergyMethods.
	/// eg etable name
	methods::EnergyMethodOptionsOP energy_method_options_;

	///////////////////////////////////////////////////////////
	// scoring methods
	/// NOTE: if you add another array here, be sure to reset it inside
	/// initialize_methods_arrays or the copy c-tor won't work

	/// these are the primary data -- others can be rederived by calling
	/// update_method_vectors
	/// filled by the private member function add_method
	CI_2B_Methods     ci_2b_methods_;
	CD_2B_Methods     cd_2b_methods_;
	CI_1B_Methods     ci_1b_methods_;
	CD_1B_Methods     cd_1b_methods_;
	CI_LR_2B_Methods  ci_lr_2b_methods_;
	CD_LR_2B_Methods  cd_lr_2b_methods_;
	LR_2B_Methods     lr_2b_methods_;
	WS_Methods        ws_methods_;

	/////////////////////////////////////////////////////////////////////////////
	/// convenience access to the methods/alternate lookups
	/// filled by add_method along with the primary methods vectors

	/// vector of COPs to all active methods
	AllMethods all_methods_;

	/// @brief  Map from ScoreType to the method for evaluating it.
	/// will be 0 for score_types without a method, ie ones with 0 weight
	AllMethods methods_by_score_type_;


	utility::vector1< ScoreTypes > score_types_by_method_type_;

	mutable bool score_function_info_current_;
	mutable ScoreFunctionInfoOP score_function_info_;

	// 2body energies that define intra-residue energies are allowed to
	// define "one body" intra-energies; the loops to evaluate intrares
	// energies are evaluated only if some intrares energies exist
	// These lists include both short and long range energy funcitons
	bool any_intrares_energies_;
	TWO_B_Methods ci_2b_intrares_;
	TWO_B_Methods cd_2b_intrares_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // ScoreFunction

inline
std::ostream &
operator <<( std::ostream & out, ScoreFunction const & sf )
{
	sf.show( out );
	return out;
}

/// @brief Utility function to locate a weights or patch file, either with a fully qualified path,
/// in the local directory, or in the database. Names may be passes either with or without the
/// optional extension.

std::string
find_weights_file(std::string name, std::string extension=".wts");

} // namespace scoring
} // namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_ScoreFunction )
#endif // SERIALIZATION


#endif // INCLUDED_core_scoring_ScoreFunction_HH
