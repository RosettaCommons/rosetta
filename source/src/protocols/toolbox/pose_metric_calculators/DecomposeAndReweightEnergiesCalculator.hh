// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Colin A. Smith


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_DecomposeAndReweightEnergiesCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_DecomposeAndReweightEnergiesCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <basic/MetricValue.fwd.hh>
#include <core/graph/UpperEdgeGraph.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

/*
class EnergiesEdge : public core::graph::Edge {
public:
core::scoring::EnergyMap energies;

virtual EnergiesEdge () {};

virtual void copy_from (Edge const *source) {
energies = static_cast< EnergiesEdge const * > ( source ).energies;
}

virtual unsigned int count_static_memory () const { return sizeof(EnergiesEdge); }
virtual unsigned int count_dynamic_memory () const { return 0; }
}

class EnergiesGraph : public core::graph::Graph {

public:

virtual ~Graph() {

}

virtual void  delete_edge(Edge *edge) {

}

virtual unsigned int count_static_memory() const { return sizeof(EnergiesGraph); }
//virtual unsigned int count_dynamic_memory() const { return 0; }
//virtual Node * create_new_node(int node_index) {}

virtual Edge * create_new_edge(int index1, int index2) {

}

virtual Edge * create_new_edge(Edge const *example_edge) {

}

private:


}

class WeightsEdge : public core::graph::Edge {

public:

core::scoring::EnergyMap weights;
bool use_original_weights;
core::Real master_weight;

virtual WeightsEdge () {};

virtual void copy_from (Edge const *source) {
energies = static_cast< EnergyMapEdge const * > ( source ).energies;
use_original_weights = static_cast< EnergyMapEdge const * > ( source ).use_original_weights;
master_weight = static_cast< EnergyMapEdge const * > ( source ).master_weight;
}

virtual unsigned int count_static_memory () const { return sizeof(WeightsEdge); }
virtual unsigned int count_dynamic_memory () const { return 0; }
}
*/

class EnergiesData {

public:

	EnergiesData():
		use_original_weights_(true),
		master_weight_(1.)
	{}

	core::scoring::EnergyMap & energy_map() { return energy_map_; }
	core::scoring::EnergyMap const & energy_map() const { return energy_map_; }
	void energy_map(core::scoring::EnergyMap const & energy_map) { energy_map_ = energy_map; }

	core::scoring::EnergyMap & weight_map() { return weight_map_; }
	core::scoring::EnergyMap const & weight_map() const { return weight_map_; }
	void weight_map(core::scoring::EnergyMap const & weight_map) { weight_map_ = weight_map; }

	bool use_original_weights() const { return use_original_weights_; }
	void use_original_weights(bool use_original_weights) { use_original_weights_ = use_original_weights; }

	core::Real master_weight() const { return master_weight_; }
	void master_weight(core::Real master_weight) { master_weight_ = master_weight; }

	core::scoring::EnergyMap weighted_energy_map() const {
		core::scoring::EnergyMap result(energy_map_);
		result *= weight_map_;
		result *= master_weight_;
		return result;
	}
	core::Real weighted_total_no_master() const { return energy_map_.dot(weight_map_); }
	core::Real weighted_total() const { return energy_map_.dot(weight_map_)*master_weight_; }

private:

	core::scoring::EnergyMap energy_map_;
	core::scoring::EnergyMap weight_map_;
	bool use_original_weights_;
	core::Real master_weight_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class EmptyVertexData {

public:

	// set to 1 because there will often only be two vertices
	static int const NUM_EDGES_TO_RESERVE = 1;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & ) const {}
	template< class Archive > void load( Archive & ) {}
#endif // SERIALIZATION

};

class DecomposeAndReweightEnergiesCalculator : public core::pose::metrics::EnergyDependentCalculator {

public:

	typedef core::graph::UpperEdgeGraph<EmptyVertexData, EnergiesData> EnergiesGraph;
	typedef core::graph::UEVertex<EmptyVertexData, EnergiesData> EnergiesVertex;
	typedef core::graph::UEEdge<EmptyVertexData, EnergiesData> EnergiesEdge;
	typedef utility::vector1<core::graph::UEEdge<EmptyVertexData, EnergiesData> >::iterator EnergiesUpperEdgeListIter;
	typedef utility::vector1<core::graph::UEEdge<EmptyVertexData, EnergiesData> >::const_iterator EnergiesUpperEdgeListConstIter;

	// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
	DecomposeAndReweightEnergiesCalculator(
		std::string const & NameOfResidueDecompositionCalculator
	);

	DecomposeAndReweightEnergiesCalculator(
		DecomposeAndReweightEnergiesCalculator const & calculator
	);

	core::pose::metrics::PoseMetricCalculatorOP clone() const
	{ return core::pose::metrics::PoseMetricCalculatorOP( new DecomposeAndReweightEnergiesCalculator( *this ) ); }

	std::string const & residue_decomposition_calculator() const { return name_of_ResidueDecompositionCalculator_; }
	core::scoring::EnergyMap const & original_weights() const { return original_weights_; }
	EnergiesData const & other_energies() const { return other_energies_; }
	utility::vector1<EnergiesData> const & onebody_energies() const { return onebody_energies_; }
	EnergiesGraph const & twobody_energies() const { return twobody_energies_; }
	utility::vector1<std::string> const & set_names() const { return set_names_; }
	core::Real weighted_total() const { return weighted_total_; }

	core::Size
	num_sets() const;

	void
	num_sets(
		core::Size num_sets
	);

	core::Size
	num_components() const;

	EnergiesData const &
	component(
		core::Size index
	) const;

	utility::vector1<core::Real>
	master_weight_vector() const;

	void
	master_weight_vector(
		utility::vector1<core::Real> const & master_weight_vector
	);

	utility::vector1<std::string>
	names_vector() const;

	utility::vector1<core::scoring::EnergyMap>
	weighted_energy_map_vector() const;

	utility::vector1<core::Real>
	weighted_total_no_master_vector() const;

	utility::vector1<core::Real>
	weighted_total_vector() const;

	utility::vector1<core::scoring::ScoreType>
	nonzero_weight_score_types() const;

	void
	show(
		std::ostream & out
	) const;

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	EnergiesData &
	component(
		core::Size index
	);

	void
	clear_energies();

	void
	update_original_weights();

	void
	update_weighted_total();

	std::string name_of_ResidueDecompositionCalculator_;

	core::scoring::EnergyMap original_weights_;
	EnergiesData other_energies_;
	utility::vector1<EnergiesData> onebody_energies_;
	EnergiesGraph twobody_energies_;
	utility::vector1<std::string> set_names_;

	core::Real weighted_total_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DecomposeAndReweightEnergiesCalculator();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::pointer::shared_ptr< DecomposeAndReweightEnergiesCalculator > DecomposeAndReweightEnergiesCalculatorOP;
typedef utility::pointer::shared_ptr< DecomposeAndReweightEnergiesCalculator const > DecomposeAndReweightEnergiesCalculatorCOP;


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_DecomposeAndReweightEnergiesCalculator )
#endif // SERIALIZATION


#endif
