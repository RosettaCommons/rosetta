// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/DecomposeAndReweightEnergiesCalculator.cc
/// @brief  DecomposeAndReweightEnergiesCalculator class
/// @author Colin A. Smith

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/DecomposeAndReweightEnergiesCalculator.hh>

// project headers
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/Methods.hh> //for long range energies
#include <core/scoring/LREnergyContainer.hh> //long range energies
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <utility/assert.hh>

#include <utility/vector1.hh>
#include <set>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace utility;
using namespace ObjexxFCL::format;

#ifdef    SERIALIZATION
// Package serialization headers
#include <utility/graph/UpperEdgeGraph.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {


// preferred constructor - use an existing InterfaceNeighborDefinitionCalculator
DecomposeAndReweightEnergiesCalculator::DecomposeAndReweightEnergiesCalculator( std::string const & NameOfResidueDecompositionCalculator ) :
	EnergyDependentCalculator(),
	name_of_ResidueDecompositionCalculator_(NameOfResidueDecompositionCalculator)
{
	if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( name_of_ResidueDecompositionCalculator_ ) ) {
		basic::Error() << "Tried to tie DecomposeAndReweightEnergiesCalculator to ResidueDecompositionCalculator " <<
			name_of_ResidueDecompositionCalculator_ << " but this calculator does not exist." << std::endl;
		utility_exit();
	}
}

DecomposeAndReweightEnergiesCalculator::DecomposeAndReweightEnergiesCalculator(
	DecomposeAndReweightEnergiesCalculator const & calculator
) :
	EnergyDependentCalculator(),
	name_of_ResidueDecompositionCalculator_(calculator.residue_decomposition_calculator()),
	original_weights_(calculator.original_weights()),
	other_energies_(calculator.other_energies()),
	onebody_energies_(calculator.onebody_energies()),
	twobody_energies_(calculator.twobody_energies()),
	set_names_(calculator.set_names()),
	weighted_total_(calculator.weighted_total())
{}

void DecomposeAndReweightEnergiesCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "weighted_total" ) {
		basic::check_cast( valptr, &weighted_total_, "weighted_total expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( weighted_total_ );

	} else if ( key == "master_weight_vector" ) {
		utility::vector1<core::Real> const master_weight_vector_copy(master_weight_vector());
		basic::check_cast( valptr, &master_weight_vector_copy, "master_weight_vector expects to return a utility::vector1< core::Real >" );
		(static_cast<basic::MetricValue<utility::vector1<core::Real> > *>(valptr))->set( master_weight_vector_copy );

	} else if ( key == "weighted_total_vector" ) {
		utility::vector1<core::Real> const weighted_total_vector_copy(weighted_total_vector());
		basic::check_cast( valptr, &weighted_total_vector_copy, "weighted_total_vector expects to return a utility::vector1< core::Real >" );
		(static_cast<basic::MetricValue<utility::vector1<core::Real> > *>(valptr))->set( weighted_total_vector_copy );

	}  else if ( key == "weighted_total_no_master_vector" ) {
		utility::vector1<core::Real> const weighted_total_no_master_vector_copy(weighted_total_no_master_vector());
		basic::check_cast( valptr, &weighted_total_no_master_vector_copy, "weighted_total_no_master_vector expects to return a utility::vector1< core::Real >" );
		(static_cast<basic::MetricValue<utility::vector1<core::Real> > *>(valptr))->set( weighted_total_no_master_vector_copy );

	} else if ( key == "summary" ) {
		std::ostringstream sstream;
		show(sstream);
		std::string const & summary(sstream.str());
		basic::check_cast( valptr, &summary, "summary expects to return a std::string" );
		(static_cast<basic::MetricValue<std::string> *>(valptr))->set( summary );

	} else {
		basic::Error() << "DecomposeAndReweightEnergiesCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}
}


std::string DecomposeAndReweightEnergiesCalculator::print( std::string const & key ) const {

	if ( key == "weighted_total" ) {
		return utility::to_string( weighted_total_ );
	} else if ( key == "master_weight_vector" ) {
		utility::vector1<core::Real> const master_weight_vector_copy(master_weight_vector());
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= master_weight_vector_copy.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << master_weight_vector_copy[i];
		}
		sstream << "]";
		return sstream.str();
	} else if ( key == "weighted_total_vector" ) {
		utility::vector1<core::Real> const weighted_total_vector_copy(weighted_total_vector());
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= weighted_total_vector_copy.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << weighted_total_vector_copy[i];
		}
		sstream << "]";
		return sstream.str();
	} else if ( key == "weighted_total_no_master_vector" ) {
		utility::vector1<core::Real> const weighted_total_no_master_vector_copy(weighted_total_no_master_vector());
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= weighted_total_no_master_vector_copy.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << weighted_total_no_master_vector_copy[i];
		}
		sstream << "]";
		return sstream.str();
	} else if ( key == "summary" ) {
		std::ostringstream sstream;
		show(sstream);
		return sstream.str();
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}

void
DecomposeAndReweightEnergiesCalculator::recompute(
	Pose const & this_pose
)
{
	runtime_assert(this_pose.energies().energies_updated());

	basic::MetricValue<utility::vector1<std::set<core::Size> > > residue_decomposition_value;
	basic::MetricValue<utility::vector1<core::Size> > residue_set_numbers_value;
	basic::MetricValue<utility::vector1<std::string > > set_names_value;
	this_pose.metric(name_of_ResidueDecompositionCalculator_, "residue_decomposition", residue_decomposition_value);
	this_pose.metric(name_of_ResidueDecompositionCalculator_, "residue_set_numbers", residue_set_numbers_value);
	this_pose.metric(name_of_ResidueDecompositionCalculator_, "set_names", set_names_value);
	utility::vector1<std::set<core::Size> > const & residue_decomposition(residue_decomposition_value.value());
	utility::vector1<core::Size> const & residue_set_numbers(residue_set_numbers_value.value());
	set_names_ = set_names_value.value();
	runtime_assert(residue_set_numbers.size() == this_pose.size());

	core::Size num_sets_from_decomposition = residue_decomposition.size();
	// if num_sets is 0, set num_sets to be the size of the decomposition
	if ( num_sets() == 0 ) num_sets(num_sets_from_decomposition);
	// make sure the weight vector and weight graph are the right size
	runtime_assert(onebody_energies_.size() == num_sets_from_decomposition &&
		twobody_energies_.num_vertices() == num_sets_from_decomposition);

	clear_energies();

	// add one body energies to the correct one body location
	for ( core::Size i = 1; i <= this_pose.size(); ++i ) {
		core::Size const set_num(residue_set_numbers[i]);
		if ( set_num ) onebody_energies_[set_num].energy_map() += this_pose.energies().onebody_energies(i);
	}

	core::scoring::EnergyGraph const & energy_graph(this_pose.energies().energy_graph());

	// add two body energies to the correct one body or two body locations
	for ( utility::graph::Graph::EdgeListConstIter iter = energy_graph.const_edge_list_begin();
			iter != energy_graph.const_edge_list_end(); ++iter ) {
		core::scoring::EnergyEdge const * const energy_edge = static_cast<core::scoring::EnergyEdge const *>(*iter);
		core::Size const first_set_num(residue_set_numbers[energy_edge->get_first_node_ind()]);
		core::Size const second_set_num(residue_set_numbers[energy_edge->get_second_node_ind()]);
		// only place the energies if both residues in a set
		if ( first_set_num && second_set_num ) {
			if ( first_set_num == second_set_num ) {
				// add two body energies between residues in the same set to the one body energy for that set
				energy_edge->add_to_energy_map( onebody_energies_[first_set_num].energy_map() );
			} else {
				// add two body energies between residues in different sets to the EnergyEdge between those sets
				energy_edge->add_to_energy_map( twobody_energies_.get_edge(first_set_num, second_set_num)->data().energy_map() );
			}
		}
	}

	// add long range two body energies to the correct onebody or twobody locations
	for ( Size lr = 1; lr <= scoring::methods::n_long_range_types; lr++ ) {
		scoring::methods::LongRangeEnergyType lr_type = scoring::methods::LongRangeEnergyType( lr );
		scoring::LREnergyContainerCOP lrec = this_pose.energies().long_range_container( lr_type );
		if ( !lrec ) continue;
		if ( lrec->empty() ) continue;

		// iterate over all residues
		for ( core::Size i = 1; i <= this_pose.size(); ++i ) {
			core::Size const first_set_num(residue_set_numbers[i]);
			if ( !first_set_num ) continue;
			// iterate over all neighbors of the residue
			for ( core::scoring::ResidueNeighborConstIteratorOP rni = lrec->const_upper_neighbor_iterator_begin(i);
					*rni != *( lrec->const_upper_neighbor_iterator_end(i) ); ++(*rni) ) {
				core::Size const second_set_num(residue_set_numbers[rni->upper_neighbor_id()]);
				// only place the energies if both residues both part of a set
				if ( !second_set_num ) continue;
				if ( first_set_num == second_set_num ) {
					// add two body energies between residues in the same set to the one body energy for that set
					rni->accumulate_energy(onebody_energies_[first_set_num].energy_map());
				} else {
					// add two body energies between residues in different sets to the EnergyEdge between those sets
					rni->accumulate_energy(twobody_energies_.get_edge(first_set_num, second_set_num)->data().energy_map());
				}
			}
		}
	}

	other_energies_.energy_map() = this_pose.energies().total_energies();

	// subtract out the one body and two body energies from other_energies_
	core::Size nc = num_components();
	for ( core::Size i = 2; i <= nc; ++i ) {
		other_energies_.energy_map() -= component(i).energy_map();
	}

	original_weights_ = this_pose.energies().weights();
	update_original_weights();

	update_weighted_total();
}

core::Size
DecomposeAndReweightEnergiesCalculator::num_sets() const
{
	return onebody_energies_.size();
}

void
DecomposeAndReweightEnergiesCalculator::num_sets(core::Size num_sets) {

	onebody_energies_.resize(num_sets);

	// avoid destruction of edges if the graphs are already the correct size
	if ( twobody_energies_.num_vertices() != num_sets ) twobody_energies_.set_num_vertices(num_sets);

	for ( core::Size i = 1; i < twobody_energies_.num_vertices(); ++i ) {
		for ( core::Size j = i+1; j <= twobody_energies_.num_vertices(); ++j ) {
			if ( !twobody_energies_.edge_exists(i, j) ) {
				twobody_energies_.add_edge(i, j, EnergiesData());
			}
		}
	}
}

core::Size
DecomposeAndReweightEnergiesCalculator::num_components() const
{
	return 1+(num_sets()+1)*num_sets()/2;
}

EnergiesData const &
DecomposeAndReweightEnergiesCalculator::component(
	core::Size index
) const
{
	runtime_assert(index >= 1 && index <= num_components());

	if ( index == 1 ) return other_energies_;
	--index;

	if ( index <= onebody_energies_.size() ) return onebody_energies_[index];
	index -= onebody_energies_.size();

	core::Size set1;
	for ( set1 = 1; index > num_sets() - set1; ++set1 ) {
		index -= num_sets() - set1;
	}
	core::Size set2 = set1+index;
	return twobody_energies_.get_vertex(set1).get_edge(set2)->data();
}

EnergiesData &
DecomposeAndReweightEnergiesCalculator::component(
	core::Size index
)
{
	runtime_assert(index >= 1 && index <= num_components());

	if ( index == 1 ) return other_energies_;
	--index;

	if ( index <= onebody_energies_.size() ) return onebody_energies_[index];
	index -= onebody_energies_.size();

	core::Size set1;
	for ( set1 = 1; index > num_sets() - set1; ++set1 ) {
		index -= num_sets() - set1;
	}
	core::Size set2 = set1+index;
	return twobody_energies_.get_vertex(set1).get_edge(set2)->data();
}

utility::vector1<core::Real>
DecomposeAndReweightEnergiesCalculator::master_weight_vector() const
{
	utility::vector1<core::Real> master_weight_vector(num_components());

	for ( core::Size i = 1; i <= master_weight_vector.size(); ++i ) {
		master_weight_vector[i] = component(i).master_weight();
	}

	return master_weight_vector;
}

void
DecomposeAndReweightEnergiesCalculator::master_weight_vector(
	utility::vector1<core::Real> const & master_weight_vector
)
{
	// determine the number sets that the vector corresponds to
	core::Size num_sets_for_vector = static_cast<core::Size> (floor(sqrt(float((master_weight_vector.size()-1)*2))));
	// apply the reverse mapping to make sure it is the correct length
	runtime_assert(master_weight_vector.size() == 1+(num_sets_for_vector+1)*num_sets_for_vector/2);
	// if num_sets is 0, set num_sets to be the size of the vector
	if ( num_sets() == 0 ) num_sets(num_sets_for_vector);
	// make sure that the vector and num_sets are are compatible
	runtime_assert(num_sets() == num_sets_for_vector);

	for ( core::Size i = 1; i <= master_weight_vector.size(); ++i ) {
		component(i).master_weight(master_weight_vector[i]);
	}
}

utility::vector1<std::string>
DecomposeAndReweightEnergiesCalculator::names_vector() const
{
	utility::vector1<std::string> names_vec;

	names_vec.push_back("Other");

	names_vec.insert(names_vec.end(), set_names_.begin(), set_names_.end());

	for ( core::Size i = 1; i < set_names_.size(); ++i ) {
		for ( core::Size j = i+1; j <= set_names_.size(); ++j ) {
			names_vec.push_back(set_names_[i]+"-"+set_names_[j]);
		}
	}

	return names_vec;
}

utility::vector1<core::scoring::EnergyMap>
DecomposeAndReweightEnergiesCalculator::weighted_energy_map_vector() const
{
	utility::vector1<core::scoring::EnergyMap> weighted_energy_map_vector_result(num_components());

	for ( core::Size i = 1; i <= weighted_energy_map_vector_result.size(); ++i ) {
		weighted_energy_map_vector_result[i] = component(i).weighted_energy_map();
	}

	return weighted_energy_map_vector_result;
}

utility::vector1<core::Real>
DecomposeAndReweightEnergiesCalculator::weighted_total_no_master_vector() const
{
	utility::vector1<core::Real> weighted_total_vector_result(num_components());

	for ( core::Size i = 1; i <= weighted_total_vector_result.size(); ++i ) {
		weighted_total_vector_result[i] = component(i).weighted_total_no_master();
	}

	return weighted_total_vector_result;
}

utility::vector1<core::Real>
DecomposeAndReweightEnergiesCalculator::weighted_total_vector() const
{
	utility::vector1<core::Real> weighted_total_vector_result(num_components());

	for ( core::Size i = 1; i <= weighted_total_vector_result.size(); ++i ) {
		weighted_total_vector_result[i] = component(i).weighted_total();
	}

	return weighted_total_vector_result;
}

utility::vector1<core::scoring::ScoreType>
DecomposeAndReweightEnergiesCalculator::nonzero_weight_score_types() const
{
	utility::vector1<core::scoring::EnergyMap> weight_map_vector(num_components());

	for ( core::Size i = 1; i <= weight_map_vector.size(); ++i ) {
		weight_map_vector[i] = component(i).weight_map();
	}

	utility::vector1<core::scoring::ScoreType> score_types;

	for ( int i = 1; i < core::scoring::n_score_types; ++i ) {
		for ( core::Size j = 1; j < weight_map_vector.size(); ++j ) {
			if ( weight_map_vector[j][core::scoring::ScoreType(i)] ) {
				score_types.push_back(core::scoring::ScoreType(i));
				break;
			}
		}
	}

	return score_types;
}

void
DecomposeAndReweightEnergiesCalculator::show(
	std::ostream & out
) const
{
	utility::vector1<core::scoring::ScoreType> const score_types(nonzero_weight_score_types());
	utility::vector1<core::scoring::EnergyMap> const weighted_energy_map_vec(weighted_energy_map_vector());
	utility::vector1<core::Real> const weighted_total_vec(weighted_total_vector());
	utility::vector1<std::string> const names_vec(names_vector());

	out << "------------------------";
	for ( core::Size i = 1; i <= names_vec.size(); ++i ) out << "---------";
	out << "---------" << std::endl;

	out << LJ(24,"Score Type");
	for ( core::Size i = 1; i <= names_vec.size(); ++i ) out << RJ(9,names_vec[i]);
	out << RJ(9,"Total") << std::endl;

	out << "------------------------";
	for ( core::Size i = 1; i <= names_vec.size(); ++i ) out << "---------";
	out << "---------" << std::endl;

	for ( core::Size i = 1; i <= score_types.size(); ++i ) {
		core::Real total = 0;
		out << LJ(24,score_types[i]);
		for ( core::Size j = 1; j <= weighted_energy_map_vec.size(); ++j ) {
			out << F(9,3, weighted_energy_map_vec[j][score_types[i]]);
			total += weighted_energy_map_vec[j][score_types[i]];
		}
		out << F(9,3, total) << std::endl;
	}

	out << "------------------------";
	for ( core::Size i = 1; i <= names_vec.size(); ++i ) out << "---------";
	out << "---------" << std::endl;

	out << LJ(24,"Total");
	for ( core::Size i = 1; i <= weighted_total_vec.size(); ++i ) out << F(9,3, weighted_total_vec[i]);
	out << F(9,3, weighted_total_) << std::endl;
}

void
DecomposeAndReweightEnergiesCalculator::clear_energies()
{
	core::Size nc = num_components();
	for ( core::Size i = 1; i <= nc; ++i ) {
		component(i).energy_map().clear();
	}

	weighted_total_ = 0;
}

void
DecomposeAndReweightEnergiesCalculator::update_original_weights()
{
	core::Size nc = num_components();
	for ( core::Size i = 1; i <= nc; ++i ) {
		EnergiesData & this_component(component(i));
		if ( this_component.use_original_weights() ) this_component.weight_map(original_weights_);
	}
}

void
DecomposeAndReweightEnergiesCalculator::update_weighted_total()
{
	weighted_total_ = 0;

	core::Size nc = num_components();
	for ( core::Size i = 1; i <= nc; ++i ) {
		weighted_total_ += component(i).weighted_total();
	}
}

#ifdef    SERIALIZATION
template < class Archive >
void
save(
	Archive & arc,
	DecomposeAndReweightEnergiesCalculator::EnergiesGraph const & g
)
{
	utility::graph::save_to_archive( arc, g );
}

template < class Archive >
void
load(
	Archive & arc,
	DecomposeAndReweightEnergiesCalculator::EnergiesGraph & g )
{
	utility::graph::load_from_archive( arc, g );
}

EXTERNAL_SAVE_AND_LOAD_SERIALIZABLE( DecomposeAndReweightEnergiesCalculator::EnergiesGraph );

#endif

} // PoseMetricCalculators
} // toolbox
} // protocols


#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::EnergiesData::save( Archive & arc ) const {
	arc( CEREAL_NVP( energy_map_ ) ); // core::scoring::EnergyMap
	arc( CEREAL_NVP( weight_map_ ) ); // core::scoring::EnergyMap
	arc( CEREAL_NVP( use_original_weights_ ) ); // _Bool
	arc( CEREAL_NVP( master_weight_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::EnergiesData::load( Archive & arc ) {
	arc( energy_map_ ); // core::scoring::EnergyMap
	arc( weight_map_ ); // core::scoring::EnergyMap
	arc( use_original_weights_ ); // _Bool
	arc( master_weight_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::EnergiesData );

/// @brief Default constructor required by cereal to deserialize this class
protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator::DecomposeAndReweightEnergiesCalculator() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( CEREAL_NVP( name_of_ResidueDecompositionCalculator_ ) ); // std::string
	arc( CEREAL_NVP( original_weights_ ) ); // core::scoring::EnergyMap
	arc( CEREAL_NVP( other_energies_ ) ); // class protocols::toolbox::pose_metric_calculators::EnergiesData
	arc( CEREAL_NVP( onebody_energies_ ) ); // utility::vector1<EnergiesData>
	arc( CEREAL_NVP( twobody_energies_ ) ); // EnergiesGraph
	arc( CEREAL_NVP( set_names_ ) ); // utility::vector1<std::string>
	arc( CEREAL_NVP( weighted_total_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::EnergyDependentCalculator >( this ) );
	arc( name_of_ResidueDecompositionCalculator_ ); // std::string
	arc( original_weights_ ); // core::scoring::EnergyMap
	arc( other_energies_ ); // class protocols::toolbox::pose_metric_calculators::EnergiesData
	arc( onebody_energies_ ); // utility::vector1<EnergiesData>
	arc( twobody_energies_ ); // EnergiesGraph
	arc( set_names_ ); // utility::vector1<std::string>
	arc( weighted_total_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::DecomposeAndReweightEnergiesCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_DecomposeAndReweightEnergiesCalculator )

#endif // SERIALIZATION
