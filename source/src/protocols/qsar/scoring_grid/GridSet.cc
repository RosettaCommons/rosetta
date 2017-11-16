// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/qsar/scoring_grid/GridSet.cc
/// @brief A set of related grids
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <protocols/qsar/scoring_grid/GridSet.hh>

#include <protocols/qsar/scoring_grid/GridFactory.hh>
#include <protocols/qsar/scoring_grid/GridBase.hh>
#include <protocols/qsar/scoring_grid/ScoreNormalization.hh>
#include <protocols/qsar/qsarMap.hh>

#include <core/conformation/UltraLightResidue.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/Pose.hh>

#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.qsar.scoring_grid.GridSet" );


namespace protocols {
namespace qsar {
namespace scoring_grid {

GridSet::GridSet() = default;

GridSet::~GridSet() = default;

GridSet::GridSet( GridSet const & other) :
	grid_weights_( other.grid_weights_ ),
	width_( other.width_ ),
	resolution_( other.resolution_ ),
	chain_( other.chain_),
	qsar_map_( other.qsar_map_ ), // COP
	norm_function_( other.norm_function_ )
{
	// Need to do a deep copy on the grids, so the reinitialize() function only modifies this object.
	grids_.clear();
	for ( MappingType::value_type const & pair: other.grids_ ) {
		grids_[ pair.first ] = pair.second->clone();
	}
}

GridSetOP
GridSet::clone() const {
	return GridSetOP( new GridSet( *this ) );
}

void GridSet::set_normalization_function(std::string norm_function_name)
{
	norm_function_ = get_score_normalization_function(norm_function_name);
}

bool GridSet::is_normalization_enabled() const
{
	if ( norm_function_ ) {
		return true;
	} else {
		return false;
	}
}

void GridSet::set_qsar_map(qsarMapCOP qsar_map)
{
	qsar_map_ = qsar_map;
}

bool GridSet::is_qsar_map_attached() const
{
	if ( qsar_map_ != 0 ) {
		return true;
	} else {
		return false;
	}
}

void
GridSet::reinitialize( core::pose::Pose const & pose,
	core::Vector const & center )
{
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		GridBaseOP current_grid( grid_set_entry.second );
		TR.Debug <<"updating grid " << grid_set_entry.first << std::endl;
		current_grid->initialize( center, width_, resolution_ );
		current_grid->set_chain( chain_ );
		current_grid->refresh( pose, center );
	}
}

GridBaseCOP
GridSet::get_grid( std::string const & name ) const {
	if ( grids_.count( name ) == 0 ) {
		utility_exit_with_message("Grid " + name + " not found in grid set.");
	}
	return grids_.at( name );
}

void
GridSet::make_new_grid(utility::tag::TagCOP tag)
{
	std::string name = tag->getOption< std::string >( "grid_name" );
	core::Real weight = tag->getOption<core::Real>("weight", 1.0);

	if ( has_grid( name ) ) {
		// Should we print an error message? The previous version did not.
	} else {
		TR.Debug << "Making new grid " << name << std::endl;
		GridBaseOP new_grid(GridFactory::get_instance()->new_grid(tag));
		add_grid( name, new_grid, weight );
	}
}

void
GridSet::add_grid( std::string const & name, GridBaseOP grid, core::Real weight ) {
	if ( grids_.count( name ) != 0 ) {
		utility_exit_with_message( "GridSet already has a grid named " + name );
	}
	grids_[ name ] = grid;
	grid_weights_[ name ] = weight;
}

bool
GridSet::has_grid( std::string const & name ) const {
	return grids_.count( name ) != 0;
}

utility::vector1<std::string>
GridSet::get_grid_names() const
{
	utility::vector1<std::string> grid_names;

	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		grid_names.push_back(grid_set_entry.first);
	}
	return grid_names;
}

core::Size GridSet::size() const
{
	return grids_.size();
}



core::Real
GridSet::ideal_score(utility::vector1<core::conformation::UltraLightResidue> & residues) const
{

	core::Real score=0.0;

	//Does not use weighted average because the maximum score of each residue already scales with atom size
	//Hence, the contribution to total_score from larger ligands is already greater

	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		score += ideal_score(residue);
	}

	score = (score)/(residues.size());
	return score;

}

core::Real
GridSet::ideal_score(core::conformation::UltraLightResidue const & residue) const
{

	core::Real score = 0;
	score= -1.0 * residue.natoms();

	return score;

}

core::Real
GridSet::average_score(utility::vector1<core::conformation::UltraLightResidue> & residues) const
{
	core::Real score=0.0;

	//Does not use weighted average because the maximum score of each residue already scales with atom size
	//Hence, the contribution to total_score from larger ligands is already greater

	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		score += total_score(residue);
	}

	score = (score)/(residues.size());
	return score;

}

core::Real
GridSet::total_score(core::conformation::UltraLightResidue const & residue) const
{
	utility::vector1< core::conformation::UltraLightResidue const *> items;
	items.push_back( & residue );
	return total_score( items );
}

core::Real
GridSet::total_score(core::conformation::Residue const & residue) const
{
	utility::vector1< core::conformation::Residue const *> items;
	items.push_back( & residue );
	return total_score( items );
}

core::Real
GridSet::total_score(core::pose::Pose const & pose, core::Size const chain_id) const
{
	utility::vector1< core::Size > residues_in_chain( core::pose::get_resnums_for_chain_id( pose, chain_id ) );
	return total_score( pose, residues_in_chain );
}

core::Real
GridSet::total_score(core::pose::Pose const & pose, utility::vector1< core::Size > const & residues) const
{
	utility::vector1< core::conformation::Residue const *> items;
	for ( core::Size ii: residues ) {
		items.push_back( & pose.residue(ii) );
	}
	return total_score( items );
}

GridSet::ScoreMap
GridSet::grid_scores( core::conformation::Residue const & residue ) const {
	utility::vector1< core::conformation::Residue const *> items;
	items.push_back( & residue );
	return grid_scores( items );
}

GridSet::ScoreMap
GridSet::grid_scores(core::pose::Pose const & pose, utility::vector1< core::Size > const & residues) const
{
	utility::vector1< core::conformation::Residue const *> items;
	for ( core::Size ii: residues ) {
		items.push_back( & pose.residue(ii) );
	}
	return grid_scores( items );
}


GridSet::ScoreMap
GridSet::atom_score(core::pose::Pose const & /*pose*/, core::conformation::Residue const & residue, core::Size atomindex ) const
{
	GridSet::ScoreMap score_map;
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		GridBaseCOP current_grid(grid_set_entry.second);
		core::Real weight = 1.0;
		if ( grid_weights_.count( grid_set_entry.first ) != 0 ) {
			weight = grid_weights_.at(grid_set_entry.first);
		}
		core::Real atom_score = current_grid->atom_score(residue,atomindex,qsar_map_);
		std::string grid_type = current_grid->get_type();
		score_map.insert(std::make_pair(grid_type,atom_score*weight));
	}
	return score_map;
}


bool
GridSet::is_in_grid(utility::vector1<core::conformation::UltraLightResidue> const & residues) const
{
	for ( core::conformation::UltraLightResidue const & residue: residues ) {
		if ( is_in_grid(residue) == false ) {
			return false;
		}
	}

	return true;

}

bool
GridSet::is_in_grid(core::conformation::UltraLightResidue const & residue) const
{
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		GridBaseCOP current_grid(grid_set_entry.second);
		if ( !current_grid->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

bool
GridSet::is_in_grid(core::conformation::Residue const & residue) const
{
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		GridBaseCOP current_grid(grid_set_entry.second);
		if ( !current_grid->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

std::string
GridSet::hash_fingerprint() const {
	std::stringstream ss;

	// Use tab characters as separators, as that's unlikely to occur in the middle of strings
	char const sep('\t');

	ss << width_ << sep << resolution_ << sep << chain_;

	if ( qsar_map_ ) {
		// Probably should actually get a fingerprint from the qsar_map
		// but we're not really using it - hash on identity (pointer) for now.
		ss << sep << qsar_map_.get();
	} else {
		ss << sep << 0; // doesn't need to be interpretable
	}

	if ( norm_function_ ) {
		ss << sep << norm_function_->get_name();
	} else {
		ss << sep << "NONORM";
	}

	ss << sep << grid_weights_.size();
	for ( auto const & pair: grid_weights_ ) {
		ss << sep << pair.first << sep << pair.second;
	}

	ss << sep << grids_.size();
	for ( auto const & pair: grids_ ) {
		ss << sep << pair.first << sep << pair.second->hash_fingerprint();
	}

	return ss.str();
}

void
GridSet::write_grids(std::string prefix) const
{
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		GridBaseCOP current_grid(grid_set_entry.second);
		current_grid->dump_BRIX(prefix);
	}
}

utility::json_spirit::Value
GridSet::serialize() const
{
	using utility::json_spirit::Value;
	std::vector<Value> gridmap_data;
	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		Value grid_name(grid_set_entry.first);
		Value grid(grid_set_entry.second->serialize());
		std::vector<Value> grid_pair_values;
		grid_pair_values.push_back(grid_name);
		grid_pair_values.push_back(grid);
		Value grid_pair(grid_pair_values);
		gridmap_data.push_back(grid_pair);
	}
	return Value(gridmap_data);
}

void GridSet::deserialize(utility::json_spirit::mArray data)
{
	grids_.clear();
	for ( utility::json_spirit::mArray::iterator it = data.begin(); it != data.end(); ++it ) {
		utility::json_spirit::mArray grid_data(it->get_array());
		std::string grid_name = grid_data[0].get_str();
		GridBaseOP grid(GridFactory::get_instance()->new_grid(grid_data[1].get_obj()));
		grids_[grid_name] = grid;
	}
}


template< class GridScorable >
core::Real
GridSet::total_score(utility::vector1< GridScorable const * > const & items) const
{
	ScoreMap score_map( grid_scores(items) );

	core::Real total_score = 0.0;

	for ( ScoreMap::value_type const & score_map_entry: score_map ) {
		std::string const & name( score_map_entry.first );
		core::Real const & component_score( score_map_entry.second );
		core::Real weight = 1.0;
		if ( grid_weights_.count( name ) != 0 ) {
			weight = grid_weights_.at(name);
		}

		total_score += component_score*weight;
	}

	if ( norm_function_ ) {
		return normalize( total_score, items );
	} else {
		return total_score;
	}
}

core::Real
GridSet::normalize( core::Real total_score, utility::vector1< core::conformation::UltraLightResidue const* > const & items ) const
{
	if ( ! norm_function_ ) { return total_score; }

	core::conformation::ResidueCOPs residue_cops;
	for ( core::conformation::UltraLightResidue const* resi : items ) {
		residue_cops.push_back( resi->residue() );
	}
	core::Real normalized_score = (*norm_function_)(total_score,residue_cops);
	TR.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
	return normalized_score;
}

core::Real
GridSet::normalize( core::Real total_score, utility::vector1< core::conformation::Residue const* > const & items ) const {
	if ( ! norm_function_ ) { return total_score; }

	core::conformation::ResidueCOPs residue_cops;
	for ( core::conformation::Residue const* resi : items ) {
		residue_cops.push_back( resi->get_self_ptr() );
	}
	core::Real normalized_score = (*norm_function_)(total_score,residue_cops);
	TR.Trace << "Score normalized from " << total_score << " to "<< normalized_score << std::endl;
	return normalized_score;
}

template< class GridScorable >
GridSet::ScoreMap
GridSet::grid_scores(utility::vector1< GridScorable const * > const & items) const
{
	ScoreMap score_map;

	const core::Real max_score = 9999.0;

	for ( MappingType::value_type const & grid_set_entry: grids_ ) {
		core::Real component_score = 0;
		GridBaseCOP current_grid(grid_set_entry.second);
		for ( GridScorable const * item : items ) {
			core::Real current_score(current_grid->score(*item,max_score,qsar_map_));
			component_score += current_score;
		}
		score_map[ grid_set_entry.first ] = component_score;
	}

	return score_map;
}

} //protocols
} //qsar
} //scoring_grid






