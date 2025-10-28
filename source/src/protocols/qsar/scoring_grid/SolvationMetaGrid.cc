// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationMetaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/SolvationMetaGrid.hh>
#include <protocols/qsar/scoring_grid/SolvationMetaGridCreator.hh>

#include <protocols/qsar/scoring_grid/SingleGrid.hh>
#include <protocols/qsar/scoring_grid/SolvationGrid.fwd.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>

#include <utility/tag/XMLSchemaGeneration.hh>

#include <map>

#include <protocols/qsar/scoring_grid/SolvationGrid.hh> // AUTO IWYU For SolvationGrid

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string SolvationMetaGridCreator::keyname() const
{
	return SolvationMetaGrid::grid_name();
}

GridBaseOP SolvationMetaGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP solvation_meta_grid( new SolvationMetaGrid() );
	solvation_meta_grid->parse_my_tag(tag);
	return solvation_meta_grid;
}

GridBaseOP SolvationMetaGridCreator::create_grid() const
{
	return utility::pointer::make_shared< SolvationMetaGrid >();
}

void SolvationMetaGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SolvationMetaGrid::provide_xml_schema( xsd );
}

SolvationMetaGrid::SolvationMetaGrid() :type_("SolvationMetaGrid")
{}

SolvationMetaGrid::SolvationMetaGrid( SolvationMetaGrid const & other):
	GridBase( other ),
	type_( other.type_ )
{
	for ( std::map<core::ShortSize,SingleGridOP>::value_type entry: grid_set_ ) {
		grid_set_[ entry.first ] = utility::pointer::dynamic_pointer_cast< SingleGrid >( entry.second->clone() );
	}
}

SolvationMetaGrid::~SolvationMetaGrid() = default;

/// @brief Make a copy of the grid, respecting the subclassing.
GridBaseOP SolvationMetaGrid::clone() const {
	return utility::pointer::make_shared< SolvationMetaGrid >( *this );
}

void SolvationMetaGrid::initialize(core::Vector const & center, core::Real width, core::Real resolution)
{
	core::chemical::AtomTypeSetCOP atom_type_set(
		core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard"));

	core::ShortSize max_atom_type = atom_type_set->n_atomtypes();
	for ( core::ShortSize atom_type = 1; atom_type <= max_atom_type; ++atom_type ) {

		SolvationGridOP new_grid( new SolvationGrid() );
		new_grid->set_probe_atom_type(atom_type);
		new_grid->initialize(center, width, resolution);
		grid_set_[atom_type] = SingleGridOP(new_grid);
	}
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center,utility::vector1<core::Size> )
{
	refresh(pose,center);
}

void SolvationMetaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center)
{
	auto it = grid_set_.begin();
	for ( ; it != grid_set_.end(); ++it ) {
		std::cout << "initializing solvation grid for atom type " <<it->first <<std::endl;
		it->second->refresh(pose,center);
	}
}

void SolvationMetaGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{
}

core::Real SolvationMetaGrid::score(
	core::conformation::UltraLightResidue const & residue,
	core::Real const /*max_score*/,
	qsarMapCOP /*qsar_map*/) const
{
	core::Real total_score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		//core::Vector const & coords = residue.xyz(atom_index);
		core::conformation::Atom current_atom(residue.residue()->atom(atom_index));
		auto grid_iterator(grid_set_.find(current_atom.type()));
		if ( grid_iterator == grid_set_.end() ) {
			utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
		}

		SingleGridOP current_grid = grid_iterator->second;
		total_score += current_grid->get_point(residue[atom_index]);
	}

	return total_score;
}

core::Real SolvationMetaGrid::atom_score(
	core::conformation::UltraLightResidue const & residue,
	core::Size atomno,
	qsarMapCOP /*qsar_map*/) const
{
	core::conformation::Atom current_atom(residue.residue()->atom(atomno));
	auto grid_iterator(grid_set_.find(current_atom.type()));
	if ( grid_iterator == grid_set_.end() ) {
		utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
	}

	SingleGridOP current_grid = grid_iterator->second;
	return  current_grid->get_point(residue[atomno]);
}

core::Real SolvationMetaGrid::score(core::conformation::Residue const & residue, core::Real const /*max_score*/, qsarMapCOP) const
{
	core::Real total_score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
		//core::Vector const & coords = residue.xyz(atom_index);
		core::conformation::Atom current_atom(residue.atom(atom_index));
		auto grid_iterator(grid_set_.find(current_atom.type()));
		if ( grid_iterator == grid_set_.end() ) {
			utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
		}

		SingleGridOP current_grid = grid_iterator->second;
		total_score += current_grid->get_point(current_atom.xyz());
	}

	return total_score;
}

core::Real SolvationMetaGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapCOP /*qsar_map*/) const
{
	core::conformation::Atom current_atom(residue.atom(atomno));
	auto grid_iterator(grid_set_.find(current_atom.type()));
	if ( grid_iterator == grid_set_.end() ) {
		utility_exit_with_message("Ligands must be parameterized with the FA_STANDARD atom type set for use in the SolvationMetaGrid");
	}

	SingleGridOP current_grid = grid_iterator->second;
	return  current_grid->get_point(current_atom.xyz());
}

std::string SolvationMetaGrid::get_type() const
{
	//return "SolvationMetaGrid";
	return grid_name();
}

void SolvationMetaGrid::set_chain(std::string const & chain)
{
	auto it = grid_set_.begin();
	for ( ; it != grid_set_.end(); ++it ) {
		it->second->set_chain(chain);
	}
}

void SolvationMetaGrid::dump_BRIX(std::string const & /*prefix*/) const
{
	utility_exit_with_message("SolvationMetaGrid is currently unable to output a BRIX grid, sorry :(");
}

utility::json_spirit::Value SolvationMetaGrid::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair type_record("type",Value(type_));
	std::vector<Value> grid_set_vector;

	for ( const auto & it : grid_set_ ) {
		std::vector<Value> current_map_data(utility::tools::make_vector(Value(core::ShortSize(it.first)),Value(it.second->serialize())));
		grid_set_vector.emplace_back(current_map_data);
	}

	Pair grid_set_data("grids",grid_set_vector);

#ifdef PYROSETTA
		Value _;  return _;
#endif

	return Value(utility::tools::make_vector(type_record,grid_set_data));

}

void SolvationMetaGrid::deserialize(utility::json_spirit::mObject data)
{
	type_ = data["type"].get_str();
	utility::json_spirit::mArray grid_set_data(data["grids"].get_array());
	for ( auto & it : grid_set_data ) {
		utility::json_spirit::mArray grid_pair(it.get_array());
		core::ShortSize atom_type = grid_pair[0].get_int();
		SingleGridOP grid( new SolvationGrid() );
		grid->deserialize(grid_pair[1].get_obj());
		grid_set_[atom_type] = grid;
	}
}

bool SolvationMetaGrid::is_in_grid(core::conformation::UltraLightResidue const & residue) const
{
	for ( const auto & it : grid_set_ ) {
		if ( !it.second->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

bool SolvationMetaGrid::is_in_grid(core::conformation::Residue const & residue) const
{
	for ( const auto & it : grid_set_ ) {
		if ( !it.second->is_in_grid(residue) ) {
			return false;
		}
	}
	return true;
}

std::string SolvationMetaGrid::grid_name()
{
	return "SolvationMetaGrid";
}

void SolvationMetaGrid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "grid_name", xs_string, "The name used to insert the scoring grid into the GridSet" );

	xsd_type_definition_w_attributes( xsd, grid_name(), "A collection of scoring grids based on the EEF1 (aka Lazaridis Karplus) solvation energy where each individual grid represents the desolvation energy of a particular atom type against the protein. This grid offers no customizable data.", attributes );
}

/// @brief Print a brief summary about this grid to the provided output stream
void SolvationMetaGrid::show( std::ostream & out ) const {
	out << "SolvationMetaGrid of type " << type_ << " with " << grid_set_.size() << " subcomponents." << std::endl;
}

std::string
SolvationMetaGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << type_;
	ss << sep << grid_set_.size();
	for ( auto const & entry: grid_set_ ) {
		ss << sep << entry.first << sep << entry.second->hash_fingerprint();
	}
	return ss.str();
}

}
}
}
