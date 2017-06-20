// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/ligand_docking/scoring_grid/SolvationGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/SolvationGrid.hh>
#include <protocols/qsar/scoring_grid/SolvationGridCreator.hh>

#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/EtableEnergy.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <list>

namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string SolvationGridCreator::keyname() const
{
	return SolvationGrid::grid_name();
}

GridBaseOP SolvationGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP solvation_grid( new SolvationGrid() );

	solvation_grid->parse_my_tag(tag);

	return solvation_grid;
}

GridBaseOP SolvationGridCreator::create_grid() const
{
	return GridBaseOP( new SolvationGrid() );
}

void SolvationGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SolvationGrid::provide_xml_schema( xsd );
}



//std::string SolvationGridCreator::grid_name()
//{
// return "SolvationGrid";
//}

SolvationGrid::SolvationGrid() :SingleGrid("SolvationGrid")
{}

SolvationGrid::~SolvationGrid()
{}

GridBaseOP SolvationGrid::clone() const {
	return GridBaseOP( new SolvationGrid( *this ) );
}

void SolvationGrid::refresh(core::pose::Pose const & pose, core::Vector const & /*center*/)
{
	core::scoring::methods::EnergyMethodOptions default_options; // initialized from the command line
	core::scoring::etable::EtableCOP etable(
		core::scoring::ScoringManager::get_instance()->etable( default_options ));

	core::scoring::etable::TableLookupEvaluator etable_evaluator(*etable);

	//put all the atoms in a list so we don't have to deal with extracting them all more than once

	std::list<core::conformation::Atom> atom_list;
	for ( core::Size resnum = 1; resnum <=pose.size(); ++resnum ) {
		core::conformation::Residue current_residue = pose.residue(resnum);

		for ( core::Size atomnum = 1; atomnum <= current_residue.natoms(); ++atomnum ) {
			atom_list.push_back(current_residue.atom(atomnum));
		}
	}

	numeric::xyzVector<core::Size> dimensions = get_dimensions();
	for ( core::Size x_index =0; x_index < dimensions.x(); ++x_index ) {
		for ( core::Size y_index = 0; y_index < dimensions.y(); ++y_index ) {
			for ( core::Size z_index = 0; z_index < dimensions.z(); ++z_index ) {
				core::Vector pdb_coords(get_pdb_coords(x_index,y_index,z_index));
				core::conformation::Atom probe(pdb_coords,probe_atom_type_,1);

				core::Real total_solvation = 0.0;


				for ( std::list<core::conformation::Atom>::iterator it = atom_list.begin(); it != atom_list.end(); ++it ) {
					//the interface for the etable evaluator gets atr,rep,distance and solvation at once
					//atr,rep and d2 are dummy variables
					core::Real atr = 0.0;
					core::Real rep = 0.0;
					core::Real sol = 0.0;
					core::Real d2 = 0.0;

					etable_evaluator.atom_pair_energy(probe, *it, 1.0, atr, rep, sol, d2);
					total_solvation += sol;
				}

				set_point(pdb_coords,total_solvation);

			}
		}
	}


}

void SolvationGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}


void SolvationGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

void SolvationGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{

}


utility::json_spirit::Value SolvationGrid::serialize() const
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair probe_type_data("probe_type",probe_atom_type_);
	Pair base_data("base_data",SingleGrid::serialize());

#ifdef PYROSETTA
		Value _;  return _;
#endif

	return Value(utility::tools::make_vector(probe_type_data,base_data));

}

void SolvationGrid::deserialize(utility::json_spirit::mObject data )
{
	probe_atom_type_ = data["probe_type"].get_int();
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void SolvationGrid::set_probe_atom_type(core::ShortSize const & atom_type)
{
	probe_atom_type_ = atom_type;
}

std::string SolvationGrid::grid_name()
{
	return "SolvationGrid";
}

void SolvationGrid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "grid_name", xs_string, "The name used to insert the scoring grid into the GridSet" );

	xsd_type_definition_w_attributes( xsd, grid_name(), "A scoring grid based on the EEF1 (aka Lazaridis Karplus) solvation energy for a probe atom for this grid. DO NOT USE! The probe atom is not correctly initialized and this will give you garbage. Contact Sam Deluca for advice.", attributes );
}

std::string
SolvationGrid::hash_fingerprint() const {
	std::stringstream ss;
	const char sep('\t');
	ss << grid_name();
	ss << sep << get_type(); // Only thing of interest from parent class
	ss << sep << probe_atom_type_;
	return ss.str();
}

}
}
}
