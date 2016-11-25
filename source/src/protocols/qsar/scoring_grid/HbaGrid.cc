// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/HbaGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/HbaGrid.hh>
#include <protocols/qsar/scoring_grid/HbaGridCreator.hh>

#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>

#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>

#include <basic/database/open.hh>

#include <numeric/interpolation/util.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string HbaGridCreator::keyname() const
{
	return HbaGrid::grid_name();
}

GridBaseOP HbaGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP hba_grid( new HbaGrid() );

	hba_grid->parse_my_tag(tag);

	return hba_grid;
}

GridBaseOP HbaGridCreator::create_grid() const
{
	return GridBaseOP( new HbaGrid() );
}

void HbaGridCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HbaGrid::provide_xml_schema( xsd );
}

//std::string HbaGridCreator::grid_name()
//{
// return "HbaGrid";
//}

HbaGrid::HbaGrid() : SingleGrid("HbaGrid")
{
	std::string lj_file(basic::database::full_name("scoring/qsar/hb_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.05).get_interpolator();
}


HbaGrid::~HbaGrid()
{

}

utility::json_spirit::Value HbaGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair spline_data("spline",lj_spline_->serialize());
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(spline_data,base_data));

}

void HbaGrid::deserialize(utility::json_spirit::mObject data)
{
	lj_spline_->deserialize(data["spline"].get_obj());
	SingleGrid::deserialize(data["base_data"].get_obj());
}

void
HbaGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/){

}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{

	this->fill_with_value(0.0);

	for ( core::Size residue_index = 1; residue_index <= pose.size(); ++residue_index ) {
		core::conformation::Residue const residue = pose.residue(residue_index);
		if ( !residue.is_protein() ) {
			continue;
		}
		for ( core::Size atom_index=1; atom_index <= residue.natoms(); ++atom_index ) {
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if ( atom_type.is_acceptor() ) {
				core::id::AtomID atom_id(atom_index,residue_index);
				core::Vector xyz(pose.xyz(atom_id));
				this->set_score_sphere_for_atom(lj_spline_,xyz,5.0);
			}
		}
	}
}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void HbaGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

core::Real HbaGrid::score(
	core::conformation::UltraLightResidue const & residue,
	core::Real const max_score,
	qsarMapOP /*qsar_map*/)
{
	core::Real score = 0.0;
	//GridBaseTracer << "map size is: " << qsar_map->size() <<std::endl;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue[atom_index]);
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::chemical::AtomType atom_type(residue.residue()->atom_type(atom_index));
			if ( atom_type.is_hydrogen() ) {
				utility::vector1<core::Size> bonded_to_hydrogen(residue.residue()->bonded_neighbor(atom_index));
				for ( core::Size index = 1; index <= bonded_to_hydrogen.size(); ++index ) {
					if ( residue.residue()->atom_type(bonded_to_hydrogen[index]).is_donor() ) {
						core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
						score += grid_value;
					}
				}
			}
		}
	}

	return score;
}

core::Real HbaGrid::atom_score(
	core::conformation::UltraLightResidue const & residue,
	core::Size atomno,
	qsarMapOP /*qsar_map*/)
{
	core::Real score = 0;
	core::Vector const & atom_coord(residue[atomno]);
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::chemical::AtomType atom_type(residue.residue()->atom_type(atomno));
		if ( atom_type.is_hydrogen() ) {
			utility::vector1<core::Size> bonded_to_hydrogen(residue.residue()->bonded_neighbor(atomno));
			for ( core::Size index = 1; index <= bonded_to_hydrogen.size(); ++index ) {
				if ( residue.residue()->atom_type(bonded_to_hydrogen[index]).is_donor() ) {
					core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
					score += grid_value;
				}
			}
		}
		return score;
	}
	return 0;
}

core::Real HbaGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP /*qsar_map*/)
{
	core::Real score = 0.0;
	//GridBaseTracer << "map size is: " << qsar_map->size() <<std::endl;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue.xyz(atom_index));
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if ( atom_type.is_hydrogen() ) {
				utility::vector1<core::Size> bonded_to_hydrogen(residue.bonded_neighbor(atom_index));
				for ( core::Size index = 1; index <= bonded_to_hydrogen.size(); ++index ) {
					if ( residue.atom_type(bonded_to_hydrogen[index]).is_donor() ) {
						core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
						score += grid_value;
					}
				}
			}
		}
	}

	return score;
}

core::Real HbaGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP /*qsar_map*/)
{
	core::Real score = 0;
	core::Vector const & atom_coord(residue.xyz(atomno));
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::chemical::AtomType atom_type(residue.atom_type(atomno));
		if ( atom_type.is_hydrogen() ) {
			utility::vector1<core::Size> bonded_to_hydrogen(residue.bonded_neighbor(atomno));
			for ( core::Size index = 1; index <= bonded_to_hydrogen.size(); ++index ) {
				if ( residue.atom_type(bonded_to_hydrogen[index]).is_donor() ) {
					core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
					score += grid_value;
				}
			}
		}
		return score;
	}
	return 0;
}

std::string HbaGrid::grid_name()
{
	return "HbaGrid";
}

void HbaGrid::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute( "grid_name", xs_string, "The name used to insert the scoring grid into the GridManager" );

	xsd_type_definition_w_attributes( xsd, grid_name(), "A scoring grid that computes the hydrogen-bonding energy as given by the location of hydrogen-bond acceptors -- donor atoms can be queried against this grid; no parameters may be customized currently", attributes );

}

}
}
}
