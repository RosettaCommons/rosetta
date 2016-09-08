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

#include <protocols/qsar/scoring_grid/HbdGrid.hh>
#include <protocols/qsar/scoring_grid/HbdGridCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/io/mpistream.hh>

#include <numeric/interpolation/util.hh>

#include <basic/database/open.hh>

namespace protocols {
namespace qsar {
namespace scoring_grid {


std::string HbdGridCreator::keyname() const
{
	return HbdGridCreator::grid_name();
}

GridBaseOP HbdGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP hbd_grid( new HbdGrid() );

	hbd_grid->parse_my_tag(tag);

	return hbd_grid;
}

GridBaseOP HbdGridCreator::create_grid() const
{
	return GridBaseOP( new HbdGrid() );
}


std::string HbdGridCreator::grid_name()
{
	return "HbdGrid";
}

HbdGrid::HbdGrid(): SingleGrid ("HbdGrid")
{
	std::string lj_file(basic::database::full_name("scoring/qsar/hb_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.05).get_interpolator();
}


HbdGrid::~HbdGrid()
{

}

utility::json_spirit::Value HbdGrid::serialize()
{

	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair spline_data("spline",lj_spline_->serialize());
	Pair base_data("base_data",SingleGrid::serialize());

	return Value(utility::tools::make_vector(spline_data,base_data));
}

void HbdGrid::deserialize(utility::json_spirit::mObject data)
{

	lj_spline_->deserialize(data["spline"].get_obj());
	SingleGrid::deserialize(data["base_data"].get_obj());
}


void
HbdGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{
	this->fill_with_value(0.0);

	for ( core::Size residue_index = 1; residue_index <= pose.size(); ++residue_index ) {
		core::conformation::Residue const residue = pose.residue(residue_index);
		if ( !residue.is_protein() ) {
			continue;
		}
		for ( core::Size atom_index=1; atom_index <= residue.natoms(); ++  atom_index ) {
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if ( atom_type.is_hydrogen() ) {
				utility::vector1<core::Size> bonded_to_hydrogen(residue.bonded_neighbor(atom_index));
				for ( core::Size index = 1; index <= bonded_to_hydrogen.size(); ++index ) {
					if ( residue.atom_type(bonded_to_hydrogen[index]).is_donor() ) {
						core::id::AtomID atom_id(atom_index,residue_index);
						core::Vector xyz(pose.xyz(atom_id));
						this->set_score_sphere_for_atom(lj_spline_,xyz,5.0);
					}
				}
			}
		}
	}
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void HbdGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

core::Real HbdGrid::score(
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
			if ( atom_type.is_acceptor() ) {
				core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
				score += grid_value;
			}

		}
	}

	return score;
}

core::Real HbdGrid::atom_score(
	core::conformation::UltraLightResidue const & residue,
	core::Size atomno,
	qsarMapOP /*qsar_map*/)
{
	core::Vector const & atom_coord(residue[atomno]);
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::chemical::AtomType atom_type(residue.residue()->atom_type(atomno));
		if ( atom_type.is_acceptor() ) {
			core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			return grid_value;
		}
	}
	return 0;
}

core::Real HbdGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP )
{
	core::Real score = 0.0;
	//GridBaseTracer << "map size is: " << qsar_map->size() <<std::endl;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue.xyz(atom_index));
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::chemical::AtomType atom_type(residue.atom_type(atom_index));
			if ( atom_type.is_acceptor() ) {
				core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
				score += grid_value;
			}

		}
	}

	return score;
}

core::Real HbdGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP /*qsar_map*/)
{
	core::Vector const & atom_coord(residue.xyz(atomno));
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::chemical::AtomType atom_type(residue.atom_type(atomno));
		if ( atom_type.is_acceptor() ) {
			core::Real grid_value = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			return grid_value;
		}
	}
	return 0;
}

}
}
}
