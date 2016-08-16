// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/VdwGrid.cc
/// @author Sam DeLuca

#include <protocols/qsar/scoring_grid/VdwGrid.hh>
#include <protocols/qsar/scoring_grid/VdwGridCreator.hh>


#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

#include <basic/database/open.hh>

#include <numeric/interpolation/util.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>

#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace qsar {
namespace scoring_grid {

std::string VdwGridCreator::keyname() const
{
	return VdwGridCreator::grid_name();
}

GridBaseOP VdwGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP vdw_grid( new VdwGrid() );

	vdw_grid->parse_my_tag(tag);

	return vdw_grid;
}

GridBaseOP VdwGridCreator::create_grid() const
{
	return GridBaseOP( new VdwGrid() );
}


std::string VdwGridCreator::grid_name()
{
	return "VdwGrid";
}

VdwGrid::VdwGrid() : SingleGrid("VdwGrid"), cutoff_(10.0)
{
	std::string lj_file(basic::database::full_name("scoring/qsar/lj_table.txt"));
	lj_spline_ = numeric::interpolation::spline_from_file(lj_file,0.01).get_interpolator();
}

void
VdwGrid::parse_my_tag(utility::tag::TagCOP  /*tag*/){

}


void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const &  )
{
	// loop through all the atoms in the pose
	// get the VDW radius of the atom
	// for each square within cutoff of the atom, update the score
	// continue

	core::Size chain_id = core::pose::get_chain_id_from_chain(get_chain(),pose);
	//core::Size chain_begin = pose.conformation().chain_begin(chain_id);
	//core::Size chain_end = pose.conformation().chain_end(chain_id);

	this->fill_with_value(cutoff_);


	for ( core::Size residue_index = 1; residue_index <= pose.n_residue(); ++residue_index ) {
		core::conformation::Residue residue = pose.residue(residue_index);
		if ( residue.chain() == chain_id ) {
			continue;
		}
		for ( core::Size atom_index = 1; atom_index <= residue.natoms(); ++atom_index ) {
			core::id::AtomID atom_id(atom_index,residue_index);
			core::Vector xyz(pose.xyz(atom_id));
			core::Real const & radius(residue.atom_type(atom_index).lj_radius());
			this->set_distance_sphere_for_atom(radius,xyz,cutoff_);
		}
	}
}

void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void VdwGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}

VdwGrid::~VdwGrid()
{

}

core::Real VdwGrid::score(
	core::conformation::UltraLightResidue const & residue,
	core::Real const max_score,
	qsarMapOP /*qsar_map*/)
{
	core::Real score = 0.0;

	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue[atom_index]);
		core::Real const & radius(residue.residue()->atom_type(atom_index).lj_radius());
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			core::Real spline_score = 0.0;
			core::Real spline_score_deriv = 0.0;

			lj_spline_->interpolate(max_radius-radius,spline_score,spline_score_deriv);

			score += spline_score;
		}
	}
	return score;
}

core::Real VdwGrid::atom_score(
	core::conformation::UltraLightResidue const & residue,
	core::Size atomno,
	qsarMapOP /*qsar_map*/)
{
	core::Vector const & atom_coord(residue[atomno]);
	core::Real const & radius(residue.residue()->atom_type(atomno).lj_radius());
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
		core::Real spline_score = 0.0;
		core::Real spline_score_deriv = 0.0;

		lj_spline_->interpolate(max_radius-radius,spline_score,spline_score_deriv);

		return spline_score;
	}
	return 0;
}

core::Real VdwGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP )
{
	core::Real score = 0.0;

	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue.xyz(atom_index));
		core::Real const & radius(residue.atom_type(atom_index).lj_radius());
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
			core::Real spline_score = 0.0;
			core::Real spline_score_deriv = 0.0;

			lj_spline_->interpolate(max_radius-radius,spline_score,spline_score_deriv);

			score += spline_score;
		}
	}
	return score;
}

core::Real VdwGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP /*qsar_map*/)
{
	core::Vector const & atom_coord(residue.xyz(atomno));
	core::Real const & radius(residue.atom_type(atomno).lj_radius());
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::Real max_radius = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
		core::Real spline_score = 0.0;
		core::Real spline_score_deriv = 0.0;

		lj_spline_->interpolate(max_radius-radius,spline_score,spline_score_deriv);

		return spline_score;
	}
	return 0;
}

utility::json_spirit::Value VdwGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair cutoff_data("cutoff",Value(cutoff_));
	Pair spline_data("spline",lj_spline_->serialize());
	Pair base_data("base_data",SingleGrid::serialize());
#ifdef PYROSETTA
		Value _;  return _;
#endif
	return Value(utility::tools::make_vector(cutoff_data,spline_data,base_data));
}

void VdwGrid::deserialize(utility::json_spirit::mObject data)
{
	cutoff_ = data["cutoff"].get_real();
	numeric::interpolation::spline::InterpolatorOP interp( new numeric::interpolation::spline::SimpleInterpolator );
	interp->deserialize(data["spline"].get_obj());
	lj_spline_ = interp;
	SingleGrid::deserialize(data["base_data"].get_obj());
}

}
}
}
