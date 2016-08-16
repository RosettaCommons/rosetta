// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/qsar/scoring_grid/ChargeGrid.cc
/// @author Sam DeLuca
/// @brief This is an implementation of equation 3 in Goodford, J. Med. Chem. 1985,28,849-857.  doi:10.1021/jm00145a002

#include <protocols/qsar/scoring_grid/ChargeGrid.hh>
#include <protocols/qsar/scoring_grid/ChargeGridCreator.hh>
#include <protocols/qsar/scoring_grid/SingleGrid.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/UltraLightResidue.hh>
#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector.hh>
#include <utility/json_spirit/json_spirit_value.h>

#include <numeric/xyz.json.hh>

#include <boost/math/constants/constants.hpp>


namespace protocols {
namespace qsar {
namespace scoring_grid {

utility::json_spirit::Value ChargeAtom::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair xyz_record("xyz",numeric::serialize(xyz));
	Pair nc_record("nc",Value(static_cast<boost::uint64_t>(neighbor_count)));
	Pair charge_record("charge",Value(charge));

	return Value(utility::tools::make_vector(xyz_record,nc_record,charge_record));

}

void ChargeAtom::deserialize(utility::json_spirit::mObject data)
{
	charge = data["charge"].get_real();
	neighbor_count = data["nc"].get_int();
	xyz = numeric::deserialize<core::Real>(data["xyz"].get_array());
}
std::string ChargeGridCreator::keyname() const
{
	return ChargeGridCreator::grid_name();
}

GridBaseOP ChargeGridCreator::create_grid(utility::tag::TagCOP tag) const
{
	GridBaseOP charge_grid( new ChargeGrid() );
	charge_grid->parse_my_tag(tag);
	return charge_grid;
}

GridBaseOP ChargeGridCreator::create_grid() const
{

	return GridBaseOP( new ChargeGrid() );
}


std::string ChargeGridCreator::grid_name()
{
	return "ChargeGrid";
}

ChargeGrid::ChargeGrid() :
	SingleGrid("ChargeGrid"),
	zeta_(4.0),
	// epsilon_(80.0),
	indirect_numerator_( (4.0 - 80.0) / (4.0 + 80.0) ),
	epsilon_0_(8.854187817E-12)
{
	//
}

ChargeGrid::ChargeGrid(core::Real /*charge*/) :
	SingleGrid("ChargeGrid"),
	zeta_(4.0),
	// epsilon_(80),
	indirect_numerator_( (4.0 - 80.0) / (4.0 + 80.0) ),
	epsilon_0_(8.854187817E-12)
{
	//
}

utility::json_spirit::Value ChargeGrid::serialize()
{
	using utility::json_spirit::Value;
	using utility::json_spirit::Pair;

	Pair base_data("base_data",SingleGrid::serialize());

	std::vector<Value> charge_atom_data;
	for ( std::list<ChargeAtom>::iterator it = charge_atom_list_.begin(); it != charge_atom_list_.end(); ++it ) {
		charge_atom_data.push_back(it->serialize());
	}

	Pair charge_atom_record("atoms",charge_atom_data);
	return Value(utility::tools::make_vector(charge_atom_record,base_data));

}

void ChargeGrid::deserialize(utility::json_spirit::mObject data)
{

	charge_atom_list_.clear();
	utility::json_spirit::mArray charge_atom_data(data["atoms"].get_array());
	for ( utility::json_spirit::mArray::iterator it = charge_atom_data.begin(); it != charge_atom_data.end(); ++it ) {
		ChargeAtom current_atom;
		current_atom.deserialize(it->get_obj());
		charge_atom_list_.push_back(current_atom);
	}

	SingleGrid::deserialize(data["base_data"].get_obj());

}

void ChargeGrid::refresh(core::pose::Pose const & pose, core::Vector const & )
{

	setup_charge_atoms(pose);

	numeric::xyzVector<core::Size> dimensions = get_dimensions();
	for ( core::Size x_index =0; x_index < dimensions.x(); ++x_index ) {
		for ( core::Size y_index = 0; y_index < dimensions.y(); ++y_index ) {
			for ( core::Size z_index = 0; z_index < dimensions.z(); ++z_index ) {

				core::Vector pdb_coords(get_pdb_coords(x_index,y_index,z_index));
				core::Size neighbor_count = 0;
				//first loop through the charge atoms to calculate neighbor count
				//Maybe we should use a b-tree to do atom count calculation but it might not be worth the effort
				std::list<ChargeAtom>::iterator protein_charge_atom_it = charge_atom_list_.begin();
				for ( ; protein_charge_atom_it != charge_atom_list_.end(); ++protein_charge_atom_it ) {
					if ( protein_charge_atom_it->xyz.distance(pdb_coords) <= 4.0 ) {
						++neighbor_count;
					}
					if ( neighbor_count >= 12 ) {
						break;
					}
				}

				//dummy charge value is ignored on the grid atom side during
				//el calculations.  fix this later
				ChargeAtom grid_charge_atom(pdb_coords,0.0,neighbor_count);

				//second loop through the charge atoms to calculate charge
				protein_charge_atom_it = charge_atom_list_.begin();
				core::Real total_el = 0.0;
				for ( ; protein_charge_atom_it != charge_atom_list_.end(); ++protein_charge_atom_it ) {
					total_el += charge_charge_electrostatic(pose,*protein_charge_atom_it,grid_charge_atom);
				}

				set_point(pdb_coords,total_el);

			}
		}
	}
}

void ChargeGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, core::Size const & )
{
	refresh(pose,center);
}

void ChargeGrid::refresh(core::pose::Pose const & pose, core::Vector const & center, utility::vector1<core::Size> )
{
	refresh(pose,center);
}


void ChargeGrid::parse_my_tag(utility::tag::TagCOP const /*tag*/)
{

}


core::Real ChargeGrid::score(core::conformation::UltraLightResidue const & residue, core::Real const max_score, qsarMapOP)
{
	core::Real score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue[atom_index]);
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::Real protein_charge = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());


			score += protein_charge*residue.residue()->atomic_charge(atom_index);
		}
	}
	return score;
}

core::Real ChargeGrid::atom_score(core::conformation::UltraLightResidue const & residue, core::Size atomno, qsarMapOP)
{
	core::Vector const & atom_coord(residue[atomno]);
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::Real protein_charge = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
		return protein_charge*residue.residue()->atomic_charge(atomno);
	} else {
		return 0;
	}
}

core::Real ChargeGrid::score(core::conformation::Residue const & residue, core::Real const max_score, qsarMapOP /*qsar_map*/)
{
	core::Real score = 0.0;
	for ( core::Size atom_index = 1; atom_index <= residue.natoms() && score < max_score; ++atom_index ) {
		core::Vector const & atom_coord(residue.xyz(atom_index));
		if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
			core::Real protein_charge = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());



			score += protein_charge*residue.atomic_charge(atom_index);
		}
	}
	return score;
}

core::Real ChargeGrid::atom_score(core::conformation::Residue const & residue, core::Size atomno, qsarMapOP /*qsar_map*/)
{
	core::Vector const & atom_coord(residue.xyz(atomno));
	if ( this->get_grid().is_in_grid(atom_coord.x(),atom_coord.y(),atom_coord.z()) ) {
		core::Real protein_charge = this->get_point(atom_coord.x(),atom_coord.y(),atom_coord.z());
		return protein_charge*residue.atomic_charge(atomno);
	} else {
		return 0;
	}
}

core::Real ChargeGrid::nominal_depth(core::Size const & n_atoms) const
{
	switch(n_atoms)
			{
			case 7 :
				return(0.4);
			case 8 :
				return(0.9);
			case 9 :
				return(1.4);
			case 10 :
				return(1.9);
			case 11 :
				return(2.6);
			default :
				if ( n_atoms <= 6 ) {
					return(0.0);
				} else if ( n_atoms >= 12 ) {
					return(4.0);
				} else {
					utility_exit_with_message("This should never have happened");
				}
				break;
			}
	return(-1.0); //Something has gone rather wrong
}

core::Real ChargeGrid::charge_charge_electrostatic(core::pose::Pose const & /*pose*/, ChargeAtom const & atom_q, ChargeAtom const & atom_p) const
{
	core::Real distance = atom_q.xyz.distance(atom_p.xyz);
	core::Real s_p(nominal_depth(atom_p.neighbor_count));
	core::Real s_q(nominal_depth(atom_q.neighbor_count));

	assert(s_p >= 0);
	assert(s_q >= 0);

	core::Real direct_term =  ((atom_q.charge)/(4*boost::math::constants::pi<core::Real>()*epsilon_0_*distance) *zeta_);
	core::Real indirect_term = (1.0/distance) + (indirect_numerator_ / std::sqrt(distance*distance + 4*s_p*s_q));
	return direct_term*indirect_term;

}


void ChargeGrid::setup_charge_atoms(core::pose::Pose const & pose)
{
	core::Size chain_id = core::pose::get_chain_id_from_chain(get_chain(),pose);
	core::Size chain_begin = pose.conformation().chain_begin(chain_id);
	core::Size chain_end = pose.conformation().chain_end(chain_id);

	for ( core::Size current_residue_index = chain_begin; current_residue_index <= chain_end; ++current_residue_index ) {
		core::conformation::Residue current_residue = pose.residue(current_residue_index);
		for ( core::Size current_atom_index = 1; current_atom_index <= current_residue.natoms(); ++current_atom_index ) {
			core::id::AtomID current_atom_id(current_atom_index,current_residue_index);
			core::Vector current_atom_coords(pose.xyz(current_atom_id));
			core::Real current_atom_charge = current_residue.atomic_charge(current_atom_index);
			core::Size current_atom_neighbor_count = 0;

			for ( core::Size neighbor_residue_index = chain_begin; neighbor_residue_index <= chain_end; ++neighbor_residue_index ) {
				//calculate neighbors within 12 A.  stop counting after 12 since thats where the lookup table cuts off.
				core::conformation::Residue neighbor_residue = pose.residue(neighbor_residue_index);
				for ( core::Size neighbor_atom_index = 1; neighbor_atom_index <= neighbor_residue.natoms(); ++neighbor_atom_index ) {
					core::id::AtomID neighbor_atom_id(neighbor_atom_index,neighbor_residue_index);
					core::Vector neighbor_atom_coords(pose.xyz(neighbor_atom_id));
					if ( neighbor_atom_coords.distance(neighbor_atom_coords) <= 4.0 ) {
						current_atom_neighbor_count++;
					}

					if ( current_atom_neighbor_count >= 12 ) {
						break;
					}
				}
				charge_atom_list_.push_back(ChargeAtom(current_atom_coords,current_atom_charge,current_atom_neighbor_count));
			}
		}
	}
}

}
}
}

