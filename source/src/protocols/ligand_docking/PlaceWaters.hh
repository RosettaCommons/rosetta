// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/PlaceWaters.hh
/// @brief Places waters in the interface
/// @author Laura Shub (laurashub@gmail.com)

#ifndef INCLUDED_protocols_ligand_docking_PlaceWaters_HH
#define INCLUDED_protocols_ligand_docking_PlaceWaters_HH

// Unit headers
#include <protocols/ligand_docking/PlaceWaters.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>
#include <core/grid/CartGrid.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/UltraLightResidue.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace ligand_docking {

struct water_info {

public:

	core::Vector centroid;
	utility::vector1<core::Vector> ligand_acceptors;
	utility::vector1<core::Vector> ligand_donors;
	utility::vector1<core::Vector> protein_acceptors;
	utility::vector1<core::Vector> protein_donors;
	utility::vector1<core::Vector> other_waters;

};

///@brief Places waters in the interface
class PlaceWaters: public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PlaceWaters();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PlaceWaters(PlaceWaters const & src);

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PlaceWaters()
	override ;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &
	) override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	//Additional methods

	/// @brief determines size of grids
	core::Real grid_size(core::conformation::UltraLightResidueOP lig);

	utility::vector1<water_info>
	check_possible_hbonds(utility::vector1<core::Vector> const & vs, core::Real inner_radius, core::Real outer_radius);

	utility::vector1<core::Vector>
	empty_spaces(core::Size cutoff);

	void
	fill_occupied(core::pose::Pose const & pose);

	void
	detect_hydrogens(core::pose::Pose const & pose);

	void
	set_sphere(
		core::grid::CartGrid<core::Real> & grid,
		core::Vector const & coords,
		core::Real radius,
		core::Real value
	);

	void
	set_ring(
		core::grid::CartGrid<core::Real> & grid,
		core::Vector const & coords,
		core::Real inner_radius,
		core::Real outer_radius,
		core::Real value
	);

	bool
	check_sphere(
		core::grid::CartGrid<core::Real> const & grid,
		core::Vector const & coords,
		core::Real radius,
		core::Real value,
		core::Vector & hit
	);

	bool
	check_ring(
		core::grid::CartGrid<core::Real> const & grid,
		core::Vector const & coords,
		core::Real inner_radius,
		core::Real outer_radius,
		core::Real value,
		core::Vector & hit
	);

	void
	find_nearby_hbonds( core::grid::CartGrid<core::Real> const & grid, core::Vector const coords, water_info & wi, core::Real inner_radius, core::Real outer_radius);

	void
	cluster_groups(utility::vector1<core::Vector> & vs, core::Real group_cutoff, core::Real merge_cutoff);

	void
	check_other_waters( core::Vector coords, utility::vector1<core::Vector> const & waters, water_info & wi);

	void
	merge_close(utility::vector1<core::Vector> & centers, core::Real merge_cutoff);

	//lots of fetch functions for all the private variables
	core::Real
	get_grid_size();

	core::Real
	get_chain();

	std::string
	get_water();

	core::Real
	get_resolution();

	core::Real
	get_cutoff();

	core::Real
	get_group_cutoff();

	core::Real
	get_merge_cutoff();

	core::Real
	get_hbonds_inner();

	core::Real
	get_hbonds_outer();

private:// methods

	void
	transform_water(
		core::conformation::Residue & water,
		water_info & wi,
		core::pose::Pose & pose);

	core::Size
	get_closest(
		utility::vector1<core::Vector> const & means,
		core::Vector const & v,
		core::Real group_cutoff);

	core::Vector
	recalculate_mean(utility::vector1<core::Vector> col);

	bool
	prot_lig_hbonds(water_info const & wi);

	core::Size
	too_close(utility::vector1<core::Vector> & centers,
		core::Vector & a,
		core::Vector & b,
		core::Real cutoff);

	void
	initialize_grids(
		core::Vector const & center,
		core::Real width,core::Real resolution);

	void
	output_water_coords_file(
		utility::vector1<core::Vector> const & waters,
		std::string outfile);

	void
	print_contacts(water_info const & wi);

	void
	update_info(
		core::Vector box_center,
		core::Real val,
		water_info & wi);

private:// data

	core::grid::CartGrid<core::Real> void_grid_;
	core::grid::CartGrid<core::Real> h_grid_;
	core::Real grid_size_;
	std::string ligand_chain_;
	char water_chain_;
	std::string water_type_;
	core::Real resolution_;
	core::Real cutoff_;
	core::Size cycles_;
	core::Size group_cutoff_;
	core::Real merge_cutoff_;
	core::Real hbonds_inner_;
	core::Real hbonds_outer_;
	bool output_water_coords_;
};

std::ostream &
operator<<(std::ostream & os, PlaceWaters const & mover);

} //protocols
} //ligand_docking

#endif //protocols_ligand_docking_PlaceWaters_HH
