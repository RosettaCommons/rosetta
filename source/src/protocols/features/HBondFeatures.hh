// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/HBondFeatures.hh
/// @brief  report HBond geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_features_HBondFeatures_hh
#define INCLUDED_protocols_features_HBondFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/HBondFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

/// @brief How a hydrogen bond should be defined See note for
/// definition_type_ below.
enum HBDefType {
	hbdef_NONE = 1,
	hbdef_ENERGY, // define what constitutes a hydrogen bond by its energy
	hbdef_AHDIST, // define what constitutes a hydrogen bond by the A-H distance
	hbdef_MAX = hbdef_AHDIST
};



class HBondFeatures : public protocols::features::FeaturesReporter {
public:
	HBondFeatures();

	HBondFeatures(
		core::scoring::ScoreFunctionOP scfxn);

	HBondFeatures( HBondFeatures const & src );

	virtual ~HBondFeatures();

	///@brief return string with class name
	std::string
	type_name() const;

	///@brief generate the table schemas and write them to the database
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	void
	write_hbond_chem_types_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_sites_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_sites_pdb_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_site_environment_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_site_atoms_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbonds_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_lennard_jones_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_geom_coords_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_dehydrons_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief get what criteria should be used to define what
	///constitutes a hydrogen bond
	// Undefined, commenting out to fix PyRosetta build  HBDefType definition_type() const;

	///@brief set what criteria should be used to define what
	///constitutes a hydrogen bond
	// Undefined, commenting out to fix PyRosetta build  void definition_type( HBDefType definition_type );


	///@brief get the definition threshold that should be
	///used to define what constitutes a hydrogen bond
	// Undefined, commenting out to fix PyRosetta build  core::Real definition_threshold() const;

	///@brief set the definition threshold that should be
	///used to define what constitutes a hydrogen bond
	// Undefined, commenting out to fix PyRosetta build  void definition_threshold( core::Real definition_threshold );

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_row(
		core::pose::Pose const & pose,
		StructureID struct_id,
		core::Size site_id,
		core::Size resNum,
		core::Size atmNum,
		bool is_donor,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_pdb_row(
		core::pose::Pose const & pose,
		core::Size resNum,
		core::Size atmNum,
		core::Size heavy_atmNum,
		StructureID struct_id,
		core::Size site_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_environment_row(
		core::pose::Pose const & pose,
		core::Size resNum,
		core::Size atmNum,
		StructureID struct_id,
		core::Size site_id,
		core::id::AtomID_Map< core::Real > const & atom_sasa_s,
		core::id::AtomID_Map< core::Real > const & atom_sasa_m,
		core::id::AtomID_Map< core::Real > const & atom_sasa_l,
		core::id::AtomID_Map< utility::vector1< core::scoring::hbonds::HBondCOP > > const & site_partners,
		core::id::AtomID_Map< core::Real > const & site_hbond_energies,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_atoms_row(
		core::pose::Pose const & pose,
		core::Size resNum,
		core::Size atmNum,
		StructureID struct_id,
		core::Size site_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_row(
		core::scoring::hbonds::HBond const & hbond,
		StructureID struct_id,
		core::Size hbond_id,
		core::id::AtomID_Map< core::Size > const & site_ids,
		core::id::AtomID_Map< utility::vector1< core::scoring::hbonds::HBondCOP > > const & site_partners,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_geom_coords(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBondOptions const & hbond_options,
		core::scoring::hbonds::HBond const & hbond,
		StructureID struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_lennard_jones_row(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBond const & hbond,
		StructureID struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_dehydron_row(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBond const & hbond,
		StructureID struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);


private:

	core::scoring::ScoreFunctionOP scfxn_;

	// The standard Rosetta definition of a hydrogen bond is when the
	// hbond energy < 0.  However, to detect the distribution around
	// this threhsold it can be useful to define hydrogen bonds
	// differently, either by a higher energy thershold or using some
	// other threshold, like the distance between the acceptor and the
	// hydrogen.
	HBDefType definition_type_;
	core::Real definition_threshold_;
};

} // namespace
} // namespace

#endif // include guard
