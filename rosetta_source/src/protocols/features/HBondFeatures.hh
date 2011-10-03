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
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_HBondFeatures_hh
#define INCLUDED_protocols_features_HBondFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/HBondFeatures.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.fwd.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Utility Headers
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

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

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_row(
		core::pose::Pose const & pose,
		core::Size struct_id,
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
		core::Size struct_id,
		core::Size site_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_site_environment_row(
		core::pose::Pose const & pose,
		core::Size resNum,
		core::Size atmNum,
		core::Size struct_id,
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
		core::Size struct_id,
		core::Size site_id,
		utility::sql_database::sessionOP db_session);

	core::Size
	insert_hbond_row(
		core::scoring::hbonds::HBond const & hbond,
		core::Size struct_id,
		core::Size hbond_id,
		core::id::AtomID_Map< core::Size > const & site_ids,
		core::id::AtomID_Map< utility::vector1< core::scoring::hbonds::HBondCOP > > const & site_partners,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_geom_coords(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBond const & hbond,
		core::Size struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_lennard_jones_row(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBond const & hbond,
		core::Size struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);

	void
	insert_hbond_dehydron_row(
		core::pose::Pose const & pose,
		core::scoring::hbonds::HBond const & hbond,
		core::Size struct_id,
		core::Size hbond_id,
		utility::sql_database::sessionOP db_session);


private:

	core::scoring::ScoreFunctionOP scfxn_;

};

} // namespace
} // namespace

#endif // include guard
