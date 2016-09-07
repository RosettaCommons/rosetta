// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   PackStatMover.hh
///
/// @brief
/// @author Monica Berrondo

#ifndef INCLUDED_protocols_analysis_PackStatMover_hh
#define INCLUDED_protocols_analysis_PackStatMover_hh


#include <core/types.hh>

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

#include <protocols/analysis/PackStatMover.fwd.hh>

namespace protocols {
namespace analysis {

class PackStatMover : public moves::Mover
{
public:

	PackStatMover();

	void apply(
		core::pose::Pose & pose
	) override;
	std::string get_name() const override;

	// void set_verbose( bool _verbose ) { verbose_ = _verbose; }
	// bool get_verbose( ) { return verbose_; }
	//
	// void set_packstat_pdb( bool _packstat_pdb ) { packstat_pdb_ = _packstat_pdb; }
	// bool get_packstat_pdb( ) { return packstat_pdb_; }
	//
	// void set_include_water( bool _include_water ) { include_water_ = _include_water; }
	// bool get_include_water( ) { return include_water_; }
	//
	// void set_surface_accessibility( bool _surface_accessibility ) { surface_accessibility_ = _surface_accessibility; }
	// bool get_surface_accessibility( ) { return surface_accessibility_; }
	//
	// void set_residue_scores( bool _residue_scores ) { residue_scores_ = _residue_scores; }
	// bool get_residue_scores( ) { return residue_scores_; }
	//
	// void set_cavity_burial_probe_radius( bool _cavity_burial_probe_radius ) { cavity_burial_probe_radius_ = _cavity_burial_probe_radius; }
	// bool get_cavity_burial_probe_radius( ) { return cavity_burial_probe_radius_; }
	//
	// void set_oversample( core::Real _oversample ) { oversample_ = _oversample; }
	// core::Real get_oversample( ) { return oversample_; }

private:

	/// information about the mode
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool verbose_;
	bool packstat_pdb_;
	bool include_water_;
	bool surface_accessibility_;
	bool residue_scores_;
	core::Real cavity_burial_probe_radius_;
	int oversample_;

};

} // analysis
} // protocols

#endif // INCLUDED_protocols_analysis_PackStatMover_hh

