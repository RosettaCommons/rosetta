// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/SampleGrid.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_magnesium_SampleGrid_HH
#define INCLUDED_protocols_magnesium_SampleGrid_HH

#include <protocols/moves/Mover.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray3D.hh>

namespace protocols {
namespace magnesium {

class SampleGrid: public utility::pointer::ReferenceCount {

public:

	//constructor
	SampleGrid( core::pose::Pose const & pose );

	//destructor
	~SampleGrid();

public:

	utility::vector1< core::Vector > get_mg_positions( core::pose::Pose const & pose );

	void set_input_scan_res( utility::vector1< core::Size > const & setting ){ input_scan_res_ = setting; }
	utility::vector1< core::Size > input_scan_res() const { return input_scan_res_; }

	void set_tether_to_closest_res( bool const & setting ){ tether_to_closest_res_ = setting; }
	bool tether_to_closest_res() const { return tether_to_closest_res_; }

	void set_xyz_step( core::Real const & setting ){ xyz_step_ = setting; }
	core::Real xyz_step() const { return xyz_step_; }

	core::Real xmin() const { return xmin_; }
	core::Real xmax() const { return xmax_; }
	core::Real ymin() const { return ymin_; }
	core::Real ymax() const { return ymax_; }
	core::Real zmin() const { return zmin_; }
	core::Real zmax() const { return zmax_; }

private:

	void
	figure_out_box_bounds( core::pose::Pose const & pose );

	void
	create_grid();

	utility::vector1< Size >
	figure_out_scan_res( utility::vector1< Size > const & input_scan_res,
		core::pose::Pose const & pose );

	void
	define_bins( core::Real const x,
		core::Real const subgrid_radius,
		core::Real const xmin,
		Size const xgridsize,
		core::Real const xyz_increment,
		Size & xbinmin,
		Size & xbinmax ) const;

	core::Real
	get_position( Size const xbin, core::Real const xmin, core::Real const xyz_increment ) const;

private:

	bool tether_to_closest_res_;
	core::Real xyz_step_;
	core::Real xmax_, xmin_, ymax_, ymin_, zmax_, zmin_;
	utility::vector1< core::Size > input_scan_res_;
	utility::vector1< core::Size > scan_res_;

	ObjexxFCL::FArray3D< core::Real > min_distance_grid_; // contains closest distance to a pose acceptor atom.

};

} //magnesium
} //protocols

#endif
