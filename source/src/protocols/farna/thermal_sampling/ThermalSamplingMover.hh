// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/thermal_sampling/ThermalSamplingMover.hh
/// @brief Use a simulated tempering simulation to refine a pose
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_thermal_sampling_ThermalSamplingMover_hh
#define INCLUDED_protocols_farna_thermal_sampling_ThermalSamplingMover_hh

// Unit headers
#include <protocols/farna/thermal_sampling/ThermalSamplingMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>

// Basic/Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace farna {
namespace thermal_sampling {

///@brief Use a simulated tempering simulation to refine a pose
class ThermalSamplingMover : public protocols::moves::Mover {

	typedef core::Size Size;
	typedef core::Real Real;

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	ThermalSamplingMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	ThermalSamplingMover( ThermalSamplingMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	virtual ~ThermalSamplingMover();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	virtual void
	apply( core::pose::Pose & pose );

	/// @brief Show the contents of the Mover
	static std::string
	class_name();

	virtual void
	show( std::ostream & output = std::cout ) const;

	/// @brief Get the name of the Mover
	virtual std::string
	get_name() const;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//ThermalSamplingMover & operator=( ThermalSamplingMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	void set_residues( utility::vector1< Size > const & residues ) { residues_ = residues; }
	utility::vector1< Size > residues() const { return residues_; }
	void set_free_rsd( utility::vector1< Size > const & free_rsd ) { free_rsd_ = free_rsd; }
	utility::vector1< Size > free_rsd() const { return free_rsd_; }

	void set_recces_turner_mode( bool const setting ) { recces_turner_mode_ = setting; }
	void set_dumping_app( bool const setting ) { dumping_app_ = setting; }
	void set_n_cycle( Size const setting ) { n_cycle_ = setting; }
	void set_dump_silent( bool const setting ) { dump_pdb_ = setting; }
	void set_dump_pdb( bool const setting ) { dump_silent_ = setting; }
	void set_temps( utility::vector1< Size > const & temps ) { temps_ = temps; }
	void set_weights( utility::vector1< Size > const & weights ) { st_weights_ = weights; }

	void set_residue_sampling_from_pose_and_movemap( core::pose::Pose const & pose, core::kinematics::MoveMap const & mm );

private: // methods

private: // data
	utility::vector1< Size > residues_;
	Size total_sampled_;
	utility::vector1< Size > free_rsd_;
	bool recces_turner_mode_ = false;
	bool dumping_app_ = false;
	Size n_cycle_;
	bool dump_pdb_ = false;
	bool dump_silent_ = false;
	Real angle_range_chi_;
	Real angle_range_bb_;
	utility::vector1< Real > angle_range_bb_varying_;
	utility::vector1< Real > temps_;
	utility::vector1< Real > st_weights_;

	core::kinematics::MoveMapOP mm_;
};

std::ostream &
operator<<( std::ostream & os, ThermalSamplingMover const & mover );

} //protocols
} //farna
} //thermal_sampling

#endif //protocols/farna/thermal_sampling_ThermalSamplingMover_hh
