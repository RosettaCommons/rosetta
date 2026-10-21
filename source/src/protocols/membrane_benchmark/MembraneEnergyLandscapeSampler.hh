// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.hh
/// @brief Sample all points on a 2D membrane energy landscape given implicit model and protein dimensions
/// @author Rebecca Alford (rfalford12@gmail.com)
///edited by Rituparna Samanta (rituparna@utexas.edu)

#ifndef INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_HH
#define INCLUDED_protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_HH

// Unit headers
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/minimization_packing/PackRotamersMover.fwd.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/util.hh>
#include <core/energy_methods/pHEnergy.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/membrane/TranslationRotationMover.fwd.hh>


// Core headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/io/ozstream.hh>

namespace protocols {
namespace membrane_benchmark {

///@brief Sample all points on a 2D membrane energy landscape given implicit model and protein dimensions
class MembraneEnergyLandscapeSampler : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	MembraneEnergyLandscapeSampler();

	/// @brief Non-defualt cnstructor for landscape sampling
	MembraneEnergyLandscapeSampler(
		std::string sfxn_weights,
		std::string rotation_type,
		bool interface,
		core::Real start_z,
		core::Real end_z,
		core::Real flag_axis,
		core::Real azimuthal_delta,
		bool repack
	);

	/// @brief Copy constructor (not needed unless you need deep copies)
	MembraneEnergyLandscapeSampler( MembraneEnergyLandscapeSampler const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~MembraneEnergyLandscapeSampler() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output ) const override;

	/// @brief Declaring an overload
	void
	show() const;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	//MembraneEnergyLandscapeSampler & operator=( MembraneEnergyLandscapeSampler const & src );

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

	core::Vector
	getaxis(
		core::pose::Pose & pose,
		core::Real flag
	);
	core::Vector
	getcenter(
		core::pose::Pose & pose,
		core::Real flag
	);



private: // methods

	void
	init_from_options();

	//protocols::minimization_packing::PackRotamersMoverOP
	//get_packer( core::pose::Pose const & pose );

	protocols::minimization_packing::PackRotamersMoverOP
	get_pH_aware_packer( core::scoring::ScoreFunctionOP sfxn ) const;

	core::Vector
	get_rotation_axis();

	protocols::membrane::RotationMoverOP
	get_rotation(
		core::Real const normal_angle,
		core::Vector const axis,
		core::Vector const rot_center,
		core::Size membrane_jump,
		std::string rotation_type
	);

	core::Real get_pH_value();

	core::Real
	count_res( std::string res,
		core::pose::Pose const & pose);

	core::Real
	count_diff( core::pose::Pose const & pose1,
		core::pose::Pose const & pose2);

private: // data

	core::scoring::ScoreFunctionOP sfxn_;
	std::string sfxn_weights_;
	std::string rotation_type_="YZ";
	// Landscape generation options
	bool interface_=false;
	core::Real start_z_=-60;
	core::Real end_z_=60;
	core::Real flag_axis_=1.0;
	core::Real azimuthal_delta_=30;
	bool repack_=false;
	bool pH_mode_=false;
	bool dump_structures_=false;
};

std::ostream &
operator<<( std::ostream & os, MembraneEnergyLandscapeSampler const & mover );

} //protocols
} //membrane_benchmark

#endif //protocols_membrane_benchmark_MembraneEnergyLandscapeSampler_HH
