// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/benchmark/SampleTiltAngles.hh
/// @brief Calculates the energy at all possible tilt angles (0->180 degrees)
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_SampleTiltAngles_hh
#define INCLUDED_protocols_membrane_benchmark_SampleTiltAngles_hh

// Unit headers
#include <protocols/membrane/benchmark/SampleTiltAngles.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers
#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <cstdlib>

namespace protocols {
namespace membrane {
namespace benchmark {

///@brief Calculates the energy at all possible tilt angles (0->180 degrees)
class SampleTiltAngles : public protocols::moves::Mover {

public:

	/// @brief Sample all tilt angles and compare energies using
	/// talaris (2014), and high-resolution membrane energy functions from
	/// 2007 and 2012
	SampleTiltAngles();

	/// @brief Sample all tilt angles and compare with user specified energy functions
	/// and a user specified output prefix for files
	SampleTiltAngles(
		std::string prefix,
		core::scoring::ScoreFunctionOP ref_sfxn1,
		core::scoring::ScoreFunctionOP ref_sfxn2,
		core::scoring::ScoreFunctionOP ref_sfxn3
	);

	// copy constructor
	SampleTiltAngles( SampleTiltAngles const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~SampleTiltAngles();

	void
	apply( core::pose::Pose & pose ) override;

public:

	void
	show( std::ostream & output=std::cout ) const override;

	// XRW TEMP  std::string
	// XRW TEMP  get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;


	/// @brief required in the context of the parser/scripting scheme
	moves::MoverOP
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


private:

	void
	write_score_to_outfiles(
		utility::vector1< core::Real > angles,
		utility::vector1< core::Real > ref_sfxn1,
		utility::vector1< core::Real > ref_sfxn2,
		utility::vector1< core::Real > ref_sfxn3
	);

private:

	// Prefix for output files
	std::string prefix_;

	// Energy Functions (might add more later, depending on project needs)
	core::scoring::ScoreFunctionOP ref_sfxn1_;
	core::scoring::ScoreFunctionOP ref_sfxn2_;
	core::scoring::ScoreFunctionOP ref_sfxn3_;

};

} // benchmark
} // membrane
} // protocols


#endif // INCLUDED_protocols_membrane_benchmark_SampleTiltAngles_hh
