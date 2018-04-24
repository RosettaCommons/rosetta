// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetric.hh
/// @brief Headers for a string metric that generates a string of PyMol commands to colour a
/// structure's polar groups based on burial.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_HH
#define INCLUDED_protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_HH

#include <protocols/analysis/burial_metrics/PolarGroupBurialPyMolStringMetric.fwd.hh>
#include <core/simple_metrics/StringMetric.hh>

// Core headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace analysis {
namespace burial_metrics {

///@brief A string metric that generates a string of PyMol commands to colour a structure's polar groups based on burial.
class PolarGroupBurialPyMolStringMetric : public core::simple_metrics::StringMetric{

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PolarGroupBurialPyMolStringMetric();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PolarGroupBurialPyMolStringMetric( PolarGroupBurialPyMolStringMetric const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PolarGroupBurialPyMolStringMetric() override;

	/////////////////////
	/// Metric Methods ///
	/////////////////////

	///Defined in RealMetric:
	///
	/// @brief Calculate the metric and add it to the pose as a score.
	///           labeled as prefix+metric+suffix.
	///
	/// @details Score is added through setExtraScorePose and is output
	///            into the score tables/file at pose output.
	//void
	//apply( pose::Pose & pose, prefix="", suffix="" ) override;

	///@brief Calculate the metric.
	std::string
	calculate( core::pose::Pose const & original_pose ) const override;

public:

	/// @brief Name of the class
	/// @details Calls name_static().
	std::string
	name() const override;

	/// @brief Name of the class for creator.
	static
	std::string
	name_static();

	/// @brief Name (descriptive) of the metric.
	std::string
	metric() const override;


public:

	/// @brief Given an XML tag, configure an instance of this class.
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data ) override;

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	core::simple_metrics::SimpleMetricOP
	clone() const override;

public:

	/// @brief Set the scorefunction.
	/// @details Used for definition of burial.  Input pointer is copied without cloning.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Set whether this metric should also output to tracer when evaluated.
	void set_verbose( bool const setting );

private:

	/// @brief A scorefunction, from which configurations will be copied.
	core::scoring::ScoreFunctionCOP sfxn_;

	/// @brief Should this metric also output to tracer when evaluated?  Default false.
	bool verbose_;

};

} //protocols
} //analysis
} //burial_metrics



#endif //protocols_analysis_burial_metrics_PolarGroupBurialPyMolStringMetric_HH





