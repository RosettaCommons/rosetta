// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/mmTFPoseOutputter.hh
/// @brief  Definition of the %mmTFPoseOutputter class
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) - PDBPoseOutputter bases of this class

#ifndef INCLUDED_protocols_jd3_pose_outputters_mmTFPoseOutputter_HH
#define INCLUDED_protocols_jd3_pose_outputters_mmTFPoseOutputter_HH

//unit headers
#include <protocols/jd3/pose_outputters/mmTFPoseOutputter.fwd.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>
//package headers
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace pose_outputters {

/// @brief The %mmTFPoseOutputter
class mmTFPoseOutputter : public PDBPoseOutputter
{
public:

	mmTFPoseOutputter();
	virtual ~mmTFPoseOutputter();

	static
	bool
	outputter_specified_by_command_line();

	/// @brief Create the PoseOutputSpecification for a particular job
	PoseOutputSpecificationOP
	create_output_specification(
		LarvalJob const & job,
		JobOutputIndex const & output_index,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) override;

	/// @brief Write a pose out to permanent storage (whatever that may be).
	void write_output(
		output::OutputSpecification const & specification,
		JobResult const & result
	) override;

	std::string
	class_key() const override;

	static
	std::string
	keyname();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

	std::string
	output_name(
		LarvalJob const & job,
		JobOutputIndex const & output_index,
		utility::options::OptionCollection const & options,
		utility::tag::TagCOP tag // possibly null-pointing tag pointer
	) const override;

private:



	// NOTE THERE IS NO PRIVATE DATA AND THERE SHOULD NOT BE:
	// The mmTFPoseOutputter should not accumulate state unless the behavior of
	// its outputter_for_job method changes. Currently, this method returns the empty string
	// signifying to the StandardJobQueen that it has no state, and therefore, it
	// is safe to use the same instance of the class in multiple threads. If that pledge
	// should have to change (because this class is given state in the form of private data)
	// then the outputter_for_job_method must be updated to respect the job-distributor
	// assigned suffix.

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace pose_outputters
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_pose_outputters_mmTFPoseOutputter )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_mmTFPoseOutputter_HH
