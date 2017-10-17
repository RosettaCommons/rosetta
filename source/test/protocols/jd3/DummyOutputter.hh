// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/jd3/StandardJobQueenOverrideOutputter.cxxtest.hh
/// @brief  test suite for the StandardJobQueen testing when derived JQs intend to use
///         a limited subset of PoseOutputters.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/standard/StandardJobQueen.hh>

// Package headers
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/standard/StandardInnerLarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/jd3/JobOutputIndex.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_inputters/PoseInputSource.hh>
#include <protocols/jd3/deallocation/InputPoseDeallocationMessage.hh>
#include <protocols/jd3/output/MultipleOutputSpecification.hh>
#include <protocols/jd3/output/MultipleOutputter.hh>
#include <protocols/jd3/output/OutputSpecification.hh>
#include <protocols/jd3/output/ResultOutputter.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/pose_outputters/PoseOutputSpecification.hh>
#include <protocols/jd3/pose_outputters/PoseOutputter.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterCreator.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/PDBPoseOutputter.hh>

#include <core/pose/Pose.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/options/keys/OptionKeyList.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


using namespace protocols::jd3;
using namespace protocols::jd3::output;
using namespace protocols::jd3::pose_outputters;
using namespace utility::tag;

namespace basic { namespace options { namespace OptionKeys {
extern basic::options::BooleanOptionKey const dummy_outputter_arg;
extern basic::options::BooleanOptionKey const dummy_outputter;
}}}


class DummyOutputSpecification : public PoseOutputSpecification
{
public:
	inline DummyOutputSpecification();
	inline DummyOutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index );

	std::string out_fname() const { return out_fname_; }
	void out_fname( std::string const & setting ) { out_fname_ = setting; }
private:
	std::string out_fname_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
DummyOutputSpecification::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( CEREAL_NVP( out_fname_ ) ); // std::string
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
DummyOutputSpecification::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::pose_outputters::PoseOutputSpecification >( this ) );
	arc( out_fname_ ); // std::string
}

SAVE_AND_LOAD_SERIALIZABLE( DummyOutputSpecification );
CEREAL_REGISTER_TYPE( DummyOutputSpecification )

// put this in only one place CEREAL_REGISTER_DYNAMIC_INIT( DummyOutputSpecification )
#endif // SERIALIZATION

typedef utility::pointer::shared_ptr< DummyOutputSpecification > DummyOutputSpecificationOP;

////
////

class DummyPoseOutputter : public PoseOutputter
{
public:

	DummyPoseOutputter() :
		PoseOutputter()
	{}

	DummyPoseOutputter( std::shared_ptr< std::map< JobResultID, JobResultOP > > results ) :
		PoseOutputter(),
		results_( results )
	{}

	void determine_job_tag(
		utility::tag::TagCOP,
		utility::options::OptionCollection const &,
		InnerLarvalJob & job
	) const override {
		job.job_tag( job.input_tag() );
	}

	std::string
	outputter_for_job(
		utility::tag::TagCOP,
		utility::options::OptionCollection const &,
		InnerLarvalJob const &
	) const override {
		return "dummy";
	}

	std::string
	outputter_for_job(
		PoseOutputSpecification const &
	) const override
	{
		return "dummy";
	}

	bool job_has_already_completed( LarvalJob const &, utility::options::OptionCollection const & ) const override { return false; }

	void mark_job_as_having_started( LarvalJob const &, utility::options::OptionCollection const & ) const override {}

	PoseOutputSpecificationOP
	create_output_specification(
		LarvalJob const & job,
		JobOutputIndex const & output_index,
		utility::options::OptionCollection const &,
		utility::tag::TagCOP // possibly null-pointing tag pointer
	) override
	{
		utility::pointer::shared_ptr< DummyOutputSpecification > spec;
		spec.reset( new DummyOutputSpecification );
		spec->out_fname( job.job_tag_with_index_suffix( output_index ) );
		return spec;
	}


	void write_output(
		OutputSpecification const & spec,
		JobResult const & result
	) override {
		using namespace core::pose;
		using standard::PoseJobResult;
		debug_assert( dynamic_cast< PoseJobResult const * > ( &result ));
		PoseJobResult const & pose_result( static_cast< PoseJobResult const & > ( result ));
		core::pose::Pose const & pose( *pose_result.pose() );

		DummyOutputSpecification const & dummy_spec( dynamic_cast< DummyOutputSpecification const & > (spec) );
		std::string name = dummy_spec.out_fname();
		//std::cout << "Saving Pose with name: " << name << " " << std::endl;
		//std::cout << output_index.secondary_output_index << " " << output_index.n_secondary_outputs << std::endl;
		named_poses_[ name ] = PoseOP( new Pose( pose ) );

		if ( results_ ) {
			(*results_)[ spec.result_id() ] = JobResultOP( new PoseJobResult(
				static_cast< PoseJobResult const & > (result) ));
		}

	}

	void flush() override {}

	static std::string keyname() { return "Dummy"; }

	std::string
	class_key() const override { return keyname(); }

public:

	std::shared_ptr< std::map< JobResultID, JobResultOP > > results_;
	std::map< std::string, core::pose::PoseOP > named_poses_;

};


typedef utility::pointer::shared_ptr< DummyPoseOutputter > DummyPoseOutputterOP;

DummyOutputSpecification::DummyOutputSpecification() {
	outputter_type( DummyPoseOutputter::keyname() );
}

DummyOutputSpecification::DummyOutputSpecification( JobResultID const & result_id, JobOutputIndex const & output_index ) :
	PoseOutputSpecification( result_id, output_index )
{
	outputter_type( DummyPoseOutputter::keyname() );
}

// An unregistered outputter
class DummyOutputterCreator : public PoseOutputterCreator
{
public:

	PoseOutputterOP create_outputter() const override { return PoseOutputterOP( new DummyPoseOutputter ); }
	std::string keyname() const override { return DummyPoseOutputter::keyname(); }
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override {
		AttributeList attrs;
		attrs + XMLSchemaAttribute( "dummy_attribute", xs_string, "A dummy attribute" );
		pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(), "testing 123", attrs );
	}
	void list_options_read( utility::options::OptionKeyList & read_options ) const override
	{
		read_options + basic::options::OptionKeys::dummy_outputter_arg;
	}

	bool outputter_specified_by_command_line() const override {
		return basic::options::option[ basic::options::OptionKeys::dummy_outputter ];
	}


};
