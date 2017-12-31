// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd3/DeNovoSilentFilePoseOutputter.cc
/// @brief  Definition of the %DeNovoSilentFilePoseOutputter class's methods
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), Andy Watkins (amw579@stanford.edu)

//unit headers
#include <protocols/jd3/pose_outputters/DeNovoSilentFilePoseOutputter.hh>
#include <protocols/jd3/pose_outputters/DeNovoSilentFilePoseOutputterCreator.hh>

//package headers
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/InnerLarvalJob.hh>
#include <protocols/jd3/pose_outputters/PoseOutputterFactory.hh>
#include <protocols/jd3/pose_outputters/SilentFilePoseOutputSpecification.hh>
#include <protocols/jd3/pose_outputters/pose_outputter_schemas.hh>
#include <protocols/jd3/chunk_library/MoverAndChunkLibraryJob.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>

// ObjexxFCL
#include <ObjexxFCL/string.functions.hh>

// basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
//#include <basic/options/keys/run.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/keys/OptionKeyList.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>


namespace protocols {
namespace jd3 {
namespace pose_outputters {

static basic::Tracer TR( "protocols.jd3.pose_outputters.DeNovoSilentFilePoseOutputter" );

using namespace core::io;
using namespace core::io::silent;

DeNovoSilentFilePoseOutputter::DeNovoSilentFilePoseOutputter() :
	buffer_limit_( 20 )
{}

DeNovoSilentFilePoseOutputter::~DeNovoSilentFilePoseOutputter() {}

bool
DeNovoSilentFilePoseOutputter::outputter_specified_by_command_line()
{
	return basic::options::option[ basic::options::OptionKeys::out::file::silent ].user();
}

void
DeNovoSilentFilePoseOutputter::determine_job_tag(
	utility::tag::TagCOP /*output_tag*/,
	utility::options::OptionCollection const & /*job_options*/,
	InnerLarvalJob & job
) const {
	job.job_tag( job.input_tag() );
}

/// @details If this is just being given an InnerLarvalJob, that means that the
/// outputter is at most a function of the input source--as there's a bijection
/// b/n ILJ and input source. So, let's just use the job_tag. What else can we do?
std::string
DeNovoSilentFilePoseOutputter::outputter_for_job(
	utility::tag::TagCOP sf_tag,
	utility::options::OptionCollection const & opts,
	InnerLarvalJob const &
) const
{
	using namespace basic::options::OptionKeys;
	if ( sf_tag ) {
		runtime_assert( sf_tag->hasOption( "filename" ));
		return sf_tag->getOption< std::string >( "filename" );
	} else {
		//return opts[ basic::options::OptionKeys::out::file::silent ];

		// For denovo, we actually want the current directory name as the 'real default'. Overwrite here.
		// This should be able to come from RNA_DeNovoProtocolOptions...
		if ( opts[ out::file::silent ].user() ) {
			return opts[ out::file::silent ]();
		} else {
			std::string tag;
			if ( opts[ rna::denovo::tag ].user() ) {
				tag = opts[ rna::denovo::tag ]();
			} else {
				tag = utility::file_basename( utility::file::cwd() );
				// Can't do this in two places or we get extraneous output. But so much safer to
				// do it in both places.
				//TR << TR.Green << "Setting silent file name based on directory: " << tag << ".out" << std::endl;
			}
			return tag  + ".out";
		}
	}
}

std::string
DeNovoSilentFilePoseOutputter::outputter_for_job(
	PoseOutputSpecification const & spec
) const {
	using SFPOS = SilentFilePoseOutputSpecification;
	debug_assert( dynamic_cast< SFPOS const * > ( &spec ) );
	auto const & sf_spec( static_cast< SFPOS const & > ( spec ) );

	return sf_spec.out_fname();
}

/// @brief Create the PoseOutputSpecification for a particular job
PoseOutputSpecificationOP
DeNovoSilentFilePoseOutputter::create_output_specification(
	LarvalJob const & job,
	JobOutputIndex const & output_index,
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP tag // possibly null-pointing tag pointer
) {
	using namespace core::io::silent;
	using namespace basic::options::OptionKeys;
	debug_assert( !tag || tag->getName() == keyname() ); // I expect this Tag to point to my data

	core::io::silent::SilentFileOptions opts( job_options );
	std::string fname_out;
	core::Size buffer_limit(0);
	if ( tag ) {
		opts.read_from_tag( tag );
		fname_out = tag->getOption< std::string >( "filename" );
		if ( tag->hasOption( "buffer_limit" ) ) {
			buffer_limit = tag->getOption< core::Size >( "buffer_limit" );
		}
	} else {
		using namespace basic::options::OptionKeys;
		// For denovo, we actually want the current directory name as the 'real default'. Overwrite here.
		// This should be able to come from RNA_DeNovoProtocolOptions...
		if ( job_options[ out::file::silent ].user() ) {
			fname_out = job_options[ out::file::silent ]();
		} else {
			std::string tag;
			if ( job_options[ rna::denovo::tag ].user() ) {
				tag = job_options[ rna::denovo::tag ]();
			} else {
				tag = utility::file_basename( utility::file::cwd() );
				TR << TR.Green << "Setting silent file name based on directory: " << tag << ".out" << std::endl;
			}
			fname_out = tag  + ".out";
		}
	}

	SilentFilePoseOutputSpecificationOP sf_pos( new SilentFilePoseOutputSpecification );
	sf_pos->sf_opts( opts );
	sf_pos->out_fname( fname_out );
	sf_pos->buffer_limit( buffer_limit );
	sf_pos->pose_tag( ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
		+ job.job_tag_with_index_suffix( output_index, 6 ) );
	dump_ = job_options[ rna::denovo::out::dump ].value();
	return sf_pos;
}

bool DeNovoSilentFilePoseOutputter::job_has_already_completed( LarvalJob const & /*job*/, utility::options::OptionCollection const & ) const
{
	return false;
}


void DeNovoSilentFilePoseOutputter::mark_job_as_having_started( LarvalJob const & /*job*/, utility::options::OptionCollection const & ) const
{
	// This is not a behavior supported by the DeNovoSilentFilePoseOutputter
}

std::string
DeNovoSilentFilePoseOutputter::class_key() const
{
	return keyname();
}

void DeNovoSilentFilePoseOutputter::write_output(
	output::OutputSpecification const & spec,
	JobResult const & result )
{
	using namespace core::pose::rna;

	using chunk_library::ChunkLibraryJobResult;
	debug_assert( dynamic_cast< ChunkLibraryJobResult const * > ( &result ));
	auto const & denovo_result( static_cast< ChunkLibraryJobResult const & > ( result ));
	core::pose::PoseOP pose( denovo_result.pose() );
	core::pose::PoseOP lores_pose( denovo_result.lores_pose() );

	// special denovo one??
	using SFPOS = SilentFilePoseOutputSpecification;
	debug_assert( dynamic_cast< SFPOS const * > ( &spec ) );
	auto const & sf_spec( static_cast< SFPOS const & > ( spec ) );

	if ( ! opts_ ) {
		// I.e., we are a fresh outputter
		opts_.reset( new core::io::silent::SilentFileOptions( sf_spec.sf_opts() ) );
		fname_out_ = sf_spec.out_fname();
		buffer_limit_ = sf_spec.buffer_limit();
	}

	/*
	void align_and_output_to_silent_file(
	core::pose::Pose & pose,
	std::string const & silent_file,
	std::string const & out_file_tag
	) const {
	rna_fragment_monte_carlo_->align_pose( pose, true verbose );
	output_to_silent_file( pose, silent_file, out_file_tag, false score_only );

	using namespace core::io::silent;
	using namespace core::scoring;

	if ( rna_params_->is_rna_and_protein() && !options_->minimize_structure() ) {
	// convert back to full atom (should give stupid coords... ok for now b/c protein doesn't move)
	// if protein sidechains have moved, then the pose should already be in full atom by now (?!)
	// if the structure is getting minimized, then it should already be converted back to full atom
	core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD, false  no sloppy match , true  only switch protein residues , true  keep energies!  );
	// but as soon as I score again, it tries to recalculate rnp scores (so they get set to 0)
	}

	// Silent file setup?
	SilentFileOptions opts;
	SilentFileData silent_file_data( opts );

	// What is all this rigamarole, making the silent struct data?
	// Why do I need to supply the damn file name? That seems silly.
	TR << "Making silent struct for " << out_file_tag << std::endl;

	SilentStructOP s = ( options_->binary_rna_output() ) ? SilentStructOP( new BinarySilentStruct( opts, pose, out_file_tag ) ) :
	SilentStructOP( new RNA_SilentStruct( opts, pose, out_file_tag ) );

	if ( options_->use_chem_shift_data() ) add_chem_shift_info( *s, pose);

	output_silent_struct( *s, silent_file_data, silent_file, pose, out_file_tag, score_only );
	}
	*/
	// What is all this rigamarole, making the silent struct data?
	// Why do I need to supply the damn file name? That seems silly.


	// opts_ should take the thing that sets options_->binary_rna_output() into account
	//core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose, *opts_ );

	//std::string output_tag = ( job.status_prefix() == "" ? "" : ( job.status_prefix() + "_" ) )
	//+ job.inner_job()->job_tag() + "_" + ObjexxFCL::lead_zero_string_of( pose_ind_of_total.first,
	//                std::max( 6, 1 + int( std::log10( pose_ind_of_total.second ))));
	//std::cout << "WORKS. Making silent struct for " << output_tag << std::endl;


	TR << "Making silent struct for " << sf_spec.pose_tag() << std::endl;
	// SilentStructOP( new RNA_SilentStruct( opts, pose, out_file_tag ) );

	core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( *pose, sf_spec.sf_opts() );
	ss->fill_struct( *pose, sf_spec.pose_tag() );
	//if ( options_->use_chem_shift_data() ) add_chem_shift_info( *s, pose);


	//if ( denovo_result.get_native_pose() ) calc_rmsds( ss, pose, out_file_tag );
	add_number_base_pairs( *pose, *ss );
	if ( denovo_result.get_native_pose() ) {
		add_number_native_base_pairs( *pose, *denovo_result.get_native_pose(), *ss );
	}


	buffered_structs_.push_back( ss );

	if ( lores_pose ) {
		core::io::silent::SilentStructOP lores_ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( *lores_pose, *opts_ );
		lores_ss->fill_struct( *lores_pose, sf_spec.pose_tag() );

		add_number_base_pairs( *lores_pose, *lores_ss );
		if ( denovo_result.get_native_pose() ) {
			add_number_native_base_pairs( *lores_pose, *denovo_result.get_native_pose(), *lores_ss );
		}

		buffered_lores_structs_.push_back( lores_ss );
	}

	if ( buffered_structs_.size() >= buffer_limit_ ) {
		flush();
	}


	if ( dump_ ) {
		std::string const out_file_name = sf_spec.pose_tag() + ".pdb";
		core::io::pdb::dump_pdb( *pose, out_file_name );
	}
}

/*
void
DeNovoSilentFilePoseOutputter::write_output_pose(
LarvalJob const & job,
std::pair< core::Size, core::Size > const & pose_ind_of_total,
utility::options::OptionCollection const & job_options,
utility::tag::TagCOP tag, // possibly null-pointing tag pointer
core::pose::Pose const & pose
)
{


}
*/

void DeNovoSilentFilePoseOutputter::flush()
{
	if ( ! opts_ ) {
		debug_assert( buffered_structs_.size() == 0 );
		return;
	}

	// Here or in SilentFileData, we also need to EXPLICITLY output OTHER_STRUCTs
	core::io::silent::SilentFileData sfd( *opts_ );
	for ( auto const & iter : buffered_structs_ ) {
		sfd.add_structure( *iter );
		// One must explicitly add these because of the structure of
		// SilentFileData. Maybe it can or should change! Because then on
		// the other side, it has to change whether the struct being
		// added is an other_struct. This perturbs the code less though.
		for ( auto const & other_struct : iter->other_struct_list() ) {
			sfd.add_structure( *other_struct );
		}
	}
	sfd.write_all( fname_out_ );
	buffered_structs_.clear();


	// Ditto for lores

	core::io::silent::SilentFileData lores_sfd( *opts_ );
	for ( auto const & iter : buffered_lores_structs_ ) {
		lores_sfd.add_structure( *iter );
		// One must explicitly add these because of the structure of
		// SilentFileData. Maybe it can or should change! Because then on
		// the other side, it has to change whether the struct being
		// added is an other_struct. This perturbs the code less though.
		for ( auto const & other_struct : iter->other_struct_list() ) {
			lores_sfd.add_structure( *other_struct );
		}
	}
	// LORES filename is LORES suffixed but before extension if one exists.
	if ( fname_out_.find( ".out" ) == std::string::npos ) {
		sfd.write_all( fname_out_ + "_LORES" );
	} else {
		lores_sfd.write_all( fname_out_.substr( 0, fname_out_.find( ".out" ) ) + "_LORES.out" );

	}

	buffered_lores_structs_.clear();
}

std::string
DeNovoSilentFilePoseOutputter::keyname() { return "DeNovoSilentFile"; }

void
DeNovoSilentFilePoseOutputter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::required_attribute( "filename", xs_string , "The name of the output silent file that should be written to." )
		+ XMLSchemaAttribute( "buffer_limit", xsct_non_negative_integer, "The number of Poses that should be held in memory between each write to disk" );
	core::io::silent::SilentFileOptions::append_attributes_for_tag_parsing( xsd, attributes );


	pose_outputter_xsd_type_definition_w_attributes( xsd, keyname(),
		"A PoseOutputter that writes structures out in a rosetta-specific format; a single "
		" silent file can hold hundreds or thousdands of output structures, lessening the load"
		" on file systems, and making output management easier.",
		attributes );
}

void
DeNovoSilentFilePoseOutputter::list_options_read(
	utility::options::OptionKeyList & read_options
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::io::silent::SilentFileOptions::list_read_options( read_options );
	read_options
		+ out::silent_gz
		+ out::file::silent
		+ rna::denovo::out::dump ;
}
/*
std::string lores_fname( std::string const & fname_out ) {
auto period_pos = fname_out.find_last_of( '.' );
std::string lores_fname = fname_out;
if ( period_pos == std::string::npos ) {
// No period. Therefore, just append "_LORES".
lores_fname.replace( fname_out.size(), "_LORES" );
} else {
// There is a final period. Therefore, preserve substring
// from it to end -- that's the desired extension.
lores_fname.insert( period_pos, "_LORES" );
}
return lores_fname;
}
*/
void
DeNovoSilentFilePoseOutputter::initialize_sf_options(
	utility::options::OptionCollection const & job_options,
	utility::tag::TagCOP tag // possibly null-pointing tag pointer
)
{
	using namespace core::io::silent;
	opts_ = SilentFileOptionsOP( new SilentFileOptions( job_options ) );
	using namespace basic::options::OptionKeys;
	dump_ = job_options[ rna::denovo::out::dump ].value();
	if ( tag ) {
		opts_->read_from_tag( tag );
		fname_out_ = tag->getOption< std::string >( "filename" );
		//lores_fname_out_ = fname_out_;

		// Turn foo.out into foo_LORES.out.
		// But don't care about extension and copy until final .
		//lores_fname_out_ = lores_fname( fname_out_ );

		if ( tag->hasOption( "buffer_limit" ) ) {
			buffer_limit_ = tag->getOption< core::Size >( "buffer_limit" );
		}
	} else {
		fname_out_ = job_options[ out::file::silent ]();
		// buffer limit?
	}

}


PoseOutputterOP DeNovoSilentFilePoseOutputterCreator::create_outputter() const
{
	return PoseOutputterOP( new DeNovoSilentFilePoseOutputter );
}

std::string DeNovoSilentFilePoseOutputterCreator::keyname() const
{
	return DeNovoSilentFilePoseOutputter::keyname();
}

void DeNovoSilentFilePoseOutputterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeNovoSilentFilePoseOutputter::provide_xml_schema( xsd );
}

void DeNovoSilentFilePoseOutputterCreator::list_options_read( utility::options::OptionKeyList & read_options ) const
{
	DeNovoSilentFilePoseOutputter::list_options_read( read_options );
}

bool
DeNovoSilentFilePoseOutputterCreator::outputter_specified_by_command_line() const
{
	return DeNovoSilentFilePoseOutputter::outputter_specified_by_command_line();
}


} // namespace pose_outputters
} // namespace jd3
} // namespace protocols
