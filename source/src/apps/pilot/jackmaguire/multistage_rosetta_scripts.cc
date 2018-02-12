// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/multistage_rosetta_scripts/MRSJobQueen.hh>
#include <protocols/multistage_rosetta_scripts/MRSJobQueenChecker.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd3/JobDistributor.hh>
#include <protocols/jd3/JobDistributorFactory.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/Jobs.hh>

#include <fstream>

#ifdef SERIALIZATION
#include <protocols/jd3/job_distributors/MPIWorkPoolJobDistributor.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <core/pose/Pose.hh>
#endif

#include <apps/pilot/jackmaguire/reverse_conversion_util.hh>

OPT_1GRP_KEY( Boolean, mrs, convert )
OPT_1GRP_KEY( Boolean, mrs, revert )
OPT_1GRP_KEY( Boolean, mrs, reverse_convert )
OPT_1GRP_KEY( Boolean, mrs, xml_template )
OPT_1GRP_KEY( Boolean, mrs, estimate_memory )
OPT_1GRP_KEY( StringVector, mrs, unarchive )

static basic::Tracer TR( "apps.pilot.jackmaguire.MultistageRosettaScripts" );

#ifdef SERIALIZATION
//We want to access the protected deserialize_larval_job_and_result method
//so let's just make a derived class that has access to this method
class JobDistributorForDeserializing : public protocols::jd3::job_distributors::MPIWorkPoolJobDistributor {
public:
	core::pose::PoseOP deserialize( std::string const & filename ){
		std::string const serilized_result = utility::file_contents( filename );
		std::pair< protocols::jd3::LarvalJobOP, protocols::jd3::JobResultOP > result = deserialize_larval_job_and_result( serilized_result );
		protocols::jd3::standard::PoseJobResult & pose_result = static_cast< protocols::jd3::standard::PoseJobResult & >( * result.second );
		return pose_result.pose();
	}
};

void unarchive( std::string filename_of_archive ){

	JobDistributorForDeserializing jd;
	core::pose::PoseOP pose = jd.deserialize( filename_of_archive );
	pose->dump_pdb( filename_of_archive + ".pdb" );

}
#endif

void insert_stage_tag( std::string & s ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::stringstream ss1;

	ss1 << std::endl << "\t\t\t<Stage num_runs_per_input_struct=\"";
	if ( option[ out::nstruct ].user() ) {
		ss1 << option[ out::nstruct ]();
	} else {
		ss1 << "TODO";
	}
	ss1 << "\" total_num_results_to_keep=\"TODO\">" << std::endl;

	core::Size protocols_open = s.find( std::string("<PROTOCOLS>") );
	if ( protocols_open < 2 || protocols_open == std::string::npos ) {
		core::Size protocols_one_line = s.find( std::string("<PROTOCOLS/>") );
		if ( protocols_one_line == std::string::npos && false ) { //&& false because this does not seem to work
			utility_exit_with_message("could not find <PROTOCOLS> : \n" + s );
		} else {
			s.replace( protocols_one_line, 12, "<PROTOCOLS>\n</PROTOCOLS>");
			protocols_open = s.find( std::string("<PROTOCOLS>") );
		}
	}
	//utility_exit_with_message( s );

	s.insert( protocols_open + 11, ss1.str() );

	//replace <Add with \t<Add
	core::Size add_count = 0;
	core::Size add = s.find("<Add", protocols_open );
	core::Size protocols_close = s.find( std::string("</PROTOCOLS>") );
	while ( add != std::string::npos && add < protocols_close ) {
		++add_count;
		s.insert( add, std::string("\t") );

		add = s.find("<Add", add+6 );
		protocols_close = s.find( std::string("</PROTOCOLS>") );
	}

	std::stringstream ss2;
	if ( add_count == 0 ) ss2 << "\t\t";
	ss2 << "\t\t<Sort filter=\"TODO\"/>" << std::endl << std::endl << "\t\t\t</Stage>" << std::endl << "\t\t";
	s.insert( protocols_close, ss2.str() );
}

void reverse_convert(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	runtime_assert( option[ in::file::job_definition_file ].user() );

	std::string const in_fname = option[ in::file::job_definition_file ]();
	std::string const in_string = utility::file_contents( in_fname );
	utility::tag::TagOP mrs_tag = utility::tag::Tag::create( in_string );

	std::string const out_string = reverse_convert( * mrs_tag );

	bool print_to_cout = true;
	std::string out_fname;

	if ( option[ parser::protocol ].user() ) {
		out_fname = option[ parser::protocol ]();
		print_to_cout = false;
	} else if ( option[ out::output ].user() ) {
		out_fname = option[ out::output ]();
		print_to_cout = false;
	}

	if ( ! print_to_cout ) {
		std::ofstream out( out_fname );
		if ( ! out.good() ) {
			print_to_cout = true;
		} else {
			out << out_string;
		}
		out.close();
	}
	if ( print_to_cout ) {
		std::cout << out_string;
	}

}

void convert( ){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::tag::TagOP rosetta_scripts_tag;
	if ( option[ parser::protocol ].user() ) {
		std::string const in_fname = option[ parser::protocol ]();

		std::string const in_string = utility::file_contents( in_fname );
		rosetta_scripts_tag = utility::tag::Tag::create( in_string );
	} else {
		std::stringstream rss_template;
		rss_template << "<ROSETTASCRIPTS>\n"
			<< "\t<SCOREFXNS>\n"
			<< "\t</SCOREFXNS>\n"
			<< "\t<RESIDUE_SELECTORS>\n"
			<< "\t</RESIDUE_SELECTORS>\n"
			<< "\t<TASKOPERATIONS>\n"
			<< "\t</TASKOPERATIONS>\n"
			<< "\t<FILTERS>\n"
			<< "\t</FILTERS>\n"
			<< "\t<MOVERS>\n"
			<< "\t</MOVERS>\n"
			<< "\t<PROTOCOLS>\n"
			<< "\t</PROTOCOLS>\n"
			<< "</ROSETTASCRIPTS>\n\n";
		rosetta_scripts_tag = utility::tag::Tag::create( rss_template.str() );
	}

	std::stringstream ss;
	ss << "<JobDefinitionFile>" << std::endl;

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();
	if ( input_jobs.size() > 0 ) {
		for ( auto const & ii : input_jobs ) {
			ss << std::endl << "\t<Job>" << std::endl;
			ss << "\t\t<Input>" << std::endl;
			ss << "\t\t\tTODO: The script converter is not smart enough to know if this file is actually a pdb file. Please correct this tag if necessary!" << std::endl;
			ss << "\t\t\t<PDB filename=\"" << ii->input_tag() << "\"/>" << std::endl;
			ss << "\t\t</Input>" << std::endl;
			ss << "\t</Job>" << std::endl;
		}
	} else {
		ss << std::endl << "\t<Job>" << std::endl;
		ss << "\t\t<Input>" << std::endl;
		ss << "\t\t\t<PDB filename=\"TODO\"/>" << std::endl;
		ss << "\t\t</Input>" << std::endl;
		ss << "\t</Job>" << std::endl;
	}

	ss << "\t<Common>" << std::endl << std::endl;
	for ( auto subtag : rosetta_scripts_tag->getTags() ) {
		if ( subtag->getName() == "APPLY_TO_POSE" ) continue;
		if ( subtag->getName() == "OUTPUT" ) continue;
		subtag->write( ss, 2 );
		ss << std::endl;
	}
	ss << "\t</Common>" << std::endl;

	ss << std::endl << "</JobDefinitionFile>" << std::endl;

	std::string out_string = ss.str();
	insert_stage_tag( out_string );

	bool print_to_cout = true;
	std::string out_fname;
	if ( option[ in::file::job_definition_file ].user() ) {
		out_fname = option[ in::file::job_definition_file ]();
		print_to_cout = false;
	} else if ( option[ out::output ].user() ) {
		out_fname = option[ out::output ]();
		print_to_cout = false;
	}

	if ( ! print_to_cout ) {
		std::ofstream out( out_fname );
		if ( ! out.good() ) {
			print_to_cout = true;
		} else {
			out << out_string;
		}
		out.close();
	}
	if ( print_to_cout ) {
		std::cout << out_string;
	}

}

int main( int argc, char* argv[] ){

	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT( mrs::convert, "Convert a traditional rosetta script to one that is compatible with this program. Use -parser:protocol to declare the original and -job_definition_file to give the filename for the new file. The converter may take advantage of other flags (such as -s and -nstruct) so it is recommended to run this with all of your production flags present.", false );
		NEW_OPT( mrs::revert, "Exactly the opposite of '-convert'. Takes the MRS scripts defined by -job_definition_file and converts it to a traditional rosetta script (will only keep the datamap elements in <Common/>). Will print to console or to file if -parser:protocol is used.", false );
		NEW_OPT( mrs::reverse_convert, "Idential to '-revert'", false );

		NEW_OPT( mrs::xml_template, "Dump template to console.", false );
		NEW_OPT( mrs::estimate_memory, "(Experimental) estimate the maximum amount of memory needed by the archive nodes at any given time. So far, this has been shown to underestimate but not by more than 2x.", false );
		NEW_OPT( mrs::unarchive, "If you previously ran with -archive_on_disk, you can pass in one of those archive files here and it will be dumped as a pdb using the naming scheme of the archive file. Please run with all of your original command line options present", utility::vector1<std::string>() );

		devel::init( argc, argv );

		if ( option[ mrs::convert ] || option[ mrs::xml_template ] ) {
			convert();
			return 0;
		} else if ( option[ mrs::revert ] || option[ mrs::reverse_convert ] ) {
			reverse_convert();
			return 0;
		}

		if ( option[ mrs::estimate_memory ] ) {
			TR.Debug << "Estimating Memory" << std::endl;
			protocols::multistage_rosetta_scripts::MRSJobQueenCheckerOP queen ( new protocols::multistage_rosetta_scripts::MRSJobQueenChecker );
			queen->initial_job_dag();
			queen->estimate_number_of_bytes_needed_for_archiving();
			return 0;
		}

		if ( option[ mrs::unarchive ].user() ) {
			TR.Debug << "Unarchiving" << std::endl;
#ifdef SERIALIZATION
			utility::vector1< std::string > filenames = option[ mrs::unarchive ];
			runtime_assert( filenames.size() > 0 );
			for ( auto filename : filenames ) {
				unarchive( filename );
			}
			//unarchive( option[ mrs::unarchive ] );
#else
			TR << "The -unarchive option only works when you compile with serialization" << std::endl;
#endif
			return 0;
		}

		if ( option[ parser::info ].user() ) {
			TR << "Printing Info" << std::endl;
			protocols::rosetta_scripts::print_information( option[ parser::info ]() );
			return 0;
		}

		protocols::multistage_rosetta_scripts::MRSJobQueenOP queen ( new protocols::multistage_rosetta_scripts::MRSJobQueen );
		protocols::jd3::JobDistributorOP job_dist = protocols::jd3::JobDistributorFactory::create_job_distributor();
		TR.Debug << "Running MRSJobQueen" << std::endl;
		job_dist->go( queen );
		return 0;

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
