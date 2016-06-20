// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file antibody.cc
///
/// @brief C++ port of Python antibody script
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/scs_multi_template.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/grafting/scs_functor.hh>
#include <protocols/antibody/grafting/grafter.hh>

// Grafting util function related includes
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/sequence/util.hh>

#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

#include <basic/report.hh>
#include <basic/execute.hh>

#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//#include <utility/CSI_Sequence.fwd.hh>
#include <utility/file/file_sys_util.hh>

#include <devel/init.hh>

#include <fstream>
#include <cassert>
#include <cstdlib>
#include <iostream>

#include <string>
#include <regex>
#include <cstdlib>

//#include <typeinfo>

using std::string;
using core::Size;

static THREAD_LOCAL basic::Tracer TR("antibody");

// Command-line options relevant only to this application
OPT_KEY( String,       heavy )
OPT_KEY( String,       light )
OPT_KEY( String,       fasta )
//OPT_KEY( Integer,      n)
OPT_KEY( StringVector, multi_template_regions)
OPT_KEY( Boolean,      no_relax)


void register_options()
{
	NEW_OPT( heavy, "Name of fasta file with antibody heavy chain", "" );
	NEW_OPT( light, "Name of fasta file with antibody light chain", "" );
	NEW_OPT( fasta, "Name of fasta file from which both antibody heavy and light chain should be read", "" );
	//NEW_OPT( n, "Number of templates to generate", 10);
	NEW_OPT( multi_template_regions, "Specify sequence regions for which multi-template is generated. Avalible regions are: h1, h2, h3, l1, l2, l3, frh, frl, orientation", "orientation");
	NEW_OPT( no_relax, "do not relax grafted model", false);
}

void relax_model(core::pose::PoseOP &pose)
{
	using namespace basic::options;

	// manually set options for user, if not set (probs not the way to do this, but BDW said it should be ok for now...)
	if ( !option[ OptionKeys::relax::constrain_relax_to_start_coords ].user() ) {
		option[ OptionKeys::relax::constrain_relax_to_start_coords ].value(true);
	}

	if ( !option[ OptionKeys::relax::coord_constrain_sidechains ].user() ) {
		option[ OptionKeys::relax::coord_constrain_sidechains ].value(true);
	}

	if ( !option[ OptionKeys::relax::ramp_constraints ].user() ) {
		option[ OptionKeys::relax::ramp_constraints ].value(false);
	}

	if ( !option[ OptionKeys::packing::ex1::ex1 ].user() ) {
		option[ OptionKeys::packing::ex1::ex1 ].value(true);
	}

	if ( !option[ OptionKeys::packing::ex2::ex2 ].user() ) {
		option[ OptionKeys::packing::ex2::ex2 ].value(true);
	}

	if ( !option[ OptionKeys::packing::use_input_sc ].user() ) {
		option[ OptionKeys::packing::use_input_sc ].value(true);
	}

	if ( !option[ OptionKeys::out::file::scorefile ].user() ) {
		option[ OptionKeys::out::file::scorefile ].value("score-relax.sf");
	}

	protocols::relax::RelaxProtocolBaseOP protocol = protocols::relax::generate_relax_from_cmd();

	//relax
	protocol->apply( *pose );
}

int antibody_main()
{
	using namespace utility;
	using namespace protocols::antibody::grafting;
	using protocols::antibody::grafting::uint;

	using std::string;

	string heavy_fasta_file, light_fasta_file, heavy_chain_sequence, light_chain_sequence;

	int n_templates = basic::options::option[ basic::options::OptionKeys::antibody::n_multi_templates ];

	string grafting_database = basic::options::option[basic::options::OptionKeys::antibody::grafting_database]();
	TR << TR.Cyan << "Using antibody grafting_database at: " << TR.Bold << grafting_database << std::endl;
	
	// add code to unpack PDBs here
	vector1<string> file_names;

	// assemble full path to antibody_database
	std::string full_ab_db_path = grafting_database + "/antibody_database/";
	
	// check if there are any files in the grafting database
	file::list_dir( full_ab_db_path, file_names );

	TR << TR.Magenta << "Unzipping files in antibody_databse (if any). This will only be done once." << TR.Reset << std::endl;
	
	// iterate over vector looking for bz2's and unzip
	for (auto it = file_names.begin(); it != file_names.end() ; ++it) {
		// check for bz2 and not already unzipped
		if ( file::file_extension( *it ).compare( "bz2" ) == 0 and ! file::file_exists( full_ab_db_path + file::file_basename( *it ) ) ) {
				// bzipped files, unzip
				basic::execute("Unzipping " + *it ,"cd " + full_ab_db_path + " && bunzip2 -k " + *it);
		}
	}
	TR << TR.Magenta << "Done unzipping." << TR.Reset << std::endl;
	

	if( basic::options::option[ basic::options::OptionKeys::heavy ].user()  and  basic::options::option[ basic::options::OptionKeys::light ].user()  and  !basic::options::option[ basic::options::OptionKeys::fasta ].user() ) {
		heavy_fasta_file = basic::options::option[ basic::options::OptionKeys::heavy ]();
		light_fasta_file = basic::options::option[ basic::options::OptionKeys::light ]();
	} else if( basic::options::option[ basic::options::OptionKeys::fasta ].user() ) {
		heavy_fasta_file = light_fasta_file = basic::options::option[ basic::options::OptionKeys::fasta ]();
	} else if( basic::options::option[ basic::options::OptionKeys::heavy ].user()  and  basic::options::option[ basic::options::OptionKeys::light ].user()  and  basic::options::option[ basic::options::OptionKeys::fasta ].user() ) {
		utility_exit_with_message( "Error: options --fasta and [--light, heavy] can't be specifed at the same time!");
	} else {
		utility_exit_with_message( "Error: no input sequences was specified!");
	}

	TR << TR.Green << "Heavy chain fasta input file: " << TR.Bold << heavy_fasta_file << TR.Reset << std::endl;
	TR << TR.Blue  << "Light chain fasta input file: " << TR.Bold << light_fasta_file << TR.Reset << std::endl;

	heavy_chain_sequence = core::sequence::read_fasta_file_section(heavy_fasta_file, "heavy");
	light_chain_sequence = core::sequence::read_fasta_file_section(light_fasta_file, "light");

	{
		TR << TR.Red << "Antibody grafting protocol" << TR.Reset << " Running with input:" << std::endl;
		if( n_templates > 1 ) {
			TR << TR.Red << "SCS_MultiTemplate protocol enabled" << TR.Reset << " multi-template regions: ";
			for(auto & s : basic::options::option[basic::options::OptionKeys::multi_template_regions]() ) TR << s << " ";
			TR << std::endl;
		}

		TR << TR.Green << "Heavy chain sequence: " << TR.Bold << heavy_chain_sequence << TR.Reset << std::endl;
		TR << TR.Blue  << "Light chain sequence: " << TR.Bold << light_chain_sequence << TR.Reset << std::endl;
	}

	string const prefix = basic::options::option[basic::options::OptionKeys::antibody::prefix]();

	// strip directory from prefix, then make it recursively
	string const prefix_path = prefix.substr( 0, prefix.find_last_of( "/\\" ) );
	
	if ( !file::is_directory( prefix_path ) ) {
		file::create_directory_recursive( prefix_path );
	}

	basic::ReportOP report = std::make_shared<basic::Report>(prefix+"report");

	AntibodySequence as(heavy_chain_sequence, light_chain_sequence);

	RegEx_based_CDR_Detector(report).detect(as);
	std::cout << std::endl;

	AntibodyFramework heavy_fr = as.heavy_framework();
	AntibodyFramework light_fr = as.light_framework();

	std::cout << "Heavy" << heavy_fr << std::endl;
	std::cout << "Light" << light_fr << std::endl;
	std::cout << std::endl;

	try {
		// AntibodyNumbering an( Chothia_Numberer().number(as, heavy_fr, light_fr) );
		// std::cout << "Heavy Numbering: " << an.heavy << std::endl << std::endl;
		// std::cout << "Light Numbering: " << an.heavy << std::endl << std::endl;

		SCS_BlastPlusOP blast = n_templates == 1 ? std::make_shared<SCS_BlastPlus>(report) : std::make_shared<SCS_BlastPlus>();
		blast->init_from_options();

		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_sequence_length) );
		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_alignment_length) );

		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_sequence_identity) );
		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_template_resolution) );
		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_outlier) );
		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_template_bfactor) );
		blast->add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_OCD) );

		blast->set_sorter( SCS_FunctorCOP(new SCS_BlastComparator_BitScore_Resolution) );

		SCS_BaseOP scs = n_templates == 1 ? blast : SCS_BaseOP(new SCS_MultiTemplate(blast, basic::options::option[basic::options::OptionKeys::multi_template_regions](), report) );

		SCS_ResultsOP scs_results { scs->select(n_templates, as) };

		report->write();

		if( scs_results->frh.size()  and  scs_results->frl.size()  and  scs_results->orientation.size() and
			scs_results->h1.size() and scs_results->h2.size() and scs_results->h3.size() and
			scs_results->l1.size() and scs_results->l2.size() and scs_results->l3.size()
			) {
			TR << "SCS frh best template: " << TR.Bold << TR.Blue  << scs_results->frh[0]->pdb << TR.Reset << std::endl;
			TR << "SCS frl best template: " << TR.Bold << TR.Green << scs_results->frl[0]->pdb << TR.Reset << std::endl;
			TR << "SCS Orientation best template: " << TR.Bold << TR.Yellow << scs_results->orientation[0]->pdb << TR.Reset << std::endl;

			for(int i=0; i<n_templates; ++i) {
				SCS_ResultSet r;
				try { r = scs_results->get_result_set(i); }
				catch(std::out_of_range e) { break; }

				string const suffix = '-'+std::to_string(i);

				// override scs_resultset->resultop->pdb if user specifies a template
				// needs to check for identical length, but how?
				// all only for manual specification if single template grafting
				if ( basic::options::option[ basic::options::OptionKeys::antibody::l1_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified L1 CDR!" << TR.Reset << std::endl;
					r.l1->pdb = basic::options::option[ basic::options::OptionKeys::antibody::l1_template ].value().substr( 0, 4 );
				}
				if ( basic::options::option[ basic::options::OptionKeys::antibody::l2_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified L2 CDR!" << TR.Reset << std::endl;
					r.l2->pdb = basic::options::option[ basic::options::OptionKeys::antibody::l2_template ].value().substr( 0, 4 );
				}
				if ( basic::options::option[ basic::options::OptionKeys::antibody::l3_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified L3 CDR!" << TR.Reset << std::endl;
					r.l3->pdb = basic::options::option[ basic::options::OptionKeys::antibody::l3_template ].value().substr( 0, 4 );
				}

				if ( basic::options::option[ basic::options::OptionKeys::antibody::h1_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified H1 CDR!" << TR.Reset << std::endl;
					r.h1->pdb = basic::options::option[ basic::options::OptionKeys::antibody::h1_template ].value().substr( 0, 4 );
				}
				if ( basic::options::option[ basic::options::OptionKeys::antibody::h2_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified H2 CDR!" << TR.Reset << std::endl;
					r.h2->pdb = basic::options::option[ basic::options::OptionKeys::antibody::h2_template ].value().substr( 0, 4 );
				}
				if ( basic::options::option[ basic::options::OptionKeys::antibody::h3_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified H3 CDR!" << TR.Reset << std::endl;
					r.h3->pdb = basic::options::option[ basic::options::OptionKeys::antibody::h3_template ].value().substr( 0, 4 );
				}

				if ( basic::options::option[ basic::options::OptionKeys::antibody::frl_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified light FR!" << TR.Reset << std::endl;
					r.frl->pdb = basic::options::option[ basic::options::OptionKeys::antibody::frl_template ].value().substr( 0, 4 );
				}
				if ( basic::options::option[ basic::options::OptionKeys::antibody::frh_template ].user() ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified heavy FR!" << TR.Reset << std::endl;
					r.frh->pdb = basic::options::option[ basic::options::OptionKeys::antibody::frh_template ].value().substr( 0, 4 );
				}

				if ( basic::options::option[ basic::options::OptionKeys::antibody::light_heavy_template ].user() and
						 n_templates == 1 ) {
					TR << TR.Bold << TR.Red << "Note: Grafting manually specified light-heavy orientation!" << TR.Reset << std::endl;
					r.orientation->pdb = basic::options::option[ basic::options::OptionKeys::antibody::light_heavy_template ].value().substr( 0, 4 );
				}
				else if ( basic::options::option[ basic::options::OptionKeys::antibody::light_heavy_template ].user() and
									n_templates > 1 )
				{
					TR << "Cannot both do multi-template grafting and specify a single light-heavy orientation template!" << std::endl;
					TR << "Defaulting to multi-template grafting." << std::endl;
					TR << "Please use -antibody:n_multi_templates 1 in conjuction with -antibody:light_heavy_template 1abc.pdb" << std::endl;
				}

				core::pose::PoseOP model = graft_cdr_loops(as, r, prefix, suffix, grafting_database );

				if( !basic::options::option[ basic::options::OptionKeys::no_relax ]() ) {
					relax_model(model);
					model->dump_pdb(prefix + "model" + suffix + ".relaxed.pdb");
				}
			}
		}
		else {
			TR << TR.Red << "Blast+ results for some of the regions are empty!!! Exiting..." << std::endl;
			return 1;
		}
	}
	catch(Grafting_Base_Exception e) {
		std::cout << e << std::endl;
	}

	return 0;
}


int main(int argc, char * argv [])
{
	try {
		register_options();

		basic::options::option.add_relevant( basic::options::OptionKeys::antibody::prefix );
		basic::options::option.add_relevant( basic::options::OptionKeys::antibody::grafting_database );
		basic::options::option.add_relevant( basic::options::OptionKeys::antibody::blastp );
		basic::options::option.add_relevant( basic::options::OptionKeys::antibody::n_multi_templates );

		devel::init(argc, argv);

		string grafting_database = basic::options::option[basic::options::OptionKeys::antibody::grafting_database]();
		if( grafting_database.size() == 0 ) {
			string location = "../../tools/antibody";

			if( std::ifstream(location+"/antibody_database/list").good() ) grafting_database = location;
			else {
				if( char * r = std::getenv("ROSETTA") ) {
					location = string(r) + "/tools/antibody";
					TR << "Trying: " << location << std::endl;
					if( std::ifstream(location+"/antibody_database/list").good() ) grafting_database = location;
				}
				if( grafting_database.size() == 0 ) utility_exit_with_message( "Could not guess location of antibody grafting_database, please use antibody::grafting_database option to specify it!");
			}
		}
		basic::options::option[basic::options::OptionKeys::antibody::grafting_database](grafting_database);

		return antibody_main();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return 1;
	}
}


#else // __ANTIBODY_GRAFTING__

#include <iostream>

int main( int /* argc */, char * /* argv */ [] )
{
	std::cerr << "Antibody protocol require to be build with full C++11 support! Please use GCC-4.9+ or Clang-3.6+ and add extras=cxx11 to your scons build command line.";
	return 1;
}

#endif // __ANTIBODY_GRAFTING__
