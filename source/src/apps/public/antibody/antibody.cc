// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

#include <protocols/antibody/AntibodyNumberingConverterMover.hh>

// Grafting util function related includes
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/sequence/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/FastRelax.hh>
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
OPT_KEY( Boolean,      optimal_graft)
OPT_KEY( Boolean,      optimize_cdrs)

void register_options()
{
	NEW_OPT( heavy, "Name of fasta file with antibody heavy chain", "" );
	NEW_OPT( light, "Name of fasta file with antibody light chain", "" );
	NEW_OPT( fasta, "Name of fasta file from which both antibody heavy and light chain should be read", "" );
	//NEW_OPT( n, "Number of templates to generate", 10);
	NEW_OPT( multi_template_regions, "Specify sequence regions for which multi-template is generated. Avalible regions are: h1, h2, h3, l1, l2, l3, frh, frl, orientation", "orientation");
	NEW_OPT( no_relax, "Do not relax grafted model", false);
	NEW_OPT( optimal_graft, "If the graft is not closed, use a new graft mover.  Before any full relax, relax the CDRs with dihedral constraints to make them fit better.   May become default after testing", false);
	NEW_OPT( optimize_cdrs, "Optimize CDRs.  If doing optimal graft, we will optimize the CDRs.", false);
}

void relax_model(core::pose::PoseOP &pose)
{
	using namespace basic::options;
	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	if ( !option[ OptionKeys::out::file::scorefile ].user() ) {
		option[ OptionKeys::out::file::scorefile ].value("score-relax.sf");
	}

	// have to specify no ramping, because options are not yet set on the command line
	protocols::relax::FastRelaxOP relax_protocol( new protocols::relax::FastRelax(core::scoring::get_score_function(), "NO CST RAMPING") );

	// set packer task related defaults on pose
	TaskFactoryOP task_factory( new TaskFactory() );

	// I'm pretty sure that this is handled by FastRelax... Yet scary things happen when it is not included.
	task_factory->push_back( operation::TaskOperationCOP( new operation::RestrictToRepacking ) ); // repack only, no design
	task_factory->push_back( TaskOperationCOP( new ExtraRotamers( 0 /*all*/, 1 /*ex1*/, 1 /*level*/ ) ) );
	task_factory->push_back( TaskOperationCOP( new ExtraRotamers( 0 /*all*/, 2 /*ex1*/, 1 /*level*/ ) ) );

	task_factory->push_back(TaskOperationCOP( new IncludeCurrent() ));

	relax_protocol->set_task_factory( task_factory );

	// set relax defaults
	relax_protocol->constrain_relax_to_start_coords( true );
	relax_protocol->constrain_coords( true );
	relax_protocol->coord_constrain_sidechains( true );

	relax_protocol->ramp_down_constraints( false );

	//relax
	relax_protocol->apply( *pose );
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

	TR << TR.Magenta << "Unzipping files in antibody_database (if any). This will only be done once." << TR.Reset << std::endl;

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

				bool optimal_graft = basic::options::option[ basic::options::OptionKeys::optimal_graft]();
				bool optimize_cdrs = basic::options::option[ basic::options::OptionKeys::optimize_cdrs]();

				core::pose::PoseOP model = graft_cdr_loops(as, r, prefix, suffix, grafting_database, optimal_graft, optimize_cdrs );

				if( !basic::options::option[ basic::options::OptionKeys::no_relax ]() ) {
					relax_model(model);
					model->dump_pdb(prefix + "model" + suffix + ".relaxed.pdb");
				}

				if ( basic::options::option[ basic::options::OptionKeys::antibody::output_ab_scheme].user()){
					protocols::antibody::AntibodyNumberingConverterMover converter = protocols::antibody::AntibodyNumberingConverterMover();
					converter.apply(*model);
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
				if( grafting_database.size() == 0 ) utility_exit_with_message( "Could not guess location of antibody grafting_database (typically $ROSETTA/tools/antibody), please use -antibody::grafting_database option to specify it!");
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
