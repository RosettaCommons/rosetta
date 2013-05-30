// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file demo/doug/rotamer_prediction_benchmark.cc
/// @brief this demo in conjuntion with a external script determines the percent of rotamers correctly predictied during a repack
/// @author P. Douglas Renfrew (renfrew@unc.edu)

// core headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>

// c++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


#include <utility/excn/Exceptions.hh>


// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// application specific options
namespace rpb {
	BooleanOptionKey const mm_score( "rpb::mm_score" );
	RealOptionKey const start_weight( "rpb::start_weight" );
	RealOptionKey const end_weight( "rpb::end_weight" );
	RealOptionKey const inc_weight( "rpb::inc_weight" );
}

void rotamer_prediction_benchmark( std::string pdb_filename, ScoreFunctionOP scorefxn, std::string id );

int
main( int argc, char * argv [] )
{
    try {
	// add application specific options to core options system
	option.add( rpb::mm_score, "Use MM Torsion Score instead of Dunbrack" );
	option.add( rpb::start_weight, "Weight to start testing MM Torsion term at" );
	option.add( rpb::end_weight, "Weight to start testing MM Torsion term at" );
	option.add( rpb::inc_weight, "How much to increment the MM Torsion term each round" );

	// init
	devel::init(argc, argv);

	// concatenate -s and -l flags together to get total list of PDB files
	// (This was taken from Ian's early job distributor, thanks Ian)
	std::vector< FileName > pdb_file_names;
	if ( option[ in::file::s ].active() )
		pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
	std::vector< FileName > list_file_names;
	if ( option[ in::file::l ].active() )
		list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)

	for(std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i) {
		std::string filename( i->name() );
		std::ifstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			pdb_file_names.push_back( FileName(line) );
		}
		data.close();
	}

	if ( option[ rpb::mm_score ] )
		{
			// get info from options
			Real start_weight =  option[ rpb::start_weight ].value();
			Real end_weight = option[ rpb::end_weight ].value();
			Real inc_weight = option[ rpb::inc_weight ].value();

			// set weight on term and run
			for( float weight = start_weight; weight <= end_weight; weight += inc_weight )
				{
					// create id based on weight
					std::stringstream weight_ss;
					weight_ss << "std_nodun_mm_" <<std::setprecision(4) << std::fixed <<weight << ".";
					std::string id = 	weight_ss.str();

					// create score function and modify it
					ScoreFunctionOP scfxn( getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS ) );
					scfxn->set_weight(fa_dun, 0.00);
					scfxn->set_weight(mm_twist, weight);

					// call RPB for each name in list
					for(std::vector< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i) {
						rotamer_prediction_benchmark( i->name(), scfxn, id );
					}
				}
		}
	else // use standard score function
		{
			// create id
			std::string id = "std_";

			// create score function
			ScoreFunctionOP scfxn( getScoreFunctionLegacy( PRE_TALARIS_2013_STANDARD_WTS ) );

			// run RPB for each name in list
			for(std::vector< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i) {
				rotamer_prediction_benchmark( i->name(), scfxn, id );
			}
		}
}

void rotamer_prediction_benchmark( std::string pdb_filename, ScoreFunctionOP scorefxn, std::string id )
{
	// create trace based on filename
	basic::Tracer GEN( "rpb." + id + pdb_filename ); // general output
	basic::Tracer RES( "rpb.result." + id +pdb_filename ); // results output

	// create 2 poses from the same pdb file
	Pose pose_orig, pose_pack;
	core::import_pose::pose_from_pdb( pose_orig, pdb_filename + ".pdb");
	core::import_pose::pose_from_pdb( pose_pack, pdb_filename + ".pdb");

	// create a packer task
	pack::task::PackerTaskOP repacktask( pack::task::TaskFactory::create_packer_task( pose_pack ));
	repacktask->initialize_from_command_line().restrict_to_repacking().or_include_current( true );
	//repacktask->set_bump_check( true );

	// run pack rotamers with first pose, scfxn, and ptsk, calc new energy
	Energy orig_score = (*scorefxn)( pose_pack );
	clock_t starttime = clock(); pack::pack_rotamers( pose_pack, *scorefxn, repacktask); clock_t stoptime = clock();
	Energy pack_score = (*scorefxn)( pose_pack );
	std::cout << "pack score: " << pack_score << " orig score: " << orig_score << " " << ((double) stoptime - starttime)/CLOCKS_PER_SEC << "seconds" << std::endl;
	io::pdb::dump_pdb( pose_pack, pdb_filename + id + "_packed.pdb" );

	// loop over all residues and calculate the chi angle differences between pose_orig and pose_pack and output the results

	int prot_total_chi1 = 0;
	int prot_total_chi2 = 0;
	int prot_chi1_correct = 0;
	int prot_chi1_chi2_correct = 0;

	// loop over all residues and calculate the chi angle differences between pose_orig and pose_pack and output the results
	Size nres = pose_orig.n_residue();
	for ( Size i = 1; i <= nres; ++i )
		{
			// get type name
			std::string name = pose_orig.residue( i ).name();

			// get mainchain torsions
			utility::vector1< Real > res_back_tor = pose_orig.residue( i ).mainchain_torsions();

			// get chi angles
			utility::vector1< Real > res_orig_chi = pose_orig.residue( i ).chi();
			utility::vector1< Real > res_pack_chi = pose_pack.residue( i ).chi();

			// print some stuff
			bool res_chi1_correct = false;
			bool res_chi2_correct = false;

			if( res_orig_chi.size() >= 1 ) prot_total_chi1++;
			if( res_orig_chi.size() >= 2 ) prot_total_chi2++;

			GEN << name << " " << i << "\t";
			for ( Size i = 1; i <= 4; ++i ){
				if( i <= res_orig_chi.size() )
					{
						Real orig_chi = numeric::nonnegative_principal_angle_degrees( res_orig_chi[ i ] );
						Real pack_chi = numeric::nonnegative_principal_angle_degrees( res_pack_chi[ i ] );
						Real diff=1e200;
						if( orig_chi >= pack_chi ) diff = orig_chi - pack_chi;
						if( pack_chi >  orig_chi ) diff = pack_chi - orig_chi;
						if( diff >= 180 ) diff = 360 - diff;
						if( diff <= 30 && i == 1 ) res_chi1_correct = true;
						if( diff <= 30 && i == 2 ) res_chi2_correct = true;
						std::cout << std::setprecision(3) << std::fixed << std::setw(9) << diff;
					}
				else
					{
						std::cout << std::setprecision(3) << std::fixed << std::setw(9) << 0.0;
					}
			}
			GEN << "\n";
			if( res_chi1_correct ) prot_chi1_correct++;
			if( res_chi1_correct && res_chi2_correct ) prot_chi1_chi2_correct++;
		}
	RES << "Chi1:  " << prot_chi1_correct << "/" << prot_total_chi1 << " = " << std::setprecision(3) << std::fixed
			<< static_cast<Real>(100)*static_cast<Real>(prot_chi1_correct)/static_cast<Real>(prot_total_chi1) << "%\n";
	RES << "Chi12: " << prot_chi1_chi2_correct << "/" << prot_total_chi2 << " = " << std::setprecision(3) << std::fixed
			<< static_cast<Real>(100)*static_cast<Real>(prot_chi1_chi2_correct)/static_cast<Real>(prot_total_chi2) << "%\n";
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
                                  }
    return 0;

}
