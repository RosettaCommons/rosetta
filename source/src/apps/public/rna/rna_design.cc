// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/chemical/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_ScoringInfo.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/pack_rotamers.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <core/init/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

// AUTO-REMOVED #include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
// AUTO-REMOVED #include <protocols/rna/RNA_DeNovoProtocol.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/excn/Exceptions.hh>

//Auto using namespaces
namespace ObjexxFCL { } using namespace ObjexxFCL; // AUTO USING NS
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end



using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, disable_o2star_rotamers )
OPT_KEY( Boolean, disable_include_current )
OPT_KEY( Boolean, sample_chi )
OPT_KEY( Boolean, ss_ds_ts_assign )
OPT_KEY( Boolean, dump )

///////////////////////////////////////////////////////////////////////////////
void
rna_sequence_recovery_metrics( pose::Pose const & reference_pose, utility::vector1< pose::PoseOP > const & pose_list, std::string const sequence_recovery_file )
{

	// make a local copy
	pose::Pose pose( reference_pose );

	// Get information on ss, ds, ts residues in native pose.
	Size const nres = pose.total_residue();
	FArray1D_int struct_type( nres, -1 );
	protocols::rna::check_base_pair( pose, struct_type );

	FArray1D_float recovery( nres, 0.0 );
	for (Size n = 1; n <= pose_list.size(); n++ ) {
		for (Size i = 1; i <= nres; i++ ) {
			if ( (pose_list[n])->residue(i).aa() == pose.residue(i).aa() ) {
				recovery( i ) += 1.0;
			}
		}
	}

	recovery /= pose_list.size();

	Size num_ss( 0 ), num_ds( 0 ), num_ts( 0 );
	Real frac_ss( 0.0 ), frac_ds( 0.0 ), frac_ts( 0.0 ), frac_overall( 0.0 );
	for (Size i = 1; i <= nres; i++ ) {
		frac_overall += recovery(i);
		switch ( struct_type(i) ) {
		case 0:
			num_ss++;
			frac_ss += recovery(i);
			break;
		case 1:
			num_ds++;
			frac_ds += recovery(i);
			break;
		case 2:
			num_ts++;
			frac_ts += recovery(i);
			break;
		}
	}

	if ( num_ss > 0.0 ) frac_ss /= num_ss;
	if ( num_ds > 0.0 ) frac_ds /= num_ds;
	if ( num_ts > 0.0 ) frac_ts /= num_ts;
	if ( nres > 0.0 ) frac_overall /= nres;

	std::map <Size, char > struct_symbol;
	struct_symbol[ 0 ] = 'S';
	struct_symbol[ 1 ] = 'D';
	struct_symbol[ 2 ] = 'T';

	utility::io::ozstream out( sequence_recovery_file );

	for (Size i = 1; i <= nres; i++ ) {
		out << pose.residue(i).name1() << I(3,i) << " " << F(8,3,recovery(i)) << " " << struct_symbol[ struct_type(i) ] << std::endl;
	}

	out << std::endl;
	out << "SINGLE_STRANDED " << I(3,num_ss) << " " << F(9,4,frac_ss) << std::endl;
	out << "DOUBLE_STRANDED " << I(3,num_ds) << " " << F(9,4,frac_ds) << std::endl;
	out << "TERTIARY_STRUCT " << I(3,num_ts) << " " << F(9,4,frac_ts) << std::endl;
	out << "OVERALL         " << I(3,nres) << " " << F(9,4,frac_overall) << std::endl;

	out.close();

	std::cout << "Wrote stats to: " << sequence_recovery_file << std::endl;

}



///////////////////////////////////////////////////////////////////////////////
void
rna_design_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	std::string pdb_file  = option[ in::file::s ][1];
	core::import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );
	protocols::rna::ensure_phosphate_nomenclature_matches_mini( pose );

	dump_pdb( pose, "start.pdb");
	pose::Pose save_pose( pose );

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line();

	if( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ){
		pack::task::parse_resfile(pose, *task);
	} else {
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			task->nonconst_residue_task( ii ).allow_aa( na_rad );
			task->nonconst_residue_task( ii ).allow_aa( na_ura );
			task->nonconst_residue_task( ii ).allow_aa( na_rgu );
			task->nonconst_residue_task( ii ).allow_aa( na_rcy );
			assert( task->design_residue(ii) );
		}
	}

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {

		//Hmmm, extras.
		task->nonconst_residue_task( ii ).and_extrachi_cutoff( 0 );

		if ( option[ disable_o2star_rotamers ]() ) task->nonconst_residue_task(ii).sample_proton_chi( false );

		if ( option[ sample_chi ]() ) task->nonconst_residue_task(ii).sample_rna_chi( true );

		// Can input this from command line:
		//		task->nonconst_residue_task( ii ).or_ex4( true );

		// Can input this from command line?
		//		task->nonconst_residue_task( ii ).or_ex1( true );

		// Screw this, can figure this out from command line.
		//		if ( !option[ disable_include_current ]() ) task->nonconst_residue_task( ii ).or_include_current( true );

	}

	ScoreFunctionOP scorefxn = getScoreFunction();

	//	scorefxn->energy_method_options().exclude_DNA_DNA( exclude_DNA_DNA );
	methods::EnergyMethodOptions options( scorefxn->energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn->set_energy_method_options( options );

	scorefxn->show( std::cout,pose );

	pose.dump_pdb( "start.pdb" );

	Size const nstruct = option[ out::nstruct ];
	utility::vector1< std::pair< Real, std::string > > results;
	utility::vector1< pose::PoseOP > pose_list;

	pack::pack_rotamers_loop( pose, *scorefxn, task, nstruct, results, pose_list);

	std::string outfile( pdb_file );
	Size pos( pdb_file.find( ".pdb" ) );
	outfile.replace( pos, 4, ".pack.txt" );
	protocols::rna::export_packer_results( results, pose_list, scorefxn, outfile, option[ dump ] );

	std::string sequence_recovery_file( pdb_file );
	sequence_recovery_file.replace( pos, 4, ".sequence_recovery.txt" );
	rna_sequence_recovery_metrics( save_pose, pose_list, sequence_recovery_file );

	//	std::string const out_file_tag = "S_"+lead_zero_string_of( n, 4 );
	//	dump_pdb( pose, out_file_tag + ".pdb" );

	pose = *( pose_list[1] );
	scorefxn->show( std::cout,pose );


}

///////////////////////////////////////////////////////////////
void
ss_ds_ts_assign_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::PoseOP pose_op( new pose::Pose );
	utility::vector1 < std::string > pdb_files  = option[ in::file::s ]();

	for (Size n = 1; n <= pdb_files.size(); n++ )  {
		std::string const & pdb_file = pdb_files[ n ] ;
		core::import_pose::pose_from_pdb( *pose_op, *rsd_set, pdb_file );
		protocols::rna::ensure_phosphate_nomenclature_matches_mini( *pose_op );

		std::string sequence_recovery_file( pdb_file );
		Size pos( pdb_file.find( ".pdb" ) );
		sequence_recovery_file.replace( pos, 4, ".ss_ds_ts.txt" );

		utility::vector1< pose::PoseOP > pose_list;
		pose_list.push_back( pose_op );
		rna_sequence_recovery_metrics( *pose_op, pose_list, sequence_recovery_file );
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( option[ ss_ds_ts_assign ] ) {
		ss_ds_ts_assign_test();
	} else {
		rna_design_test();
	}
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
	using namespace basic::options;

	std::cout << std::endl << "Basic usage:  " << argv[0] << "  -s <pdb file> [ -resfile <resfile>] " << std::endl;
	std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

	//Uh, options? MOVE THESE TO OPTIONS NAMESPACE INSIDE CORE/OPTIONS.
	NEW_OPT( disable_o2star_rotamers, "In designing, don't sample 2'-OH",false);
	NEW_OPT( disable_include_current, "In designing, don't include current",false);
	NEW_OPT( sample_chi,  "In designing RNA, chi torsion sample",false);
	NEW_OPT( ss_ds_ts_assign, "Figure out assignment of residues to single-stranded, double-stranded, tertiary contact categories",false);
	NEW_OPT( dump, "Dump pdb", false );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
}
