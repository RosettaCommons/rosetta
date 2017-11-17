// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file bur_nodon_tyr.cc
/// @brief
/// @author Kevin Houlihan (khouli@unc.edu)

// Unit headers
#include <devel/init.hh>

//project Headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/pose/PDBInfo.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

// Option keys
#include <basic/options/keys/run.OptionKeys.gen.hh>

// Numeric Headers

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>

// C++ headers
#include <sstream>

//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/vector0.hh>


static basic::Tracer TR( "bur_nodon_tyr" );

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace ObjexxFCL::fmt;

namespace bur_nodon_tyr {
FileOptionKey const pdb_list( "bur_nodon_tyr:pdb_list" );
FileOptionKey const parse_taskops_file( "bur_nodon_tyr:parse_taskops_file" );
RealOptionKey const probe_radius( "bur_nodon_tyr:probe_radius" );
BooleanOptionKey const layered_sasa( "bur_nodon_tyr:layered_sasa" );
RealOptionKey const sasa_bur_cutoff( "bur_nodon_tyr:sasa_bur_cutoff" );
}


/// @brief load custom TaskOperations according to an xml-like utility::tag file
core::pack::task::TaskFactoryOP setup_tf( core::pack::task::TaskFactoryOP task_factory_ ) {

	using namespace core::pack::task::operation;

	if ( option[ bur_nodon_tyr::parse_taskops_file ].user() ) {
		std::string tagfile_name( option[ bur_nodon_tyr::parse_taskops_file ]() );
		TaskOperationFactory::TaskOperationOPs tops;
		TaskOperationFactory::get_instance()->newTaskOperations( tops, tagfile_name );
		for ( TaskOperationFactory::TaskOperationOPs::iterator it( tops.begin() ), itend( tops.end() ); it != itend; ++it ) {
			task_factory_->push_back( *it );
		}
	} else {
		task_factory_->push_back( new pack::task::operation::InitializeFromCommandline );
	}

	return task_factory_;
}

/// @brief return the set of residues that are designable based given pose
std::set< Size > fill_designable_set( pose::Pose & pose, pack::task::TaskFactoryOP & tf ) {

	//we need to score the pose for many of the task operations passed from cmd line
	if ( option[ bur_nodon_tyr::parse_taskops_file ].user() ) {
		scoring::ScoreFunctionOP scorefxn = scoring::get_score_function();
		(*scorefxn)( pose );
	}
	std::set< Size > designable_set;
	core::pack::task::PackerTaskOP design_task( tf->create_task_and_apply_taskoperations( pose ) );

#ifndef NDEBUG
	//TR<< "Task for " << pose.pdb_info()->name() << " is: \n" << *(design_task)  << std::endl;
#endif

	// iterate over all residues
	for ( Size ii = 1; ii<= design_task->size(); ++ii ) {
		if ( design_task->being_designed( ii ) ) {
			designable_set.insert( ii );
		}
	}

	return designable_set;

}


/// @brief iterates over all designed positions and determines identity to native. outputs recoveries to file.
void measure_sequence_recovery( core::pose::Pose & pose ) {

	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	scorefxn->score( pose );

	// figure out the task & neighbor info
	core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
	std::set< Size > design_set;

	std::string sasa_calc_name("sasa_calc_name");
	basic::MetricValue< id::AtomID_Map< Real > > atom_sasa;

	if ( basic::options::option[bur_nodon_tyr::layered_sasa] ) {
		devel::vardist_solaccess::VarSolDistSasaCalculator varsol_sasa_calc;
		varsol_sasa_calc.get("atom_sasa", atom_sasa, pose);
	} else {
		pose::metrics::simple_calculators::SasaCalculator basic_sasa_calc;
		basic_sasa_calc.get("atom_sasa", atom_sasa, pose);
	}

	// setup what residues we are going to look at...
	setup_tf( task_factory );
	design_set = fill_designable_set( pose, task_factory );

	// record native sequence
	// native_sequence vector is sized for the WHOLE pose not just those being designed
	// it doesn't matter because we only iterate over the number of designed positions
	Size const nres( pose.size() );
	utility::vector1< chemical::AA > sequence( nres );

	Size bad_tyrs = 0;
	std::stringstream bad_tyr_ss;

	// iterate over designable positions
	for ( std::set< core::Size >::const_iterator it = design_set.begin(), end = design_set.end(); it != end; ++it ) {

		Size resi = *it;
		if ( pose.residue(resi).name3() != "TYR" ) continue;


		chemical::ResidueType rsd_type = pose.residue_type(resi);
		Size at_HH = rsd_type.atom_index("HH");
		Size at_OH = rsd_type.atom_index("OH");

		core::id::AtomID atid_HH( at_HH, resi);
		core::id::AtomID atid_OH( at_OH, resi);

		Real cursasa = ( atom_sasa.value()[ atid_HH ] + atom_sasa.value()[ atid_OH ] );

		if ( cursasa < basic::options::option[bur_nodon_tyr::sasa_bur_cutoff] ) {
			//TR << "Buried tyrosine OH at resi " << resi << ", checking for hbond donation." << std::endl;

			core::scoring::hbonds::HBondSet hbond_set;
			hbond_set.setup_for_residue_pair_energies( pose, false, false );

			core::scoring::hbonds::fill_hbond_set(pose, false, hbond_set, false, false, false, false);

			std::string res_debug = pose.pdb_info()->pose2pdb(resi);
			std::string atm_debug = pose.residue(resi).atom_name(static_cast< int >(at_HH));

			if ( hbond_set.nhbonds( atid_HH, false) > 0 ) {
				//TR << "Tyrosine HH donating at " << res_debug << atm_debug << std::endl;
			} else {
				bad_tyr_ss << "Buried tyrosine hydroxyl NOT donating at " << res_debug << std::endl;
				bad_tyrs++;
			}

		}

	} // end finding native seq

	if ( bad_tyrs > 0 ) {
		TR << "Bad tyrosines detected." << std::endl;
		TR << bad_tyr_ss.str();
		TR.flush();
	} else {
		TR << "No bad tyrosines detected." << std::endl;
	}

}


//@brief main method for the sequence recovery protocol
int main( int argc, char* argv[] ) {

	try {

		using utility::file::file_exists;
		using utility::file::FileName;

		option.add( bur_nodon_tyr::pdb_list, "List of pdb files." );
		option.add( bur_nodon_tyr::parse_taskops_file, "XML file which contains task operations to apply before measuring recovery (optional)" );
		option.add( bur_nodon_tyr::probe_radius, "Probe radius for SASA calculation." ).def( 1.2 );
		option.add( bur_nodon_tyr::layered_sasa, "Use vardist sasa calc").def( 'true' );
		option.add( bur_nodon_tyr::sasa_bur_cutoff, "Maximum SASA to be considered 'buried'." ).def( 0.01 );

		devel::init( argc, argv );

		// changing this so that native_pdb_list and redesign_pdb_list do not have default values. giving these options can lead
		// to users measuring recovery against the wrong set of PDBs.
		if ( !option[ bur_nodon_tyr::pdb_list ].user() ) {
			utility_exit_with_message_status( "no pdb_list given", 1 );
		}

		// read list file. open the file specified by the flag 'pdb_list' and read in all the lines in it
		std::vector< FileName > pdb_file_names;
		std::string pdb_list_file_name( option[ bur_nodon_tyr::pdb_list ].value() );
		std::ifstream pdb_data( pdb_list_file_name.c_str() );
		std::string pdb_line;
		if ( !pdb_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_file_name + '\n' );
		}
		while ( getline( pdb_data, pdb_line ) ) {
			pdb_file_names.push_back( FileName( pdb_line ) );
		}

		pdb_data.close();

		// iterate over both FileName vector and read in the PDB files
		utility::vector1< pose::Pose > poses;

		std::vector< FileName >::iterator pdb( pdb_file_names.begin() ), last_pdb(pdb_file_names.end());

		while ( pdb != last_pdb ) {

			// check to make sure the file exists
			if ( !file_exists( *pdb ) ) {
				utility_exit_with_message( "Pdb " + std::string(*pdb) + " not found! skipping" );
			}

			TR << "Reading in pose " << *pdb << std::endl;
			core::pose::Pose pose;
			core::import_pose::pose_from_file( pose, *pdb , core::import_pose::PDB_file);
			//TR << "Looking for bad tyrosines." << std::endl;
			measure_sequence_recovery( pose );

			pdb++;
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}

