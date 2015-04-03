// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RedesignCrystalStructure - send pdbs to silent file and filter decoys during run to give
/// @brief Score RMSD PackStats SeqSim SeqDiv and counts of each AA
/// @author Grant


// Unit headers
#include <devel/init.hh>
#include <core/pack/task/PackerTask_.hh>

//project Headers
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <core/pose/util.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/rms_util.hh>

#include <protocols/evaluation/RmsdEvaluator.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <protocols/analysis/PackStatMover.hh>

#include <core/chemical/AA.hh>


#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax_protocols.hh>

#include <utility/io/util.hh>
#include <core/io/silent/PDBSilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <basic/options/option.hh>
// Utility Headers
#include <utility/file/file_sys_util.hh>

// Numeric Headers

// ObjexxFCL Headers

// C++ headers
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <utility/assert.hh> //ASSERT_ONLY makes release build happy
#include <core/chemical/AA.hh>
#include <devel/RedesignCrystalStructure/RedesignCrystalStructureMover.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/JobDistributors.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{

	try {

devel::init(argc, argv);

using namespace core::scoring::packstat;
// set up scorefxn
core::scoring::ScoreFunctionOP scorefxn( ( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) ));

// set up calculators


 core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy;
 core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

 core::pose::metrics::PoseMetricCalculatorOP num_hbonds_calculator = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
 core::pose::metrics::CalculatorFactory::Instance().register_calculator( "num_hbonds", num_hbonds_calculator );

 core::pose::metrics::PoseMetricCalculatorOP unsat_calculator = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator("sasa", "num_hbonds");
 core::pose::metrics::CalculatorFactory::Instance().register_calculator( "unsat", unsat_calculator );

  // set up some containers to hold data
 std::map< core::Real, std::string > SequenceMap;
 SequenceMap.insert( std::pair<core::Real, std::string>(  1, "A") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  2, "C") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  3, "D") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  4, "E") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  5, "F") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  6, "G") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  7, "H") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  8, "I") );
 SequenceMap.insert( std::pair<core::Real, std::string>(  9, "K") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 10, "L") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 11, "M") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 12, "N") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 13, "P") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 14, "Q") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 15, "R") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 16, "S") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 17, "T") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 18, "V") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 19, "W") );
 SequenceMap.insert( std::pair<core::Real, std::string>( 20, "Y") );

 std::map< std::string, core::Real > nativeRealmap;
 std::map< std::string, std::map< std::string , core::Real > > nativeParentMap;
 std::map<std::string, std::map< std::string, core::Real> >::iterator nativeParentMapIterator;

 std::map< std::string, core::Real > Realmap;
 std::map< std::string, std::map< std::string , core::Real > > ParentMap;
 std::map<std::string, std::map< std::string, core::Real> >::iterator ParentMapIterator;

 std::map< std::string, std::map< std::string , core::Real > > DecoysBetterThanNativeMap;


 std::string outname( basic::options::option[ basic::options::OptionKeys::in::file::native ]().name() + "_data.txt" );
 std::ofstream out(outname.c_str(), std::ios::out | std::ios::binary);

  std::string tag;

  core::io::silent::SilentFileDataOP silent_score_file_;
  silent_score_file_ = new core::io::silent::SilentFileData();
  silent_score_file_->set_filename( basic::options::option [ basic::options::OptionKeys::out::file::silent ]() );


  // do some stuff to the native

  core::pose::Pose nativepose;
  core::import_pose::pose_from_pdb( nativepose, basic::options::option[ core::options::OptionKeys::in::file::native ]().name()  );

  core::Real seqlength = nativepose.n_residue();
  std::string nativesequence = nativepose.sequence();
  utility::vector1< core::Real  >nativeseqdis;
  nativeseqdis.resize(20,0);
  // get counts of the number of each amino acid - this a vector1 of size 20 so element one is A and element 20 is Y
  for( core::Size ii = 1; ii<=seqlength; ++ii){
  core::Size aa(nativepose.residue_type(ii).aa());
		 ++nativeseqdis[aa];
	 }
  for( core::Size ii = 1; ii <=nativeseqdis.size(); ++ii){
  nativeRealmap.insert( std::pair< std::string, core::Real >( SequenceMap[ii], nativeseqdis[ii] ) );
}

  // get score
core::Energy nativescore = (*scorefxn)( nativepose );
/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( nativepose );
// get packstats
PosePackData nativepd = pose_to_pack_data( nativepose, false );
core::Real nativepackstats = compute_packing_score( nativepd, 0 );
// get unsatisfied hbonds
std::string nativeunsathbonds_string = nativepose.print_metric("unsat", "all_bur_unsat_polars");
core::Real nativeunsathbonds;

std::stringstream convert_ss( nativeunsathbonds_string );
convert_ss >> nativeunsathbonds;

nativeRealmap.insert( std::pair< std::string, core::Real >( "UNSAT HBONDS", nativeunsathbonds ) );
nativeRealmap.insert( std::pair< std::string, core::Real >( "RMSD", 0.000 ) );
nativeRealmap.insert( std::pair< std::string, core::Real >( "SCORE", nativescore ) );
nativeRealmap.insert( std::pair< std::string, core::Real >( "PACKSTATS", nativepackstats ) );
nativeRealmap.insert( std::pair< std::string, core::Real >( "SEQSIM", 1.0 ) );
nativeRealmap.insert( std::pair< std::string, core::Real >( "SEQDIVERISTY", 0 ) );

std::pair< std::string, std::map< std::string, core::Real> > nativerealpair( basic::options::option[ basic::options::OptionKeys::in::file::native ]().name(), nativeRealmap);
nativeParentMap.insert( nativerealpair );


	out << " Native " << std::endl;
  for ( ParentMapIterator = nativeParentMap.begin(); ParentMapIterator != nativeParentMap.end(); ++ParentMapIterator){
		 out << (*ParentMapIterator).first << " "
               << nativeParentMap[(*ParentMapIterator).first]["SCORE"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["RMSD"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["PACKSTATS"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["UNSAT HBONDS"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["SEQSIM"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["SEQDIVERSITY"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["A"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["C"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["D"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["E"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["F"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["G"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["H"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["I"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["K"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["L"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["M"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["N"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["P"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["Q"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["R"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["S"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["T"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["V"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["W"] << " "
               << nativeParentMap[(*ParentMapIterator).first]["Y"] << " "
 << std::endl;
	}

  out << "------------------------------------------------------" << std::endl;

  core::pose::Pose pose;

  // do the design steps - printing to silent file after each decoy is done

  core::Size nstructs = basic::options::option [ basic::options::OptionKeys::out::nstruct]();

  for( core::Size ii = 1; ii <= nstructs; ++ii){

 pose = nativepose;

 // this is where the magic happens
 devel::RedesignCrystalStructure::RedesignCrystalStructureMover myCrystalStructureRedesign;
 myCrystalStructureRedesign.apply( pose ); // not using job distributor because I am doing my own silent file output

 tag = lead_zero_string_of( ii, 6 );

 core::io::silent::ProteinSilentStruct ss;
 ss.fill_struct( pose, tag );
 // write all decoys to the silent file incase we want to check stats after run
 silent_score_file_->write_silent_struct( ss,  silent_score_file_->filename(), false );

 // if the decoy is comparable or better than wild type print its score info to a seperate out file
 std::string currentsequence = pose.sequence();

 core::Energy score = (*scorefxn)( pose );
 /// Now handled automatically.  scorefxn->accumulate_residue_total_energies( pose );


 std::string unsathbonds_string = pose.print_metric("unsat", "all_bur_unsat_polars");
 core::Real unsathbonds;
 std::stringstream convert_ss( unsathbonds_string );
 convert_ss >> unsathbonds;

 core::Real rmsd = core::scoring::CA_rmsd(pose, nativepose);

 core::Real SeqSim = 0;

 utility::vector1< core::Real > currentseqdis;
 currentseqdis.resize(20,0);

 for( core::Size ii = 1; ii<=seqlength; ++ii){
 core::Size aa(pose.residue_type(ii).aa());
 ++currentseqdis[aa];
 }

 for( core::Size ii = 1; ii <= currentsequence.size(); ++ii ){

 if( currentsequence[ii-1] == nativesequence[ii-1]){ ++SeqSim; }

 }

 for( core::Size ii = 1; ii <=currentseqdis.size(); ++ii){
 Realmap.insert( std::pair< std::string, core::Real >( SequenceMap[ii], currentseqdis[ii] ) );
}

 PosePackData pd = pose_to_pack_data( pose, false );
 core::Real packstats = compute_packing_score( pd, 0 );
 core::Real diversity = 0;

 for( core::Size ii = 1; ii <=currentseqdis.size(); ++ii){
 diversity += std::abs(currentseqdis[ii] - nativeseqdis[ii]);
 }


 Realmap.insert( std::pair< std::string, core::Real >( "UNSAT HBONDS", unsathbonds ) );
 Realmap.insert( std::pair< std::string, core::Real >( "RMSD", rmsd ) );
 Realmap.insert( std::pair< std::string, core::Real >( "SCORE", score ) );
 Realmap.insert( std::pair< std::string, core::Real >( "PACKSTATS", packstats ) );
 Realmap.insert( std::pair< std::string, core::Real >( "SEQSIM", SeqSim/seqlength ) );
 Realmap.insert( std::pair< std::string, core::Real >( "SEQDIVERSITY", diversity ) );


 std::pair< std::string, std::map< std::string, core::Real> > realpair( tag, Realmap);
 // contains all data for all decoys
 ParentMap.insert( realpair );

 // now lets keep a map of only the best decoys
 core::Real sig = 0.3;

	 if( rmsd < 6 ){
	 if( packstats > nativepackstats){
	 if( unsathbonds < nativeunsathbonds ){
	   // lets check all small amino acids A,G,S,T to make sure they are not over represented

if( currentseqdis[1] <= nativeseqdis[1]+sig*nativeseqdis[1] and currentseqdis[6] <= nativeseqdis[6]+sig*nativeseqdis[6] and currentseqdis[16] <= nativeseqdis[16]+sig*nativeseqdis[16] and currentseqdis[17] <= nativeseqdis[17]+sig*nativeseqdis[17]){
 DecoysBetterThanNativeMap.insert( realpair );
 }
	 }
	 }
	 }

  for ( ParentMapIterator = DecoysBetterThanNativeMap.begin(); ParentMapIterator != DecoysBetterThanNativeMap.end(); ++ParentMapIterator){
		 out << (*ParentMapIterator).first << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["SCORE"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["RMSD"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["PACKSTATS"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["UNSAT HBONDS"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["SEQSIM"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["SEQDIVERSITY"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["A"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["C"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["D"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["E"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["F"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["G"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["H"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["I"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["K"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["L"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["M"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["N"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["P"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["Q"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["R"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["S"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["T"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["V"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["W"] << " "
               << DecoysBetterThanNativeMap[(*ParentMapIterator).first]["Y"] << " "
 << std::endl;
	}

	 Realmap.clear();
	 currentseqdis.empty();

  } // end nstruct for loop


  // when we are done will everything we can do some more stuff

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}
