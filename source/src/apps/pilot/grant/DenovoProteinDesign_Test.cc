// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Grant


// Unit headers
#include <devel/init.hh>
#include <core/pack/task/PackerTask_.hh>

//project Headers
#include <core/io/pdb/pdb_writer.hh>
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
#include <basic/MetricValue.hh>
#include <protocols/loops/loops_main.hh> //for getting ss from dssp
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax_protocols.hh>


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
#include <devel/denovo_protein_design/CreateStartingStructureMover.hh>
#include <devel/denovo_protein_design/DesignRelaxMover.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/MoverContainer.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using utility::file::FileName;

	devel::init(argc, argv);
	std::cout << " Test Application for creating starting structures for de novo protein design " << std::endl;
	std::cout << " 2 styles of making starting structures - from a template pdb or from a secondary structure profile " << std::endl;
	std::cout << " This mover explores sequence and conformation space for the starting structure being generated " << std::endl;
	std::cout << " Running options " << std::endl;
	std::cout << " -use_template_sequence (doesn't search sequence space), -use_template_topology (limits conformational search) " << std::endl;
	std::cout << " -create_from_secondary_structure -secondary_structure_file SS.ss (file contains H/E/L to describe desired decoy secondary structure) " << std::endl;
	std::cout << " you can specify fragments  - default is to read from the VALL and choose fragments on the fly ( working in progress )" << std::endl;

	core::pose::Pose pose;

	core::scoring::ScoreFunctionOP fullfxn( ( core::scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) ));

	devel::denovo_protein_design::CreateStartingStructureMoverOP myCreateStartingStructureMover ( new devel::denovo_protein_design::CreateStartingStructureMover() );
	devel::denovo_protein_design::DesignRelaxMoverOP myDesignRelaxMover ( new devel::denovo_protein_design::DesignRelaxMover() );


	/*
	if( option [ in::file::s ].active() or  option[ in::file::l ].active() ){

	  std::vector< FileName > pdb_file_names;
	  if ( option[ in::file::s ].active() ) pdb_file_names = option[ in::file::s ]().vector();
	  std::vector< FileName > list_file_names;
	  if ( option[ in::file::l ].active() ) list_file_names = option[ in::file::l ]().vector();

	  for(std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i)
	    {
	      std::string filename( i->name() );

	      std::ifstream data( filename.c_str() );

	      if ( !data.good() ) utility_exit_with_message( "Unable to open file: " + filename + '\n' );

	      std::string line;
	      while( getline(data, line) )
		{
		  pdb_file_names.push_back( FileName(line) );
		}
	      data.close();
	    }

	  for( std::vector< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i ){

	    core::import_pose::pose_from_file( pose, *i, core::import_pose::PDB_file);

	  protocols::moves::SequenceMover DenovoDesignProtocol_template( myCreateStartingStructureMover, myDesignRelaxMover );
	  protocols::jobdist::main_plain_pdb_mover(DenovoDesignProtocol_template, fullfxn);
	}

	} else {
	*/

        protocols::moves::SequenceMover DenovoDesignProtocol( myCreateStartingStructureMover, myDesignRelaxMover );
	//	DenovoDesignProtocol.apply( pose );
        protocols::jobdist::main_plain_pdb_mover(DenovoDesignProtocol, fullfxn);

	//	}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

  return 0;
}
