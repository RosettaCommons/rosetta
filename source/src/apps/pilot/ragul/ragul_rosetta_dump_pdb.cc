// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author Ragul Gowthaman

#include <cfloat>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>

#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>

#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/pdb_dynamic_reader.hh>

#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/access_ptr.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;
using namespace basic::options::OptionKeys;

static thread_local basic::Tracer TR( "apps.pilot.dump_pdb.main" );


int main( int argc, char * argv [] )
{
	try{
		devel::init(argc, argv);

		// create native pose from pdb
		pose::Pose pose_init;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_pdb( pose_init, input_pdb_name );
		std::string out_pdb_name = "rosetta_" + input_pdb_name;
		pose_init.dump_pdb(out_pdb_name);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

