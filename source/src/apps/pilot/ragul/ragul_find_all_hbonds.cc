// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief
/// @author jk

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
#include <core/pose/PDB_Info.hh>
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

static basic::Tracer TR( "apps.pilot.hbonds.main" );


int main( int argc, char * argv [] )
{
	try{

	devel::init(argc, argv);

	// create native pose from pdb
	pose::Pose pose_init;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_pdb( pose_init, input_pdb_name );

	// fill hbond set
	core::scoring::hbonds::HBondSet hbs1;
	scoring::ScoreFunctionOP scorefxn(get_score_function());
	(*scorefxn)(pose_init);
	core::scoring::hbonds::fill_hbond_set(pose_init, false, hbs1);

	TR << "\n\nHBOND SET:\n";
	for ( Size i=1; i<= Size(hbs1.nhbonds()); ++i ) {
		core::scoring::hbonds::HBond const & hb( hbs1.hbond(i) );
		core::conformation::Residue const & donor =	pose_init.residue( hb.don_res() );
		core::conformation::Residue const & accep =	pose_init.residue( hb.acc_res() );
		TR << i << ":" <<
			chemical::oneletter_code_from_aa(donor.aa()) <<
			pose_init.pdb_info()->number(donor.seqpos()) << pose_init.pdb_info()->chain(donor.seqpos()) << ' ' <<
			'(' << donor.seqpos() << ')' <<
			donor.atom_name( hb.don_hatm()) << " --- " <<
			chemical::oneletter_code_from_aa(accep.aa()) <<
			pose_init.pdb_info()->number(accep.seqpos()) << pose_init.pdb_info()->chain(accep.seqpos()) << ' ' <<
			'(' << accep.seqpos() << ')' <<
			accep.atom_name( hb.acc_atm()) << "\n";
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
}

