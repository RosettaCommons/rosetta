// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/krishna/mtransferase.cc
/// @Build loops for a zinc finger - mtase complex?

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/pHEnergy.hh>
#include <core/scoring/Energies.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <basic/options/option.hh>
#include <core/init.hh>
#include <core/conformation/Residue.hh>

#include <protocols/relax/ClassicRelax.hh>
#include <protocols/moves/PackRotamersMover.hh>
#include <protocols/moves/SidechainMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <basic/Tracer.hh>
#include <string>
#include <utility/vector1.hh>

//neighbors & hbonds
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>

//temp includes krishna
#include <iostream>
#include <fstream>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/string_util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>

//JD headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>


// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

static basic::Tracer TR( "apps.pilot.krishna.mtransferase" );

class mtransferase : public protocols::moves::Mover {

public:
	mtransferase()
	{
	}

	virtual ~mtransferase(){};

	virtual
	void
	apply( core::pose::Pose & pose ){

		using namespace core;
		using namespace chemical;
		using namespace conformation;
		using namespace scoring;
		using namespace io::pdb;

		core::chemical::ResidueTypeSetCAP fa_standard(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));

		//accumulate required residues for loop 1
		for ( Size i=0; i<15 ; ++i ) {
			if ( i%5 == 0 ) pose.append_polymer_residue_after_seqpos( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("SER")) ), i+85, false );
			else pose.append_polymer_residue_after_seqpos( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("GLY")) ), i+85, false );
		}

		//accumulate required residues for loop 2
		for ( Size j=0; j<15 ; ++j ) {
			if ( (j+1)%5 == 0 ) pose.append_polymer_residue_after_seqpos( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("SER")) ), j+427, false );
			else pose.append_polymer_residue_after_seqpos( *(core::conformation::ResidueFactory::create_residue(fa_standard->name_map("GLY")) ), j+427, false );
		}

		//  pose.dump_pdb("loop1_complex.pdb");
	}

	std::string get_name() const { return "mtransferase"; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const {
		return new mtransferase;
	}

	virtual
	bool
	reinitialize_for_each_job() const { return false; }

	virtual
	bool
	reinitialize_for_new_input() const { return false; }

	//private:

};

typedef utility::pointer::owning_ptr< mtransferase > mtransferaseOP;

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	core::init(argc, argv);
	mtransferaseOP mtase_test(new mtransferase);
	protocols::jd2::JobDistributor::get_instance()->go(mtase_test);
}


