// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @details
///

#include <devel/init.hh>
#include <core/types.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/ResidueKinWriter.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/Atom.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <cstdlib>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace basic::options;
using namespace core::scoring;

// Tracer will output messages to user
basic::Tracer TR("virtuals");


int
main( int argc, char * argv [] )
{
		
	devel::init( argc, argv );
	// read in pose
	std::string fname = option[ OptionKeys::in::file::s ]()[1];
	core::pose::Pose p;
    core::import_pose::pose_from_file(p,fname, core::import_pose::PDB_file);
	
	//work on a copy of the restype. Must be done in drug design. Later we will add the restype to the pose
	core::chemical::ResidueTypeOP newRT( new core::chemical::ResidueType ( p.residue_type(3) ) );
	
	/*// Change to virtuals	
	std::string VIRT = "VIRT";
	for(Size i=1; i<=newRT->natoms(); ++i){
		newRT->set_atom_type( (newRT->atom_name(i) ), VIRT);
		newRT->atom(i).charge(0.0);
		newRT->atom(i).is_virtual( true );
	}
	// VERY IMPORTANT to get the carbohydrate_info() for newRT
	newRT->finalize();
	*/
	//newRT->fa_to_virt();

	//core::pose::replace_pose_residue_copying_existing_coordinates(p,3,*newRT);
	ScoreFunctionOP score_fxn = get_score_function();

	score_fxn->show(TR, p);

	TR << "Should be FA " << std::endl;
	for (core::Size i = 1; i <= p.residue( 3 ).natoms(); ++i ){
		TR << i << " " <<p.residue_type( 3 ).is_virtual( i ) << std::endl;
	}


	p.real_to_virtual(3);
	
	TR << "Should be virt " << std::endl;
	for (core::Size i = 1; i <= p.residue_type( 3 ).natoms(); ++i ){
		TR << i << " " << p.residue_type( 3 ).is_virtual( i ) << std::endl;
	}

	(*score_fxn)(p);
	score_fxn->show(TR, p);

	// Reverse
	core::chemical::ResidueTypeOP oldRT( new core::chemical::ResidueType ( p.residue_type(3) ) );
	std::string const base_name( residue_type_base_name( *oldRT ) );
    chemical::ResidueTypeSetCOP base_RTset( oldRT->residue_type_set() );
	core::chemical::ResidueTypeCOP RT = core::chemical::ResidueTypeFinder(*base_RTset).residue_base_name( base_name ).get_representative_type(); 
	// swap in new (old) residue type
	core::pose::replace_pose_residue_copying_existing_coordinates(p,3,*RT);

	score_fxn->show(TR, p);

	TR << "Should be FA " << std::endl;
	for (core::Size i = 1; i <= p.residue( 3 ).natoms(); ++i ){
		TR << i << " " << p.residue( 3 ).is_virtual( i ) << std::endl;
	}

	p.dump_pdb("testing.pdb");

}

