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


// libRosetta headers

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/scoring/sasa.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/excn/Exceptions.hh>

#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
//#include <core/chemical/residue_io.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.cc>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <protocols/cartesian/md.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>

//#include <core/mm/MMTorsionLibrary.hh>
//#include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


#include <basic/Tracer.hh>


using namespace core;
using namespace protocols;

using utility::vector1;




class Protocol_MolecularDynamics: public moves::Mover {
public:
	Protocol_MolecularDynamics(
		core::scoring::ScoreFunctionOP scorefxn_in
	) : Mover(),
			scorefxn_(scorefxn_in),
			native_set_(false)
	{

	}

	virtual std::string get_name() const { return "MolecularDynamics"; }
	virtual moves::Mover* clone() { return new Protocol_MolecularDynamics( *this );  }

	virtual void apply( core::pose::Pose & inpose ){

		pose::PoseOP pose = new pose::Pose ( inpose );

		kinematics::MoveMap mm;
		mm.set_bb (  true );
		mm.set_chi(  true );

		using namespace optimization;
		using namespace protocols::cartesian;

		(*scorefxn_)(*pose);
		scorefxn_->show(std::cout, *pose);

		// setup the options
		MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.000010, true ,
				false , false );
		AtomTreeMinimizer minimizer;
		std::cout << "MINTEST: p_aa_pp" << "\n";
		std::cout << "start score: " << (*scorefxn_)( *pose ) << "\n";
		minimizer.run( *pose, mm, *scorefxn_, options );
		pose->dump_pdb( "min_intrares.pdb" );
		std::cout << "end score: " << (*scorefxn_)( *pose ) << "\n";

		pose->energies().clear();
		(*scorefxn_)(*pose);
		scorefxn_->show(std::cout, *pose);

		MolecularDynamics md( pose, *scorefxn_ );
		//md.testCartesianDerivatives( scorefxn_ );
		md.doMD( *scorefxn_, 10000, 150, 1 );
		md.doMD( *scorefxn_, 1000, 1, 0.00001 );

		pose->energies().reset_nblist();
		(*scorefxn_)(*pose);
		scorefxn_->show(std::cout, *pose);
		pose->energies().clear(  );
		(*scorefxn_)(*pose);
		scorefxn_->show(std::cout, *pose);

	}

	void set_native_pose( core::pose::PoseOP & n_pose ){ native_pose_ = n_pose; native_set_ = true; }
protected:
	// protocol stuff
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pose::PoseOP native_pose_;
	bool native_set_;

	float rms_to_native(  core::pose::Pose & pose );
};


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
	devel::init(argc, argv);

	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::scoring::ScoreFunctionOP scorefxn;
	scorefxn = get_score_function();

	// Build overall docking protocol Mover
	Protocol_MolecularDynamics *md;
	md = new Protocol_MolecularDynamics( scorefxn );

	if(  option[ OptionKeys::in::file::native ].active() ){
		pose::PoseOP npose = new pose::Pose;
		core::import_pose::pose_from_file( *npose, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file); // default is standard fullatom residue_set
		md->set_native_pose( npose );
	}

	MoverOP protocol = md;

	try{
    protocols::jd2::JobDistributor::get_instance()->go( protocol );
  } catch ( utility::excn::EXCN_Base& excn ) {
    std::cerr << "Exception: " << std::endl;
    excn.show( std::cerr );
    std::cout << "Exception: " << std::endl;
    excn.show( std::cout ); //so its also seen in a >LOG file
  }

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
        return -1;
    }
    return 0;
}
