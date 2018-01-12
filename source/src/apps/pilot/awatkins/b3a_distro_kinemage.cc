// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// awatkins: based heavily on kdrew/oop_creator.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>

//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

//kdrew: this app adds hbs patches to the given pdb strucure

// tracer - used to replace cout
static basic::Tracer TR( "B3A_distro" );


// application specific options
namespace b3a_distro {
// pert options
StringOptionKey const b3aa ( "b3a_distro::b3aa" );
RealOptionKey const bin_size ( "b3a_distro::bin_size" );
RealOptionKey const print_cutoff ( "b3a_distro::print_cutoff" );
StringOptionKey const scoremethod ( "b3a_distro::scoremethod" );

}

class HbsCreatorMover : public Mover {

public:

	//default ctor
	HbsCreatorMover(): Mover("HbsCreatorMover"){}

	//default dtor
	~HbsCreatorMover() override= default;

	void apply( core::pose::Pose & pose ) override;
	std::string get_name() const override { return "HbsCreatorMover"; }

};

using HbsCreatorMoverOP = utility::pointer::shared_ptr<HbsCreatorMover>;
using HbsCreatorMoverCOP = utility::pointer::shared_ptr<const HbsCreatorMover>;


int
main( int argc, char* argv[] )
{
	try {
		utility::vector1< core::Size > empty_vector(0);

		option.add( b3a_distro::b3aa, "Chain from PDB to be mimicked. Default 'A'. Use letters." ).def("B3A");
		option.add( b3a_distro::bin_size, "Residue number of the final residue for mimicry. Default 10." ).def(10);
		option.add( b3a_distro::print_cutoff, "Residue number of the final residue for mimicry. Default 0." ).def(0);
		option.add( b3a_distro::scoremethod, "Method: repack_min, min, score." ).def("score");

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		ScoreFunctionOP scorefxn = get_score_function();
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		Pose pose;
		Real bin_size = option[ b3a_distro::bin_size ]();
		Real score_cutoff = option[ b3a_distro::print_cutoff ]();
		std::string aa_name = option[ b3a_distro::b3aa ]();
		std::string method = option[ b3a_distro::scoremethod ]();

		core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory );
		task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
		//need these to keep pack_rotamers from redesigning the residue.
		operation::RestrictToRepackingOP rtrop ( new operation::RestrictToRepacking );
		task_factory->push_back( rtrop );

		ResidueType const & restype_first = residue_set_cap->name_map( "ACE" );
		ResidueType const & restype_cent = residue_set_cap->name_map( aa_name );
		ResidueType const & restype_last = residue_set_cap->name_map( "NME" );
		Residue res_first( restype_first, true );
		Residue res_cent( restype_cent, true );
		Residue res_last( restype_last, true );
		pose.append_residue_by_jump( res_first, 1 );
		pose.append_residue_by_bond( res_cent, true );
		pose.append_residue_by_bond( res_last, true );

		id::TorsionID eps( 1, id::BB, 2 ); //epsilon
		id::TorsionID bb2( 2, id::BB, 1 ); //phi
		id::TorsionID bb3( 2, id::BB, 2 ); //theta
		id::TorsionID bb4( 2, id::BB, 3 ); //psi
		id::TorsionID omg( 2, id::BB, 4 ); //omg

		pose.set_torsion( eps, -180.0 );
		pose.set_torsion( omg, -180.0 );

		kinematics::MoveMapOP movemap_sc( new kinematics::MoveMap );
		movemap_sc->set_chi(2, true);
		movemap_sc->set(eps, false);
		movemap_sc->set(bb2, false);
		movemap_sc->set(bb3, false);
		movemap_sc->set(bb4, false);
		movemap_sc->set(omg, false);

		protocols::minimization_packing::MinMover minmover_sc(movemap_sc, scorefxn, "lbfgs_armijo_nonmonotone", 0.0001, true);

		ObjexxFCL::FArray3D< Real > energies( 36, 36, 36 );
		Real minscore = 9999999999;
		Size i = 1, j = 1, k = 1;
		for ( Real phi = -180; phi < 180; phi += bin_size )  {
			//std::cout << "phi = " << phi << std::endl;
			for ( Real tht = -180; tht < 180; tht += bin_size )  {
				for ( Real psi = -180; psi < 180; psi += bin_size )  {
					pose.set_torsion( bb2, phi );
					pose.set_torsion( bb3, tht );
					pose.set_torsion( bb4, psi );

					if ( method == "repack_min" ) {
						PackerTaskOP packer_task = task_factory->create_task_and_apply_taskoperations ( pose );
						utility::vector1< bool > repackable ( 3, false );
						repackable[2] = true;
						packer_task->restrict_to_residues( repackable );
						core::pack::pack_rotamers ( pose, *scorefxn, packer_task );
					}
					if ( method == "repack_min" || method == "min" ) {
						minmover_sc.apply ( pose );
					}

					Real score = ( * scorefxn ) ( pose );
					energies( i, j, k++ ) = score;
					if ( score < minscore ) minscore = score;
					if ( score < score_cutoff ) {

						//std::cout << score << "\t";
						// cout point; for first draft don't print to kinemage
						//std::cout << "pt " << phi << " " << tht << " " << psi << std::endl;
					} else {
						//std::cout << "99999999" << "\t";
					}

					if ( k > 36 ) {
						k = 1;
						j++;
					}
					if ( j > 36 ) {
						j = 1;
						i++;
					}
				} //std::cout << std::endl;
			} //std::cout << std::endl << std::endl;
		}

		std::cout << aa_name << " minscore is " << minscore << std::endl << std::endl;
		// what we print is an energy range from min
		for ( Real phi = -180; phi < 180; phi += bin_size )  {
			std::cout << "phi = " << phi << std::endl;
			for ( Real tht = -180; tht < 180; tht += bin_size )  {
				for ( Real psi = -180; psi < 180; psi += bin_size )  {
					std::cout << ( energies( i, j, k++ ) - minscore ) << "\t";

					if ( k > 36 ) {
						k = 1;
						j++;
					}
					if ( j > 36 ) {
						j = 1;
						i++;
					}
				} std::cout << std::endl;
			} std::cout << std::endl << std::endl;
		}

		//call job distributor
	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
