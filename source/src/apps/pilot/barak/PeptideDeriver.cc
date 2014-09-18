// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University
// (C) 199x-2008 Hebrew University, Jerusalem
//
/// @file   PeptideDeriver.cc
//
/// @brief Application that reads in dimer (composed of two chains) and outputs the peptide which contibutes most to the interface.

/// @author Nir London
/// @date Nov. 15, 2009

//#define GL_GRAPHICS

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/util/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/options/util.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/disulfide_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// C++ headers

using core::util::T;
using core::util::Error;
using core::util::Warning;


static core::util::Tracer TR( "PeptideDeriver" );

using namespace core;
using namespace core::options;
using namespace OptionKeys;

core::pose::PoseOP get_subPose(core::pose::PoseCOP, core::pose::PoseCOP,
														 core::Size,core::Size,char);
Real calculateInterfaceScore (core::pose::Pose&, core::scoring::ScoreFunctionCOP );

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace pose;
		using namespace scoring;
		using namespace conformation;
		using namespace core::chemical;
		using namespace protocols;
		using namespace core::scoring::constraints;

		//setup random numbers and options
		devel::init(argc, argv);
		//setup scorefxn
		core::scoring::ScoreFunctionOP scorefxn(get_score_function());
		//create a pose
		pose::Pose pose;
		io::pdb::pose_from_pdb(pose, options::start_file());

		//////////////// set native constraints prior to minimization
		ConstraintSetOP cst_set( new ConstraintSet() );
		HarmonicFuncOP spring = new HarmonicFunc( 0 /*mean*/, 1 /*std-dev*/);
		conformation::Conformation const & conformation( pose.conformation() );
		for (Size i=1; i <= pose.total_residue(); ++i) {
				Residue const  & reside = pose.residue( i );
				id::AtomID CAi ( pose.residue(i).atom_index( " CA " ), i );
				cst_set->add_constraint
						(  new CoordinateConstraint
							 ( CAi, CAi, conformation.xyz( CAi ), spring )
							 );
		}
		pose.constraint_set( cst_set );
		scorefxn->set_weight(coordinate_constraint, 1.0);
		/////////////////////////////////////////////////////////////////////

		//minimize pose
		TR << "minimizing" << std::endl;
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(false);
		protocols::simple_moves::MinMover minimizer(movemap, scorefxn, "dfpmin_armijo_atol", 0.01, true /*nb_list*/);

		//////////////////////!!!!!!!!!!!!1111
		/////////////!!!!!!!!!!!!!!!!
		/////////////////********   FOR DEBUG THIS IS COMMENTED OUT
		//		minimizer.apply(pose);

		//pose.dump_pdb("./after_minimize.pdb");

		//remove constraint weight from scorefxn
		scorefxn->set_weight(coordinate_constraint, 0.0);

		//spilt into two chains
		utility::vector1< PoseOP > chains = pose.split_by_chain();
		pose::PoseOP chainA = chains[1];
		pose::PoseOP chainB = chains[2];
		TR << "INITIAL: let's score chain-A two times:" << std::endl;
		(*scorefxn)(*chainA);
		(*scorefxn)(*chainA);


		conformation::Conformation& a = chainA->conformation();
		conformation::Conformation& b = chainB->conformation();
		core::Size a_start = a.chain_begin(1);
		core::Size a_end = a.chain_end(1);
		core::Size b_start = b.chain_begin(1);
		core::Size b_end = b.chain_end(1);

		Real total_isc = calculateInterfaceScore(pose,scorefxn);

		//////////////////////////// pep from chain B ///////////////////////////////
		Real maxIsc = 1000;
		pose::Pose best_pose;
		Size best_lenb = 0;
		Size best_posb = 0;
		Real isc;
		pose::PoseOP curr_pose;

		for (Size pep_length=5; pep_length<=10; ++pep_length) {
				if (pep_length > b_end) break;
				maxIsc = 1000.0;
				for (Size i=b_start; i<=b_end-pep_length+1; ++i) {
						curr_pose = get_subPose(chainA,chainB,i,i+pep_length-1,'B');
						TR << "let's score chainA first time:" << std::endl;
						(*scorefxn)(*chainA);
						TR << "let's score chainA second time:" << std::endl;
						(*scorefxn)(*chainA);
						TR << "let's score chainB first time:" << std::endl;
						(*scorefxn)(*chainB);
						TR << "let's score chainB second time:" << std::endl;
						(*scorefxn)(*chainB);
						TR << "let's score first time:" << std::endl;
						(*scorefxn)(*curr_pose);
						TR << "let's score second time:" << std::endl;
						(*scorefxn)(*curr_pose);
						TR << "let's calc intrf score:" << std::endl;
						isc = calculateInterfaceScore(*curr_pose,scorefxn);
						//debug
						TR << "after Isc calc .. " << std::endl;
						if (isc < maxIsc ) {
								maxIsc = isc;
								best_pose = *new pose::Pose(*curr_pose);
								best_posb=i;
						}
				}
				TR << "Best B len: " << pep_length << " pos: " << best_posb <<" isc: "
					 << maxIsc << " %tot " << maxIsc/total_isc << std::endl;
		}

		TR << "Best pep from B at pos: "<<best_posb<<" length: "<<best_lenb<<" w/ isc: "<<maxIsc<< std::endl;
		best_pose.dump_pdb("./bestPepB.pdb");

		//exit(0);
		//////////////////////////// pep from chain A ///////////////////////////////

		Size best_lena = 0;
		Size best_posa = 0;

		for (Size pep_length=5; pep_length<=10; ++pep_length) {
				if (pep_length > a_end) break;
				maxIsc = 1000.0;
				for (Size i=a_start; i<=a_end-pep_length+1; ++i) {
						curr_pose = get_subPose(chainA,chainB,i,i+pep_length-1,'A');
						isc = calculateInterfaceScore(*curr_pose,scorefxn);
						if (isc < maxIsc /*normalize by length - pep_length + 5 */) {
								maxIsc = isc;
								best_pose = *new pose::Pose(*curr_pose);
								best_posa=i;
						}
				}
				TR << "Best A len: " << pep_length << " pos: " << best_posa <<" isc: "
					 << maxIsc << " %tot " << maxIsc/total_isc << std::endl;
		}

		TR << "Best pep from A at pos: "<<best_posa<<" length: "<<best_lena<<" w/ isc: "<<maxIsc<< std::endl;
		best_pose.dump_pdb("./bestPepA.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


// Given a pose and a score function - calculate the interface score by seperating the monomers
// and calculating the scores of the complex and unbound monomers.
Real calculateInterfaceScore ( core::pose::Pose & pose, core::scoring::ScoreFunctionCOP scorefxn ){
		//debug
		pose.dump_pdb("./beforeIsc.pdb");
		TR << "calculating Isc " << std::endl;

    //make a new pose in which the monomers are far apart from one another
    core::pose::Pose unbound_pose = pose;
    float trans_magnitude = 1000; // 1000A
    Size rb_jump = 1;     //rigid-body jump number is assumed to be 1 - this will work for two monomers.
		protocols::rigid::RigidBodyTransMoverOP translate_away (
		                new protocols::rigid::RigidBodyTransMover( unbound_pose, rb_jump ) );
		translate_away->step_size( trans_magnitude );

		//calculate scores
    float bound_energy, unbound_energy;
		bound_energy = (*scorefxn)( unbound_pose ); // before the move - as a complex.
		translate_away->apply( unbound_pose );
		unbound_energy = (*scorefxn)( unbound_pose ); // after the move - unbound monomers.
		TR << "after " << bound_energy << " " << unbound_energy << std::endl;
    //answer
    return (bound_energy - unbound_energy);
}


core::pose::PoseOP get_subPose(
				 core::pose::PoseCOP chainA,
  			 core::pose::PoseCOP chainB,
	  		 core::Size pep_start,
				 core::Size pep_end,
				 char pep_chain )
{
		using namespace core;
		using namespace core::chemical;

		pose::PoseOP pose = NULL;
		pose::PoseCOP other;

		if (pep_chain == 'A') {
				pose = new pose::Pose(*chainB);
				other = chainA;
		}
		else {//pep_chain == B
				pose = new pose::Pose(*chainA);
				other = chainB;
		}

		pose->append_residue_by_jump(other->residue(pep_start), 1, "", "", true);
		add_variant_type_to_pose_residue( *pose , LOWER_TERMINUS, pose->total_residue() );

		for (core::Size i=pep_start+1; i<=pep_end; ++i){

				pose->append_residue_by_bond(other->residue(i),false,0,0,0,false);
		}
		add_variant_type_to_pose_residue( *pose , UPPER_TERMINUS, pose->total_residue() );
		pose->conformation().detect_disulfides();

		return pose;
}





