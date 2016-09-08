// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
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
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/options/util.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/util/disulfide_util.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <protocols/scoring/Interface.hh> //core/conformation
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/file_data.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers

//Auto using namespaces
namespace ObjexxFCL { namespace fmt { } } using namespace ObjexxFCL::fmt; // AUTO USING NS
//Auto using namespaces end

using basic::T;
using basic::Error;
using basic::Warning;


namespace peptide_deriver {
	basic::options::StringOptionKey chain_to_derive_from("peptide_deriver:chain_to_derive_from");
	basic::options::RealOptionKey   length_to_derive("peptide_deriver:length_to_derive");
	basic::options::BooleanOptionKey optimize("peptide_deriver:optimize");
	basic::options::StringOptionKey dump_peptide_pdb("peptide_deriver:dump_peptide_pdb");
}

static thread_local basic::Tracer TR( "PeptideDeriver" );

using namespace core;
using namespace basic::options;
using namespace OptionKeys;


//utility functions
core::pose::PoseOP get_subPose(core::pose::PoseOP, core::pose::PoseOP,
														 core::Size,core::Size,char);
Real calculateInterfaceScore (core::pose::Pose, core::scoring::ScoreFunctionOP );
core::pose::PoseOP twoChainPose (core::pose::Pose const & src, Size chainA, Size chainB);
core::pose::PoseOP oneChainPose (core::pose::Pose const & src, Size chainA);
void derive(core::pose::Pose const & pose);

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
	using namespace core::scoring::func;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add( peptide_deriver::chain_to_derive_from,
		"Chain from which to derive the peptide, options are FIRST/SECOND/BOTH").def("BOTH");
	option.add( peptide_deriver::length_to_derive,"Lengths of peptides to check" ).def(10);
	option.add( peptide_deriver::optimize, "Makes deriver go faster by skipping peptides with 0 I_sc" ).def(true);
	option.add( peptide_deriver::dump_peptide_pdb, "Output a pair of PDBs for best peptides in each pair of chains, using this name as a prefix");


	//setup random numbers and options
	devel::init(argc, argv);

	//create a pose
	pose::Pose orig_pose;
	core::import_pose::pose_from_pdb(orig_pose, basic::options::start_file());

	//setup scorefxn
	//TODO: consider using talaris instead - Oriel
	core::scoring::ScoreFunctionOP scorefxn(
		ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH ));
	
		//check number of chains
		utility::vector1< PoseOP > chains = orig_pose.split_by_chain();
		if (chains.size() < 2) {
			TR << "ERROR: pdb contains only 1 chain - can not run Peptideriver" << std::endl;
			exit(1);
		}
		//if (chains.size() > 2) {
		//	TR << "ERROR: pdb contains more than 2 chains - can not run Peptideriver" << std::endl;
		//	exit(1);
		//}
		

	//////////////// set native constraints prior to minimization
	ConstraintSetOP cst_set( new ConstraintSet() );
	HarmonicFuncOP spring( new HarmonicFunc( 0 /*mean*/, 1 /*std-dev*/) );
	conformation::Conformation const & conformation( orig_pose.conformation() );
	for (Size i=1; i <= orig_pose.size(); ++i) {
//		Residue const  & reside = orig_pose.residue( i ); // commented out to fix compilation error: unused variable 'reside' [-Werror=unused-variable] - Oriel
		id::AtomID CAi ( orig_pose.residue(i).atom_index( " CA " ), i );
		cst_set->add_constraint
		(  ConstraintCOP( ConstraintOP( new CoordinateConstraint
		 ( CAi, CAi, conformation.xyz( CAi ), spring ) ) )
		 );
	}
	orig_pose.constraint_set( cst_set );
	scorefxn->set_weight(coordinate_constraint, 1.0);
	/////////////////////////////////////////////////////////////////////

	////////////// minimize bb and sc
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_bb(true);
	movemap->set_chi(true);
	protocols::simple_moves::MinMover minimizer(movemap, scorefxn, "lbfgs_armijo_atol", 0.01, true /*nb_list*/);
	/////////////////////////////////////////////////////////////////////

	//minimize pose
	TR << "minimizing" << std::endl;
	minimizer.apply(orig_pose);
	TR << "minimized" << std::endl;
	//detect SS
	orig_pose.conformation().detect_disulfides();
	//remove constraint weight from scorefxn
	scorefxn->set_weight(coordinate_constraint, 0.0);

	for (Size i=1; i < chains.size(); ++i) {
		for (Size j=i+1; j <=chains.size(); ++j) { 
			pose::Pose pose = *(twoChainPose (orig_pose, i, j));
			derive(pose);
		}
	}

//	pose::Pose pose = *(oneChainPose (orig_pose, 1));
//	derive(pose);
//	pose =  *(oneChainPose (orig_pose, 2));
//	derive(pose);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

//the main derivation function
void derive(core::pose::Pose const & pose) {
	using namespace core::pose;
	using namespace scoring;
	using namespace conformation;
	using namespace core::chemical;
	using namespace protocols;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

		//setup scorefxn
		core::scoring::ScoreFunctionOP scorefxn(
			ScoreFunctionFactory::create_score_function( PRE_TALARIS_2013_STANDARD_WTS, SCORE12_PATCH ));

		utility::vector1< PoseOP > chains = pose.split_by_chain();
		TR << "Number of chains " << chains.size() << std::endl;
		TR << "ft " << pose.fold_tree() << std::endl;
		pose::PoseOP chainA = chains[1];
		pose::PoseOP chainB = chains[2];

		conformation::Conformation & a = chainA->conformation();
		conformation::Conformation & b = chainB->conformation();
		core::Size a_start = a.chain_begin(1);
		core::Size a_end = a.chain_end(1);
		core::Size b_start = b.chain_begin(1);
		core::Size b_end = b.chain_end(1);

		Real total_isc = calculateInterfaceScore(pose,scorefxn);
		TR << "Total_isc: " << total_isc << std::endl;
		//set length of peptide to traverse
		//if user provides 0 - scans from length 5 to length 10
		Size from, to;
		if (option[ peptide_deriver::length_to_derive ]() == 0) {
				from = 1;
				to = 10;
		} else {
				from = (Size)(option[ peptide_deriver::length_to_derive ]());
				to = from;
		}

		Real maxIsc = 1000;
		pose::PoseOP best_pose;
//		Size best_lenb = 0;  // commented out to fix compilation error: unused variable 'best_lenb' [-Werror=unused-variable] - Oriel
		Size best_posb = 0;
		Real isc;
		pose::PoseOP curr_pose;

		//////////////////////////// pep from chain B ///////////////////////////////
		if ((option[ peptide_deriver::chain_to_derive_from ]() == "BOTH") ||
			(option[ peptide_deriver::chain_to_derive_from ]() == "SECOND")) {

		for (Size pep_length=from; pep_length<=to; ++pep_length) {
				if (pep_length > b_end) break;
				maxIsc = 1000.0;
				for (Size i=b_start; i<=b_end-pep_length+1; ++i) {
						curr_pose = get_subPose(chainA,chainB,i,i+pep_length-1,'B');
						isc = calculateInterfaceScore(*curr_pose,scorefxn);
						//optimization
						if (option[ peptide_deriver::optimize ]()){
								if (isc==0) { i=i+pep_length-1; }
						}
						//end optimiztaion
						TR << "ALL: Recp: "<< chainA->pdb_info()->chain(a_start) << " Pep: " <<  chainB->pdb_info()->chain(b_start) << " , pos: " << i << " len: " << pep_length << " sc: " << isc << " %tot: " << isc/total_isc << " ";
						//detect putative SS bridge positions
						for (Size sA=curr_pose->conformation().chain_begin(2); sA<curr_pose->conformation().chain_end(2); sA++) {
							for (Size sB=sA+1; sB<=curr_pose->conformation().chain_end(2); sB++) {

//								if ((curr_pose->aa(sA) != chemical::aa_gly) && (curr_pose->aa(sB) != chemical::aa_gly))
								if ((curr_pose->residue(sA).name1() != 'G') && (curr_pose->residue(sB).name1() != 'G')) {
									float dist = curr_pose->residue(sA).xyz("CB").distance(curr_pose->residue(sB).xyz("CB"));
									if ((dist<5) && (dist>3)) { 
									TR << sA << "-" << sB << "(" << i+sA-a_end-1 << "-" << i+sB-a_end-1 << "):" << dist << " "; 
									}
								}
							}
						}
						TR <<  std::endl;
						if (isc < maxIsc ) {
								maxIsc = isc;
								best_pose = pose::PoseOP( new pose::Pose(*curr_pose) );
								best_posb=i;
						}
				}
				TR << "Best Recp: "<< chainA->pdb_info()->chain(a_start) << " Pep: " <<  chainB->pdb_info()->chain(b_start) << " len: " << pep_length << " pos: " << best_posb <<" isc: "
					 << maxIsc << " %tot " << maxIsc/total_isc << std::endl;
		}

		//TR << "Best pep Recp: "<< chainA->pdb_info()->chain(a_start) << " Pep: " <<  chainB->pdb_info()->chain(b_start) << " " << best_posb <<" length: "<<best_lenb<<" w/ isc: "<<maxIsc<< std::endl;
		if( option[ peptide_deriver::dump_peptide_pdb ].user() ) {
			std::string const dump_peptide_pdb_prefix = option[ peptide_deriver::dump_peptide_pdb ]();
			char recp_chain = chainA->pdb_info()->chain(a_start);
			char pep_chain = chainB->pdb_info()->chain(b_start);
			best_pose->dump_scored_pdb(dump_peptide_pdb_prefix + ".recp" + recp_chain + ".pep" + pep_chain + ".pdb", *scorefxn);
		}

		}// end if derive second

		//////////////////////////// pep from chain A ///////////////////////////////
		if ((option[ peptide_deriver::chain_to_derive_from ]() == "BOTH") ||
				(option[ peptide_deriver::chain_to_derive_from ]() == "FIRST")) {
//		Size best_lena = 0; // commented out to fix compilation error: unused variable 'best_lena' [-Werror=unused-variable] - Oriel
		Size best_posa = 0;

		for (Size pep_length=from; pep_length<=to; ++pep_length) {
				if (pep_length > a_end) break;
				maxIsc = 1000.0;
				for (Size i=a_start; i<=a_end-pep_length+1; ++i) {
						curr_pose = get_subPose(chainA,chainB,i,i+pep_length-1,'A');
						isc = calculateInterfaceScore(*curr_pose,scorefxn);
						//optimization
						if (option[ peptide_deriver::optimize ]()){
								if (isc==0) { i=i+pep_length-1; }
						}
						//end optimiztaion
					TR << "ALL: Recp: "<< chainB->pdb_info()->chain(b_start) << " Pep: " <<  chainA->pdb_info()->chain(a_start) << " , pos: " << i << " len: " << pep_length << " sc: " << isc << " %tot: " << isc/total_isc << " ";
					TR <<  std::endl;	
					if (isc < maxIsc /*normalize by length - pep_length + 5 */) {
								maxIsc = isc;
								best_pose = pose::PoseOP( new pose::Pose(*curr_pose) );
								best_posa=i;
						}
				}
				TR << "Best Recp: "<< chainB->pdb_info()->chain(b_start) << " Pep: " <<  chainA->pdb_info()->chain(a_start) << " len: " << pep_length << " pos: " << best_posa <<" isc: "
				<< maxIsc << " %tot " << maxIsc/total_isc << std::endl;
		}

		if( option[ peptide_deriver::dump_peptide_pdb ].user() ) {
                        std::string const dump_peptide_pdb_prefix = option[ peptide_deriver::dump_peptide_pdb ]();
                        char recp_chain = chainB->pdb_info()->chain(b_start);
                        char pep_chain = chainA->pdb_info()->chain(a_start);
                        best_pose->dump_scored_pdb(dump_peptide_pdb_prefix + ".recp" + recp_chain + ".pep" + pep_chain + ".pdb", *scorefxn);
                }

		//TR << "Best pep from A at pos: "<<best_posa<<" length: "<<best_lena<<" w/ isc: "<<maxIsc<< std::endl;
		//best_pose->dump_pdb("./bestPepA.pdb");
		} // end if derive first
}


// Given a pose and a score function - calculate the interface score by seperating the monomers
// and calculating the scores of the complex and unbound monomers.
Real calculateInterfaceScore ( core::pose::Pose pose, core::scoring::ScoreFunctionOP scorefxn ){

    //make a new pose in which the monomers are far apart from one another
    core::pose::PoseOP unbound_pose( new core::pose::Pose(pose) );
    //translation size = 1000A.
    float trans_magnitude = 1000;
    //rb jump number is assumed to be 1 - this will work for two monomers.
    Size rb_jump = 1;
    //create a Rigid body mover - the mover works on a rigid body jump
    protocols::rigid::RigidBodyTransMoverOP translate_away( new protocols::rigid::RigidBodyTransMover( *unbound_pose, rb_jump ) );

    //set the size of the move
    translate_away->step_size( trans_magnitude );

		//calculate scores
    float bound_energy = (*scorefxn)( *unbound_pose ); // before the move - as a complex.
    translate_away->apply( *unbound_pose );
		float unbound_energy = (*scorefxn)( *unbound_pose ); // after the move - unbound monomers.

    //answer
    return (bound_energy - unbound_energy);
}


//create a sub-pose of the receptor peptide only.
core::pose::PoseOP get_subPose(
			 core::pose::PoseOP chainA,
  			 core::pose::PoseOP chainB,
	  		 core::Size pep_start,
			 core::Size pep_end,
				 char pep_chain )
{
		using namespace core;
		using namespace core::chemical;

		pose::PoseOP pose;
		pose::PoseOP other;

		if (pep_chain == 'A') {
				pose = pose::PoseOP( new pose::Pose(*chainB) );
				other = chainA;
		}
		else {//pep_chain == B
				pose = pose::PoseOP( new pose::Pose(*chainA) );
				other = chainB;
		}

		pose->append_residue_by_jump(other->residue(pep_start), 1, "", "", true);
		core::pose::add_variant_type_to_pose_residue( *pose , LOWER_TERMINUS_VARIANT, pose->size() );

		for (core::Size i=pep_start+1; i<=pep_end; ++i){

				pose->append_residue_by_bond(other->residue(i),false,0,0,0,false);
		}
		core::pose::add_variant_type_to_pose_residue( *pose , UPPER_TERMINUS_VARIANT, pose->size() );

		pose->conformation().detect_disulfides();
		return pose;
}

//create a pose from two given chains of the original pose
core::pose::PoseOP twoChainPose (core::pose::Pose const & src, Size chainA, Size chainB) {

	Size startA = src.conformation().chain_begin(chainA);
	Size startB = src.conformation().chain_begin(chainB);
	Size endA = src.conformation().chain_end(chainA);
	Size endB = src.conformation().chain_end(chainB);	
	utility::vector1< core::Size > residue_indices;
	
	for(Size i=startA ; i <= endA; ++i){
		residue_indices.push_back(i);
	}
	for(Size i=startB ; i <= endB; ++i){
		residue_indices.push_back(i);
	}
	
	pose::PoseOP pose( new pose::Pose() );
	core::io::pdb::pose_from_pose(*pose, src, residue_indices);
	return pose;
}

//create a pose from two given chains of the original pose
core::pose::PoseOP oneChainPose (core::pose::Pose const & src, Size chainA) {
	
	Size startA = src.conformation().chain_begin(chainA);
	Size endA = src.conformation().chain_end(chainA);
	utility::vector1< core::Size > residue_indices;
	
	for(Size i=startA ; i <= endA; ++i){
		residue_indices.push_back(i);
	}
	for(Size i=1 ; i <= src.size(); ++i){
		if ((i<startA) || (i>endA))
			residue_indices.push_back(i);
	}
	
	pose::PoseOP pose( new pose::Pose() );
	core::io::pdb::pose_from_pose(*pose, src, residue_indices);
	return pose;
}
