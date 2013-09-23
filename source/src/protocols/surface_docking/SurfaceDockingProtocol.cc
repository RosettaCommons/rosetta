// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;
//     rm-trailing-spaces:t -*-
//     vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
//     under license.
// (c) The Rosetta software is developed by the contributing members of the
//     Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about
//     this can be
// (c) addressed to University of Washington UW TechTransfer,
//     email: license@u.washington.edu.

/// @file SurfaceDockingProtocol.cc
/// @author Robin A Thottungal (rathottungal@gmail.com)

// Unit Headers
#include <protocols/surface_docking/SurfaceDockingProtocol.hh>

// Package Headers
#include <protocols/surface_docking/CentroidRelaxMover.hh>
#include <protocols/surface_docking/FullatomRelaxMover.hh>
#include <protocols/surface_docking/SurfaceOrientMover.hh>
#include <protocols/surface_docking/SlideIntoSurface.hh>


// Project headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/surface_docking/SurfaceParameters.fwd.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <protocols/relax/ClassicRelax.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

// for adding data to pose
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/datacache/CacheableDataType.hh>


#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <core/scoring/solid_surface/SurfaceEnergies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>

//Abinitio
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.hh>

//option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>


//Utility Headers
#include <utility/exit.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/FoldTree.hh>
#include <string>


using basic::T;
using basic::Error;
using basic::Warning;

using namespace protocols::surface_docking;
using namespace protocols;
using namespace protocols::moves;
//using namespace ObjexxFCL::format;

//using core::pose::datacache::CacheableDataType::SURFACE_PARAMS;
static basic::Tracer TR("protocols.SurfaceDocking.SurfaceDockingProtocol");

namespace protocols {
namespace surface_docking {

using namespace core;
using protocols::jd2::JobDistributor;

//constructor
SurfaceDockingProtocol::SurfaceDockingProtocol() : Mover(){
	Mover::type( "SurfaceDockingProtocol");
	score_sidechain_pack_ = scoring::getScoreFunction();
	// setting weighs for score12
	TR << "Setting Weights for score12" << std::endl;
	score_sidechain_pack_->set_weight( core::scoring::fa_elec, 1.0 );
	//score_sidechain_pack_->set_weight(core::scoring::fa_atr,0.97);
	//score_sidechain_pack_->set_weight(core::scoring::fa_sol,0.92);
	//score_sidechain_pack_->set_weight(core::scoring::fa_rep,0.64);
	}

//destructor
SurfaceDockingProtocol::~SurfaceDockingProtocol() {}

void SurfaceDockingProtocol::setupFoldTree(pose::Pose & pose){
	// Make sure the HETATM is first
	TR<<"Total Residue:"<<pose.total_residue()<<std::endl;
	kinematics::FoldTree ft=pose.fold_tree();
	//TR<<"Old FoldTree"<<ft<<std::endl;
	TR<<"FoldTree"<<ft<<std::endl;
	TR<<"Number of Jumps:"<<pose.num_jump()<<std::endl;
}



void SurfaceDockingProtocol::apply(pose::Pose & pose) {

	Size const first_protein_residue( pose.num_jump() + 1 );

	// Initialize the SurfaceEnergies object; this protocol assumes that the block of residues from num_jump to total_residue
	// represents the single range of non-surface residues
	core::scoring::solid_surface::SurfaceEnergiesOP surfEs = new core::scoring::solid_surface::SurfaceEnergies;
	surfEs->set_total_residue( pose.total_residue() );
	surfEs->set_residue_range_not_surface( first_protein_residue, pose.total_residue() );
	pose.set_new_energies_object( surfEs );

	// Attaching SurfaceParameters to pose
	utility::vector0< std::string > SurfVectors(3);
	TR<<"Remark Info attached to pose SurfaceDockingProtocol:"<<std::endl;
	TR<<"Size of Remarks:"<<pose.pdb_info()->remarks().size()<<std::endl;
	for(Size i=0; i<pose.pdb_info()->remarks().size(); i++) {
		TR << pose.pdb_info()->remarks().at(i).num << " " << pose.pdb_info()->remarks().at(i).value << std::endl;
		if ( i > 2 ) { TR << "TOO MANY REMARKS!" << std::endl; continue; }
		SurfVectors[i]=pose.pdb_info()->remarks().at(i).value;
	}


	surface_docking::SurfaceParametersOP surfaceParameters =
		new surface_docking::SurfaceParameters(SurfVectors[0],SurfVectors[1],SurfVectors[2]);
	pose.data().set(core::pose::datacache::CacheableDataType::SURFACE_PARAMS,surfaceParameters);

	// grabing the information attached to the pose to save temporarily
	surface_docking::SurfaceParameters & surfaceVectors_tmp=
		*( static_cast< surface_docking::SurfaceParameters * >
		( pose.data().get_ptr( core::pose::datacache::CacheableDataType::SURFACE_PARAMS )() ));

	// Creating a job compatible with JD2
	protocols::jd2::JobOP job = jd2::JobDistributor::get_instance()->current_job();
	std::string job_name (JobDistributor::get_instance()->
		job_outputter()->output_name( job ) );

	//Testing
	TR<<"Total Residue"<<pose.total_residue()<<std::endl;

	// Setting up the foldTree
	setupFoldTree(pose);
	// Copying the pdb_info((
	core::pose::PDBInfoOP pose_pdbinfo = pose.pdb_info();
	// getting info about the surface vectors from remarks

	simple_moves::SwitchResidueTypeSetMover to_centroid("centroid");
	simple_moves::SwitchResidueTypeSetMover to_all_atom( "fa_standard" );

	TR<<"Applying surfaceOrient Mover"<<std::endl;
	surface_docking::SurfaceOrientMoverOP sf=
		new surface_docking::SurfaceOrientMover();
	sf->apply(pose);

	//split the pose into separate chains
	utility::vector1< pose::PoseOP > singlechain_poses;
	singlechain_poses = pose.split_by_chain();
	assert( singlechain_poses.size() == 2);
	singlechain_poses[ 2 ]->set_new_energies_object( new core::scoring::Energies ); // replace the SurfaceEnergies object

	// singlechain_poses[2] contains the protein; can be made generic
	// for multiple proteins too
	TR<<"CentroidRelaxed Started................"<<std::endl;
	to_centroid.apply( *singlechain_poses[2] );
	// Starting abinitio protocol
	abinitio(*singlechain_poses[2]);
	CalcSecondayStruct(*singlechain_poses[2]);
	TR<<"SS After Abinitio()"<<SecStruct_<<std::endl;
	surface_docking::CentroidRelaxMoverOP centroidrelax =
		new surface_docking::CentroidRelaxMover();
	centroidrelax->set_nmoves(10);
	centroidrelax->apply(*singlechain_poses[2]);
	CalcSecondayStruct(*singlechain_poses[2]);
	TR<<"SS After CentroidRelax:"<<SecStruct_<<std::endl;
	TR<<"CentroidRelaxed Done................"<<std::endl;
	to_all_atom.apply(*singlechain_poses[2]);
	CalcSecondayStruct(*singlechain_poses[2]);
	TR<<"SS After to_all_atom:"<<SecStruct_<<std::endl;
	pose.copy_segment((*singlechain_poses[2]).total_residue(),
		(*singlechain_poses[2]),
		(*singlechain_poses[1]).total_residue()+1,1);
	setupFoldTree(pose);
	// re-attach the information to the pose
	pose.data().set(core::pose::datacache::CacheableDataType::SURFACE_PARAMS,&surfaceVectors_tmp);
	// re-attach pdb_info to pose
	pose.pdb_info(pose_pdbinfo);
	// setting the secondary structure info
	SetSecondayStruct(pose);
	TR<<"Secondary Structure:"<<pose.secstruct()<<std::endl;
	// added on Jun 22, to see if to_full atom is introducing high ene structures
	pack::task::PackerTaskOP my_task_fullatom
		( pack::task::TaskFactory::create_packer_task( pose ));

	// allow repacking for only the protein side chains!
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if (ii < first_protein_residue) {
			my_task_fullatom->nonconst_residue_task( ii ).prevent_repacking();
		} else {
			my_task_fullatom->nonconst_residue_task( ii ).restrict_to_repacking();
		}
	}
	protocols::simple_moves::PackRotamersMoverOP side_chain_pack_fullatom
		( new protocols::simple_moves::PackRotamersMover( score_sidechain_pack_,  my_task_fullatom ) );
	(*score_sidechain_pack_)(pose);
	TR << " Score Before Sidechain Packing:" << (*score_sidechain_pack_)(pose) << std::endl;
	side_chain_pack_fullatom->apply(pose);
	(*score_sidechain_pack_)(pose);
	TR << " Score After Sidechain Packing:" << (*score_sidechain_pack_)(pose) << std::endl;
	TR<<"Packing Sidechain Done................"<<std::endl;
	//Fullatom Relax starts!
	TR<<"FullatomeRelaxed Started................"<<std::endl;
	sf->apply(pose);
	// Moving the protein away from the surface
	FaSlideAwayFromSurface MoveProteinAway( pose.num_jump() );
	MoveProteinAway.apply(pose);

	surface_docking::FullatomRelaxMoverOP allatomrelax = new surface_docking::FullatomRelaxMover();
	if ( basic::options::option[ basic::options::OptionKeys::run::benchmark ] ) {
		allatomrelax->set_nmoves(3);
	} else {
		allatomrelax->set_nmoves(6);
	}

	core::Size lj_ramp_cycle=5; // inside the constructor
	// set the random number
	core::Real lj_increment = ( 1.0 - 0.02 )/ lj_ramp_cycle;
	for (Size i=1;i<=lj_ramp_cycle;++i) {
		TR<<"FullatomRelax Execution Number:"<<i<<std::endl;
		if ( basic::options::option[ basic::options::OptionKeys::run::benchmark ] ) {
			allatomrelax->set_smallmovesize(3);
		} else {
			allatomrelax->set_smallmovesize(30/i);
		}

		allatomrelax->set_ljrepulsion_weight((0.02+i*lj_increment));
		allatomrelax->set_ecounter(i); // when value of i matches the random number
        // in the allatomrelax, protein is slide into the surface
		allatomrelax->apply(pose);
	}
	// Final Side Chain re-packing using rtmin; adopted from Dave's protocol
	TR<<"FullatomeRelaxed Done................"<<std::endl;
	sf->apply(pose);
	TR<<"Packing Sidechains Started................"<<std::endl;
	pack::task::PackerTaskOP my_task( pack::task::TaskFactory::create_packer_task( pose ));

	// allow repacking for only the protein side chains!
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if (ii < first_protein_residue ){
			my_task->nonconst_residue_task( ii ).prevent_repacking();
		} else {
			my_task->nonconst_residue_task( ii ).restrict_to_repacking();
		}

	}
	protocols::simple_moves::PackRotamersMoverOP side_chain_pack
		( new protocols::simple_moves::PackRotamersMover( score_sidechain_pack_, my_task ) );
	(*score_sidechain_pack_)(pose);
	TR << " Score Before Sidechain Packing:" << (*score_sidechain_pack_)(pose) << std::endl;
	side_chain_pack->apply(pose);
	(*score_sidechain_pack_)(pose);
	TR << " Score After Sidechain Packing:" << (*score_sidechain_pack_)(pose) << std::endl;
	TR<<"Packing Sidechain Done................"<<std::endl;

	// Side Chain Packing done
	CalcSecondayStruct_withSurface(pose);
	TR<<"Secondary Structure in Solution:"<<allatomrelax->getSecondayStruct()<<std::endl;
	TR<<"Secondary Structure on Surface:"<<SecStruct_<<std::endl;
	job->add_string_string_pair("SolState_SecondaryStructure:",allatomrelax->getSecondayStruct());
	job->add_string_string_pair("Adsorbed_SecondaryStructure:",SecStruct_);
	JobDistributor::get_instance()->job_outputter()->other_pose( job,pose, "Surface_");
}

void SurfaceDockingProtocol::CalcSecondayStruct(core::pose::Pose & pose){
	SecStruct_="";
	core::scoring::dssp::Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		SecStruct_+=pose.secstruct(ii);
	}
}

void SurfaceDockingProtocol::CalcSecondayStruct_withSurface(core::pose::Pose & pose){
	SecStruct_="";
	// creating a new pose and splitting!
	//split the pose into separate chains
	pose::PoseOP pose_tmp=new pose::Pose(pose);
	utility::vector1< pose::PoseOP > singlechain_poses;
	singlechain_poses = pose_tmp->split_by_chain();
	core::scoring::dssp::Dssp dssp( *singlechain_poses[2] );
	dssp.insert_ss_into_pose( *singlechain_poses[2] );
	for ( Size ii = 1; ii <= singlechain_poses[2]->total_residue(); ++ii ) {
		SecStruct_+=singlechain_poses[2]->secstruct(ii);
	}
}

void SurfaceDockingProtocol::SetSecondayStruct(core::pose::Pose & pose){
	Size index=0;
	for ( Size ii = pose.num_jump()+1; ii <= pose.total_residue(); ++ii ) {
		pose.set_secstruct(ii,SecStruct_[index]);
		index=index+1;
	}
}

void SurfaceDockingProtocol::abinitio (core::pose::Pose & pose){
	using namespace core;
	using namespace basic::options;
	//ClassicAbinitio
	TR << "starting ClassicAbinitio...." << std::endl;
	protocols::abinitio::ClassicAbinitio::register_options();
	std::string frag_large_file = option[ OptionKeys::in::file::frag9 ]();
	std::string frag_small_file = option[ OptionKeys::in::file::frag3 ]();
	TR<<"Loading Fragment 9:"<<std::endl;
	core::fragment::ConstantLengthFragSetOP fragset_large = new core::fragment::ConstantLengthFragSet;
	fragset_large->read_fragment_file(frag_large_file);
	TR<<"Loading Fragment 3:"<<std::endl;
	core::fragment::ConstantLengthFragSetOP fragset_small = new core::fragment::ConstantLengthFragSet;
	fragset_small->read_fragment_file(frag_small_file);
	core::kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	movemap->set_bb_true_range(1,pose.total_residue());
	protocols::abinitio::ClassicAbinitio abinitio(fragset_small,fragset_large,movemap);
	abinitio.set_score_weight(scoring::rg,0);
	abinitio.init(pose);
	abinitio.apply(pose);

}

std::string SurfaceDockingProtocol::get_name() const {
	return "SurfaceDockingProtocol";
}

}	//surfaceDockingProtocol

}	//protocol


