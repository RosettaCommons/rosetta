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
/// @author Michael Pacella (mpacella88@gmail.com)

// Unit Headers
#include <protocols/surface_docking/SurfaceDockingProtocol.hh>

// Package Headers
#include <protocols/surface_docking/CentroidRelaxMover.hh>
#include <protocols/surface_docking/FullatomRelaxMover.hh>
#include <protocols/surface_docking/SurfaceOrientMover.hh>
#include <protocols/surface_docking/SurfaceParameters.hh>
#include <core/scoring/solid_surface/SurfaceEnergies.hh>

// Project headers
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
//Abinitio
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.hh>

//Basic Headers
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

//Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

//C++ Headers
#include <string>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace protocols::surface_docking;
using namespace protocols;
using namespace protocols::moves;

static basic::Tracer TR("protocols.SurfaceDocking.SurfaceDockingProtocol");

namespace protocols {
namespace surface_docking {

using namespace core;

//Default constructor
SurfaceDockingProtocol::SurfaceDockingProtocol() : Mover()
	{
		init();
	}

SurfaceDockingProtocol::SurfaceDockingProtocol(SurfaceDockingProtocol const & src) : Mover(src)
	{
		copy_data(*this, src);
	}
	
/// @brief
//Destructor
SurfaceDockingProtocol::~SurfaceDockingProtocol() {}
	
protocols::moves::MoverOP SurfaceDockingProtocol::clone() const
	{
		return new SurfaceDockingProtocol(*this);
	}
	
protocols::moves::MoverOP SurfaceDockingProtocol::fresh_instance() const
	{
		return new SurfaceDockingProtocol();
	}
	
void SurfaceDockingProtocol::show(std::ostream & output) const
	{
		using namespace std;
		Mover::show(output);  // name, type, tag
	}
	
void SurfaceDockingProtocol::apply(pose::Pose & pose)
	{
	TR<<"Starting Surface Docking Protocol..."<<std::endl;
	Size const first_protein_residue( pose.num_jump() + 1 );
		
	initialize_surface_energies( pose, first_protein_residue );
	set_surface_parameters( pose );
	setup_movers ( pose, first_protein_residue );

	// Grabbing the job from JD2 for appending job specific information to the pdb
	protocols::jd2::JobOP job = jd2::JobDistributor::get_instance()->current_job();
	
	TR<<"Reorienting protein above the surface..."<<std::endl;
	surface_orient_->apply(pose);
	
	TR<<"Abinitio started..."<<std::endl;
	core::pose::Pose surface, protein;
	split_protein_surface_poses(pose, surface, protein);
	to_centroid_->apply( protein );
	abinitio_->init(protein);
	abinitio_->apply(protein);

	TR<<"Centroid mode refinement started..."<<std::endl;
	centroid_relax_->apply(protein);
	
	TR<<"Adding sidechains and packing..."<<std::endl;
	to_full_atom_->apply(protein);
	merge_protein_surface_poses(pose, surface, protein);
	pack_rotamers_fullatom_->apply( pose );
	
	TR<<"Repositioning protein above the surface..."<<std::endl;
	position_above_surface_->apply(pose);
	 pose.dump_pdb("/Users/mpacella/Rosetta_Surface_Test/after_position_above.pdb");
	
	TR<<"Fullatom refinement started..."<<std::endl;
	fullatom_relax_->apply(pose);
		
	TR<<"Outputting decoy..."<<std::endl;
	surface_orient_->apply(pose);
	pack_rotamers_fullatom_->apply(pose);
	job->add_string_real_pair("Total weighted score: ", pose.energies().total_energy());
	job->add_string_string_pair("SolState_SecondaryStructure:",fullatom_relax_->get_sol_secondary_struct());
	job->add_string_string_pair("Adsorbed_SecondaryStructure:",fullatom_relax_->get_ads_secondary_struct());
	jd2::JobDistributor::get_instance()->job_outputter()->other_pose( job,pose, "Surface_");
}

std::string SurfaceDockingProtocol::get_name() const
	{
		return "SurfaceDockingProtocol";
	}

//Private methods////////////////////////////////////////
	
void SurfaceDockingProtocol::copy_data(SurfaceDockingProtocol object_to_copy_to, SurfaceDockingProtocol object_to_copy_from)
	{
		object_to_copy_to.score_sidechain_pack_ = object_to_copy_from.score_sidechain_pack_;
		object_to_copy_to.sec_struct_ = object_to_copy_from.sec_struct_;
		object_to_copy_to.surface_parameters_ = object_to_copy_from.surface_parameters_;
		object_to_copy_to.to_centroid_ = object_to_copy_from.to_centroid_;
		object_to_copy_to.to_full_atom_ = object_to_copy_from.to_full_atom_;
		object_to_copy_to.surface_orient_ = object_to_copy_from.surface_orient_;
		object_to_copy_to.abinitio_ = object_to_copy_from.abinitio_;
		object_to_copy_to.centroid_relax_ = object_to_copy_from.centroid_relax_;
		object_to_copy_to.pack_rotamers_fullatom_ = object_to_copy_from.pack_rotamers_fullatom_;
		object_to_copy_to.slide_away_from_surface_ = object_to_copy_from.slide_away_from_surface_;
		object_to_copy_to.slide_into_surface_ = object_to_copy_from.slide_into_surface_;
		object_to_copy_to.fullatom_relax_ = object_to_copy_from.fullatom_relax_;
	}
	
void SurfaceDockingProtocol::init()
	{
		Mover::type( "SurfaceDockingProtocol");
		score_sidechain_pack_ = scoring::getScoreFunction();
		TR << "Setting Weights for talaris" << std::endl;
		score_sidechain_pack_->set_weight( core::scoring::fa_elec, 1.0 );
		surface_parameters_ = NULL;
		sec_struct_ = "";
		to_centroid_ = NULL;
		to_full_atom_ = NULL;
		surface_orient_ = NULL;
		abinitio_ = NULL;
		centroid_relax_ = NULL;
		pack_rotamers_fullatom_ = NULL;
		slide_away_from_surface_ = NULL;
		fullatom_relax_ = NULL;		
	}
	
bool
SurfaceDockingProtocol::valid_surface_pose(core::pose::Pose const & pose)
	{
		if (pose.conformation().chain_endings().size() == 1)
		{
			return true;
		}
		return false;
	}
	
void SurfaceDockingProtocol::calc_secondary_structure(core::pose::Pose & pose)
	{
		sec_struct_="";
		core::scoring::dssp::Dssp dssp( pose );
		dssp.insert_ss_into_pose( pose );
		
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii )
		{
			sec_struct_+=pose.secstruct(ii);
		}
	}
	
void SurfaceDockingProtocol::calc_secondary_structure_with_surface(core::pose::Pose const & pose)
	{
		sec_struct_="";
		// creating a new pose and splitting!
		//split the pose into separate chains
		pose::Pose pose_tmp= pose::Pose(pose);
		utility::vector1< pose::PoseOP > singlechain_poses;
		singlechain_poses = pose_tmp.split_by_chain();
		core::scoring::dssp::Dssp dssp( *singlechain_poses[2] );
		dssp.insert_ss_into_pose( *singlechain_poses[2] );
		for ( Size ii = 1; ii <= singlechain_poses[2]->total_residue(); ++ii )
		{
			sec_struct_+=singlechain_poses[2]->secstruct(ii);
		}
	}
	
void SurfaceDockingProtocol::set_secondary_structure(core::pose::Pose & pose)
	{
		Size index=0;
		for ( Size ii = pose.num_jump()+1; ii <= pose.total_residue(); ++ii )
		{
			pose.set_secstruct(ii,sec_struct_[index]);
			index=index+1;
		}
	}
	
void SurfaceDockingProtocol::initialize_surface_energies(core::pose::Pose & pose, Size first_protein_residue)
	{
		// Initialize the SurfaceEnergies object; this protocol assumes that the block of residues from num_jump to total_residue
		// represents the single range of non-surface residues
		core::scoring::solid_surface::SurfaceEnergiesOP surfEs = new core::scoring::solid_surface::SurfaceEnergies;
		surfEs->set_total_residue( pose.total_residue() );
		surfEs->set_residue_range_not_surface( first_protein_residue, pose.total_residue() );
		pose.set_new_energies_object( surfEs );
	}
	
void SurfaceDockingProtocol::set_surface_parameters ( core::pose::Pose &)
{
	if (! basic::resource_manager::ResourceManager::get_instance()->has_resource_with_description("surface_vectors"))
	{
		throw utility::excn::EXCN_Msg_Exception(" Either a resource definition file or the command line option "\
				"-in:file:surface_vectors must be specified for surface docking");
	}
	
	surface_parameters_ = basic::resource_manager::get_resource<SurfaceParameters>("surface_vectors");
}

void SurfaceDockingProtocol::setup_movers ( core::pose::Pose const & pose, Size const first_protein_residue )
	{
		
		to_centroid_ = new simple_moves::SwitchResidueTypeSetMover("centroid");
		to_full_atom_ = new simple_moves::SwitchResidueTypeSetMover( "fa_standard" );
		surface_orient_ = new surface_docking::SurfaceOrientMover();
		surface_orient_->set_surface_parameters(surface_parameters_);
		setup_abinitio();
		centroid_relax_ = new surface_docking::CentroidRelaxMover();
		core::Size protein_length = pose.total_residue()-first_protein_residue+1;
		
		if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
			centroid_relax_->set_nmoves(2);
			centroid_relax_->set_outer_loop_cycles(1);
			centroid_relax_->set_inner_loop_cycles(1);
		}
		else {
			centroid_relax_->set_nmoves(10);
			centroid_relax_->set_outer_loop_cycles(6);
			centroid_relax_->set_inner_loop_cycles(protein_length);
		}
		setup_slide_movers(pose);
		pack::task::PackerTaskOP my_task_fullatom = create_surface_packer_task(pose, first_protein_residue);
		pack_rotamers_fullatom_ = new protocols::simple_moves::PackRotamersMover( score_sidechain_pack_,  my_task_fullatom );
		
		fullatom_relax_ = new protocols::surface_docking::FullatomRelaxMover();
		fullatom_relax_->set_surface_contact_mover(slide_into_surface_);
		fullatom_relax_->set_surface_orient_mover(surface_orient_);
		fullatom_relax_->set_surface_parameters(surface_parameters_);
		
		if ( basic::options::option[ basic::options::OptionKeys::run::test_cycles ] ) {
			fullatom_relax_->set_nmoves(3);
		} else {
			fullatom_relax_->set_nmoves(6);
		}
	}
	
void SurfaceDockingProtocol::setup_abinitio()
	{
		using namespace core;
		using namespace basic::options;
		//ClassicAbinitio
		TR << "Setting up classic abinitio" << std::endl;
		protocols::abinitio::ClassicAbinitio::register_options();
		
		//if option[ OptionKeys::in::file::frag9 ]
		std::string frag_large_file = option[ OptionKeys::in::file::frag9 ]();
		std::string frag_small_file = option[ OptionKeys::in::file::frag3 ]();
		TR<<"Loading Fragment 9:"<<std::endl;
		core::fragment::ConstantLengthFragSetOP fragset_large = new core::fragment::ConstantLengthFragSet;
		fragset_large->read_fragment_file(frag_large_file);
		TR<<"Loading Fragment 3:"<<std::endl;
		core::fragment::ConstantLengthFragSetOP fragset_small = new core::fragment::ConstantLengthFragSet;
		fragset_small->read_fragment_file(frag_small_file);
		core::kinematics::MoveMapOP movemap = new kinematics::MoveMap;
		movemap->set_bb( true );
		abinitio_ = new protocols::abinitio::ClassicAbinitio(fragset_small,fragset_large,movemap);
		abinitio_->set_score_weight(scoring::rg,0);
	}
	
void SurfaceDockingProtocol::setup_slide_movers( core::pose::Pose const & pose )
	{
		Vector protein_centroid, surf_centroid;
		Size const rb_jump=pose.num_jump();
			
		Vector const slide_axis = surface_parameters_->slide_axis();
		Vector slide_into, slide_away;
			
		protocols::geometry::centroids_by_jump (pose, rb_jump, surf_centroid, protein_centroid);
			
		//asses which side of the surface the protein is and pick the correct direction for the sliding vectors accordingly
		if (surf_centroid.distance_squared(protein_centroid + slide_axis) < surf_centroid.distance_squared(protein_centroid))
		{
			slide_into = surface_parameters_->slide_axis();
		}
		else
		{
			slide_into = surface_parameters_->slide_axis().negated();
		}
			
		slide_away = slide_into.negated();
		
		//setting the slide axis in the shared surface parameters so everyone can agree which side of the surface the protein should be!
		surface_parameters_->set_slide_axis(slide_into);
		slide_away_from_surface_ = new protocols::rigid::RigidBodyTransMover( slide_away, pose.num_jump());
		slide_away_from_surface_->step_size(20);
		
		slide_into_surface_ = new protocols::docking::FaDockingSlideIntoContact( pose.num_jump(), slide_into.negated());
		
		//getting a point 30 angstroms above surface centroid
		Vector point_above = surf_centroid+slide_away.normalized()*100;
		//vector needed to move protein centroid to point above
		Vector position_above = point_above-protein_centroid;
		position_above_surface_ = new protocols::rigid::RigidBodyTransMover(position_above, pose.num_jump());
		position_above_surface_->step_size(position_above.magnitude());
		//TODO: going to manually bring the centroids of the protein and surface within a reasonable distance of eachother before beginning slide-into
			
	}
	
void SurfaceDockingProtocol::split_protein_surface_poses (core::pose::Pose const & pose, core::pose::Pose & surface, core::pose::Pose & protein )
	{
		//split the pose into separate chains
		utility::vector1< pose::PoseOP > singlechain_poses;
		singlechain_poses = pose.split_by_chain();
		assert( valid_surface_pose( pose ));
		singlechain_poses[ 2 ]->set_new_energies_object( new core::scoring::Energies ); // replace the SurfaceEnergies object
		surface = *singlechain_poses[ 1 ];
		protein = *singlechain_poses[ 2 ];
	}
	
void SurfaceDockingProtocol::merge_protein_surface_poses(core::pose::Pose & pose, const core::pose::Pose & surface, const core::pose::Pose & protein)
	{
		pose.copy_segment((protein).total_residue(),
								(protein),
								(surface).total_residue()+1,1);
		
		set_secondary_structure(pose);		
	}
	
core::pack::task::PackerTaskOP
SurfaceDockingProtocol::create_surface_packer_task ( core::pose::Pose const & pose, Size const first_protein_residue )
	{
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
		return my_task_fullatom;
	}
	


}	//surface_docking

}	//protocols


//TODO: Add check for surface pose correctness, fix #includes, make tracers neater, implement register_options(), implement show(), remove surface params from list of cacheable data types, add 2 additional adsorbed refinement cycles to match emily
