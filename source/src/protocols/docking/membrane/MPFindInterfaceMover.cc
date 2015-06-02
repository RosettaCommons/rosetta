// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/docking/membrane/MPFindInterfaceMoverCreator.hh
/// @brief      Sample protein-protein interface in the membrane
/// @details	The foldtree after the mover is reset to the original one as it
///				was before - so it should work with both fixed and movable membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_cc
#define INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_cc

// Unit Headers
#include <protocols/docking/membrane/MPFindInterfaceMover.hh>
#include <protocols/docking/membrane/MPFindInterfaceMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/FlipMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/MPQuickRelaxMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/SpinAroundPartnerMover.hh>
#include <protocols/membrane/TiltMover.hh>
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/util.hh>
#include <protocols/docking/util.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh> 
#include <core/conformation/membrane/SpanningTopology.hh> 
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/pose/util.hh>
#include <protocols/docking/metrics.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.MPFindInterfaceMover" );

namespace protocols {
namespace docking {
namespace membrane {
	
using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
using namespace protocols::docking;
using namespace protocols::membrane;
using namespace protocols::membrane::geometry;
		
/////////////////////
/// Constructors  ///
/////////////////////

	
/// @brief Default Constructor
/// @details Takes topology from MembraneInfo, gets the jump from docking partners
MPFindInterfaceMover::MPFindInterfaceMover() :
	protocols::moves::Mover(),
	topo_( new SpanningTopology() ),
	topo_up_( new SpanningTopology() ),
	topo_down_( new SpanningTopology() )
{}
	
/// @brief Copy Constructor
MPFindInterfaceMover::MPFindInterfaceMover( MPFindInterfaceMover const & src ) :
	protocols::moves::Mover( src ),
	native_( src.native_ ),
	partners_( src.partners_ ),
	jump_( src.jump_ ),
	jumps_( src.jumps_ ),
	topo_( src.topo_ ),
	topo_up_( src.topo_up_ ),
	topo_down_( src.topo_down_ ),
	sfxn_( src.sfxn_ ),
	flips_( src.flips_ ),
	flexible_bb_( src.flexible_bb_ )
{}
	
/// @brief Destructor
MPFindInterfaceMover::~MPFindInterfaceMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPFindInterfaceMover::clone() const {
	return ( protocols::moves::MoverOP( new MPFindInterfaceMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPFindInterfaceMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPFindInterfaceMover() );
}

/// @brief Parse Rosetta Scripts Options for this Mover
void
MPFindInterfaceMover::parse_my_tag(
	 utility::tag::TagCOP tag,
	 basic::datacache::DataMap &,
	 protocols::filters::Filters_map const &,
	 protocols::moves::Movers_map const &,
	 core::pose::Pose const &
	 ) {

	// Read in docking partners
	if ( tag->hasOption( "partners" ) ) {
		partners_ = tag->getOption< std::string >( "partners" );
	}
	
	// TODO: to implement
	
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPFindInterfaceMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPFindInterfaceMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPFindInterfaceMoverCreator::keyname() const {
	return MPFindInterfaceMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPFindInterfaceMoverCreator::mover_name() {
	return "MPFindInterfaceMover";
}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPFindInterfaceMover)
std::string
MPFindInterfaceMover::get_name() const {
	return "MPFindInterfaceMover";
}

/// @brief Sampling protein-protein interface in the membrane
void
MPFindInterfaceMover::apply( Pose & pose ) {
	
	using namespace numeric;
	using namespace basic::options;
	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;
	using namespace protocols::rigid;
	using namespace protocols::docking;
	using namespace protocols::scoring;
	using namespace protocols::simple_moves;

	TR << "Sampling protein-protein interface in the membrane..." << std::endl;

	// initializations
	register_options();
	init_from_cmd();
	finalize_setup( pose );

	// compute embedding for partners (compute structure-based embedding with split topologies)
	EmbeddingDef emb_up, emb_down;
	update_partner_embeddings( pose, jump_, emb_up, emb_down );
	
	core::Vector center = pose.conformation().membrane_info()->membrane_center();
	core::Vector normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 1: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

//	// get initial membrane position
//	MembranePositionFromTopologyMoverOP mempos( new MembranePositionFromTopologyMover() );
//	mempos->apply( pose );
//	
//	// optimize membrane position using smooth scorefunction
//	OptimizeMembranePositionMoverOP optmem( new OptimizeMembranePositionMover() );
//	optmem->apply( pose );
	
	// superimpose upstream partner with the native
	superimpose_upstream_partner( pose );

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 3a: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	////////////////////// MOVE PARTNERS APART //////////////////////////
	
	TR << "Moving partners apart..." << std::endl;
	
	// get membrane axis from docking metrics function
	core::Vector slide_axis = membrane_axis( pose, jump_ );
//	slide_axis.z( 0 );
	bool vary_stepsize = false;
	
	// get distance between points
	core::Real dist1 = ( emb_down.center() - emb_up.center() ).length();
	TR << "distance between partners: " << dist1 << std::endl;
	
	slide_axis.negate();
	slide_axis.normalize();
	TR << "slide axis: " << slide_axis.to_string() << std::endl;
	
	// move apart
	TR << "Moving apart" << std::endl;
	RigidBodyTransMoverOP mover( new rigid::RigidBodyTransMover( slide_axis, jump_, vary_stepsize ) );
	mover->step_size( 100 );
	mover->apply( pose );

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 2: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;
	
	///////////////////////// RUN RELAX ////////////////////////////
	
	// run quick relax
	MPQuickRelaxMoverOP relax( new MPQuickRelaxMover() );
	if ( flexible_bb_ == true ) {
		relax->add_membrane_again( false );
		relax->membrane_from_topology( false );
		relax->optimize_membrane( false );
		relax->apply( pose );
	}

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 3: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;
	pose.dump_pdb( "after_relax.pdb" );

	// superimpose upstream partner with the native
	superimpose_upstream_partner( pose );

	// set membrane to fixed one again
	SetMembranePositionMoverOP setmem( new SetMembranePositionMover() );
	setmem->apply( pose );

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 3a: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

//	//////////////// MOVE INTO FIXED MEMBRANE /////////////////////
//
//	// translate and rotate pose into membrane
//	EmbeddingDefOP embed( compute_structure_based_embedding( pose ) );
//	core::Size memjump = pose.conformation().membrane_info()->membrane_jump();
//	TranslationRotationMoverOP rt1( new TranslationRotationMover( embed->center(), embed->normal(), mem_center, mem_normal, memjump ) );
//	rt1->apply( pose );

	//////////////// CENTROID OR FULL-ATOM /////////////////////
	
	if ( option[ OptionKeys::mp::dock::lowres ].user() ) {
		SwitchResidueTypeSetMoverOP centroid( new SwitchResidueTypeSetMover( "centroid" ) );
		centroid->apply( pose );
	}

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 4: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	TR << "foltree before protocol: " << std::endl;
	pose.fold_tree().show( TR );

	//////////////// RUN REST OF PROTOCOL ////////////////////////
	
	// create MC object
	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *sfxn_, 1.0 ) );

	// do this for a certain number of iterations
	for ( Size i = 1; i <= 10; ++i ){
		
		center = pose.conformation().membrane_info()->membrane_center();
		normal = pose.conformation().membrane_info()->membrane_normal();
		TR << "membrane center 5: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

		// SPIN MOVER
		// get a random spin angle between 0 and 360 degrees
		int spin_angle( numeric::random::random_range( 0, 360 ) );
		TR << "=========================SPIN MOVER==============================" << std::endl;
		TR << "random spin angle: " << spin_angle << std::endl;
		
		// spin downstream partner around spin angle
		update_partner_embeddings( pose, jump_, emb_up, emb_down );
		RigidBodyDeterministicSpinMoverOP spin( new RigidBodyDeterministicSpinMover(
					jump_, emb_down.normal(), emb_down.center(), spin_angle ) );
		spin->apply( pose );
		
		// slide into contact
		TR << "=========================SLIDE-INTO-CONTACT======================" << std::endl;
		DockingSlideIntoContactOP slide( new DockingSlideIntoContact( jump_ ) );
		slide->apply( pose );
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
		
		// SPIN-AROUND-PARTNER-MOVER!!!
		TR << "=========================SPIN-AROUND-PARTNER MOVER===============" << std::endl;
		SpinAroundPartnerMoverOP around( new SpinAroundPartnerMover( jump_, 100 ) );
		around->apply( pose );
		
		// slide into contact
		TR << "=========================SLIDE-INTO-CONTACT======================" << std::endl;
		slide->apply( pose );
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
		
		// TILT MOVER
		// get axis perpendicular to downstream partner embedding normal and
		// perpendicular to vector between embedding centers of two partners
		TR << "=========================TILT MOVER==============================" << std::endl;
		TiltMoverOP tilt( new TiltMover( jump_ ) );
		tilt->apply( pose );
		
		// slide into contact
		TR << "=========================SLIDE-INTO-CONTACT======================" << std::endl;
		slide->apply( pose );
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
		
		// FLIP MOVER, if allowed
		// FlipMover calculates embedding automatically
		if ( flips_ == true ) {
			TR << "=========================FLIP MOVER==============================" << std::endl;
			protocols::membrane::FlipMoverOP flip( new FlipMover( jump_ ) );
			flip->set_random_membrane_flip_angle();
			flip->apply( pose );
		}
		else {
			TR << "====================SMALL TILTS WITH FLIP MOVER==================" << std::endl;
			update_partner_embeddings( pose, jump_, emb_up, emb_down );
			core::Vector mem_norm = pose.conformation().membrane_info()->membrane_normal();
			core::Real angle_mem_emb = numeric::conversions::degrees( angle_of( mem_norm, emb_down.normal() ) );
			
			// if angle between membrane normal and embedding is smaller 90
			int max( 45 );
			core::Real angle;
			if ( angle_mem_emb <= 90 ) {
				angle = numeric::random::random_range( max - nearest_int( angle_mem_emb ), max + nearest_int( angle_mem_emb ) );
			}
			else {
				angle = numeric::random::random_range( max + 180 - nearest_int( angle_mem_emb ), max - 180 + nearest_int( angle_mem_emb ) );
			}
			protocols::membrane::FlipMoverOP flip( new FlipMover( jump_, angle ) );
			flip->apply( pose );
		}
		
		// slide into contact
		TR << "=========================SLIDE-INTO-CONTACT======================" << std::endl;
		slide->apply( pose );
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
		
		// SPIN-AROUND-PARTNER-MOVER!!!
		TR << "=========================SPIN-AROUND-PARTNER MOVER===============" << std::endl;
		SpinAroundPartnerMoverOP around1( new SpinAroundPartnerMover( jump_, 100 ) );
		around1->apply( pose );
		
		// slide into contact
		TR << "=========================SLIDE-INTO-CONTACT======================" << std::endl;
		slide->apply( pose );
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;
		
	} // number of iterations for search

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 6: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	// switch back to full-atom
	if ( pose.is_centroid() ) {
		SwitchResidueTypeSetMoverOP fa( new SwitchResidueTypeSetMover( "fa_standard" ) );
		fa->apply( pose );
	}

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 7: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	//////////////////////// RELAX //////////////////////

	// do another round of relax
	if ( flexible_bb_ == true ) {
		relax->add_membrane_again( false );
		relax->membrane_from_topology( false );
		relax->optimize_membrane( false );
		relax->apply( pose );
	}

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 8: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

//	// optimize membrane position using smooth scorefunction
//	optmem->apply( pose );
//	
//	center = pose.conformation().membrane_info()->membrane_center();
//	normal = pose.conformation().membrane_info()->membrane_normal();
//	TR << "membrane center 8a: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;
//	
//	//////////////// MOVE INTO FIXED MEMBRANE /////////////////////
//	
//	// translate and rotate pose into membrane
//	embed = compute_structure_based_embedding( pose );
//	memjump = pose.conformation().membrane_info()->membrane_jump();
//	TranslationRotationMoverOP rt2( new TranslationRotationMover( embed->center(), embed->normal(), mem_center, mem_normal, memjump ) );
//	rt2->apply( pose );
//
//	center = pose.conformation().membrane_info()->membrane_center();
//	normal = pose.conformation().membrane_info()->membrane_normal();
//	TR << "membrane center 8b " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	// superimpose upstream partner with the native
	superimpose_upstream_partner( pose );

	// set membrane to fixed one again
	SetMembranePositionMoverOP setmem2( new SetMembranePositionMover() );
	setmem2->apply( pose );

	center = pose.conformation().membrane_info()->membrane_center();
	normal = pose.conformation().membrane_info()->membrane_normal();
	TR << "membrane center 9: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

	/////////////////////// WRITE QUALITY MEASURES TO SCORE FILE ///////////////

	// score the pose - this hopefully writes this into the scorefile
	core::scoring::ScoreFunctionOP highres_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_docking_fa_2015.wts" );
	( *highres_sfxn )( pose );

	// get job
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// calculate and store the rmsd in the score file
	job->add_string_real_pair( "Lrms", protocols::docking::calc_Lrmsd( pose, native_, jumps_ ) );
	job->add_string_real_pair( "P1rms", protocols::docking::calc_P1rmsd( pose, native_, jumps_ ) );
	job->add_string_real_pair( "P2rms", protocols::docking::calc_P2rmsd( pose, native_, jumps_ ) );
	job->add_string_real_pair( "Irms", protocols::docking::calc_Irmsd( pose, native_, sfxn_, jumps_ ) );

	// get interface
	Interface interface = Interface( jump_ );
	interface.calculate( pose );

	job->add_string_real_pair( "Ncntct", interface.interface_nres() );
//	job->add_string_real_pair( "SASA", sasa->compute( pose ) );
//	job->add_string_real_pair( "intSASA", calculate_interface_SASA( pose, interface ) );

	job->add_string_real_pair( "Fnat", protocols::docking::calc_Fnat( pose, native_, sfxn_, jumps_ ) );
	job->add_string_real_pair( "Fnonnat", protocols::docking::calc_Fnonnat( pose, native_, sfxn_, jumps_ ) );

	if ( pose.residue( nres_protein( pose ) ).atom("CA").xyz().z() > 0.0 ){
		job->add_string_real_pair( "in/out", 0.0 );
	}
	else {
		job->add_string_real_pair( "in/out", 1.0 );
	}
	
	// surface complementarity
	core::scoring::sc::ShapeComplementarityCalculator shape;
	shape.Init();
	job->add_string_real_pair( "Shape", shape.CalcSc( pose, jump_, 1 ) );
	
	// difference of fractions of small residues in the interface
	core::Real small_frct_diff = fractions_small_residues( pose, interface).first - fractions_small_residues( pose, interface).second;
	job->add_string_real_pair( "small_res", small_frct_diff );

}// apply

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options with JD2
void
MPFindInterfaceMover::register_options() {
	
	using namespace basic::options;

	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::docking::partners );
	option.add_relevant( OptionKeys::mp::dock::lowres );
	option.add_relevant( OptionKeys::mp::dock::highres );
	option.add_relevant( OptionKeys::mp::dock::allow_flips );
	option.add_relevant( OptionKeys::mp::dock::flexible_bb );
	
}

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize Mover options from the comandline
void
MPFindInterfaceMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// docking partners
	if( option[ OptionKeys::docking::partners ].user() ) {
		partners_ = option[ OptionKeys::docking::partners ]();
	}
	
	// native
	if( option[ OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_, option[ OptionKeys::in::file::native ]() );
	}

	// setup scorefunction
	if ( option[ OptionKeys::mp::dock::lowres ].user() ) {
		sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_docking_cen_2015.wts" );
	}
	else if ( option[ OptionKeys::mp::dock::highres ].user() ) {
		sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_docking_fa_2015.wts" );
	}
	else {
		utility_exit_with_message( "Have to specify score function: either use -mp::dock:lowres or -mp::dock:highres flag!" );
	}
	
	// allow flips in the membrane?
	flips_ = true;
	if ( option[ OptionKeys::mp::dock::allow_flips ].user() ) {
		flips_ = option[ OptionKeys::mp::dock::allow_flips ]();
	}

	// allow flips in the membrane?
	flexible_bb_ = false;
	if ( option[ OptionKeys::mp::dock::flexible_bb ].user() ) {
		flexible_bb_ = option[ OptionKeys::mp::dock::flexible_bb ]();
	}

}// init from cmd

////////////////////////////////////////////////////////////////////////////////

/// @brief Finalize setup
void MPFindInterfaceMover::finalize_setup( Pose & pose ) {

	using namespace utility;
	using namespace basic::options;
	using namespace protocols::simple_moves;
	TR << "Finalizing setup... " << std::endl;

	////// SUPERIMPOSE THE FIRST PARTNER OF THE POSE WITH THE NATIVE ///////

	// call AddMembraneMover on native for RMSD calculation
	AddMembraneMoverOP addmem( new AddMembraneMover() );
	addmem->apply( native_ );
	
	// superimpose upstream partner with the native
	superimpose_upstream_partner( pose );
	
	// for debugging
	pose.dump_pdb("pose_superimposed_to_native.pdb");
	TR << "dumped superimposed pose" << std::endl;

	///////////////////// FOLDTREE STUFF ///////////////////////

	// get foldtree from partners (setup_foldtree) and movable jump
	// we are using this function to add the membrane add the anchor point
	// closest to COM
	setup_foldtree( pose, partners_, jumps_);
	jump_ = jumps_[1];
	TR << "jump_ from foldtree: " << jump_ << std::endl;
	
	// get anchor point for membrane from jump; I think this is residue closest to
	// COM of the upstream partner
	int anchor = pose.fold_tree().upstream_jump_residue( jump_ );
	TR << "anchor point: " << anchor << std::endl;

	// Add Membrane, appends MEM as jump1
	AddMembraneMoverOP add_memb( new AddMembraneMover( anchor, 0 ) );
	add_memb->apply( pose );
	
	// reorder foldtree to have membrane at root
	core::kinematics::FoldTree ft = pose.fold_tree();
	ft.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( ft );
	TR << "reordered foltree: " << std::endl;
	pose.fold_tree().show( TR );
	TR << "jump: " << jump_ << std::endl;
	
	//////////////// TOPOLOGY //////////////////////
	
	// get topology
	topo_ = pose.conformation().membrane_info()->spanning_topology();

	// splitting topology by jump into upstream and downstream topology
	split_topology_by_jump_noshift( pose, jump_, topo_, topo_up_, topo_down_ );

} // finalize setup

////////////////////////////////////////////////////////////////////////////////

/// @brief Superimpose upstream partner
/// @details Superimpose upstream partner of the pose with the native
void MPFindInterfaceMover::superimpose_upstream_partner( Pose & pose ) {
	
	using namespace protocols::simple_moves;
	
	// get partner 1: ALL CHAINS MUST BE IN PDB CONSECUTIVELY!!!
	utility::vector1< std::string > partners( utility::string_split( partners_, '_' ) );
	utility::vector1< std::string > partner1( utility::split( partners[1] ) );
	
	// get residue range for superposition: get start residue
	Size start(0);
	for ( Size i = 1; i <= pose.total_residue(); ++i ){
		if ( start == 0 &&
			partner1[1] == utility::to_string( pose.pdb_info()->chain( i ) ) ){
			start = i;
		}
	}
	
	// get end residue
	Size end(0);
	for ( Size j = pose.total_residue(); j >= 1; --j ){
		if ( end == 0 &&
			partner1[partner1.size()] == utility::to_string( pose.pdb_info()->chain( j ) ) ){
			end = j;
		}
	}
	TR << "range start: " << start << ", end: " << end << std::endl;
	
	// superimpose partner 1 with starting pose for easier rmsd calculation
	SuperimposeMoverOP super( new SuperimposeMover( native_, start, end, start, end, true ) );
	super->apply( pose );
	
} // superimpose_upstream_partner

////////////////////////////////////////////////////////////////////////////////

/// @brief Calculate interface SASA
/// @details Calculate SASA buried in the interface
core::Real MPFindInterfaceMover::calculate_interface_SASA( Pose & pose, Interface & interface ) {

	using namespace numeric;
	using namespace core::pose;
	using namespace core::scoring::sasa;
	using namespace protocols::scoring;
//	using namespace protocols::simple_filters;
	
	// get per-residue SASA
	SasaCalc calc = SasaCalc();
	utility::vector1< core::Real > sasa = calc.get_residue_sasa();

	// initialize SASA of interface
	core::Real interface_sasa(0);

	// go through residues and add SASA of that residue to total SASA
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		interface_sasa += interface.is_interface( i ) * sasa[ i ];
	}
	
	// get SASA
	//	InterfaceSasaFilterOP sasa( new InterfaceSasaFilter() );
	//	sasa->jump( jump_ );

	return interface_sasa;

} // calculate interface SASA

////////////////////////////////////////////////////////////////////////////////
/// @brief Calculate fraction of small residues in interface
/// details See Andrew Bordner, 2009, BMC Bioinformatics
std::pair< core::Real, core::Real > MPFindInterfaceMover::fractions_small_residues( Pose & pose, Interface & interface ) {
	
	core::Size cnt_intf_res( 0 );
	core::Size cnt_nonintf_res( 0 );
	core::Size small_intf_res( 0 );
	core::Size small_nonintf_res( 0 );
	
	// iterate over residues in protein
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		// check if residue is in interface
		bool intf( false );
		intf = interface.is_interface( i );
		
		// check if residue is small
		bool small( false );
		std::string res = pose.residue( i ).name3();
		if ( res == "ALA" || res == "GLY" || res == "LEU" || res == "VAL" ){
			small = true;
		}

		// count small and in interface
		if ( small == true && intf == true ) {
			++small_intf_res;
		}
		// count small and not in interface
		else if ( small == true && intf == false ) {
			++small_nonintf_res;
		}
		
		// increment total counters
		if ( intf == true ) {
			++cnt_intf_res;
		}
		else {
			++cnt_nonintf_res;
		}
	}
	
	// get fractions
	core::Real frct_small_intf = static_cast< core::Real > ( small_intf_res );
	frct_small_intf /= cnt_intf_res;
	core::Real frct_small_nonintf = static_cast< core::Real > ( small_nonintf_res );
	frct_small_nonintf /= cnt_nonintf_res;
	
	TR << "frct_small_intf " << frct_small_intf << " frct_small_nonintf " << frct_small_nonintf << std::endl;
	
	// return pair
	std::pair< core::Real, core::Real > fractions( frct_small_intf, frct_small_nonintf );
	return fractions;
	
}// fractions small residues

	
} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPFindInterfaceMover_cc
