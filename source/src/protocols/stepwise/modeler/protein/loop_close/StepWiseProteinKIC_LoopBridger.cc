// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseProteinKIC_LoopBridger
/// @brief Makes a list of (phi, psi, omega) at moving_residues that
///              could be useful for full-atom packing
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinKIC_LoopBridger.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/VariantType.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/Ramachandran.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/loops/Loop.hh>
#include <numeric/kinematic_closure/bridgeObjects.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <numeric/conversions.hh>

#include <utility/exit.hh>

#ifdef WIN32
#include <time.h>
#endif

#include <string>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Good ol' kinematic loop closure. Applied to loops whose N-terminal
// and C-terminal segment have been modeled in, e.g., different silent
// files.
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.protein.loop_close.StepWiseProteinKIC_LoopBridger" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace loop_close {


//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseProteinKIC_LoopBridger::StepWiseProteinKIC_LoopBridger( sampler::StepWiseSamplerSizedOP sampler,
	working_parameters::StepWiseWorkingParametersCOP working_parameters ):
	sampler_( sampler ),
	working_bridge_res_( working_parameters->working_bridge_res() ),
	is_pre_proline_( working_parameters->is_pre_proline() ),
	num_perturb_steps_( 0 ), // perturbations of 'takeoff' psi and phi -- currently disabled.
	perturb_torsion_( 20.0 ),
	idl_CA_C_N_(116.2), // taken from hard-coded numbers in KIC code.
	idl_C_N_CA_(121.7), // taken from hard-coded numbers in KIC code.
	idl_C_N_(1.32869), // taken from hard-coded numbers in KIC code.
	OMEGA_MEAN_(179.8), // taken from hard-coded numbers in KIC code.
	use_icoor_geometry_( false ), //this would supercede above vals with ICOOR values (actually they're very similar)
	verbose_( false )
{

	if ( working_bridge_res_.size() != 3 ) utility_exit_with_message( "Must supply three bridge residues that are covered by -sample_res, -input_res1, and -input_res2!");
	initialize_is_fixed_res( working_parameters->working_fixed_res(), working_parameters->working_sequence() );

}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseProteinKIC_LoopBridger::~StepWiseProteinKIC_LoopBridger()
{}
/////////////////////
std::string
StepWiseProteinKIC_LoopBridger::get_name() const {
	return "StepWiseProteinKIC_LoopBridger";
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::apply( core::pose::Pose & pose )
{
	clock_t const time_start( clock() );

	Pose pose_save = pose;

	setup_torsions( pose );
	figure_out_loop( pose );
	pose_count_ = 0;

	for ( sampler_->reset(); sampler_->not_end(); ++( *sampler_ ) ) {
		sampler_->apply( pose );
		KIC_loop_close_with_perturbations( pose );
	}

	// Kind of silly -- should have at least one output.
	if ( main_chain_torsion_sets_for_moving_residues_.size() == 0  ) grab_main_chain_torsion_set_list( pose );

	std::cout << "Total time in StepWiseProteinKIC_LoopBridger: " <<
		static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

	pose = pose_save;
}


///////////////////////////////////////////////////////////////////////////
// Should just do the KIC loop closure, but added the option of perturbing
// the takeoff psi (at the N-terminal end of the loop) and landing phi (at the
// the C-terminal end of the loop ).
void
StepWiseProteinKIC_LoopBridger::KIC_loop_close_with_perturbations( core::pose::Pose & pose ) {

	if ( num_perturb_steps_ == 0 ) {
		std::cout << "Loop combination " << pose_count_++ << ". " ;
		KIC_loop_close( pose );
		return;
	}

	Size const pre_loop_res  = working_bridge_res_[1] - 1;
	Size const post_loop_res = working_bridge_res_[3] + 1 ;
	Real const psi_start = pose.psi( pre_loop_res  );
	Real const phi_start = pose.phi( post_loop_res );

	for ( int offset1 = -1 * num_perturb_steps_; offset1 <= num_perturb_steps_; offset1++ ) {
		Real const psi_perturb = psi_start + Real( offset1 ) * perturb_torsion_;
		pose.set_psi( pre_loop_res, psi_perturb );

		for ( int offset2 = -1 * num_perturb_steps_; offset2 <= num_perturb_steps_; offset2++ ) {

			Real const phi_perturb = phi_start + Real(offset2) * perturb_torsion_ ;
			pose.set_phi( post_loop_res, phi_perturb );

			//std::cout << "Psi of pre_loop_res  " << pre_loop_res << " is " << psi_perturb << std::endl;
			//std::cout << "Phi of post_loop_res " << post_loop_res << " is " << phi_perturb << std::endl;

			std::cout << "Loop combination " << pose_count_++ << ". " ;
			//pose.dump_pdb( "S_" + ObjexxFCL::string_of( pose_count_ ) + ".pdb" );

			// Do it!
			KIC_loop_close( pose );
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::figure_out_loop( core::pose::Pose const & pose ){

	using namespace protocols::loops;
	using namespace core::chemical;

	Size const & middle_bridge_res = working_bridge_res_[ 2 ];
	assert( working_bridge_res_[ 1 ] = middle_bridge_res - 1 );
	assert( working_bridge_res_[ 3 ] = middle_bridge_res + 1 );

	// consistency checks.
	if ( middle_bridge_res == 0 ) {
		utility_exit_with_message( "Is there a bridge residue that is not an input residue?" );
	}

	// Check that there is a cutpoint in here.
	std::cout << "Found bridge residue: " << middle_bridge_res << std::endl;
	std::cout << pose.fold_tree() << std::endl;
	std::cout << pose.annotated_sequence( true ) << std::endl;
	Size cutpoint_ = 0;
	for ( int offset = -2; offset < 1; offset++ ) {
		Size const test_res = static_cast<int>( middle_bridge_res ) + offset;
		if ( test_res >= 1 && test_res <= pose.size() ) {
			if ( pose.fold_tree().is_cutpoint( test_res ) ) {
				cutpoint_ = test_res;
				if ( !pose.residue_type( cutpoint_ ).has_variant_type( CUTPOINT_LOWER ) ||
						!pose.residue_type( cutpoint_+1 ).has_variant_type( CUTPOINT_UPPER ) ) {
					std::cout << " cutpoint res? " << cutpoint_ << std::endl;
					utility_exit_with_message( "cutpoints not set up properly at cutpoint residue" );
				}
			}
		}
	}
	if ( cutpoint_ == 0 ) utility_exit_with_message( "could not find cutpoint!" );
	std::cout << "Found cutpoint residue: " << cutpoint_ << std::endl;

	loop_ = Loop( middle_bridge_res-1 /*start*/, middle_bridge_res+1 /*end*/, cutpoint_ /*cutpoint?*/ );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::setup_torsions( pose::Pose const & pose ){

	// Need to fill torsions for moving and bridge residues.
	using namespace core::id;
	which_torsions_.clear();
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !is_fixed_res_[ n ] ) {
			// loop residues.
			for ( Size k = 1; k <= 3; k++ ) which_torsions_.push_back( TorsionID( n, BB, k ) );
		} else if ( n < pose.size() && !is_fixed_res_[ n+1 ] ) {
			// psi,omega of 'takeoff' residues
			for ( Size k = 1; k <= 3; k++ ) which_torsions_.push_back( TorsionID( n, BB, k ) );
		} else if ( n > 1 && !is_fixed_res_[ n-1 ] ) {
			// phi of 'landing' residues
			for ( Size k = 1; k <= 3; k++ ) which_torsions_.push_back( TorsionID( n, BB, k ) );
		}
	}
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::grab_main_chain_torsion_set_list( pose::Pose const & pose ){

	utility::vector1< Real > main_chain_torsion_set_for_moving_residues;
	for ( Size n = 1; n <= which_torsions_.size(); n++ )  {
		main_chain_torsion_set_for_moving_residues.push_back(  pose.torsion( which_torsions_[ n ] ) );
		std::cout << ' ' << pose.torsion( which_torsions_[ n ] );
	}
	std::cout << ' ' << std::endl;

	main_chain_torsion_sets_for_moving_residues_.push_back( main_chain_torsion_set_for_moving_residues );
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< utility::vector1< core::Real > > const &
StepWiseProteinKIC_LoopBridger::main_chain_torsion_set_lists() const
{
	return main_chain_torsion_sets_for_moving_residues_;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::TorsionID > const &
StepWiseProteinKIC_LoopBridger::which_torsions() const
{
	return which_torsions_;
}


///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::output_chainTORS( utility::vector1< core::Real > const & dt_ang,
	utility::vector1< core::Real > const & db_ang,
	utility::vector1< core::Real > const & db_len ) const {

	std::cout << "------  chainTORS output ---- " << std::endl;
	for ( Size i = 1; i <= ( dt_ang.size()/3) ; i++ ) {

		std::cout << "TORSIONS: ";
		for ( Size j = 1; j <= 3; j++ ) std::cout << ObjexxFCL::format::F(8,3,dt_ang[ 3*(i-1)+ j ]) << " ";

		std::cout << "   BOND_ANGLES: ";
		for ( Size j = 1; j <= 3; j++ ) std::cout << ObjexxFCL::format::F(8,3,db_ang[ 3*(i-1)+ j ]) << " ";

		std::cout << "   BOND_LENGTHS: ";
		for ( Size j = 1; j <= 3; j++ ) std::cout << ObjexxFCL::format::F(8,3,db_len[ 3*(i-1)+ j ]) << " ";

		std::cout << std::endl;
	}
}


/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::fill_chainTORS_info( pose::Pose const & pose,
	utility::vector1<utility::vector1<Real> > & atoms,
	utility::vector1<Real> & dt_ang,
	utility::vector1<Real> & db_ang,
	utility::vector1<Real> & db_len,
	Size const & start_res_ ,
	Size const & end_res_ ) const {

	using namespace numeric::kinematic_closure;

	if ( verbose_ ) std::cout << "About to run chainTORS" << std::endl;
	Size ind = 1;
	for ( Size i =  start_res_ - 1;  i <= end_res_ + 1;   i++ ) {
		if ( verbose_ ) std::cout << "Filling residue " << i << std::endl;
		conformation::Residue res = pose.residue(i);
		for ( Size j=1; j<=3; j++ ) { // DJM: just keeping N, CA, C atoms. We assume these are always the first 3.  BAD -- PROTEIN ONLY ASSUMPTION -- How about metal ions with only 1 atom?
			atoms[ind].resize(3);
			atoms[ind][1] = static_cast<Real> (res.xyz(j).x());
			atoms[ind][2] = static_cast<Real> (res.xyz(j).y());
			atoms[ind][3] = static_cast<Real> (res.xyz(j).z());
			ind++;
		}
	}

	utility::vector1<utility::vector1<Real> > Q0 (3);
	utility::vector1<Real> R0 (3);

	chainTORS(atoms.size(), atoms, dt_ang, db_ang, db_len, R0, Q0);

	if ( verbose_ ) output_chainTORS( dt_ang, db_ang, db_len );
}

///////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::KIC_loop_close( pose::Pose & pose ){

	using namespace core::kinematics;
	using namespace protocols::loops;
	using namespace numeric::kinematic_closure;

	///// kinematic loop close.
	// Following copied from, e.g., KinematicMover.cc.  Need to elaborate for terminal residues!
	// inputs to loop closure
	utility::vector1<utility::vector1<Real> > atoms;
	utility::vector1<Size> pivots (3), order (3);
	// for eliminating identical solutions
	utility::vector1<Real> dt_ang, db_len, db_ang, save_t_ang, save_b_len, save_b_ang;
	utility::vector1<Real> dummy_t_ang, dummy_b_ang, dummy_b_len;
	utility::vector1<conformation::ResidueOP> save_residues;

	start_res_ = loop_.start();
	middle_res_ = loop_.start() + 1;
	end_res_ = loop_.stop();
	middle_offset_ = middle_res_ - start_res_; // is used to set central pivot atom
	seg_len_ = end_res_ - start_res_ + 1;
	atoms.resize( (seg_len_ + 2) * 3); // one extra residue on each side to establish the geometric frame

	fill_chainTORS_info( pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_ );

	order[1]=1;
	order[2]=2;
	order[3]=3;

	// Set the pivot atoms
	Size pvatom1=5; // second C-alpha
	Size pvatom2=5 + (3 * middle_offset_); // middle res C-alpha
	Size pvatom3=(3 * (seg_len_+1)) - 1; // second-to-last C-alpha
	pivots[1]=pvatom1;
	pivots[2]=pvatom2;
	pivots[3]=pvatom3;

	// Need to fix bond lengths and angles at cutpoint
	Size const cut_offset_ = loop_.cut() - start_res_;
	dt_ang[ 3 + 3*cut_offset_ + 3 ] = OMEGA_MEAN_;
	db_len[ 3 + 3*cut_offset_ + 3 ] = idl_C_N_;
	db_ang[ 3 + 3*cut_offset_ + 3 ] = idl_CA_C_N_;
	db_ang[ 3 + 3*cut_offset_ + 4 ] = idl_C_N_CA_;

	if ( use_icoor_geometry_ ) {
		Size cutpoint = loop_.cut();
		Real const bond_angle1( pose.residue( cutpoint ).upper_connect().icoor().theta() );// CA-C=N bond angle
		Real const bond_angle2( pose.residue( cutpoint+1 ).lower_connect().icoor().theta() ); // C=N-CA bond angle
		Real const bond_length( pose.residue( cutpoint+1 ).lower_connect().icoor().d() ); // C=N distance
		dt_ang[ 3 + 3*cut_offset_ + 3 ] = 180.0;
		db_len[ 3 + 3*cut_offset_ + 3 ] = bond_length;
		db_ang[ 3 + 3*cut_offset_ + 3 ] = 180.0 - numeric::conversions::degrees(bond_angle1);
		db_ang[ 3 + 3*cut_offset_ + 4 ] = 180.0 - numeric::conversions::degrees(bond_angle2);
	}


	if ( verbose_ ) {
		std::cout << "After setting desired geometry at cutpoint, relative to loop_start " << cut_offset_ << std::endl;
		output_chainTORS( dt_ang, db_ang, db_len );
	}

	///////////////////////////////////
	// Perform loop closure
	///////////////////////////////////
	// outputs from loop closure
	sample_omega_recursively( pose, -1, atoms, dt_ang, db_ang, db_len, pivots, order );

	// just for output.
	if ( verbose_ ) fill_chainTORS_info( pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_ );
}


/////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::sample_omega_recursively(
	pose::Pose & pose,
	int const offset,
	utility::vector1<utility::vector1<Real> > & atoms,
	utility::vector1<Real> & dt_ang,
	utility::vector1<Real> & db_ang,
	utility::vector1<Real> & db_len,
	utility::vector1< Size > const & pivots,
	utility::vector1< Size > const & order ){

	using namespace numeric::kinematic_closure;
	if ( offset == 2 ) {
		/* hardwired!! -- check cis-omega at residue before bridge residue (offset = 0 ),
		and at bridge residue ( offset = 1 ), and that's it */
		utility::vector1<utility::vector1<Real> > t_ang, b_ang, b_len;
		int nsol=0;

		if ( verbose_ ) std::cout << "About to run bridgeObjects" << std::endl;
		bridgeObjects(atoms, dt_ang, db_ang, db_len, pivots, order, t_ang, b_ang, b_len, nsol);
		if ( verbose_ ) std::cout << "Finished bridgeObjects" << std::endl;

		Size const num_solutions =  t_ang.size();
		std::cout << "Kinematic loop closure found this many solutions: " << num_solutions << std::endl;

		for ( Size i = 1; i <= num_solutions; i++ ) {

			for ( core::Size res = 0; res < seg_len_; res++ ) {
				pose.set_phi ( start_res_ + res, t_ang[ i ][ (3*(res+1)) + 1 ] );
				pose.set_psi ( start_res_ + res, t_ang[ i ][ (3*(res+1)) + 2 ] );
				pose.set_omega( start_res_ + res, t_ang[ i ][ (3*(res+1)) + 3 ] );
			}

			grab_main_chain_torsion_set_list( pose );

			if ( verbose_ ) pose.dump_pdb( "KIC_"+ ObjexxFCL::string_of( i )+".pdb" );

			// readout for checking.
			fill_chainTORS_info( pose, atoms, dt_ang, db_ang, db_len, start_res_, end_res_ );
		}
	} else {
		dt_ang[ 3 + 3*offset + 3 ] = OMEGA_MEAN_;
		sample_omega_recursively( pose, offset + 1, atoms, dt_ang, db_ang, db_len, pivots, order );

		//   std::cout << "RES " << pose.sequence()[ start_res_ + offset - 1 ] << " " << pose.sequence()[ start_res_ + offset ] << " " << working_parameters_->is_pre_proline( start_res_+ offset ) << std::endl;

		if ( is_pre_proline_[ start_res_ + offset ] ) {
			//    std::cout << " sample cis omega: " << start_res_ + offset;
			dt_ang[ 3 + 3*offset + 3 ] = 0.0 /*OMEGA_MEAN_*/;
			sample_omega_recursively( pose, offset + 1, atoms, dt_ang, db_ang, db_len, pivots, order );
		}
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinKIC_LoopBridger::initialize_is_fixed_res( utility::vector1< core::Size > const & fixed_res, std::string const & working_sequence ){

	is_fixed_res_.clear();
	Size const nres = core::pose::rna::remove_bracketed( working_sequence ).size();
	for ( Size n = 1; n <= nres; ++n ) {
		is_fixed_res_.push_back( false );
	}

	for ( Size i = 1; i <= fixed_res.size(); i++ ) {
		is_fixed_res_[ fixed_res[i] ] = true;
	}
}


} //loop_close
} //protein
} //modeler
} //stepwise
} //protocols
