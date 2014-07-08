// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RigidBodySampler
/// @brief Not particularly fancy, just filters a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/legacy/sampling/RigidBodySampler.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/sampling/rna/util.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Stub.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <protocols/farna/RNA_BasePairClassifier.hh>

#include <ObjexxFCL/string.functions.hh>

#include <list>
#include <iostream>

//Auto Headers
#include <core/id/AtomID.hh>

using namespace core;
using core::Real;

static basic::Tracer TR( "protocols.stepwise.legacy.sampling.RigidBodySampler" ) ;
using numeric::conversions::radians;
using numeric::conversions::degrees;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace sampling {

	////////////////////////////////////////////////////////////////////////////////////////////
	// Rigid body sample. Keep track of total number of states so that we can extract a Kd
	// Use input parameters to define fineness of sampling -- will look for convergence.
	// Save lowest energy states.
	//
	//  This only really works if we start with both the "reference residue" and the "moving residue" at the origin,
	//  with axes pointing along x, y, and z. That's pretty artificial. Another option might be to
	//  specify stubs of these fixed and moving residues?
	//
	////////////////////////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////
	//constructor!
	RigidBodySampler::RigidBodySampler( utility::vector1< Size > const & fixed_res,
																			utility::vector1< Size > const & moving_res ):
		fixed_res_( fixed_res ),
		moving_res_( moving_res ),
		reference_axes_( Matrix::identity() ),
		reference_centroid_( Vector( 0.0, 0.0, 0.0 ) ),
		scorefxn_( core::scoring::get_score_function() ),
		o2prime_trials_( false ),
		ignore_o2prime_hbonds_in_filter_( false ),
		assign_WC_edges_( false ),
		reference_energy_( 0.0 ),
		min_hbonds_( 0 ),
		fa_rep_cutoff_( 0.0 ),
		CONTACT_CUTOFF_squared_( 4.5 * 4.5 ),
		STERIC_DIST_CUTOFF_squared_( 2.5 * 2.5 ),
		MIN_NUM_CONTACTS_( 1 ),
		rmsd_cutoff_( -1.0 )
	{
		initialize_counters();
	}

RigidBodySampler::~RigidBodySampler() {}

	////////////////////////////////////////////////////////////
	void
	RigidBodySampler::initialize_counters(){

		all_rigid_body_settings_save_.clear();

		best_energy_ = reference_energy_;

		count_total_ = 0;
		count_good_ = 0;
		count_no_contact_ = 0;
		count_clash_ = 0;
	}

	////////////////////////////////////////////////////////////
	// This is the only RNA centric thing, I think...
	// also this basically doesn't work for what I want anyway.
	////////////////////////////////////////////////////////////
	void
	RigidBodySampler::initialize_reference_axes_and_centroid( conformation::Residue const & rsd ){
		using namespace scoring::rna;
		using namespace kinematics;
		static RNA_CentroidInfo rna_centroid_info;

		reference_centroid_ = rna_centroid_info.get_base_centroid( rsd );
		Stub s = rna_centroid_info.get_base_coordinate_system( rsd, reference_centroid_ );
		reference_axes_ = s.M;

		//std::cout << "REFERENCE_CENTROID" << reference_centroid_(1) << ' ' << reference_centroid_(1) << ' ' << reference_centroid_(3) <<  std::endl;
	}

	////////////////////////////////////////////////////////////
	void
	RigidBodySampler::do_the_sampling( pose::Pose & pose ){

		//clock_t const time_start( clock() );

		pose::Pose pose_start = pose;

		figure_out_reference_energy( pose );

		search_rotations_and_translations( pose );

	}

	/////////////////////////////////////////////////////////////
	// Sample Euler angles.
	////////////////////////////////////////////////////////////
	void
	RigidBodySampler::search_rotations_and_translations( pose::Pose & pose ){

		pose::Pose pose_start = pose;
		Matrix M;
		Vector const & axis1 = reference_axes_.col_x();
		Vector const & axis2 = reference_axes_.col_y();
		Vector const & axis3 = reference_axes_.col_z();

		Size const N_SAMPLE_ALPHA( static_cast<Size>( ( alpha_max_ - alpha_min_) / alpha_increment_ + 1 ) );
		Size i( 1 );

		for ( alpha_ = alpha_min_; alpha_ <= alpha_max_;  alpha_ += alpha_increment_ ){

			std::cout << i++ << " out of " << N_SAMPLE_ALPHA << ". Current count: " << count_total_ <<
				". num poses that pass cuts: " << count_good_ << std::endl;

			for ( Real cosbeta = cosbeta_min_; cosbeta <= cosbeta_max_;  cosbeta += cosbeta_increment_ ){
				if ( cosbeta < -1.0 ){
					beta_ = -1.0 * degrees( std::acos( -2.0 - cosbeta ) );
				} else if ( cosbeta > 1.0 ){
					beta_ = -1.0 * degrees( std::acos( 2.0 - cosbeta ) );
				} else {
					beta_ = degrees( std::acos( cosbeta ) );
				}

				std::cout << "BETA: " << beta_ << std::endl;

				// Try to avoid singularity at pole.
				Real gamma_min_local = gamma_min_;
				Real gamma_max_local = gamma_max_;
				Real gamma_increment_local = gamma_increment_;
				if ( (beta_<-179.999 || beta_>179.999) ){
					gamma_min_local = 0.0;
					gamma_max_local = 0.0;
					gamma_increment_local = 1.0;
				}

				for ( gamma_ = gamma_min_local; gamma_ <= gamma_max_local;  gamma_ += gamma_increment_local ){

					create_euler_rotation( M, alpha_, beta_, gamma_, axis1, axis2, axis3 );

					rotate( pose, M, pose_start, moving_res_, reference_centroid_ );

					pose::Pose pose_to_translate = pose;

					search_translations( pose,
															 pose_to_translate );

				} // gamma
			} // beta
		}// alpha

	}

	/////////////////////////////////////////////////////////////
	void
	virtualize_o2prime( pose::Pose & pose ){
		using namespace core::chemical;
		for( Size i = 1; i <= pose.total_residue(); i++ ){
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2PRIME_HYDROGEN", i );
		}
	}

	/////////////////////////////////////////////////////////////
	// xyz
	////////////////////////////////////////////////////////////
	void
	RigidBodySampler::search_translations( pose::Pose & pose,
																				 pose::Pose const & pose_to_translate )
	{
		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace protocols::stepwise;
		using namespace pose;

		//Size const N_SAMPLE_TRANSLATE  = 2 * static_cast< Size >( (box_size_ / xyz_increment_) + 0.5 ) + 1;
		//Real const beta = 1.0/ temperature_;
		//	bool const positive_Z =  option[ only_positive_Z ]();

		utility::vector1< Vector > moving_atoms, fixed_atoms;

		setup_heavy_atoms( pose_to_translate, moving_atoms, moving_res_ );
		setup_heavy_atoms( pose_to_translate, fixed_atoms, fixed_res_  );

		// virual o2prime pose -- needed for fa_rep checks.
		Pose pose_to_translate_virtual_o2prime = pose_to_translate;
		Pose pose_virtual_o2prime = pose;
		virtualize_o2prime( pose_to_translate_virtual_o2prime );
		virtualize_o2prime( pose_virtual_o2prime );

		//		std::cout << "X " << x_min_ << ' ' << x_max_ << ' ' << x_increment_ << std::endl;
		//		std::cout << "Y " << y_min_ << ' ' << y_max_ << ' ' << y_increment_ << std::endl;
		//		std::cout << "Z " << z_min_ << ' ' << z_max_ << ' ' << z_increment_ << std::endl;

		for ( delx_ = x_min_; delx_ <= x_max_; delx_ += x_increment_ ){
			for ( dely_ = y_min_; dely_ <= y_max_; dely_ += y_increment_ ){
				for ( delz_ = z_min_; delz_ <= z_max_; delz_ += z_increment_ ){

					Vector translation( delx_, dely_, delz_ );

					count_total_++;
					std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count_total_, 6 );
					// These could be sped up using grid indexing.
					if ( !check_contact(        translation, moving_atoms, fixed_atoms ) ) {
						count_no_contact_++;
						continue;
					}

					if ( !check_steric_overlap( translation, moving_atoms, fixed_atoms ) ) {
						count_clash_++;
						continue;
					}

					////////////////////////////////////////////////////////////////////

					if ( fa_rep_cutoff_ > 0.0 ){
						translate( pose_virtual_o2prime, translation, pose_to_translate_virtual_o2prime, moving_res_ );
						if ( !check_fa_rep( pose_virtual_o2prime ) ) continue;
					}

					translate( pose, translation, pose_to_translate, moving_res_ );

					if ( ignore_o2prime_hbonds_in_filter_ && min_hbonds_ > 0 && !check_num_hbonds( pose ) ) continue;
					if ( o2prime_trials_ && check_o2prime_needs_optimization( pose ) ) protocols::stepwise::sampling::rna::o2prime_trials( pose, o2prime_pack_scorefxn_ );
					if ( min_hbonds_ > 0 && !check_num_hbonds( pose ) ) continue;

					count_good_++;

					Real const energy = (*scorefxn_)( pose );
					save_rigid_body_settings( energy );

					////////////////////////////////////////////////////////////////////
					if ( energy < (best_energy_ + score_cutoff_) && sfd_ ){
						save_silent_struct( pose, tag );
					}
					if ( energy < best_energy_ ) best_energy_ = energy;


				}
			}
		}

	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::save_silent_struct( pose::Pose & pose, std::string const tag ){
		using namespace core::io::silent;
		using namespace core::scoring;

		(*scorefxn_)( pose );

		BinarySilentStruct s( pose, tag ); // this is RNA-centric -- could make it OK for proteins.

		if ( assign_WC_edges_ ) assign_WC_edges_to_base_pair12( pose, s );

		s.add_energy( "alphaRB", alpha_ );
		s.add_energy( "betaRB" , beta_ );
		s.add_energy( "gammaRB", gamma_ );
		s.add_energy( "x", delx_ );
		s.add_energy( "y", dely_ );
		s.add_energy( "z", delz_ );
		s.add_energy( "log_vol", log( x_increment_ * y_increment_ * z_increment_ * radians(alpha_increment_) * cosbeta_increment_ * radians(gamma_increment_) ) );

		if (native_pose_) {
			Real const rmsd = all_atom_rmsd( pose, *native_pose_ );
			if ( (rmsd_cutoff_ > 0.0) && (rmsd > rmsd_cutoff_) ) return;
			s.add_energy( "all_rms", rmsd );
		}

		sfd_->add_structure( s );
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::assign_WC_edges_to_base_pair12( pose::Pose & pose, core::io::silent::SilentStruct & s ){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace protocols::farna;

		utility::vector1< core::scoring::rna::Base_pair > base_pair_list;
		utility::vector1< bool > is_bulged;

		classify_base_pairs( pose, base_pair_list, is_bulged );

		Size edge1( 0 ), edge2( 0 );
		for ( Size n = 1; n <= base_pair_list.size(); n++ ) {
			core::scoring::rna::Base_pair const base_pair = base_pair_list[ n ];

			if( base_pair.res1 == 1 && base_pair.res2 == 2 ){
				edge1 = base_pair.edge1; edge2 = base_pair.edge2;
			}
			if( base_pair.res2 == 1 && base_pair.res1 == 2 ){
				edge1 = base_pair.edge2; edge2 = base_pair.edge1;
				}
		}

		s.add_energy( "edge1", edge1 );
		s.add_energy( "edge2", edge2 );

	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::apply_input_samples(
				pose::Pose & pose,
				std::string const rigid_body_sample_file ){

		using namespace core::io::silent;
		using namespace utility::io;
		izstream input( rigid_body_sample_file );

		if ( !input ) {
			std::cerr << "No file: " << rigid_body_sample_file << std::endl;
			utility_exit_with_message( "No file" );
		}

		pose::Pose pose_start = pose;

		Real alpha, beta, gamma, x, y, z;
		count_total_ = 0;
		while ( input >> alpha ) {

			input >> beta >> gamma >> x >> y >> z >> skip;
			apply_rigid_body_settings( pose, pose_start,
																 alpha, beta, gamma, x, y, z );

			// Real const energy = // Unused variable causes warning.
			(*scorefxn_)( pose );

			count_total_++;
			std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count_total_, 6 );

			BinarySilentStruct s( pose, tag ); // this is RNA-centric -- could make it OK for proteins.
			sfd_->add_structure( s );

		}

		pose = pose_start;
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::apply_rigid_body_settings( pose::Pose & pose, pose::Pose const & pose_start,
																							 Real const alpha,
																							 Real const beta,
																							 Real const gamma,
																							 Real const x,
																							 Real const y,
																							 Real const z ){

		static Matrix M;
		Vector const & axis1 = reference_axes_.col_x();
		Vector const & axis2 = reference_axes_.col_y();
		Vector const & axis3 = reference_axes_.col_z();

		create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );
		rotate( pose, M, pose_start, moving_res_, reference_centroid_ );
		translate( pose, Vector( x,y,z), pose, moving_res_ );

		// This is useful for saving silent structures...
		alpha_ = alpha;
		beta_  = beta;
		gamma_ = gamma;
		delx_  = x;
		dely_  = y;
		delz_  = z;
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::setup_heavy_atoms( pose::Pose const & pose,
																			 utility::vector1< Vector > & pose_atoms,
																			 utility::vector1< Size > const & subset_res ){

		pose_atoms.clear();

		for ( Size n = 1; n <= subset_res.size(); n++ ) {
			Size const i = subset_res[ n ];

			for ( Size j = 1; j <= pose.residue_type( i ).nheavyatoms(); j++ ){
				if ( pose.residue(i).is_virtual( j ) ) continue;
				pose_atoms.push_back( pose.xyz( core::id::AtomID(j,i) ) );
			}

		}

	}


	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_contact_cutoff( Real const setting ){
		Real const CONTACT_CUTOFF_ = setting;
		CONTACT_CUTOFF_squared_ = CONTACT_CUTOFF_ * CONTACT_CUTOFF_;
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_min_num_contacts( Size const setting ){
		MIN_NUM_CONTACTS_ = setting;
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_silent_file_data(  core::io::silent::SilentFileDataOP sfd ){
		sfd_ = sfd;
	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_o2prime_trials(  bool const setting ){

		using namespace core::scoring;

		o2prime_trials_ = setting;

		if ( o2prime_trials_ ){
			o2prime_pack_scorefxn_ = new ScoreFunction;
			o2prime_pack_scorefxn_->set_weight( fa_atr, scorefxn_->get_weight( fa_atr ) );
			o2prime_pack_scorefxn_->set_weight( fa_rep, scorefxn_->get_weight( fa_rep ) );
			o2prime_pack_scorefxn_->set_weight( hbond_lr_bb_sc, scorefxn_->get_weight( hbond_lr_bb_sc ) );
			o2prime_pack_scorefxn_->set_weight( hbond_sr_bb_sc, scorefxn_->get_weight( hbond_sr_bb_sc ) );
			o2prime_pack_scorefxn_->set_weight( hbond_sc, scorefxn_->get_weight( hbond_sc ) );
			o2prime_pack_scorefxn_->set_energy_method_options( scorefxn_->energy_method_options() );
			// note that geom_sol is not optimized well --> replace with lk_sol for now.
			o2prime_pack_scorefxn_->set_weight( fa_sol, scorefxn_->get_weight( lk_nonpolar ) );
		}
	}

	///////////////////////////////////////////////////////////////////////
	bool
	RigidBodySampler::check_contact( Vector const & translation,
																	 utility::vector1< Vector > const & moving_atoms,
																	 utility::vector1< Vector > const & partner_atoms
																	 ){

		Size num_contacts( 0 );

		for ( Size i = 1; i <= moving_atoms.size(); i++ ) {
			Vector const test_cbeta = moving_atoms[ i ] + translation;

			for ( Size j = 1; j <= partner_atoms.size(); j++ ) {
				Vector const & partner_cbeta = partner_atoms[ j ];

				if ( ( test_cbeta - partner_cbeta ).length_squared() < CONTACT_CUTOFF_squared_ ) {
					num_contacts ++;
					if ( num_contacts >= MIN_NUM_CONTACTS_ ) return true;
				}

			}
		}

		return false;

	}

	///////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_steric_dist_cutoff( Real const setting ){
		Real const STERIC_DIST_CUTOFF_ = setting;
		STERIC_DIST_CUTOFF_squared_ = STERIC_DIST_CUTOFF_ * STERIC_DIST_CUTOFF_;
	}

	///////////////////////////////////////////////////////////////////////
	bool
	RigidBodySampler::check_steric_overlap( Vector const & translation,
																					utility::vector1< Vector > const & moving_atoms,
																					utility::vector1< Vector > const & partner_atoms
																					){

		for ( Size i = 1; i <= moving_atoms.size(); i++ ) {
			Vector const test_atom = moving_atoms[ i ] + translation;

			for ( Size j = 1; j <= partner_atoms.size(); j++ ) {
				Vector const & partner_atom = partner_atoms[ j ];
				if ( ( test_atom - partner_atom ).length_squared() < STERIC_DIST_CUTOFF_squared_ ) return false;
			}
		}

		return true;

	}

	/////////////////////////////////////////////////////////////
	bool in_vector( Size const & i, utility::vector1< Size > const & vec ){
		for ( Size n = 1; n <= vec.size(); n++ ) {
			if ( i == vec[n] ) return true;
		}
		return false;
	}

	/////////////////////////////////////////////////////////////
	bool
	RigidBodySampler::check_num_hbonds( pose::Pose & pose ){

		using namespace core::scoring;

		// figure out H-bonds.
		hbonds::HBondOptionsOP hbond_options( new hbonds::HBondOptions() );
		hbond_options->use_hb_env_dep( false );
		hbonds::HBondSetOP hbond_set( new hbonds::HBondSet( hbond_options ) );
		hbonds::fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

		Size num_cross_hbonds( 0 );
		Real const HBOND_ENERGY_CUTOFF( -0.5 );

		for (Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
			hbonds::HBond const & hbond( hbond_set->hbond( i ) );

			Size const don_res_num = hbond.don_res();
			Size const don_hatm = hbond.don_hatm();
			Size const acc_res_num = hbond.acc_res();
			Size const acc_atm = hbond.acc_atm();

			if ( hbond.energy() > HBOND_ENERGY_CUTOFF ) continue;

			// could be sped up by indexing 'in_moving_res' and 'in_fixed_res'
			if ( !( ( in_vector( don_res_num, moving_res_ ) && in_vector( acc_res_num, fixed_res_ ) ) ||
							( in_vector( acc_res_num, moving_res_ ) && in_vector( don_res_num, fixed_res_ ) ) ) ) continue;

			if ( ignore_o2prime_hbonds_in_filter_ && ( pose.residue_type( don_res_num ).atom_name( don_hatm ) == "HO2'" ) ) continue;
			if ( ignore_o2prime_hbonds_in_filter_ && ( pose.residue_type( acc_res_num ).atom_name( acc_atm  ) == " O2'" ) ) continue;

			num_cross_hbonds++;

			if ( num_cross_hbonds >= min_hbonds_ ) return true;

		}

		return false;

	}

	/////////////////////////////////////////////////////////////
	bool
	RigidBodySampler::check_o2prime_needs_optimization( pose::Pose const & pose ){
		Real const DIST_CUTOFF = 4.0;

		for ( Size n = 1; n <= moving_res_.size(); n++ ){
			Vector const & o2prime_xyz = pose.residue( moving_res_[n] ).xyz( " O2'" );
			for ( Size m = 1; m <= fixed_res_.size(); m++ ){
				for ( Size k = 1; k <= pose.residue_type( fixed_res_[m] ).nheavyatoms(); k++ ) {
					if ( ( o2prime_xyz - pose.residue( fixed_res_[m] ).xyz( k ) ).length() < DIST_CUTOFF )	return true;
				}
			}
		}

		for ( Size n = 1; n <= fixed_res_.size(); n++ ){
			Vector const & o2prime_xyz = pose.residue( fixed_res_[n] ).xyz( " O2'" );
			for ( Size m = 1; m <= moving_res_.size(); m++ ){
				for ( Size k = 1; k <= pose.residue_type( moving_res_[m] ).nheavyatoms(); k++ ) {
					if ( ( o2prime_xyz - pose.residue( moving_res_[m]).xyz( k ) ).length() < DIST_CUTOFF )	return true;
				}
			}
		}

		return false;
	}

	/////////////////////////////////////////////////////////////
	bool
	RigidBodySampler::check_fa_rep( pose::Pose & pose ){

		using namespace core::scoring;

		static bool init( false );
		static ScoreFunctionOP rep_scorefxn( new ScoreFunction);
		if (!init ){
			rep_scorefxn->set_weight( fa_rep, 1.0 );
			init = true;
		}

		Real const fa_rep_score = (*rep_scorefxn)( pose );
		return ( fa_rep_score <= fa_rep_cutoff_ );

	}


	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::save_rigid_body_settings( Real const energy ){

		utility::vector1< Real > rigid_body_settings;
		rigid_body_settings.push_back( alpha_ );
		rigid_body_settings.push_back( beta_ );
		rigid_body_settings.push_back( gamma_ );
		rigid_body_settings.push_back( delx_ );
		rigid_body_settings.push_back( dely_ );
		rigid_body_settings.push_back( delz_ );
		rigid_body_settings.push_back( energy );

		all_rigid_body_settings_save_.push_back( rigid_body_settings );

	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_n_sample_alpha_full_range( Size const n_sample ){
		alpha_increment_ = ( 360.0 / static_cast< Real >( n_sample ) );
		alpha_min_ = alpha_increment_;
		alpha_max_ = n_sample * alpha_increment_;
	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_n_sample_gamma_full_range( Size const n_sample ){
		gamma_increment_ = ( 360.0 / static_cast< Real >( n_sample ) );
		gamma_min_ = gamma_increment_;
		gamma_max_ = n_sample * gamma_increment_;
	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_n_sample_cosbeta_full_range( Size const n_sample ){
		cosbeta_increment_ = 2.0 / n_sample;

		//cosbeta_min_ = -1.0 + 0.5 * cosbeta_increment_;
		//cosbeta_max_ = -1.0 + ( n_sample - 0.5 ) * cosbeta_increment_;

		// Put in a condition to not sample gamma if we're at north or south pole
		cosbeta_min_ = -1.0;
		cosbeta_max_ =  1.0;

	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::set_translation_sample(  Real const box_size, Real const xyz_increment ){

		x_increment_ = xyz_increment;
		y_increment_ = xyz_increment;
		z_increment_ = xyz_increment;

		x_min_ = -box_size;
		x_max_ = box_size;

		y_min_ = -box_size;
		y_max_ = box_size;

		z_min_ = -box_size;
		z_max_ = box_size;

	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::force_coplanar(){
		z_min_ = 0.0;
		z_max_ = 0.0;
		z_increment_ = 0.1; //will only sample z = 0.0.
	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::force_antiparallel(){
		cosbeta_min_ = -1.0;
		cosbeta_max_ = -1.0;
		cosbeta_increment_ = 0.1;
	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::force_parallel(){
		cosbeta_min_ = 1.0;
		cosbeta_max_ = 1.0;
		cosbeta_increment_ = 0.1;
	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::figure_out_reference_energy( pose::Pose & pose ){

		pose::Pose pose_perturb = pose;
		translate( pose_perturb, Vector( 1000.0, 1000.0, 1000.0 ), pose, moving_res_ );
		reference_energy_ = (*scorefxn_)( pose );

	}

	/////////////////////////////////////////////////////////////
	void
	RigidBodySampler::output_results( utility::io::ozstream & out ){

		Real const volume_element = x_increment_ * y_increment_ * z_increment_ * radians(alpha_increment_) * cosbeta_increment_ * radians(gamma_increment_);

		out << "0.0 0.0 0.0 0.0 0.0 0.0 " << ' ' << 0.0 << ' ' << count_no_contact_ * volume_element << std::endl;
		out << "0.0 0.0 0.0 0.0 0.0 0.0 " << ' ' << 999.99 << ' ' << count_clash_ * volume_element << std::endl;

		for ( Size i = 1; i <= all_rigid_body_settings_save_.size(); i++ ) {

			for ( Size n = 1; n <= all_rigid_body_settings_save_[ i ].size(); n++ ){
				out << ' ' << all_rigid_body_settings_save_[ i ][ n ];
			}
			out << ' ' << volume_element << std::endl;

		}

	}

	/////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::output_silent_file( std::string const silent_file, bool const write_score_only ){

		using namespace core::io::silent;

		SilentFileData sfd;

		SilentFileData::iterator iter = sfd_->begin();
		SilentStructOP s = *iter;
		sfd.write_silent_struct( *s, silent_file, write_score_only );

		utility::io::ozstream out;
		out.open_append( silent_file );


		for ( ++iter;  iter != sfd_->end(); ++iter ){
			SilentStructOP s = *iter;
			sfd.write_silent_struct( *s, out, write_score_only );
		}

	}


	/////////////////////////////////////////////////////////////////////////////////////////////
	void
	RigidBodySampler::output_histogram( utility::io::ozstream & out ){

		Real const volume_element = x_increment_ * y_increment_ * z_increment_ * radians(alpha_increment_) * cosbeta_increment_ * radians(gamma_increment_);

		// need to define histogram bins and energies.
		Real const energy_hist_min_ = -5.0;
		Real const energy_hist_max_ = 10.0;
		Real const energy_increment_ = 0.01;
		utility::vector1< Real > histogram_energy, volumes;
		for ( Real energy = energy_hist_min_; energy < energy_hist_max_; energy += energy_increment_ ){
			histogram_energy.push_back( energy );
			volumes.push_back( 0.0 );
		}

		/// Go through data. assign as either clash, special low energy, or in the histogram
		std::map< Size, utility::vector1< Real > > histogram_rigid_body_settings;
		utility::vector1< utility::vector1< Real > > good_energy_rigid_body_settings;
		Size count_clash_local = count_clash_;
		for ( Size i = 1; i <= all_rigid_body_settings_save_.size(); i++ ) {

			utility::vector1< Real > const & rigid_body_setting = all_rigid_body_settings_save_[ i ];

			Real const energy =  rigid_body_setting[ 7 ];

			if ( energy < energy_hist_min_ ) {
				good_energy_rigid_body_settings.push_back( rigid_body_setting );
			} else if ( energy > energy_hist_max_ ){
				count_clash_local += 1;
			} else {
				Size const bin = static_cast< Size >( ( energy - energy_hist_min_) / energy_increment_ ) + 1;
				volumes[ bin ] += volume_element;
				histogram_rigid_body_settings[ bin ] = rigid_body_setting;
			}

		}

		////////////////////////////////////////////////////////
		// output
		////////////////////////////////////////////////////////
		out << "0.0 0.0 0.0 0.0 0.0 0.0 " << ' ' << 0.00 << ' ' << count_no_contact_ * volume_element << std::endl;
		out << "0.0 0.0 0.0 0.0 0.0 0.0 " << ' ' << 999.99 << ' ' << count_clash_local * volume_element << std::endl;

		for ( Size i = 1; i <= histogram_energy.size(); i++ ){
			if ( volumes[ i ] > 0.0 ){
				utility::vector1< Real > const & rigid_body_setting = histogram_rigid_body_settings[ i ];

				for ( Size n = 1; n <= 6; n++ ){
					out << ' ' << rigid_body_setting[ n ];
				}
				out << ' ' << histogram_energy[ i ];
				out << ' ' << volumes[ i ] << std::endl;
			}
		}


		for ( Size i = 1; i <= good_energy_rigid_body_settings.size(); i++ ){
			utility::vector1< Real > const & rigid_body_setting = good_energy_rigid_body_settings[ i ];

			for ( Size n = 1; n <= rigid_body_setting.size(); n++ ){
				out << ' ' << rigid_body_setting[ n ];
			}
			out << ' ' << volume_element << std::endl;
		}

	}

	////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionOP  RigidBodySampler::score_function(){ return scorefxn_; }



} //sampling
} //legacy
} //stepwise
} //protocols
