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
/// @author

// type headers
#include <core/types.hh>


// unit headers
#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>
#include <utility/excn/Exceptions.hh>

// package headers
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/jd2/util.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
// #include "utility/sys_util.h"
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// C++ headers
#include <map>
#include <string>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <basic/prof.hh>

static basic::Tracer tr( "protocols.canonical_sampling.mc_convergence_check.Pool", basic::t_info );

// Forward declarations

namespace protocols {
namespace canonical_sampling {
namespace mc_convergence_checks {

// @brief Auto-generated virtual destructor
Pool_Evaluator::~Pool_Evaluator() = default;

// @brief Auto-generated virtual destructor
Pool_ConvergenceCheck::~Pool_ConvergenceCheck() = default;

// @brief Auto-generated virtual destructor
Pool_RMSD::~Pool_RMSD() = default;

using namespace ObjexxFCL;
using namespace core;
using namespace basic;

void Pool_RMSD::get( core::Size index , FArray2D<double>& result){
	tr.Debug << "getting " << index << " coords of " << pool_.n_decoys() << std::endl;
	FArray2P_double temp = pool_.coords( index );
	result.dimension(temp.u1(),temp.u2(),0.0);
	for ( int i = 1; i <= result.u1(); i++ ) {
		for ( int j = 1; j <= result.u2(); j++ ) {
			result( i, j ) = temp( i, j );
		}
	}
}

std::string& Pool_RMSD::get_tag( core::Size index ){
	tr.Debug << "getting tag " << index << " of " << tags_.size() << " which has value " << tags_[index] << std::endl;
	return tags_[ index ];
}

void Pool_RMSD::pop_back(){
	tags_.pop_back( );
	tr.Debug << "checking size before popping back " << pool_.n_decoys() << std::endl;
	pool_.pop_back_CA_xyz();
	pool_.center_structure( pool_.n_decoys(), weights_ );
	tr.Debug << "finished popping back. size is now " << pool_.n_decoys() << std::endl;
	tr.Debug << "just checking that size of arrays match up " << tags_.size() << " " << pool_.n_decoys() << std::endl;

	runtime_assert( tags_.size() == pool_.n_decoys() );

}

void Pool_RMSD::clear(){
	tags_.clear();
	pool_.clear();
}


void Pool_RMSD::set_excluded_residues( utility::vector1< core::Size > const& excluded_residues ) {
	excluded_residues_ = excluded_residues;
}

void Pool_RMSD::fill_pool( std::string const& silent_file ) {
	//tags_.clear();
	//pool_.clear();
	clear();

	io::silent::SilentFileOptions opts;
	io::silent::SilentFileData sfd( opts );
	try {
		sfd.read_file( silent_file );
	} catch( utility::excn::BadInput const& excn ) {
		excn.show( tr.Warning );
		tr.Warning << "Pool_RMSD did not find any structures, will output n/a 10000 for converged_tag and converged_rmsd" << std::endl;
		return;
	}
	pool_.push_back_CA_xyz_from_silent_file( sfd, false /*don't save energies*/ );
	tags_.reserve( sfd.size() );
	for ( io::silent::SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
		tags_.push_back( it->decoy_tag() );
	}
	weights_.dimension( pool_.n_atoms(), 1.0 );
	if ( pool_.n_decoys() ) pool_.center_all( weights_ ); //this centers all structures
	reserve_size_ = 100;
}


void Pool_RMSD::add( core::io::silent::SilentStruct const& pss, std::string new_tag ) {
	if ( new_tag == "keep_silent_tag" ) {
		new_tag = pss.decoy_tag();
	}
	add( pss.get_CA_xyz(), pss.nres(), new_tag );
}

void Pool_RMSD::add( core::pose::Pose const& pose, std::string new_tag ) {
	FArray2D_double coords( 3, pool_.n_atoms(), 0.0 );
	runtime_assert( pool_.n_atoms() == pose.size() );
	toolbox::fill_CA_coords( pose, pool_.n_atoms(), coords );
	add( coords, pose.size(), new_tag );
}

void Pool_RMSD::add( ObjexxFCL::FArray2D_double const& xyz, core::Size nres, std::string new_tag ) {

	if ( weights_.size() == 0 ) {
		weights_.dimension( nres, 1.0 );
		reserve_size_ = 100;
	}

	PROF_START( basic::POOL_RMSD_ADD_STRUCTURE );
	//runtime_assert( pool_.n_atoms() == nres );
	if ( pool_.n_decoys_max() >= pool_.n_decoys() ) {
		//pool_.reserve( pool_.n_decoys() + 100 );
		pool_.reserve( pool_.n_decoys() + reserve_size_ ); //ek
		//  pool_.reserve( pool_.n_decoys() + reserve_size_ ); //ek
		//  tr.Debug << "pool can now hold a max of " << pool_.n_decoys() + reserve_size_;
	}
	tags_.push_back( new_tag );
	pool_.push_back_CA_xyz( xyz, nres );
	pool_.center_structure( pool_.n_decoys(), weights_ );

	runtime_assert( tags_.size() == pool_.n_decoys() );
	PROF_STOP( basic::POOL_RMSD_ADD_STRUCTURE );
}


void Pool_RMSD::set_reserve_size( int max_size ){ //ek
	reserve_size_ = max_size;
}


core::Size Pool_RMSD::evaluate_and_add(pose::Pose const& pose, std::string& cluster_center, core::Real& rms_to_cluster, core::Real transition_threshold) {
	tr.Debug << "using Pool_RMSD::evaluate_and_add " << std::endl;
	ObjexxFCL::FArray2D_double coords;
	coords.redimension( 3, pose.size() );
	toolbox::fill_CA_coords( pose, pool_.n_atoms(), coords );
	return evaluate_and_add( coords, cluster_center, rms_to_cluster, transition_threshold );
}

core::Size Pool_RMSD::evaluate_and_add( ObjexxFCL::FArray2D_double& coords, std::string& best_decoy, core::Real& best_rmsd, core::Real transition_threshold  ){
	core::Size best_index = evaluate( coords, best_decoy, best_rmsd );
	if ( best_rmsd > transition_threshold ) {
		tr.Debug  << "adding structure to pool " << std::endl;
		using namespace ObjexxFCL;
		//form tags of type: c.<cluster_number>.<decoy_in_group_nr>
		std::string jobname = protocols::jd2::current_output_name();
		std::string new_cluster_tag = "new."+lead_zero_string_of( size(), 8 )+".0"+"_"+jobname;
		add( coords, coords.u2(), new_cluster_tag);
	}
	return best_index;
}

core::Size Pool_RMSD::evaluate( core::io::silent::SilentStruct const& pss, std::string& best_decoy, core::Real& best_rmsd ) const {
	if ( pool_.n_decoys() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}
	FArray2D_double coords( pss.get_CA_xyz() );
	runtime_assert( pss.nres() == pool_.n_atoms() );
	return evaluate( coords, best_decoy, best_rmsd );
}

core::Size Pool_RMSD::evaluate( core::pose::Pose const& fit_pose, std::string& best_decoy, core::Real& best_rmsd ) const {
	if ( pool_.n_decoys() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}
	FArray2D_double coords( 3, pool_.n_atoms(), 0.0 );
	runtime_assert( pool_.n_atoms() <= fit_pose.size() );
	toolbox::fill_CA_coords( fit_pose, pool_.n_atoms(), coords );
	return evaluate( coords, best_decoy, best_rmsd );
}

core::Size Pool_RMSD::evaluate( FArray2D_double& coords, std::string& best_decoy, core::Real& best_rmsd ) const {
	return evaluate( coords, best_decoy, best_rmsd, 1);
	//return 0; //ek commented out, i need the index of the best-decoy!!
}

//index indicates what index you start evaluation from
core::Size Pool_RMSD::evaluate( FArray2D_double& coords, std::string& best_decoy, core::Real& best_rmsd , core::Size index ) const {
	PROF_START( basic::POOL_RMSD_EVALUATE );
	if ( pool_.n_decoys() == 0 ) {
		best_decoy = "n/a";
		best_rmsd = 10000;
		return 0;
	}
	tr.Debug << "evaluating: start from index " << index << " of " << pool_.n_decoys() << std::endl;

	FArray1D_double weights( pool_.n_atoms(), 1.0 );
	for ( unsigned long excluded_residue : excluded_residues_ ) {
		weights( excluded_residue ) = 0.0;
	}

	FArray1D_double transvec( 3 );
	toolbox::reset_x( pool_.n_atoms(), coords, weights, transvec );//center coordinates
	//n_atoms is simply # CA atoms in each "pose"
	Real invn( 1.0 / pool_.n_atoms() );
	best_rmsd = 1000000;
	Size best_index( 0 );
	for ( Size i = index; i <= pool_.n_decoys(); i++ ) {
		toolbox::Matrix R;
		ObjexxFCL::FArray2P_double xx2( pool_.coords( i ) );
		toolbox::fit_centered_coords( pool_.n_atoms(), weights, xx2, coords, R );
		Real rmsd( 0 );
		for ( Size n = 1; n <= pool_.n_atoms(); n++ ) {
			for ( Size d = 1; d<=3; ++d ) {
				rmsd += ( coords( d, n ) -  xx2( d, n ) ) * ( coords( d, n ) - xx2( d, n ) ) * invn;
			}
		}
		rmsd = sqrt( rmsd );
		if ( rmsd <= best_rmsd ) {
			best_index = i;
			best_rmsd = rmsd;
		}
	}


	best_decoy = tags_[ best_index ];
	PROF_STOP( basic::POOL_RMSD_EVALUATE );
	return best_index;
}

bool Pool_ConvergenceCheck::operator() ( core::pose::Pose const & fit_pose, moves::MonteCarlo const&, bool reject ) {
	if ( reject ) return true; //dont bother for these

	core::Real best_rmsd;
	std::string best_decoy;

	rmsd_pool_->evaluate( fit_pose, best_decoy, best_rmsd );

	//store in Job-Object:
	protocols::jd2::add_string_string_pair_to_current_job( "pool_converged_tag", best_decoy );
	protocols::jd2::add_string_real_pair_to_current_job( "pool_converged_rmsd", best_rmsd );
	if ( best_rmsd <= threshold_ ) throw CREATE_EXCEPTION(EXCN_Pool_Converged, "");
	return best_rmsd >= threshold_;
}

/// @brief evaluate pose and store values in Silent_Struct
void Pool_Evaluator::apply( core::pose::Pose& fit_pose, std::string, core::io::silent::SilentStruct &pss) const {
	core::Real best_rmsd;
	std::string best_decoy;

	rmsd_pool_->evaluate( fit_pose, best_decoy, best_rmsd );

	pss.add_string_value( name( 1 ), best_decoy );
	pss.add_energy( name( 2 ), best_rmsd );
	// //store in Job-Object:
	//  protocols::jd2::add_string_string_pair_to_current_job( "pool_converged_tag", best_decoy );
	//  protocols::jd2::add_string_real_pair_to_current_job( "pool_converged_rmsd", best_rmsd );

}

/// @brief evaluate pose and store values in Silent_Struct
void Pool_Evaluator::apply( core::io::silent::SilentStruct &pss) const {
	core::Real best_rmsd;
	std::string best_decoy;

	rmsd_pool_->evaluate( pss, best_decoy, best_rmsd );

	pss.add_string_value( name( 1 ), best_decoy );
	pss.add_energy( name( 2 ), best_rmsd );
	// //store in Job-Object:
	//  protocols::jd2::add_string_string_pair_to_current_job( "pool_converged_tag", best_decoy );
	//  protocols::jd2::add_string_real_pair_to_current_job( "pool_converged_rmsd", best_rmsd );

}

} // mc_convergence_check
} // moves
} // rosetta
