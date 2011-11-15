// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

// type headers
#include <core/types.hh>

//
// unit headers
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
#include <utility/excn/Exceptions.hh>

// package headers
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <core/io/silent/SilentFileData.hh>
#include <protocols/toolbox/superimpose.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>
// #include "utility/basic_sys_util.h"
#include <basic/Tracer.hh>
#include <utility/exit.hh>
// C++ headers
#include <map>
#include <string>

#include <platform/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/io/silent/SharedSilentData.hh>
#include <core/io/silent/SilentEnergy.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/jd2/InnerJob.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
#include <protocols/jd2/JobDistributor.fwd.hh>
#include <protocols/jd2/JobInputter.fwd.hh>
#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Parser.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.fwd.hh>
#include <protocols/moves/mc_convergence_checks/ConvergenceCheck.hh>
#include <protocols/toolbox/DecoySetEvaluation.fwd.hh>
#include <utility/Bound.fwd.hh>
#include <utility/Bound.hh>
#include <utility/down_cast.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.fwd.hh>
#include <utility/file/PathName.hh>
#include <utility/keys/AutoKey.fwd.hh>
#include <utility/keys/AutoKey.hh>
#include <utility/keys/Key.fwd.hh>
#include <utility/keys/Key.hh>
#include <utility/keys/KeyLess.fwd.hh>
#include <utility/keys/KeyLookup.fwd.hh>
#include <utility/keys/KeyLookup.hh>
#include <utility/keys/NoClient.fwd.hh>
#include <utility/keys/NoClient.hh>
#include <utility/keys/SmallKeyVector.fwd.hh>
#include <utility/keys/SmallKeyVector.hh>
#include <utility/keys/UserKey.fwd.hh>
#include <utility/keys/VariantKey.fwd.hh>
#include <utility/keys/VariantKey.hh>
#include <utility/options/AnyOption.fwd.hh>
#include <utility/options/AnyOption.hh>
#include <utility/options/AnyVectorOption.fwd.hh>
#include <utility/options/AnyVectorOption.hh>
#include <utility/options/BooleanOption.fwd.hh>
#include <utility/options/BooleanOption.hh>
#include <utility/options/BooleanVectorOption.fwd.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <utility/options/FileOption.fwd.hh>
#include <utility/options/FileOption.hh>
#include <utility/options/FileVectorOption.fwd.hh>
#include <utility/options/FileVectorOption.hh>
#include <utility/options/IntegerOption.fwd.hh>
#include <utility/options/IntegerOption.hh>
#include <utility/options/IntegerVectorOption.fwd.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/options/Option.fwd.hh>
#include <utility/options/Option.hh>
#include <utility/options/OptionCollection.fwd.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/options/PathOption.fwd.hh>
#include <utility/options/PathOption.hh>
#include <utility/options/PathVectorOption.fwd.hh>
#include <utility/options/PathVectorOption.hh>
#include <utility/options/RealOption.fwd.hh>
#include <utility/options/RealOption.hh>
#include <utility/options/RealVectorOption.fwd.hh>
#include <utility/options/RealVectorOption.hh>
#include <utility/options/ScalarOption.fwd.hh>
#include <utility/options/ScalarOption.hh>
#include <utility/options/ScalarOption_T_.fwd.hh>
#include <utility/options/ScalarOption_T_.hh>
#include <utility/options/StringOption.fwd.hh>
#include <utility/options/StringOption.hh>
#include <utility/options/StringVectorOption.fwd.hh>
#include <utility/options/StringVectorOption.hh>
#include <utility/options/VariantOption.fwd.hh>
#include <utility/options/VariantOption.hh>
#include <utility/options/VectorOption.fwd.hh>
#include <utility/options/VectorOption.hh>
#include <utility/options/VectorOption_T_.fwd.hh>
#include <utility/options/VectorOption_T_.hh>
#include <utility/options/mpi_stderr.hh>
#include <utility/options/keys/AnyOptionKey.fwd.hh>
#include <utility/options/keys/AnyOptionKey.hh>
#include <utility/options/keys/AnyVectorOptionKey.fwd.hh>
#include <utility/options/keys/AnyVectorOptionKey.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/BooleanVectorOptionKey.fwd.hh>
#include <utility/options/keys/BooleanVectorOptionKey.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.fwd.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/IntegerVectorOptionKey.fwd.hh>
#include <utility/options/keys/IntegerVectorOptionKey.hh>
#include <utility/options/keys/OptionKey.fwd.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/keys/OptionKeys.hh>
#include <utility/options/keys/PathOptionKey.fwd.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/PathVectorOptionKey.fwd.hh>
#include <utility/options/keys/PathVectorOptionKey.hh>
#include <utility/options/keys/RealOptionKey.fwd.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/RealVectorOptionKey.fwd.hh>
#include <utility/options/keys/RealVectorOptionKey.hh>
#include <utility/options/keys/ScalarOptionKey.fwd.hh>
#include <utility/options/keys/ScalarOptionKey.hh>
#include <utility/options/keys/StringOptionKey.fwd.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/options/keys/StringVectorOptionKey.fwd.hh>
#include <utility/options/keys/StringVectorOptionKey.hh>
#include <utility/options/keys/VectorOptionKey.fwd.hh>
#include <utility/options/keys/VectorOptionKey.hh>
#include <utility/options/keys/all.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2P.fwd.hh>
#include <ObjexxFCL/FArray2P.hh>
#include <ObjexxFCL/FArray3.fwd.hh>
#include <ObjexxFCL/FArray3.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/Fstring.fwd.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/TypeTraits.hh>
#include <ObjexxFCL/byte.fwd.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/ubyte.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <istream>
#include <limits>
#include <list>
#include <ostream>
#include <set>
#include <sstream>
#include <time.h>
#include <utility>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/prof.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>

static basic::Tracer tr("protocols.moves.mc_convergence_check.Pool",basic::t_info);

// Forward declarations

namespace protocols {
namespace moves {
namespace mc_convergence_checks {

using namespace ObjexxFCL;
using namespace core;
using namespace basic;

void Pool_RMSD::get( core::Size index , FArray2D<double>& result){
	 tr.Debug << "getting " << index << " coords of " << pool_.n_decoys() << std::endl;
	 FArray2P_double temp = pool_.coords( index );
	 result.dimension(temp.u1(),temp.u2(),0.0);
	 for(int i = 1; i <= result.u1(); i++){
		 for(int j = 1; j <= result.u2(); j++){
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

void Pool_RMSD::fill_pool( std::string const& silent_file ) {
	//tags_.clear();
	//pool_.clear();
	clear();

	io::silent::SilentFileData sfd;
	try {
		sfd.read_file( silent_file );
	} catch( utility::excn::EXCN_BadInput const& excn ) {
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
	runtime_assert( pool_.n_atoms() == pose.total_residue() );
	toolbox::fill_CA_coords( pose, pool_.n_atoms(), coords );
	add( coords, pose.total_residue(), new_tag );
}

void Pool_RMSD::add( ObjexxFCL::FArray2D_double const& xyz, core::Size nres, std::string new_tag ) {

	if( weights_.size() == 0 ){
		weights_.dimension( nres, 1.0 );
		reserve_size_ = 100;
	}

	PROF_START( basic::POOL_RMSD_ADD_STRUCTURE );
	//runtime_assert( pool_.n_atoms() == nres );
	if ( pool_.n_decoys_max() >= pool_.n_decoys() ) {
		//pool_.reserve( pool_.n_decoys() + 100 );
		pool_.reserve( pool_.n_decoys() + reserve_size_ ); //ek
		//		pool_.reserve( pool_.n_decoys() + reserve_size_ ); //ek
		//		tr.Debug << "pool can now hold a max of " << pool_.n_decoys() + reserve_size_;
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
	coords.redimension( 3, pose.total_residue() );
	toolbox::fill_CA_coords( pose, pool_.n_atoms(), coords );
	return evaluate_and_add( coords, cluster_center, rms_to_cluster, transition_threshold );
}

core::Size Pool_RMSD::evaluate_and_add( ObjexxFCL::FArray2D_double& coords, std::string& best_decoy, core::Real& best_rmsd, core::Real transition_threshold  ){
	core::Size best_index = evaluate( coords, best_decoy, best_rmsd );
	if( best_rmsd > transition_threshold ){
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
	runtime_assert( pool_.n_atoms() <= fit_pose.total_residue() );
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

	FArray1D_double const weights( pool_.n_atoms(), 1.0 );
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
		for ( Size n = 1; n <= pool_.n_atoms(); n++) {
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
	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair( "pool_converged_tag", best_decoy );
	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair( "pool_converged_rmsd", best_rmsd );
	if ( best_rmsd <= threshold_ ) throw EXCN_Pool_Converged();
	return best_rmsd >= threshold_;
}

///@brief evaluate pose and store values in Silent_Struct
void Pool_Evaluator::apply( core::pose::Pose& fit_pose, std::string, core::io::silent::SilentStruct &pss) const {
	core::Real best_rmsd;
	std::string best_decoy;

	rmsd_pool_->evaluate( fit_pose, best_decoy, best_rmsd );

	pss.add_string_value( name( 1 ), best_decoy );
	pss.add_energy( name( 2 ), best_rmsd );
// //store in Job-Object:
// 	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair( "pool_converged_tag", best_decoy );
// 	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair( "pool_converged_rmsd", best_rmsd );

}

///@brief evaluate pose and store values in Silent_Struct
void Pool_Evaluator::apply( core::io::silent::SilentStruct &pss) const {
	core::Real best_rmsd;
	std::string best_decoy;

	rmsd_pool_->evaluate( pss, best_decoy, best_rmsd );

	pss.add_string_value( name( 1 ), best_decoy );
	pss.add_energy( name( 2 ), best_rmsd );
// //store in Job-Object:
// 	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair( "pool_converged_tag", best_decoy );
// 	protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_real_pair( "pool_converged_rmsd", best_rmsd );

}

} // mc_convergence_check
} // moves
} // rosetta
