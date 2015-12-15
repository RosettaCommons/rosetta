// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file .cc file for enzdes cacheable observer
/// @brief
/// @author Florian Richter, floric@u.washington.edu, september 09


//unit headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesSeqRecoveryCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesLoopsFile.hh>

//project headers
#include <basic/Tracer.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
//MSA
#include <core/sequence/SequenceProfile.hh>

#include <core/scoring/constraints/SequenceProfileConstraint.hh> //msa

//utility headers
#include <utility/signals/Link.hh>

//option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.EnzdesCacheableObserver" );


EnzdesCacheableObserverOP
get_enzdes_observer(
	core::pose::Pose & pose )
{
	//using namespace core::pose::datacache::CacheableObserverType;
	using namespace core::pose::datacache;

	if ( !pose.observer_cache().has( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER ) ) {

		EnzdesCacheableObserverOP enz_obs( new EnzdesCacheableObserver() );
		pose.observer_cache().set( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER, enz_obs, true );
		//std::cout << " setting new cacheable observer for pose " << std::endl;
	}
	CacheableObserverOP enz_obs = pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER );
	//return (pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER ) );
	//std::cout << " returning nonconst enzdes observer " << std::endl;
	return utility::pointer::dynamic_pointer_cast< EnzdesCacheableObserver >( enz_obs );
	//return static_cast< EnzdesCacheableObserverOP > (enz_obs);
}

EnzdesCacheableObserverCOP
get_enzdes_observer(
	core::pose::Pose const & pose )
{
	//using namespace core::pose::datacache::CacheableObserverType;
	using namespace core::pose::datacache;

	//const access: if cacheable observer hasn't been set, return NULL pointer
	if ( !pose.observer_cache().has( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER ) ) return NULL;

	CacheableObserverCOP enz_obs = pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER );
	//std::cout << " returning const enzdes observer with val " << enz_obs << std::endl;
	return utility::pointer::static_pointer_cast< EnzdesCacheableObserver const >( enz_obs );
	//return static_cast< EnzdesCacheableObserverCOP > (pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::ENZDES_OBSERVER ) );
}

EnzdesCacheableObserver::EnzdesCacheableObserver()
: CacheableObserver(),
	cst_cache_(/* NULL */),
	seq_recovery_cache_(/* NULL */),
	enz_loops_file_(/* NULL */)
{
	favor_native_constraints_.clear();
	lig_rigid_body_confs_.clear();
}

EnzdesCacheableObserver::EnzdesCacheableObserver( EnzdesCacheableObserver const & other )
: CacheableObserver( other ),
	cst_cache_( /* NULL */ ),
	seq_recovery_cache_(/* NULL */),
	favor_native_constraints_(other.favor_native_constraints_),
	enz_loops_file_(other.enz_loops_file_ ),
	lig_rigid_body_confs_( other.lig_rigid_body_confs_ )
{
	if ( other.cst_cache_ ) cst_cache_ = toolbox::match_enzdes_util::EnzdesCstCacheOP( new toolbox::match_enzdes_util::EnzdesCstCache( *(other.cst_cache_) ) );
	if ( other.seq_recovery_cache_ ) seq_recovery_cache_ = EnzdesSeqRecoveryCacheOP( new EnzdesSeqRecoveryCache( *(other.seq_recovery_cache_) ) );
}

EnzdesCacheableObserver::~EnzdesCacheableObserver(){}

core::pose::datacache::CacheableObserverOP
EnzdesCacheableObserver::clone()
{
	return core::pose::datacache::CacheableObserverOP( new EnzdesCacheableObserver( *this ) );
}

core::pose::datacache::CacheableObserverOP
EnzdesCacheableObserver::create()
{
	return core::pose::datacache::CacheableObserverOP( new EnzdesCacheableObserver() );
}

bool
EnzdesCacheableObserver::is_attached() const {
	return length_event_link_.valid(); }

void
EnzdesCacheableObserver::attach_impl( core::pose::Pose & pose ){

	length_event_link_ = pose.conformation().attach_length_obs( &EnzdesCacheableObserver::on_length_change, this );

}

void
EnzdesCacheableObserver::detach_impl(){
	length_event_link_.invalidate();
}

void
EnzdesCacheableObserver::on_length_change( core::conformation::signals::LengthEvent const & event ){

	if ( event.tag == core::conformation::signals::LengthEvent::INVALIDATE ) {
		//don't know what the best behaviour is in this case
		//probably nothing, because pose destruction is imminent
		return;
	}
	core::id::SequenceMapping smap( event );
	if ( cst_cache_ ) cst_cache_->remap_resid( smap );
	if ( seq_recovery_cache_ ) seq_recovery_cache_ -> remap_residues( smap );

	//remap favor native constraints (if they exist)
	//note: in case residues were deleted, the constraints
	//will be null pointers, so the vector might change size
	if ( favor_native_constraints_.size() > 0 ) {
		utility::vector1< core::scoring::constraints::ConstraintCOP > new_favor_native_csts;
		for ( core::scoring::constraints::ConstraintCOPs::iterator cst_it = favor_native_constraints_.begin();
				cst_it != favor_native_constraints_.end(); ++cst_it ) {
			core::scoring::constraints::ConstraintCOP remappedcst = (*cst_it)->remap_resid( smap );
			if ( remappedcst ) new_favor_native_csts.push_back( remappedcst );
		}
		favor_native_constraints_.clear();
		favor_native_constraints_ = new_favor_native_csts;
	}
}


void
EnzdesCacheableObserver::set_cst_cache(
	toolbox::match_enzdes_util::EnzdesCstCacheOP cst_cache )
{
	cst_cache_ = cst_cache;
}

toolbox::match_enzdes_util::EnzdesCstCacheOP
EnzdesCacheableObserver::cst_cache()
{
	return cst_cache_;
}

toolbox::match_enzdes_util::EnzdesCstCacheCOP
EnzdesCacheableObserver::cst_cache() const
{
	return cst_cache_;
}

std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > > const &
EnzdesCacheableObserver::lig_rigid_body_confs() const{
	return lig_rigid_body_confs_;
}


void
EnzdesCacheableObserver::set_rigid_body_confs_for_lig(
	core::Size seqpos,
	utility::vector1< core::conformation::ResidueCOP > const & rg_confs
)
{
	erase_rigid_body_confs_for_lig( seqpos );
	utility::vector1< core::conformation::ResidueCOP > resvec;
	for ( core::Size i = 1; i <= rg_confs.size(); ++i ) {
		core::conformation::ResidueOP res( rg_confs[i]->clone());
		res->seqpos( seqpos ); //security measure
		resvec.push_back( res );
	}
	lig_rigid_body_confs_.insert( std::pair< core::Size, utility::vector1< core::conformation::ResidueCOP > >( seqpos, resvec ) );

	//conf_it = lig_rigid_body_confs_.find( seqpos );
	//std::cout << "enzdes cache observer rb confs for position " << seqpos << " are being set. " << conf_it->second.size() << " rb positions were set" << std::endl;
}

void
EnzdesCacheableObserver::erase_rigid_body_confs_for_lig(
	core::Size seqpos )
{
	std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > >::iterator conf_it( lig_rigid_body_confs_.find( seqpos ) );
	if ( conf_it != lig_rigid_body_confs_.end() ) lig_rigid_body_confs_.erase( conf_it );
}

EnzdesSeqRecoveryCacheOP EnzdesCacheableObserver::get_seq_recovery_cache(){
	return seq_recovery_cache_;
}

EnzdesSeqRecoveryCacheCOP EnzdesCacheableObserver::get_seq_recovery_cache() const {
	return seq_recovery_cache_;
}

void EnzdesCacheableObserver::set_seq_recovery_cache(
	EnzdesSeqRecoveryCacheOP seq_recovery_cache
){
	seq_recovery_cache_ = seq_recovery_cache;
}

void
EnzdesCacheableObserver::set_enzdes_loops_file(
	EnzdesLoopsFileCOP loopfile_in
){
	enz_loops_file_ = loopfile_in;
}

EnzdesLoopsFileCOP
EnzdesCacheableObserver::enzdes_loops_file() const
{
	return enz_loops_file_;
}


void
EnzdesCacheableObserver::setup_favor_native_constraints(
	core::pose::Pose & pose,
	core::pack::task::PackerTaskCOP task,
	core::pose::Pose const & native_pose
)
{
	using namespace basic::options;

	if ( option[OptionKeys::enzdes::favor_native_res].user() ) {
		using namespace core::scoring::constraints;

		core::Real bonus = option[OptionKeys::enzdes::favor_native_res].value();

		tr.Info << "favor_native_res: adding a bonus of " << bonus << " for native residues to pose." << std::endl;

		//safety check first
		if ( favor_native_constraints_.size() != 0 ) {
			tr.Info << "Warning: when setting up favor native constraints, there might already be some previously generated favor_native constraints in the pose, trying to remove these first." << std::endl;
			remove_favor_native_constraints( pose );

		}

		favor_native_constraints_.clear();
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {

			if ( task->design_residue(i) ) {

				ConstraintOP resconstraint( new ResidueTypeConstraint( native_pose, i, bonus ) );
				favor_native_constraints_.push_back( resconstraint );

			}
		}
		favor_native_constraints_ = pose.add_constraints( favor_native_constraints_ );
	} else if ( option[ OptionKeys::in::file::pssm ].user() ) {
		//multiple sequence aligniment (adapted from MSA app)

		using namespace core;
		using namespace scoring;
		using namespace constraints;
		using namespace sequence;

		using namespace protocols;

		tr << " Starting MSA design " << std::endl;

		// register SequenceProfileConstraint with the ConstraintFactory so that it can be constructed from a constraint file
		//ConstraintIO::get_cst_factory().add_type(
		//new core::scoring::constraints::SequenceProfileConstraint( Size(), utility::vector1< id::AtomID >(), NULL ) );

		// add constraints to bias design toward a sequence profile
		SequenceProfileOP profile( new SequenceProfile );
		utility::file::FileName filename( option[ OptionKeys::in::file::pssm ]().front() );

		profile->read_from_file( filename );
		profile->convert_profile_to_probs( 1 ); // was previously implicit in read_from_file()

		tr << *profile << std::endl;

		for ( Size seqpos(1), end( pose.total_residue() ); seqpos <= end; ++seqpos ) {
			// add individual profile constraint for each residue position
			// because of the underlying constraint implementation, this enures that the constraint is
			// a context-independent 1-body energy, or (intra)residue constraint
			pose.add_constraint( scoring::constraints::ConstraintCOP( scoring::constraints::ConstraintOP( new core::scoring::constraints::SequenceProfileConstraint( pose, seqpos, profile ) ) ) );
		}
	}// else if ( option[ OptionKeys::constraints::in::file::pssm ].user() ){
	//}

}//setup_favor_native_constraints

//migrate in second step
void
EnzdesCacheableObserver::remove_favor_native_constraints(
	core::pose::Pose & pose
)
{
	if ( !( pose.remove_constraints( favor_native_constraints_ ) ) ) {
		tr.Info << "Warning: some of the favor native constraints that were previously added to the pose are not there anymore, something's a little unclean somewhere." << std::endl;
	}
	favor_native_constraints_.clear();
}


} //match_enzdes_util
} //toolbox
} //protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesCacheableObserver::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::datacache::CacheableObserver >( this ) );
	arc( CEREAL_NVP( cst_cache_ ) ); // toolbox::match_enzdes_util::EnzdesCstCacheOP
	arc( CEREAL_NVP( seq_recovery_cache_ ) ); // EnzdesSeqRecoveryCacheOP
	arc( CEREAL_NVP( favor_native_constraints_ ) ); // core::scoring::constraints::ConstraintCOPs
	arc( CEREAL_NVP( enz_loops_file_ ) ); // EnzdesLoopsFileCOP
	arc( CEREAL_NVP( lig_rigid_body_confs_ ) ); // std::map<core::Size, utility::vector1<core::conformation::ResidueCOP> >
	// EXEMPT length_event_link_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::match_enzdes_util::EnzdesCacheableObserver::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::datacache::CacheableObserver >( this ) );
	arc( cst_cache_ ); // toolbox::match_enzdes_util::EnzdesCstCacheOP
	arc( seq_recovery_cache_ ); // EnzdesSeqRecoveryCacheOP
	utility::vector1< std::shared_ptr< core::scoring::constraints::Constraint > > local_favor_native_constraints;
	arc( local_favor_native_constraints ); // core::scoring::constraints::ConstraintCOPs
	favor_native_constraints_ = local_favor_native_constraints; // copy the non-const pointer(s) into the const pointer(s)

	std::shared_ptr< protocols::toolbox::match_enzdes_util::EnzdesLoopsFile > local_enz_loops_file;
	arc( local_enz_loops_file ); // EnzdesLoopsFileCOP
	enz_loops_file_ = local_enz_loops_file; // copy the non-const pointer(s) into the const pointer(s)

	std::map< core::Size, utility::vector1< core::conformation::ResidueOP > > local_lig_rigid_body_confs;
	arc( local_lig_rigid_body_confs ); // std::map<core::Size, utility::vector1<core::conformation::ResidueCOP> >
	// lig_rigid_body_confs_ = local_lig_rigid_body_confs; // copy the non-const pointer(s) into the const pointer(s)
	for ( std::map< core::Size, utility::vector1< core::conformation::ResidueOP > >::const_iterator
			iter = local_lig_rigid_body_confs.begin(), iter_end = local_lig_rigid_body_confs.end();
			iter != iter_end; ++iter ) {
		lig_rigid_body_confs_[ iter->first ] = iter->second;
	}
	// EXEMPT length_event_link_

}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::match_enzdes_util::EnzdesCacheableObserver );
CEREAL_REGISTER_TYPE( protocols::toolbox::match_enzdes_util::EnzdesCacheableObserver )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_match_enzdes_util_EnzdesCacheableObserver )
#endif // SERIALIZATION
