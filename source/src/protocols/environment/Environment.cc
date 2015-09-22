// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/Environment.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/Environment.hh>

// Package headers
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MoverApplyingMover.hh>

#include <core/environment/DofPassport.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <set>
#include <boost/foreach.hpp>

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.Environment", basic::t_info );

namespace protocols {
namespace environment {

Environment::Environment( std::string name ):
	Parent( name ),
	broker_( /* NULL */ ),
	bAutoCut_( false ),
	bInheritCuts_( true ),
	bAllowPureMovers_( false )
{}

Environment::~Environment() {
	BOOST_FOREACH ( Conformation * conf, pconfs_ ) {
		ProtectedConformation * pconf = dynamic_cast< ProtectedConformation * >( conf );
		if ( pconf ) {
			pconf->env_destruction();
		}
	}

	cancel_passports();
}

void Environment::register_mover( moves::MoverOP mover ){
	// I think this is kind of ugly, but I don't think this happens too often, so it's probably ok.
	ClientMoverOP claiming_mover = utility::pointer::dynamic_pointer_cast< protocols::environment::ClientMover > ( mover );
	moves::MoverContainerOP mover_container = utility::pointer::dynamic_pointer_cast< moves::MoverContainer > ( mover );
	moves::MoverApplyingMoverOP mover_applier = utility::pointer::dynamic_pointer_cast< moves::MoverApplyingMover > ( mover );

	if ( claiming_mover ) {
		if ( !is_registered( claiming_mover ) ) {
			tr.Trace << "Registered mover '" << claiming_mover->get_name() << "'" << std::endl;
			registered_movers_.insert( claiming_mover );

			std::set< ClientMoverOP > submovers;
			claiming_mover->yield_submovers( submovers );
			register_movers( submovers.begin(), submovers.end() );
		}
	} else if ( mover_container ) {
		for ( Size i = 0; i < mover_container->nr_moves(); ++i ) {
			moves::MoverOP const& mover = mover_container->movers()[ i ];
			register_mover( mover );
		}
	} else if ( mover_applier ) {
		register_mover( mover_applier->mover() );
	} else if ( !allow_pure_movers() ) {
		std::ostringstream err;
		err << "The mover '" << mover->get_name()
			<< "' is not a ClientMover or a MoverContainer, and thus cannot be used inside an BrokeredEnvironment. "
			<< "If you're sure you want to include this mover, set the option 'allow_pure_movers' to true. "
			<< "This will cause the Environment to ignore Environment-incompatible movers during registration.";
		throw utility::excn::EXCN_BadInput( err.str() );
	}
}

core::pose::Pose Environment::start( core::pose::Pose const& in_pose ){

	this->input_pose_ = in_pose;

	// This cast of reference to OP is ok because we know that Conformations are stored as pointers anyway.
	// ^^ not OK with std::shared_ptr; changing to const*.
	Conformation const* in_conf = &in_pose.conformation();

	//figure out if we've got a superenvironment... also initialize Annotations object
	ProtectedConformation const* conf_ptr = dynamic_cast<ProtectedConformation const*>( in_conf );
	if ( conf_ptr ) {
		set_superenv( conf_ptr->environment() );

		// Keep old annotations by copy-constructing a new annotations object.
		ann_ = SequenceAnnotationOP( new SequenceAnnotation( *( conf_ptr->annotations() ) ) );
	} else {
		set_superenv( core::environment::EnvCoreCAP() );

		// Unprotected Conformation objects don't have annotations (yet?). We then build the annotation object.
		ann_ = SequenceAnnotationOP( new SequenceAnnotation( in_conf->size() ) );
	}

	tr.Debug << "Start environment: '" << name() << "'" << std::endl;

	core::pose::Pose broker_result = this->broker( in_pose );

	assert( dynamic_cast< basic::datacache::WriteableCacheableMap* >( broker_result.data().get_raw_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ) );

	return broker_result;
}

core::pose::Pose Environment::end( core::pose::Pose const& pose ){
	ProtectedConformationCOP conf = utility::pointer::dynamic_pointer_cast< ProtectedConformation const >( pose.conformation_ptr() );
	if ( !conf ) {
		tr.Error << "[ERROR] Environment::end recieved a pose that contains an unprotcted Conformation."
			<< std::endl;
		throw utility::excn::EXCN_BadInput( "Nonprotected pose came in to Environment::end" );
	}

	core::pose::Pose new_pose = pose;

	// Setting a new conformation clears Pose::datacache, and we want to protect the WRITABLE_DATA.
	if ( new_pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ) {

		if ( ! broker_->result().pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ) {
			tr.Warning << "[WARNING] The brokered pose contained WriteableCacheable data in the Pose DataCache, but this "
				<< "data is no longer present in the DataCache at environment closing. This is probably "
				<< "because a call was made to pose::set_new_conformation somewhere (anywhere), which clears the "
				<< "cache. Certain features of one or more of the ClientMovers registered to the Environment '"
				<< name() << "' are nonfunctional." << std::endl;
		}

		basic::datacache::WriteableCacheableMapOP data_map = utility::pointer::dynamic_pointer_cast< basic::datacache::WriteableCacheableMap >( new_pose.data().get_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) );
		new_pose.set_new_conformation( end( conf ) );
		new_pose.data().set( core::pose::datacache::CacheableDataType::WRITEABLE_DATA, data_map );
	} else {
		new_pose.set_new_conformation( end(conf ) );
	}

	tr.Debug << "Finish environment " << name() << " with sequence: " << new_pose.annotated_sequence()
		<< " and fold tree: " << new_pose.fold_tree() << std::endl;

	if ( input_pose_.pdb_info() ) {
		tr.Trace << "  Applying old PDBInfo object with " << input_pose_.pdb_info()->nres() << " residues to pose with "
			<< new_pose.total_residue() << " residues." << std::endl;

		new_pose.pdb_info( core::pose::PDBInfoOP( new core::pose::PDBInfo( *( input_pose_.pdb_info() ) ) ) );
	} else {
		tr.Trace << "  No PDBInfo object being built in to Environment's output pose, "
			<< "as it appears the input pose didn't have one." << std::endl;
	}

	return new_pose;
}

core::conformation::ConformationOP Environment::end( ProtectedConformationCOP conf ){
	tr.Debug << "End environment: '" << name() << "'" << std::endl;

	core::pose::Pose pose;

	ConformationOP ret_conf( new Conformation( *conf ) );
	pose.set_new_conformation( ret_conf );

	remove_nonpermenant_features( pose );

	try {
		tr.Debug << "Applying original fold tree " << input_pose_.fold_tree() << std::endl;
		pose.fold_tree( input_pose_.fold_tree() );
	} catch ( ... ){
		tr.Error << "Environment " << name() << " failed to apply the input fold tree at env closing." << std::endl
			<< "Brokered fold tree was: " << conf->core::conformation::Conformation::fold_tree()
			<< "Attempted closed fold tree was: " << input_pose_.fold_tree();
		throw;
	}

// Reprotect Conformation if there's a superenvironment.
	if ( ! superenv().expired() ) {
		ret_conf = ConformationOP( new ProtectedConformation( superenv(), pose.conformation() ) );
	} else {
		ret_conf = ConformationOP( new Conformation( pose.conformation() ) );
	}

	cancel_passports();

	return ret_conf;
}

void Environment::assign_passport( ClientMoverOP mover, core::environment::DofPassportCOP passport ){
	mover->push_passport( get_self_weak_ptr() , passport );
}

void Environment::cancel_passports(){
	BOOST_FOREACH ( ClientMoverOP mover, registered_movers_ ) {
		tr.Info << "stripping passport from " << mover->get_name() << std::endl;
		if ( mover->has_passport() &&
				mover->passport()->env_id() == this->id() ) {
			mover->pop_passport( *this );
		}
	}
}

/// ENV OPEN/CLOSING
void Environment::remove_nonpermenant_features( core::pose::Pose& pose ){

	//reset chainbreak varaints
	tr.Trace << "Removing chainbreak variants from sequence "
		<< pose.annotated_sequence() << std::endl;
	BOOST_FOREACH ( std::set< core::Size >::key_type cut_res,
			broker()->result().auto_cuts ) {
		remove_chainbreak_variants( pose, cut_res, cut_res+1 );
		tr.Trace << "Chainbreak variants @ " << cut_res << "/" << cut_res+1 << " removed:"
			<< pose.annotated_sequence() << std::endl;
	}

	// remove virtual residues added in this environment.
	// jump edges must have atom names removed because the jump
	// gets reapplied to whichever residue is convenient.
	if ( broker()->result().new_vrts.size() > 0 ) {
		core::kinematics::FoldTree ft_clean( pose.fold_tree() );
		for ( core::Size i = 1; i <= ft_clean.num_jump(); ++i ) {
			ft_clean.set_jump_atoms( (int) i, "", "" );
		}
		pose.fold_tree( ft_clean );
		pose.conformation().delete_residue_range_slow( input_pose_.fold_tree().nres()+1,
			pose.total_residue() );
	}

	//Strip out pose datacache elements
	//  using namespace basic::datacache;
	//  using namespace core::pose::datacache;
	//  typedef std::map< std::string, std::set< WriteableCacheableDataOP > > DataMap;
	//  WriteableCacheableMapOP map = dynamic_cast< WriteableCacheableMap* >( pose.data().get_raw_ptr( CacheableDataType::WRITEABLE_DATA ) );
	//
	//  for( DataMap::const_iterator it = broker()->result().cached_data->begin();
	//       it != broker()->result().cached_data->end(); ++it ){
	//    if( map->find( it->first ) != map->end() ){
	//      BOOST_FOREACH( WriteableCacheableDataOP d, it->second ){
	//        map->erase( d );
	//      }
	//    }
	//  }
}

void Environment::auto_cut( bool setting ){
	if ( broker() ) {
		std::ostringstream ss;
		ss << "The Environment '" << name() << "' was asked to set auto_cut to "
			<< setting << ", but broking was already completed." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( ss.str() );
	}
	bAutoCut_ = setting;
}

void Environment::inherit_cuts( bool setting ) {
	if ( broker() ) {
		std::ostringstream ss;
		ss << "The Environment '" << name() << "' was asked to set inherit_cuts to "
			<< setting << ", but broking was already completed." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( ss.str() );
	}
	bInheritCuts_ = setting;
}


void Environment::allow_pure_movers( bool setting ) {
	if ( broker() ) {
		std::ostringstream ss;
		ss << "The Environment '" << name() << "' was asked to set allow_pure_movers to "
			<< setting << ", but broking was already completed." << std::endl;
		throw utility::excn::EXCN_Msg_Exception( ss.str() );
	}
	bAllowPureMovers_ = setting;
}

core::pose::Pose Environment::broker( core::pose::Pose const& in_pose ){
	tr.Debug << "Beginning Broking for environment " << this->name() << " with "
		<< registered_movers_.size() << " registered movers." << std::endl;
	tr.Debug << "  Registered movers are:" << std::endl;
	BOOST_FOREACH ( ClientMoverOP mover, registered_movers_ ) { tr.Debug << "    " << mover->get_name() << std::endl; }

	std::map< ClientMoverOP, core::environment::DofPassportOP > mover_passports;

	BOOST_FOREACH ( ClientMoverOP mover, registered_movers_ ) {
		assert( mover_passports.find( mover ) == mover_passports.end() );
		mover_passports[ mover ] = Parent::issue_passport( mover->get_name() );
		assign_passport( mover, mover_passports[ mover ] );
	}

	core::pose::Pose out_pose;
	try{
		broker_ = EnvClaimBrokerOP( new EnvClaimBroker( get_self_weak_ptr(), mover_passports, in_pose, ann_ ) );
		out_pose = broker_->result().pose;
	} catch ( ... ){
		// For exception safety. Revoke (possibly unfinished and/or now invalid) passports
		cancel_passports();
		throw;
	}

	return out_pose;
}

void Environment::remove_chainbreak_variants( core::pose::Pose& pose, core::Size down_res, core::Size up_res ) const {
	using namespace core::chemical;
	using namespace core::conformation;

	Conformation& conf = pose.conformation();

	Residue const& rsd_lower( conf.residue( down_res ) );
	Residue const& rsd_upper( conf.residue( up_res ) );
	ResidueTypeSet const& rsd_set( rsd_lower.residue_type_set() );

	ResidueType const& new_type_lower( rsd_set.get_residue_type_with_variant_removed( rsd_lower.type(),
		CUTPOINT_LOWER ) );
	ResidueType const& new_type_upper( rsd_set.get_residue_type_with_variant_removed( rsd_upper.type(),
		CUTPOINT_UPPER ) );

	ResidueOP new_lower( ResidueFactory::create_residue( new_type_lower, rsd_lower, conf ) );
	ResidueOP new_upper( ResidueFactory::create_residue( new_type_upper, rsd_upper, conf ) );

	copy_residue_coordinates_and_rebuild_missing_atoms( rsd_lower, *new_lower, conf );
	copy_residue_coordinates_and_rebuild_missing_atoms( rsd_upper, *new_upper, conf );

	conf.replace_residue( rsd_lower.seqpos(), *new_lower, true );
	conf.replace_residue( rsd_upper.seqpos(), *new_upper, true );
}

bool Environment::is_registered( ClientMoverOP mover ) const {
	return ( std::find( registered_movers_.begin(), registered_movers_.end(), mover ) != registered_movers_.end() );
}

EnvironmentCAP Environment::superenv() const{
	if ( ! Parent::superenv().expired() ) {
		EnvironmentCOP env = utility::pointer::dynamic_pointer_cast< Environment const >( Parent::superenv().lock() );
		assert( env != 0 );
		return EnvironmentCAP( env );
	} else {
		return EnvironmentCAP();
	}
}

} // environment
} // protocols
