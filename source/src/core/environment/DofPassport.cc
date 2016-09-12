// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/DofPassport.cc
/// @author Justin R. Porter

// Unit Headers
#include <core/environment/DofPassport.hh>

// Package headers
#include <core/environment/EnvCore.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Conformation.hh>

#include <utility/excn/Exceptions.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static THREAD_LOCAL basic::Tracer tr( "protocols.environment.DofPassport", basic::t_info );

namespace core {
namespace environment {

using core::kinematics::MoveMapOP;

DofPassport::DofPassport( std::string const & mover,
	core::Size env ) :
	mover_( mover ),
	env_id_( env )
{}

DofPassport::~DofPassport() = default;

std::set< core::id::DOF_ID >::const_iterator DofPassport::begin() const {
	if ( accessible_dofs_.size() == 0 ) {
		return accessible_dofs_.end();
	} else if ( accessible_dofs_.begin()->rsd() == 0 ) {
		return ++accessible_dofs_.begin();
	} else {
		return accessible_dofs_.begin();
	}

}


//Movemap Config

MoveMapOP DofPassport::render_movemap() const {
	kinematics::MoveMapOP mm( new kinematics::MoveMap() );
	this->render_movemap( mm );
	return mm;
}

void DofPassport::render_movemap( core::kinematics::MoveMapOP mm ) const {
	for ( auto const & accessible_dof : accessible_dofs_ ) {
		mm->set( accessible_dof, true);
	}

	for ( core::Size seqpos = 1; seqpos <= conf_->size(); ++seqpos ) {
		//    if( conf_->residue( seqpos ).is_protein() ){

		// false if no bb torsions
		Size const n_bbs = conf_->residue( seqpos ).mainchain_torsions().size();
		bool seqpos_access = (bool) n_bbs;

		for ( Size i = 1; i <= n_bbs; ++i ) {
			id::TorsionID t_id = id::TorsionID( seqpos, id::BB, i );
			id::DOF_ID d_id = conf_->dof_id_from_torsion_id( t_id );

			if ( d_id.valid() ) {
				seqpos_access &= accessible_dofs_.find( d_id ) != accessible_dofs_.end();
			}
		}
		mm->set_bb( seqpos, seqpos_access );

		// false if no chi torsions
		seqpos_access = conf_->residue( seqpos ).nchi();
		for ( core::Size chi_i = 1; chi_i <= conf_->residue( seqpos ).nchi(); ++chi_i ) {
			id::TorsionID t_id = id::TorsionID( seqpos, id::CHI, chi_i );
			id::DOF_ID d_id = conf_->dof_id_from_torsion_id( t_id );

			if ( d_id.valid() ) {
				seqpos_access &= accessible_dofs_.find( d_id ) != accessible_dofs_.end();
			}
		}
		mm->set_chi( seqpos, seqpos_access );
	}

	for ( int jump_i = 1; jump_i <= (int) conf_->fold_tree().num_jump(); ++jump_i ) {
		bool allow = has_jump_access( jump_i );
		mm->set_jump( jump_i, allow );
		mm->set_jump( conf_->fold_tree().upstream_jump_residue( jump_i ),
			conf_->fold_tree().downstream_jump_residue( jump_i ),
			allow );
	}
}

bool DofPassport::has_jump_access( int jump_num ) const {
	bool allow = true;
	for ( core::Size rb_i = id::RB1; rb_i <= id::RB6; ++rb_i ) {
		id::AtomID a_id = conf_->jump_atom_id( jump_num );
		allow &= accessible_dofs_.find( id::DOF_ID( a_id, id::DOF_Type(rb_i) ) ) != accessible_dofs_.end();
	}
	return allow;
}

utility::vector1< int > DofPassport::active_jumps() const {
	utility::vector1< int > active_jumps;
	for ( int i = 1; i <= (int) conf_->fold_tree().num_jump(); ++i ) {
		if ( has_jump_access( i ) ) {
			active_jumps.push_back( i );
		}
	}
	return active_jumps;
}

std::set< core::id::DOF_ID > const& DofPassport::active_dofs() const {
	return accessible_dofs_;
}

void DofPassport::add_dof_access( id::DOF_ID const& dof_id ){
	accessible_dofs_.insert( dof_id );
}

void DofPassport::revoke_all_access(){
	accessible_dofs_.clear();
}

//Security Checks
bool DofPassport::access_check( EnvCore const& env, bool type_specific_check ) const {
	if ( env.id() == env_id_ ) {
		if ( type_specific_check ) {
			if ( tr.Trace.visible() ) tr.Trace << "allowed" << std::endl;
			return true;
		} else {
			if ( tr.Trace.visible() ) tr.Trace << "denied (right not given)" << std::endl;
			return false;
		}
	} else {
		if ( tr.Trace.visible() ) tr.Trace << "denied (environmental conflict)" << std::endl;
		return false;
	}
}

bool DofPassport::dof_access( id::DOF_ID const& id ) const {
	return ( accessible_dofs_.find( id ) != accessible_dofs_.end() );
}

bool DofPassport::dof_access( EnvCore const& env, id::DOF_ID const& id ) const {
	if ( tr.Trace.visible() ) tr.Trace << "Verifying access to DOF_ID " << id << ": ";
	return access_check( env, dof_access( id ) );
}

void DofPassport::reference_conformation( conformation::ConformationCOP conf ) {
	conf_ = conf;
}

conformation::ConformationCOP DofPassport::reference_conformation() const {
	return conf_;
}

std::string const& DofPassport::mover() const {
	return mover_;
}

void DofPassport::show( std::ostream& str ) const {
	str << "DofPassport: " << mover_ << ", " << env_id_;
}

std::ostream& operator<<( std::ostream& str, DofPassport const& dof_passport ) {
	dof_passport.show( str );
	return str;
}

} // environment
} // core
