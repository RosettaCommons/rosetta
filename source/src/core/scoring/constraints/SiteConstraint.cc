// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/SiteConstraint.cc
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of all residues in another chain.
/// @author Brian Weitzner (brian.weitzner@jhu.edu, May 2011)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - residue subsetting

#include <core/scoring/constraints/SiteConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

static basic::Tracer TR( "core.scoring.constraints.SiteConstraint" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraint::SiteConstraint():
	AmbiguousConstraint()
{}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraint::SiteConstraint( ConstraintCOPs const & cst_in ) :
	AmbiguousConstraint( cst_in )
{}

///
ConstraintOP SiteConstraint::clone() const {
	return ConstraintOP( new SiteConstraint( *this ) );
}

bool SiteConstraint::operator == ( Constraint const & rhs ) const {
	return AmbiguousConstraint::operator == ( rhs );
}

bool SiteConstraint::same_type_as_me( Constraint const & other ) const {
	return dynamic_cast< SiteConstraint const * > (&other);
}


std::string SiteConstraint::type() const {
	return "SiteConstraint";
}


void
SiteConstraint::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	//out << "AmbiguousConstraint Active constraint:" << std::endl;
	//active_constraint()->show(out);
	out << "SiteConstraint is an AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for ( auto const & cst_it : member_constraints() ) {
		cst_it->show( out );
	}

	out << " ...all member constraints of this SiteConstraint shown." << std::endl;
}


void
SiteConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory
){
	TR.Debug << "ConstraintIO::read_site_cst" << std::endl;
	Size res;
	std::string tempres;
	std::string atm_name;
	std::string chain;
	std::string func_type;
	data >> atm_name >> tempres >> chain >> func_type;

	ConstraintIO::parse_residue( pose, tempres, res );
	TR.Info << "read: " << atm_name << " " << res << " " << "constrain to chain: " << chain <<
		" func: " << func_type << std::endl;


	func::FuncOP aFunc = func_factory.new_func( func_type );
	aFunc->read_data( data );

	if ( TR.Debug.visible() ) {
		aFunc->show_definition( TR.Debug ); TR.Debug<<std::endl;
	}
	setup_csts( res, atm_name, chain, pose, aFunc );

	if ( data.good() ) {
		//chu skip the rest of line since this is a single line defintion.
		while ( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) { data.setstate( std::ios_base::eofbit ); }
	}

	if ( TR.Debug.visible() ) {
		aFunc->show_definition( TR.Debug );
		TR.Debug << std::endl;
	}

} // read_def

void
SiteConstraint::setup_csts(
	Size res,
	std::string atm_name,
	std::string chain,
	core::pose::Pose const & pose,
	func::FuncOP const & func
) {
	//Size target_chain = pose.chain( res );
	Size constraint_chain = pose::get_chain_id_from_chain( chain, pose );

	Size start_res = pose.conformation().chain_begin( constraint_chain );
	Size end_res = pose.conformation().chain_end( constraint_chain );
	utility::vector1<bool> residues(pose.size(), false);

	for ( Size j = start_res ; j <= end_res ; ++j ) {
		residues[j] = true;
	}
	setup_csts(res, atm_name, residues, pose, func);

} // setup_csts

void
SiteConstraint::setup_csts(
	Size res,
	std::string atm_name,
	utility::vector1<bool> const & residues,
	core::pose::Pose const & pose,
	func::FuncOP const & func
) {
	debug_assert( pose.size() == residues.size() );
	id::AtomID target_atom( pose.residue_type( res ).atom_index( atm_name ), res );

	for ( Size i = 1 ; i <= pose.size() ; ++i ) {

		if ( residues[ i ] ) {
			id::AtomID atom2;
			if ( pose.residue( i ).is_protein() ) {
				atom2 = id::AtomID( pose.residue_type( i ).atom_index( "CA" ), i );
			} else if ( pose.residue( i ).is_carbohydrate() ) {
				atom2 = id::AtomID( pose.residue( i ).carbohydrate_info()->anomeric_carbon_index(), i );
			} else {
				TR.Warning << "Residue " << i << " is not a protein or a carbohydrate; ";
				TR.Warning << "Setting the nbr_atom as the alternative for 'CA' and 'C1'" << std::endl;
				atom2 = id::AtomID( pose.residue( i ).nbr_atom(), i );
			}
			runtime_assert( target_atom.valid() && atom2.valid() );
			add_individual_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( target_atom, atom2, func ) ) ) );
		}
	}

} // setup_csts

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::SiteConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< AmbiguousConstraint >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::SiteConstraint::load( Archive & arc ) {
	arc( cereal::base_class< AmbiguousConstraint >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::SiteConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::SiteConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_SiteConstraint )
#endif // SERIALIZATION
