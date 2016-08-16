// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/constraints/SiteConstraintResidues.cc
/// @brief This class is an AmbiguousConstraint in which the set is comprised of AtomPairConstraints
/// @brief of an atom of interest in one chain versus the CA of of another sets of residues
/// @author Lei Shi (shilei@uw.edu)

#include <core/scoring/constraints/SiteConstraintResidues.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <core/id/AtomID.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.constraints.SiteConstraintResidues" );

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
SiteConstraintResidues::SiteConstraintResidues() :
	AmbiguousConstraint()
{}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
SiteConstraintResidues::SiteConstraintResidues( ConstraintCOPs const & cst_in ) :
	AmbiguousConstraint( cst_in )
{}

///
ConstraintOP SiteConstraintResidues::clone() const {
	return ConstraintOP( new SiteConstraintResidues( *this ));
}

bool SiteConstraintResidues::operator == ( Constraint const & rhs ) const
{
	return AmbiguousConstraint::operator== ( rhs );
}

bool SiteConstraintResidues::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< SiteConstraintResidues const * > ( &other );
}

std::string SiteConstraintResidues::type() const {
	return "SiteConstraintResidues";
}


void
SiteConstraintResidues::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	//out << "AmbiguousConstraint Active constraint:" << std::endl;
	//active_constraint()->show(out);
	out << "SiteConstraintResidues is an AmbiguousConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
	for ( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(), end = member_constraints().end(); cst_it != end; ++cst_it ) {
		( *cst_it )->show( out );
	}

	out << " ...all member constraints of this SiteConstraintResidues shown." << std::endl;
}


void
SiteConstraintResidues::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory
) {
	TR.Debug << "read_site_cst" << std::endl;
	Size res1;
	std::string name;
	Size res2;
	Size res3;
	std::string func_type;
	data >> res1 >> name >> res2 >> res3 >> func_type;
	TR.Info << "read: " << res1 << " "<< name << " constrain to residues " << res2 << ":" << res3  << " func: " << func_type << std::endl;
	func::FuncOP aFunc = func_factory.new_func( func_type );
	aFunc->read_data( data );

	if ( TR.Debug.visible() ) {
		aFunc->show_definition( TR.Debug ); TR.Debug<<std::endl;
	}
	setup_csts( res1, name, res2,res3, pose, aFunc );

	if ( data.good() ) {
		//chu skip the rest of line since this is a single line defintion.
		while ( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( TR.Debug.visible() ) {
		aFunc->show_definition( TR.Debug );
		TR.Debug << std::endl;
	}

} // read_def

void
SiteConstraintResidues::setup_csts(
	Size res1,
	std::string name,
	Size res2,
	Size res3,
	core::pose::Pose const & pose,
	func::FuncOP const & func
) {
	id::AtomID target_atom( pose.residue_type( res1 ).atom_index( name ), res1 );
	//Size target_chain = pose.chain( res );
	for ( Size j = res2 ; j < res3 ; ++j ) {
		id::AtomID atom2( pose.residue_type( j ).atom_index( "CA" ), j );
		runtime_assert( target_atom.valid() && atom2.valid() );
		add_individual_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( target_atom, atom2, func ) ) ) );
	}
} // setup_csts

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::SiteConstraintResidues::save( Archive & arc ) const {
	arc( cereal::base_class< AmbiguousConstraint >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::SiteConstraintResidues::load( Archive & arc ) {
	arc( cereal::base_class< AmbiguousConstraint >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::SiteConstraintResidues );
CEREAL_REGISTER_TYPE( core::scoring::constraints::SiteConstraintResidues )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_SiteConstraintResidues )
#endif // SERIALIZATION
