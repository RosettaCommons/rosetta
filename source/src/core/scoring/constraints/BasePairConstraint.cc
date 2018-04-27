// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief A wrapper for a very particular AmbiguousConstraint of MultiConstraints
/// @author Andrew Watkins (amw579@stanford.edu, October 2016)

// Unit headers
#include <core/scoring/constraints/BasePairConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// Package headers
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/RNA_Info.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <core/pose/Pose.hh>
#include <core/id/SequenceMapping.hh>


// Project headers
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer tr( "core.scoring.constraints" );

namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
BasePairConstraint::BasePairConstraint():
	Constraint( base_pair_constraint )
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
BasePairConstraint::BasePairConstraint( Size const res1, Size const res2,
	core::chemical::rna::BaseEdge const edge1,
	core::chemical::rna::BaseEdge const edge2,
	core::chemical::rna::LW_BaseDoubletOrientation const orientation ):
	Constraint( base_pair_constraint ),
	res1_( res1 ),
	res2_( res2 ),
	edge1_( edge1 ),
	edge2_( edge2 ),
	orientation_( orientation )
{
}

ConstraintOP
BasePairConstraint::clone() const {
	auto bpc = BasePairConstraintOP( new BasePairConstraint( res1_, res2_, edge1_, edge2_, orientation_ ) );
	core::scoring::constraints::ConstraintCOPs constraints;
	for ( auto const & cst : constraints_ ) {
		constraints.push_back( cst->clone() );
	}
	bpc->set_constraints( constraints );
	return ConstraintOP( bpc ); //new BasePairConstraint( res1_, res2_, edge1_, edge2_, orientation_ ) );
}

Real
BasePairConstraint::dist( core::scoring::func::XYZ_Func const & /*xyz*/ ) const {
	return 0;
}

void BasePairConstraint::show_def( std::ostream& out, pose::Pose const& /*pose*/ ) const {
	out << "Base pair between bases " << res1_ << " and " << res2_ << std::endl;
}

ConstraintOP BasePairConstraint::remapped_clone( pose::Pose const& /*src*/, pose::Pose const& /*dest*/, id::SequenceMappingCOP smap ) const {

	if ( !smap ) {
		return clone();
	}
	Size const res1 = (*smap)[ res1_ ];
	Size const res2 = (*smap)[ res2_ ];

	if ( res1 == 0 || res2 == 0 ) return nullptr;
	return ConstraintOP( new BasePairConstraint( res1, res2, edge1_, edge2_, orientation_ ) );
}

ConstraintOP
BasePairConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	Size const res1 = seqmap[ res1_ ];
	Size const res2 = seqmap[ res2_ ];
	if ( seqmap[ res1_ ] != 0 && seqmap[ res2_ ] != 0 ) {
		return ConstraintOP( new BasePairConstraint( res1, res2, edge1_, edge2_, orientation_ ) );
	} else {
		return nullptr;
	}
}

std::string BasePairConstraint::type() const {
	return "BasePairConstraint";
}

bool
BasePairConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me(     *this ) ) return false;

	auto const & other( static_cast< BasePairConstraint const & > (other_cst) );
	if ( res1_ != other.res1_ ) return false;
	if ( res2_ != other.res2_ ) return false;
	if ( edge1_ != other.edge1_ ) return false;
	if ( edge2_ != other.edge2_ ) return false;
	if ( orientation_ != other.orientation_ ) return false;

	return true;
}

bool
BasePairConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< BasePairConstraint const * > (&other);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints but only reports the lowest one
void
BasePairConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	//std::cout << "Hello! You have made it to BasePairConstraint::score. There are " << constraints_.size() << " constraints." << std::endl;
	for ( auto const & cst : constraints_ ) {
		cst->score( xyz_func, weights, emap );
	}
} //score

/// @brief function to minimize lowest scoring member constraint
void
BasePairConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {
	// for each contained
	for ( auto const & cst : constraints_ ) {
		//std::cout << ", type of cst is " << (*cst_it)->score_type() << ",  ";
		cst->fill_f1_f2(atom, xyz, F1, F2, weights);
	}
}

void
BasePairConstraint::show( std::ostream& out) const
{
	/// APL -- you cannot show the active constraint in the absence of a Pose.
	out << "BasePairConstraint between " << res1_ << " and " << res2_ << std::endl;
	out << "There are " << constraints_.size() << " constraints." << std::endl;
}

Size
BasePairConstraint::show_violations( std::ostream& /*out*/, pose::Pose const& /*pose*/, Size /*verbose_level*/, Real /*threshold*/ ) const
{
	return 0; // STUBBED OUT!
}

inline
core::chemical::rna::BaseEdge get_edge_from_char( char const c ) {
	if ( c == 'W' ) return core::chemical::rna::WATSON_CRICK;
	if ( c == 'H' ) return core::chemical::rna::HOOGSTEEN;
	if ( c == 'S' ) return core::chemical::rna::SUGAR;
	utility_exit_with_message( "Impossible character specified for edge in BasePairConstraint definition: \'"+std::string(1,c)+"\'" );
}

void
BasePairConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & //func_factory
) {
	std::string tempres1, tempres2;
	char edgechar1, edgechar2, orientationchar;
	data >> tempres1 >> tempres2 >> edgechar1 >> edgechar2 >> orientationchar;

	edge1_ = get_edge_from_char( edgechar1 );
	edge2_ = get_edge_from_char( edgechar2 );
	if ( orientationchar == 'C' ) {
		orientation_ = core::chemical::rna::CIS;
	} else if ( orientationchar == 'T' ) {
		orientation_ = core::chemical::rna::TRANS;
	} else if ( orientationchar == 'A' ) {
		orientation_ = core::pose::rna::get_LW_orientation( edge1_, edge2_, core::chemical::rna::ANTIPARALLEL );
	} else if ( orientationchar == 'P' ) {
		orientation_ = core::pose::rna::get_LW_orientation( edge1_, edge2_, core::chemical::rna::PARALLEL );
	} else {
		utility_exit_with_message( "Impossible character specified for orientation in BasePairConstraint definition: \'"+std::string(1,orientationchar)+"\'" );
	}

	ConstraintIO::parse_residue( pose, tempres1, res1_ );
	ConstraintIO::parse_residue( pose, tempres2, res2_ );

	tr.Debug << "read: " << res1_ << " " << res2_ << std::endl;

	if ( res1_ > pose.size() || res2_ > pose.size() || res1_ == 0 || res2_ == 0 ) {
		tr.Warning  << "ignored constraint (requested residue numbers exceed numbers of residues in pose): " << "Total in Pose: " << pose.size() << " "
			<< res1_ << " " << res2_ << std::endl;
		tr.Warning  << "Also possible that you provided resnum-chain information,"
			<< " and the residue isn't found in the pose so the mapping returned zero." << std::endl;
		// We can't set this. It nukes parent constraints that might just not be 'right yet'
		//data.setstate( std::ios_base::failbit );
		return;
	}
	// At one time this errored if the bases weren't RNA. Now, we will error only if they do not
	// have a base_analogue.
	using namespace core::chemical;
	if ( ( !pose.residue_type( res1_ ).is_RNA() || !pose.residue_type( res2_ ).is_RNA() )
			&& ( pose.residue_type( res1_ ).base_analogue() == aa_unp || pose.residue_type( res2_ ).base_analogue() == aa_unp ) ) {
		tr.Warning  << "ignored constraint (requested residue numbers may not be RNA): " << "Total in Pose: " << pose.size() << " " << res1_ << " " << res2_ << std::endl;
		// We can't set this. It nukes parent constraints that might just not be 'right yet'
		//data.setstate( std::ios_base::failbit );
		return;
	}

	init_subsidiary_constraints( pose );
}

void
BasePairConstraint::init_subsidiary_constraints( core::pose::Pose const & pose ) {
	// OK, now set up all the AtomPairConstraints and everything...
	// NOTE: currently we use the same method rna_helix does for defining the
	// constraints, but we could easily accept functions from the constraint definition
	// too.


	// AMW TODO: enable non-WC specifications. Right now it'll just do as much with the
	// WC donors and acceptors as it can.
	// How do we handle chemically modified residues? Great question. We don't add
	// constraints on atoms that don't exist, but we don't totally rule out interaction.
	// So UR3 (3-methyl uridine) can still make a "Watson-Crick base pair" with
	// one of its carbonyls. It's just lame-o.

	using namespace core::scoring::func;
	Real const WC_distance( 1.9 );
	Real const distance_stddev( 0.25 );
	FuncOP const distance_func( new HarmonicFunc( WC_distance, distance_stddev ) );

	// Need to force base pairing -- not base stacking!
	Real const C1prime_distance( 10.5 );
	Real const C1prime_distance_stddev( 1.0 ); //Hmm. Maybe try linear instead?
	FuncOP const C1prime_distance_func( new HarmonicFunc( C1prime_distance, C1prime_distance_stddev ) );

	// Since we now support non-NA residues that have base_analogues, we can't
	// assume RNA_info() -- so condition on that.
	Size const atom1 = pose.residue_type( res1_ ).is_RNA() ?
		pose.residue_type( res1_ ).RNA_info().c1prime_atom_index() :
		pose.residue_type( res1_ ).is_TNA() ?
		pose.residue_type( res1_ ).atom_index( "C1'" ) :
		pose.residue_type( res1_ ).atom_index( "CE2" ); // only PNA supported
	Size const atom2 = pose.residue_type( res2_ ).is_RNA() ?
		pose.residue_type( res2_ ).RNA_info().c1prime_atom_index() :
		pose.residue_type( res2_ ).is_TNA() ?
		pose.residue_type( res2_ ).atom_index( "C1'" ) :
		pose.residue_type( res2_ ).atom_index( "CE2" ); // ditto
	constraints_.emplace_back( new AtomPairConstraint( id::AtomID( atom1, res1_ ), id::AtomID( atom2, res2_ ), C1prime_distance_func, base_pair_constraint ) );

	utility::vector1< std::string > atom_ids1, atom_ids2;
	using namespace chemical::rna;
	get_base_pair_atoms( pose.residue_type( res1_ ), pose.residue_type( res2_ ), atom_ids1, atom_ids2, edge1_, edge2_, orientation_ );

	for ( Size p = 1; p <= atom_ids1.size(); p++ ) {

		Size const atom1 = pose.residue_type( res1_ ).atom_index( atom_ids1[ p ] ) ;
		Size const atom2 = pose.residue_type( res2_ ).atom_index( atom_ids2[ p ] ) ;

		// AMW TODO: check to see each RT has the atom in question (only important for chemical
		// modifications) and accordingly skip. After the loop, if no constraints have been added,
		// warn as a result that you might be adding a Dumb Constraint. (Or maybe not! After all,
		// C1' distances have been restrained.)

		tr.Debug << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue_type(res1_).name1() << res1_ << " <-->  " <<
			pose.residue_type(res2_).name1() << res2_ << "   " <<
			atom_ids1[ p ] << " <--> " <<
			atom_ids2[ p ] << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;

		constraints_.emplace_back( new AtomPairConstraint( id::AtomID(atom1,res1_), id::AtomID(atom2,res2_), distance_func, base_pair_constraint ) );
	}

}

BasePairConstraint::BasePairConstraint( BasePairConstraint const & src ) :
	Constraint( src ),
	res1_( src.res1_ ),
	res2_( src.res2_ )
{}

} //constraints
} //scoring
} //core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::BasePairConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( constraints_ ) ); // ConstraintCOPs
	arc( CEREAL_NVP( res1_ ) ); // Size
	arc( CEREAL_NVP( res2_ ) ); // Size
	arc( CEREAL_NVP( edge1_ ) ); // enum
	arc( CEREAL_NVP( edge2_ ) ); // enum
	arc( CEREAL_NVP( orientation_ ) ); // enum
	arc( CEREAL_NVP( member_atoms_ ) ); // utility::vector1<AtomID>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::BasePairConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( constraints_ ); // ConstraintCOPs
	arc( res1_ ); // Size
	arc( res2_ ); // Size
	arc( edge1_ ); // enum
	arc( edge2_ ); // enum
	arc( orientation_ ); // enum
	arc( member_atoms_ ); // utility::vector1<AtomID>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::BasePairConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::BasePairConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_BasePairConstraint )
#endif // SERIALIZATION
