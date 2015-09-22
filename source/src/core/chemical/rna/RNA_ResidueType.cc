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
/// @author Parin Sripakdeevong (sripakpa@stanford.edu)

// Unit headers
#include <core/chemical/rna/RNA_ResidueType.hh>
#include <core/chemical/ResidueType.hh>

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>


namespace core {
namespace chemical {
namespace rna {

using namespace ObjexxFCL;

static THREAD_LOCAL basic::Tracer TR( "core.chemical.rna.RNA_ResidueType" );

////////////////////////////////////////////////////////////

//Could perhaps find an atom index in term of its chemical position along the chain.
//For example could call chi_atoms( chi_torsion_# ) to get the atom_index of the atoms that define this torsion and then deduce the atom_index of the desired atom.
//Problem is that there is not analogous function for the backbone atoms! Parin S. June 25, 2011.
//Can get main_chain atoms from mainchain_atoms_?
//I think: mainchain_atoms_[1] is p_atom_index_
//         mainchain_atoms_[2] is o5prime_atom_index_
//         mainchain_atoms_[3] is c5prime_atom_index_
//         mainchain_atoms_[4] is c4prime_atom_index_
//         mainchain_atoms_[5] is o3prime_atom_index_
//         mainchain_atoms_[6] is o4prime_atom_index_

////////////////////////////////////////////////////////////////////////
RNA_ResidueType::RNA_ResidueType():
	o2prime_index_( 0 ),
	ho2prime_index_( 0 ),
	p_atom_index_( 0 ),
	op2_atom_index_( 0 ),
	op1_atom_index_( 0 ),
	o5prime_index_( 0 ),
	o3prime_index_( 0 ),
	o4prime_index_( 0 ),
	c1prime_index_( 0 ),
	c2prime_index_( 0 ),
	c4prime_index_( 0 )
{
	base_atom_list_.clear();
	is_RNA_base_.clear();
	is_phosphate_.clear();
}


////////////////////////////////////////////////////////////////////////
RNA_ResidueType::~RNA_ResidueType(){}


////////////////////////////////////////////////////////////////////////

void
RNA_ResidueType::update_derived_rna_data( ResidueTypeCAP residue_type_in ){

	residue_type_ = residue_type_in;
	ResidueTypeCOP residue_type( residue_type_ );

	is_virtual_.clear();
	for ( Size i = 1; i <= residue_type->natoms(); i++ ) {
		if ( residue_type->is_virtual( i ) ) {
			is_virtual_.push_back( true );
		} else {
			is_virtual_.push_back( false );
		}
	}
	runtime_assert( is_virtual_.size() == residue_type->natoms() );

	o2prime_index_ = residue_type->atom_index( " O2'" );
	ho2prime_index_ = residue_type->atom_index( "HO2'" );

	p_atom_index_  = residue_type->atom_index( " P  " );
	op2_atom_index_ = residue_type->atom_index( " OP2" );
	op1_atom_index_ = residue_type->atom_index( " OP1" );
	o5prime_index_  = residue_type->atom_index( " O5'" );
	o3prime_index_  = residue_type->atom_index( " O3'" );

	o4prime_index_ = residue_type->atom_index( " O4'" );
	c1prime_index_ = residue_type->atom_index( " C1'" );
	c2prime_index_ = residue_type->atom_index( " C2'" );
	c4prime_index_ = residue_type->atom_index( " C4'" );

	base_atom_list_.clear();

	for ( Size i = 1; i <= residue_type->natoms(); i++ ) {
		/*assume that chi # 1 is the base chi. chi_atoms_in[ 1 ][ 3 ] is the first base atom (either N1 or N9)*/
		if ( residue_type->last_controlling_chi( i ) == 1 || ( i == residue_type->chi_atoms( 1 )[ 3 ] ) ) {
			base_atom_list_.push_back( i );
		}
	}

	is_RNA_base_.clear();
	is_phosphate_.clear();

	for ( Size i = 1; i <= residue_type->natoms(); i++ ) {

		bool is_phosphate_atom = false;
		if ( i == p_atom_index_  ) is_phosphate_atom = true;
		if ( i == op2_atom_index_ ) is_phosphate_atom = true;
		if ( i == op1_atom_index_ ) is_phosphate_atom = true;
		if ( i == o5prime_index_  ) is_phosphate_atom = true;
		if ( i == o3prime_index_  ) is_phosphate_atom = true;

		is_phosphate_.push_back( is_phosphate_atom );

		/*assume that chi # 1 is the base chi, chi_atoms_in[ 1 ][3] is the index of the first RNA_base atom */
		if ( ( residue_type->last_controlling_chi( i ) ) == 1 || ( i == residue_type->chi_atoms( 1 )[3] ) ) {
			is_RNA_base_.push_back( true );
		} else {
			is_RNA_base_.push_back( false );
		}

	}

	runtime_assert( is_RNA_base_.size() == residue_type->natoms() );
	runtime_assert( is_phosphate_.size() == residue_type->natoms() );

	if ( ( residue_type->atoms_last_controlled_by_chi( 1 ).size() + 1 ) != base_atom_list_.size() ) {
		std::cout << "residue_type_->atoms_last_controlled_by_chi( 1 ).size() = " << residue_type->atoms_last_controlled_by_chi( 1 ).size() << std::endl;
		std::cout << "base_atom_list_.size() = " << base_atom_list_.size() << std::endl;
		utility_exit_with_message( "( residue_type_->atoms_last_controlled_by_chi( 1 ).size() + 1 ) != base_atom_list_.size()" );
	}

	chi_number_pseudoalpha_ = 0;
	chi_number_pseudobeta_ = 0;
	chi_number_pseudogamma_ = 0;
	chi_number_pseudoepsilon_ = 0;
	chi_number_pseudozeta_ = 0;

	for ( Size ii = 1; ii <= residue_type->nchi(); ii++ ) {
		std::string const chi_atom_name = residue_type->atom_name( residue_type->chi_atoms( ii )[ 4 ] );
		if ( chi_atom_name == "XO5'" ) {
			chi_number_pseudogamma_ = ii;
		} else if ( chi_atom_name == "XP  " ) {
			chi_number_pseudobeta_ = ii;
		} else if ( chi_atom_name == "XO3'" ) {
			chi_number_pseudoalpha_ = ii;
		} else if ( chi_atom_name == "YP  " ) {
			chi_number_pseudoepsilon_ = ii;
		} else if ( chi_atom_name == "YO5'" ) {
			chi_number_pseudozeta_ = ii;
		}
	}

}


////////////////////////////////////////////////////////////
// ugh, sorry this was hard-coded by Parin.
// should not be that hard to fill in automatically, though.
utility::vector1< Size > const
RNA_ResidueType::figure_out_chi_order() const {

	utility::vector1< Size > chi_order;
	chi_order.push_back( 4 );  //chi_1 is furthest [nucleosidic chi]
	chi_order.push_back( 1 );  //chi_2 is nearest  [nu2]
	chi_order.push_back( 3 );  //chi_3 is 3rd furthest (the choice between chi_4 and chi 3 is somewhat arbitrary) [nu1]
	chi_order.push_back( 2 );  //chi_4 is 2nd furthest (the choice between chi_4 and chi 3 is somewhat arbitrary) [2'-OH]

	chi_order.push_back( 1 );  //chi_5 pseudo-gamma   for packable 5' phosphate
	chi_order.push_back( 2 );  //chi_6 pseudo-beta    for packable 5' phosphate
	chi_order.push_back( 3 );  //chi_7 pseudo-alpha   for packable 5' phosphate
	chi_order.push_back( 4 );  //chi_8 pseudo-epsilon for packable 3' phosphate [if 5' phosphate exists]
	chi_order.push_back( 5 );  //chi_9 pseudo-zeta    for packable 3' phosphate [if 3' phosphate exists]

	return chi_order;
}

////////////////////////////////////////////////////////////
///WARNING THIS FUNCTION SHOULD NOT ACCESS ANY DATA of the RNA_ResidueType object itself since at this point it is not yet updated!
///ALSO SHOULD MAKE THIS FUNCTION A CONST FUNCTION!
void
RNA_ResidueType::rna_note_chi_controls_atom( Size const chi, Size const atomno,
	utility::vector1< core::Size >  & last_controlling_chi,
	utility::vector1< Size > const & chi_order ){

	//RNA require special treatment: Parin Sripakdeevong, June 26, 2011
	//1) there are 4 chi torsions.
	//2)The chi are not ordered from nearest to furthest from mainchain like in Protein.
	//   The sidechain also contain 2 branch (RNA base and 2'OH sidechains)
	//3)chi_2 is NEAREST. Then branch to chi_3 and chi_4. then chi_3 connects to chi_1 which IS FURTHEST!
	ResidueTypeCOP residue_type( residue_type_ );
	runtime_assert( residue_type->nchi() >= 4 );
	runtime_assert ( chi <= residue_type->nchi() );
	runtime_assert ( chi >= 1 );

	if ( last_controlling_chi[ atomno ] != 0 ) {
		if ( chi_order[ last_controlling_chi[ atomno ] ] > chi_order[ chi ] ) return;
	}

	//Stuck in a infinite recursion? But what about ring molecule like RNA base, can't a atom have two different base_atom is that case?
	//[In effect this condition ensure that the recursion is called for each atom in residue at most once for each chi] since
	//if called the second time then last_controlling_chi[ atomno ] is already set to chi in the prior call!
	runtime_assert ( last_controlling_chi[ atomno ] != chi );

	last_controlling_chi[ atomno ] = chi;

	AtomIndices const & nbrs( residue_type->bonded_neighbor( atomno ) );
	for ( Size ii = 1; ii <= nbrs.size(); ++ii ) {
		if ( residue_type->atom_base( nbrs[ ii ] ) == atomno ) {
			rna_note_chi_controls_atom( chi, nbrs[ ii ],
				last_controlling_chi, chi_order );
		}
	}

}


////////////////////////////////////////////////////////////
///WARNING THIS FUNCTION SHOULD NOT ACCESS ANY DATA of the RNA_ResidueType object itself since at this point it is not yet updated!
///ALSO SHOULD MAKE THIS FUNCTION A CONST FUNCTION!
void
RNA_ResidueType::rna_update_last_controlling_chi( ResidueTypeCAP residue_type_in,
	utility::vector1< core::Size >  & last_controlling_chi,
	utility::vector1< AtomIndices > & atoms_last_controlled_by_chi ){
	residue_type_ = residue_type_in;
	ResidueTypeCOP residue_type( residue_type_ );

	last_controlling_chi.clear(); //figure this out
	atoms_last_controlled_by_chi.clear(); //figure this out!

	//RNA requires special treatment: Parin Sripakdeevong, Juen 26, 2011
	Size const nchi = residue_type->nchi();
	last_controlling_chi.resize( residue_type->natoms() );
	std::fill( last_controlling_chi.begin(), last_controlling_chi.end(), 0 );

	utility::vector1< Size > const & chi_order = figure_out_chi_order();

	for ( Size ii = nchi; ii >= 1; --ii ) {
		/// Note children of atom 3 of chi_ii as being controlled by chi ii.
		if ( residue_type->chi_atoms( ii ).size() == 0 ) continue;
		Size const iiat3 = residue_type->chi_atoms( ii )[ 3 ]; // atom 3
		Size const iiat2 = residue_type->chi_atoms( ii )[ 2 ]; // atom 2

		Size const iiat3base = residue_type->atom_base( iiat3 ); // don't go back to base of atom.
		AtomIndices const & ii_nbrs( residue_type->bonded_neighbor( iiat3 ) );

		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( residue_type->atom_base( jj_atom ) == iiat3 &&
					iiat3base != jj_atom  &&
					iiat2 != jj_atom ) {
				rna_note_chi_controls_atom( ii, jj_atom,
					last_controlling_chi,
					chi_order );
			}
		}

	}

	/// Now compute the atoms_last_controlled_by_chi_ arrays.
	/// get ready to allocate space in the atoms_last_controlled_by_chi_ arrays
	utility::vector1< Size > natoms_for_chi( nchi, 0 );
	for ( Size ii = 1; ii <= residue_type->natoms(); ++ii ) {
		if ( last_controlling_chi[ ii ] != 0 ) {
			++natoms_for_chi[ last_controlling_chi[ ii ] ];
		}
	}

	/// allocate space
	atoms_last_controlled_by_chi.resize( nchi );
	for ( Size ii = 1; ii <= nchi; ++ii ) {
		atoms_last_controlled_by_chi[ ii ].clear();
		atoms_last_controlled_by_chi[ ii ].reserve( natoms_for_chi[ ii ] );
	}

	/// fill the arrays
	for ( Size ii = 1; ii <= residue_type->natoms(); ++ii ) {
		if ( last_controlling_chi[ ii ] != 0 ) {
			atoms_last_controlled_by_chi[ last_controlling_chi[ ii ]].push_back( ii );
		}
	}

	// if ( nchi > 8 ){
	//  TR << residue_type_->name() << std::endl;
	//  for ( Size ii = 1; ii <= nchi; ++ii ) {
	//   TR << "chi " << ii << " controls: ";
	//   for ( Size jj = 1; jj <= atoms_last_controlled_by_chi[ ii ].size(); jj++ ){
	//    TR << ' ' << residue_type_->atom_name( atoms_last_controlled_by_chi[ ii ][ jj ] );
	//   }
	//   TR << std::endl;
	//  }
	//  TR << std::endl;
	// }

}

////////////////////////////////////////////////////////////

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
utility::vector1< bool > const &
RNA_ResidueType::is_virtual() const
{
	return is_virtual_;

}

bool
RNA_ResidueType::atom_is_virtual( Size const atomno ) const
{
	return ( is_virtual_[ atomno ] );
}


utility::vector1< bool > const &
RNA_ResidueType::is_phosphate() const
{
	return is_phosphate_;
}

bool
RNA_ResidueType::atom_is_phosphate( Size const atomno ) const
{
	return ( is_phosphate_[ atomno ] );

}

utility::vector1< bool > const &
RNA_ResidueType::is_RNA_base() const
{
	return is_RNA_base_;

}

//True for heavy, hydrogen and virtual atoms as long as it is in the RNA_base!
bool
RNA_ResidueType::is_RNA_base_atom( Size const atomno ) const
{
	return ( is_RNA_base_[ atomno ] );
}

//Note that this implement return both heavy and hydrogen atoms in the RNA_base.
//Virtual atoms in the RNA_base atoms are also returned.
AtomIndices const &
RNA_ResidueType::RNA_base_atoms() const
{
	return base_atom_list_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::ho2prime_index() const
{
	return ho2prime_index_;
}


//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o2prime_index() const
{
	return o2prime_index_;
}


//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::p_atom_index() const
{
	return p_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::op2_atom_index() const
{
	return op2_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::op1_atom_index() const
{
	return op1_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o5prime_atom_index() const
{
	return o5prime_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o3prime_atom_index() const
{
	return o3prime_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o4prime_atom_index() const
{
	return o4prime_index_;
}

Size
RNA_ResidueType::c1prime_atom_index() const
{
	return c1prime_index_;
}

Size
RNA_ResidueType::c2prime_atom_index() const
{
	return c2prime_index_;
}

Size
RNA_ResidueType::c4prime_atom_index() const
{
	return c4prime_index_;
}
////////////////////////////////RNA specific stuff...maybe move function into its own class?///////////

} // rna
} // chemical
} // core
