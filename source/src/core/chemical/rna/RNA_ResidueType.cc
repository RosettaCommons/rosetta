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

// ObjexxFCL headers
#include <core/chemical/ResidueSupport.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>


namespace core {
namespace chemical {
namespace rna {

using namespace ObjexxFCL;

static basic::Tracer tr("core.chemical.rna.RNA_ResidueType");

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
	o1p_atom_index_( 0 ),
	o2p_atom_index_( 0 ),
	o5prime_index_( 0 ),
	o3prime_index_( 0 ),
	o4prime_index_( 0 ),
	c1prime_index_( 0 ),
	c2prime_index_( 0 ),
	c4prime_index_( 0 )
{
	base_atom_list_.clear();
	Is_RNA_base_atom_list_.clear();
	Is_phosphate_atom_list_.clear();
}


////////////////////////////////////////////////////////////////////////
RNA_ResidueType::~RNA_ResidueType(){}


////////////////////////////////////////////////////////////////////////

void
RNA_ResidueType::update_derived_rna_data(ResidueTypeCOP const residue_type_in){

	residue_type_=residue_type_in;

	//std::cout << "finalizing rna_residue_type_" << std::endl;


	Is_virtual_atom_list_.clear();

	for ( Size i=1; i<= residue_type_->natoms(); i++ ) {
		if(residue_type_->is_virtual( i ) ){
			Is_virtual_atom_list_.push_back(true);
		}else{
			Is_virtual_atom_list_.push_back(false);	
		}
	}

	if(Is_virtual_atom_list_.size()!=residue_type_->natoms()) utility_exit_with_message("Is_virtual_atom_list_.size()!=residue_type_->natoms()");


	o2prime_index_=residue_type_->atom_index( " O2'" );
	ho2prime_index_=residue_type_->atom_index( "HO2'" );

	p_atom_index_  =residue_type_->atom_index( " P  " );
	o1p_atom_index_=residue_type_->atom_index( " OP2" );
	o2p_atom_index_=residue_type_->atom_index( " OP1" );
	o5prime_index_  =residue_type_->atom_index( " O5'" );
	o3prime_index_  =residue_type_->atom_index( " O3'" );

	o4prime_index_=residue_type_->atom_index( " O4'" );
	c1prime_index_=residue_type_->atom_index( " C1'" );
	c2prime_index_=residue_type_->atom_index( " C2'" );
	c4prime_index_=residue_type_->atom_index( " C4'" );

	base_atom_list_.clear(); 

	for ( Size i=1; i<= residue_type_->natoms(); i++ ) {
		/*assume that chi # 1 is the base chi. chi_atoms_in[ 1 ][ 3 ] is the first base atom (either N1 or N9)*/
		if(residue_type_->last_controlling_chi( i )==1 || ( i == residue_type_->chi_atoms( 1 )[ 3 ] ) ){ 	
			base_atom_list_.push_back( i );
		} 
	}

	Is_RNA_base_atom_list_.clear();
	Is_phosphate_atom_list_.clear();

	for ( Size i=1; i<= residue_type_->natoms(); i++ ) {

		bool Is_phosphate_atom=false;
		if(i == p_atom_index_  ) Is_phosphate_atom=true;
		if(i == o1p_atom_index_) Is_phosphate_atom=true;
		if(i == o2p_atom_index_) Is_phosphate_atom=true;
		if(i == o5prime_index_  ) Is_phosphate_atom=true;
		if(i == o3prime_index_  ) Is_phosphate_atom=true;

		Is_phosphate_atom_list_.push_back(Is_phosphate_atom);

		/*assume that chi # 1 is the base chi, chi_atoms_in[ 1 ][3] is the index of the first RNA_base atom */
		if( ( residue_type_->last_controlling_chi( i ) )==1 || ( i == residue_type_->chi_atoms( 1 )[3] ) ){
			Is_RNA_base_atom_list_.push_back(true);
		}else{
			Is_RNA_base_atom_list_.push_back(false);
		}

	}

	if(Is_RNA_base_atom_list_.size()!=residue_type_->natoms()) utility_exit_with_message("Is_RNA_base_atom_list_.size()!=residue_type_->natoms()");

	if(Is_phosphate_atom_list_.size()!=residue_type_->natoms()) utility_exit_with_message("Is_phosphate_atom_list_.size()!=residue_type_->natoms()");

			
	if( (residue_type_->atoms_last_controlled_by_chi(1).size()+1)!= base_atom_list_.size() ){
		std::cout << "residue_type_->atoms_last_controlled_by_chi(1).size()=" << residue_type_->atoms_last_controlled_by_chi(1).size() << std::endl;
		std::cout << "base_atom_list_.size()=" << base_atom_list_.size() << std::endl;
		utility_exit_with_message( "(residue_type_->atoms_last_controlled_by_chi(1).size()+1)!= base_atom_list_.size()" );
	}

}


////////////////////////////////////////////////////////////
///WARNING THIS FUNCTION SHOULD NOT ACCESS ANY DATA of the RNA_ResidueType object itself since at this point it is not yet updated!
///ALSO SHOULD MAKE THIS FUNCTION A CONST FUNCTION!
void
RNA_ResidueType::rna_note_chi_controls_atom( Size const chi, Size const atomno, utility::vector1< core::Size >  & last_controlling_chi){


	//RNA require special treatment: Parin Sripakdeevong, Juen 26, 2011
	//1) there are 4 chi torsions.
	//2)The chi are not ordered from nearest to furthest from mainchain like in Protien. The sidechain also contain 2 branch (RNA base and 2'OH sidechains)
	//3)chi_2 is NEAREST. Then branch to chi_3 and chi_4. then chi_3 connects to chi_1 which IS FURTHEST!

	if(4!= residue_type_->nchi()) utility_exit_with_message("is_RNA_==true but nchi()=" +ObjexxFCL::string_of(residue_type_->nchi())+ "!=4");
	if(chi>residue_type_->nchi()) utility_exit_with_message("chi>residue_type_->nchi()");
	if(chi<1)  utility_exit_with_message("chi<1");

	utility::vector1< Size > chi_order;
	chi_order.push_back(1000);  //chi_1 is furthest
	chi_order.push_back(1);     //chi_2 is nearest
	chi_order.push_back(100);   //chi_3 is 3rd furthest (the choice between chi_4 and chi 3 is somewhat arbitrary)
	chi_order.push_back(10);    //chi_4 is 2nd furthest (the choice between chi_4 and chi 3 is somewhat arbitrary)

	if(last_controlling_chi[ atomno ]!=0){
		if( chi_order[ last_controlling_chi[ atomno ] ] > chi_order[ chi ] ) return;
	}

	//Stuck in a infinite recursion? But what about ring molecule like RNA base, can't a atom have two different base_atom is that case?
	//[In effect this condition ensure that the recursion is called for each atom in residue at most once for each chi] since
	//if called the second time then last_controlling_chi[ atomno ] is already set to chi in the prior call!
	if( last_controlling_chi[ atomno ] == chi ) utility_exit_with_message("last_controlling_chi[ atomno ] == chi");

	/*
	// This condition doesn't apply to RNA case since the chi are not ordered from nearest to furthest from mainchain like in Protien
	if( last_controlling_chi[ atomno ] != 0 ){
		std::cout << "last_controlling_chi[ atomno ] != 0" << std::endl;
		std::cout << "atomno=" << atomno << std::endl;
		std::cout << "chi=" << chi << std::endl;
		std::cout << "last_controlling_chi[ atomno ]=" << last_controlling_chi[ atomno ] << std::endl;
		utility_exit_with_message("last_controlling_chi[ atomno ] != 0 ");
	}
	*/

	last_controlling_chi[ atomno ] = chi;

	AtomIndices const & nbrs( residue_type_->bonded_neighbor( atomno ) );
	for ( Size ii = 1; ii <= nbrs.size(); ++ii ) {
		if ( residue_type_->atom_base( nbrs[ ii ] ) == atomno ) {
			rna_note_chi_controls_atom( chi, nbrs[ ii ], last_controlling_chi );
		}
	}

}

////////////////////////////////////////////////////////////
///WARNING THIS FUNCTION SHOULD NOT ACCESS ANY DATA of the RNA_ResidueType object itself since at this point it is not yet updated!
///ALSO SHOULD MAKE THIS FUNCTION A CONST FUNCTION!
void
RNA_ResidueType::rna_update_last_controlling_chi(ResidueTypeCOP const residue_type_in,
																						  utility::vector1< core::Size >  & last_controlling_chi,
																							utility::vector1< AtomIndices > & atoms_last_controlled_by_chi){

	residue_type_=residue_type_in;

	last_controlling_chi.clear(); //figure this out
	atoms_last_controlled_by_chi.clear(); //figure this out!

	//RNA require special treatment: Parin Sripakdeevong, Juen 26, 2011

	Size const nchi=residue_type_->nchi();

	last_controlling_chi.resize( residue_type_->natoms() );

	std::fill( last_controlling_chi.begin(), last_controlling_chi.end(), 0 );

	for ( Size ii = nchi; ii >= 1; --ii ) {
		/// Note children of atom 3 of chi_ii as being controlled by chi ii.
		Size const iiat3 = residue_type_->chi_atoms( ii )[ 3 ];
		// This may be unnecessary; I believe two atoms pair as each other's bases only at the mainchain.
		Size const iiat3base = residue_type_->atom_base( iiat3 );
		AtomIndices const & ii_nbrs( residue_type_->bonded_neighbor( iiat3 ) );
		for ( Size jj = 1; jj <= ii_nbrs.size(); ++jj ) {
			Size const jj_atom = ii_nbrs[ jj ];
			if ( residue_type_->atom_base( jj_atom ) == iiat3 && iiat3base != jj_atom ) {
				rna_note_chi_controls_atom( ii, jj_atom, last_controlling_chi );
			}
		}
	}

	/// Now compute the atoms_last_controlled_by_chi_ arrays.

	/// get ready to allocate space in the atoms_last_controlled_by_chi_ arrays
	utility::vector1< Size > natoms_for_chi( nchi, 0 );
	for ( Size ii = 1; ii <= residue_type_->natoms(); ++ii ) {
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
	for ( Size ii = 1; ii <= residue_type_->natoms(); ++ii ) {
		if ( last_controlling_chi[ ii ] != 0 ) {
			atoms_last_controlled_by_chi[ last_controlling_chi[ ii ]].push_back( ii );
		}
	}

}

////////////////////////////////////////////////////////////

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
utility::vector1< bool > const &
RNA_ResidueType::Is_virtual_atom_list() const
{
	//runtime_assert( finalized_ );

	return Is_virtual_atom_list_;

}

bool
RNA_ResidueType::atom_is_virtual( Size const atomno ) const
{
	//runtime_assert( finalized_ ); 
	//runtime_assert(Is_virtual_atom_list_.size()==natoms_);

	return (Is_virtual_atom_list_[atomno]);
}


utility::vector1< bool > const &
RNA_ResidueType::Is_phosphate_atom_list() const
{
	//runtime_assert( finalized_ );

	return Is_phosphate_atom_list_;

}

bool
RNA_ResidueType::atom_is_phosphate( Size const atomno ) const
{
	//runtime_assert( finalized_ );
	//runtime_assert( atomno <= natoms_ );

	return (Is_phosphate_atom_list_[atomno]);

}

utility::vector1< bool > const &
RNA_ResidueType::Is_RNA_base_atom_list() const
{
	//runtime_assert( finalized_ );

	return Is_RNA_base_atom_list_;

}

//True for heavy, hydrogen and virtual atoms as long as it is in the RNA_base!
bool
RNA_ResidueType::is_RNA_base_atom( Size const atomno ) const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//runtime_assert( atomno <= natoms_ );

	/*assume that chi # 1 is the base chi, chi_atoms_[ 1 ][3] is the index of the first RNA_base atom */

	return (Is_RNA_base_atom_list_[ atomno ]);
}

//Note that this implement return both heavy and hydrogen atoms in the RNA_base. 
//Virtual atoms in the RNA_base atoms are also returned.
AtomIndices const &
RNA_ResidueType::RNA_base_atoms() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );

	/*
	if( (atoms_last_controlled_by_chi_[1].size()+1)!= base_atom_list_.size() ){
		std::cout << "atoms_last_controlled_by_chi_[1].size()=" << atoms_last_controlled_by_chi_[1].size() << std::endl;
		std::cout << "base_atom_list_.size()=" << base_atom_list_.size() << std::endl;
		utility_exit_with_message( "(atoms_last_controlled_by_chi_[1].size()+1)!= base_atom_list_.size() " );
	}
	*/

	return base_atom_list_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::ho2prime_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(ho2prime_index_==atom_index( "HO2'" )); //remove after testing!

	return ho2prime_index_;
}


//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o2prime_index() const
{

	//////////////////////////////
	//residue_type_->require_final();
	//runtime_assert(o2prime_index_==residue_type_->atom_index( " O2'" )); //remove after testing!
	//	std::cout << "--------------------------------------" << std::endl;
	//	std::cout << "Inside RNA_ResidueType::o2prime_index()" << std::endl;
	//	std::cout << "name_from_aa(residue_type_->aa())=" << name_from_aa(residue_type_->aa())<< std::endl;
	//	std::cout << "o2prime_index_=" << o2prime_index_ << std::endl;
	//	std::cout << "residue_type_->atom_name(o2prime_index_)=" << residue_type_->atom_name(o2prime_index_) << std::endl;
	//	std::cout << "residue_type_->first_sidechain_atom() =" << residue_type_->first_sidechain_atom()  << std::endl;
	//	std::cout << "residue_type_->atom_name(residue_type_->first_sidechain_atom())=" << residue_type_->atom_name(residue_type_->first_sidechain_atom()) << std::endl;
	//	std::cout << "--------------------------------------" << std::endl;
	//////////////////////////////

	return o2prime_index_;
}


//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::p_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(p_atom_index_==atom_index( " P  " )); //remove after testing!

	return p_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o1p_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(o1p_atom_index_==atom_index( " OP2" )); //remove after testing!

	return o1p_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o2p_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(o2p_atom_index_==atom_index( " OP1" )); //remove after testing!
	return o2p_atom_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o5prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(o5prime_index_==atom_index( " O5'" )); //remove after testing!
	return o5prime_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o3prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ ) ;
	//if(rna_check_) runtime_assert(o3prime_index_==atom_index( " O3'" )); //remove after testing!
	return o3prime_index_;
}

//For fast lookup! Parin Sripakdeevong, June 25th, 2011
Size
RNA_ResidueType::o4prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(o4prime_index_==atom_index( " O4'" )); //remove after testing!
	return o4prime_index_;
}

Size
RNA_ResidueType::c1prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(c1prime_index_==atom_index( " C1'" )); //remove after testing!
	return c1prime_index_;
}

Size
RNA_ResidueType::c2prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(c2prime_index_==atom_index( " C2'" )); //remove after testing!
	return c2prime_index_;
}

Size
RNA_ResidueType::c4prime_atom_index() const
{
	//if(is_RNA_==false) utility_exit_with_message("is_RNA_==false");
	//runtime_assert( finalized_ );
	//if(rna_check_) runtime_assert(c4prime_index_==atom_index( " C4'" )); //remove after testing!
	return c4prime_index_;
}
////////////////////////////////RNA specific stuff...maybe move function into its own class?///////////

} // rna
} // chemical
} // core
