// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/PoseMetricCalculators/NcontactsCalculator.cc
/// @brief calculate #atom contacts divided by pose.size() in a given pose
/// @details If two atom, each of which belongs to different residues, are within a
/// given distance, condist_( default 6.0A ), the two atoms are defined as contact pair.
/// The residues the atoms belong to are close ( isep_residue_ = 4 ) with each other along sequence,
/// the contact pair is ignored.
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit Headers
#include <protocols/fldsgn/NcontactsCalculator.hh>

//package headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>

// Project Headers
#include <basic/MetricValue.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

//// C++ headers
#include <cmath>

#include <core/chemical/AtomType.hh>
#include <utility/options/BooleanVectorOption.hh>


static basic::Tracer tr( "protocols.fldsgn.NcontactsCalculator" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace fldsgn {

class MyAtom {
public:
	typedef core::Size Size;
	MyAtom(){
		for ( Size j=1; j<= 50 ; ++j ) {
			atom_hydrophobic_.push_back(false);
		}
		atom_hydrophobic_[ 3 ] = true;  // CH1
		atom_hydrophobic_[ 4 ] = true;  // CH2
		atom_hydrophobic_[ 5 ] = true;  // CH3
		atom_hydrophobic_[ 6 ] = true;  // aroC
		atom_hydrophobic_[ 19 ] = true; // CAbb
	}
	bool is_hydrophobic( Size const index ){
		return atom_hydrophobic_[ index ];
	}
private:
	utility::vector1< bool > atom_hydrophobic_;
}; // MyAtom

/// @brief default constructor
NcontactsCalculator::NcontactsCalculator():
	condist_( 6.0 ),
	isep_residue_( 4 ),
	ignore_loops_( false ),
	ignore_same_sselement_( false ),
	ignore_same_sheet_( false ),
	use_only_calpha_( false )
{}

/// @brief default constructor
NcontactsCalculator::NcontactsCalculator(
	Real const cdist,
	Size const sep ):
	condist_( cdist ),
	isep_residue_( sep ),
	ignore_loops_( false ),
	ignore_same_sselement_( false ),
	ignore_same_sheet_( false ),
	use_only_calpha_( false )
{}

/// @brief copy constructor
NcontactsCalculator::NcontactsCalculator( NcontactsCalculator const & rval ):
	Super(),
	condist_( rval.condist_ ),
	isep_residue_( rval.isep_residue_ ),
	ignore_loops_( rval.ignore_loops_ ),
	ignore_same_sselement_( rval.ignore_same_sselement_ ),
	ignore_same_sheet_( rval.ignore_same_sheet_ ),
	use_only_calpha_( rval.use_only_calpha_ )
{}

/// @brief destructor
NcontactsCalculator::~NcontactsCalculator()= default;

/// @brief
void
NcontactsCalculator::lookup(
	String const & key,
	MetricValueBase * valptr
) const
{
	using namespace core;
	if ( key == "all_heavy_atm" ) {
		basic::check_cast( valptr, &nc_allatm_, "total_nlcontacts expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( nc_allatm_ );

	} else if ( key == "sidechain_heavy_apolar_atm" ) {
		basic::check_cast( valptr, &nc_hpatm_, "special_region_nlcontacts expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( nc_hpatm_ );

	} else if ( key == "sidechain_heavy_atm_apolar_aa" ) {
		basic::check_cast( valptr, &nc_hpres_, "special_region2_nlcontacts expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( nc_hpres_ );

	} else if ( key == "ss_entropy" ) {
		basic::check_cast( valptr, &ss_entrpy_, "special_region1_intra_nlcontacts expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( ss_entrpy_ );

	} else {
		tr.Fatal << "NcontactsCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}

} //lookup


// @brief
std::string
NcontactsCalculator::print( String const & key ) const
{
	String result;

	// #atom-contacts among all heavy atoms
	if ( key == "all_heavy_atm" ) {
		result = utility::to_string( nc_allatm_ );
		// #atom-contacts among hydrophobic heavy atoms
	} else if ( key == "sidechain_heavy_apolar_atm" ) {
		result = utility::to_string( nc_hpatm_ );
		//#atom-contacts among heavy atoms of sidechains of hydrophobic residues
	} else if ( key == "sidechain_heavy_atm_apolar_aa" ) {
		result = utility::to_string( nc_hpres_ );
	} else if ( key == "ss_entrpy" ) {
		result = utility::to_string( ss_entrpy_ );
	} else {
		basic::Error() << "NcontactsCalculator cannot compute metric " << key << std::endl;
	}
	return result;

} // apply

/// @brief recomute ncontacts
void
NcontactsCalculator::recompute( Pose const & pose )
{
	using protocols::fldsgn::topology::SheetSet;
	using protocols::fldsgn::topology::SheetSetOP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::StrandPairingSetOP;

	// intialize
	Size nres( pose.size() );
	Real condist2( numeric::square(condist_) );
	Size nc_allatm( 0 );
	Size nc_hpatm( 0 );
	Size nc_hpres( 0 );
	topology::SS_Info2_OP ssinfo( new topology::SS_Info2( pose.secstruct() ) );

	StrandPairingSetOP spairset( new StrandPairingSet( protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo ) ) );
	SheetSetOP sheet_set( new SheetSet( ssinfo, spairset ) );
	// std::cout << *sheet_set;

	Size max_ssele = ssinfo->ss_element_id( nres );
	utility::vector1< utility::vector1< Size > > ncon_sselements( max_ssele, (utility::vector1< Size >(max_ssele, 1)));
	utility::vector1< utility::vector1< bool > > calc_sselements( max_ssele, (utility::vector1< bool >(max_ssele, false)));

	// calc number of contacts
	MyAtom myatom;
	for ( Size iaa=1; iaa<=nres-isep_residue_; ++iaa ) {

		// skip calc if ss element is loop
		if ( ignore_loops_ && ssinfo->secstruct( iaa ) == 'L' ) continue;

		Size iaa_ssele( ssinfo->ss_element_id( iaa ) );

		for ( Size jaa=iaa+isep_residue_; jaa<=nres; ++jaa ) {

			// skip calc if ss element is loop
			if ( ignore_loops_ && ssinfo->secstruct( jaa ) == 'L' ) continue;

			Size jaa_ssele( ssinfo->ss_element_id( jaa ) );

			// skip calc pair if residues of pair belong to same ss elements
			if ( ignore_same_sselement_ && iaa_ssele == jaa_ssele ) continue;

			// skip calc pair if residues of pair belong to same ss elements
			if ( ignore_same_sheet_ &&
					ssinfo->secstruct( iaa ) == 'E' && ssinfo->secstruct( jaa ) == 'E' &&
					sheet_set->which_sheet( ssinfo->strand_id( iaa ) ) ==
					sheet_set->which_sheet( ssinfo->strand_id( jaa ) ) ) continue;

			Residue const & ires( pose.residue( iaa ) );
			for ( Size iatm=1, iatm_end=ires.natoms(); iatm<=iatm_end; ++iatm ) {

				if ( ires.atom_type( int(iatm) ).is_heavyatom() ) {
					Residue const & jres( pose.residue( jaa ) );

					for ( Size jatm=1, jatm_end=jres.natoms(); jatm<=jatm_end; ++jatm ) {

						if ( jres.atom_type(int(jatm)).is_heavyatom() ) {

							Atom const & iatom( ires.atom( iatm ) );
							Atom const & jatom( jres.atom( jatm ) );
							Size const iadex ( ires.atom_type_index( int(iatm) ) );
							Size const jadex ( jres.atom_type_index( int(jatm) ) );

							if ( (use_only_calpha_ &&
									ires.atom_type( iatm ).name() != "CAbb") || jres.atom_type( jatm ).name() != "CAbb" ) continue;

							calc_sselements[ iaa_ssele ][ jaa_ssele ] = true;

							Real const dsq( iatom.xyz().distance_squared( jatom.xyz() ));

							if ( dsq <= condist2 ) {
								nc_allatm ++;

								if ( myatom.is_hydrophobic(iadex) && myatom.is_hydrophobic(jadex) ) {
									nc_hpatm ++;

									ncon_sselements[iaa_ssele][jaa_ssele] ++ ;

									if ( !ires.atom_is_backbone( int(iatm) ) && !jres.atom_is_backbone( int(jatm) ) &&
											!ires.type().is_polar() && !jres.type().is_polar() ) {
										nc_hpres ++;
									}

									//std::cout << iaa  << " "  << jaa  << " " << iatm << " "  << jatm << " "
									//    << ires.name() << " " << jres.name() << " "
									//    << ssinfo->secstruct( iaa ) << " " << ssinfo->secstruct( jaa ) << " "
									//    << ires.atom_type( iatm ).name() << " " << jres.atom_type( jatm ).name() << " "
									//    << sheet_set->which_sheet( ssinfo->strand_id( iaa ) ) << " "
									//    << sheet_set->which_sheet( ssinfo->strand_id( jaa ) ) << " "
									//    << dsq << std::endl;

								}
							} // dsq

						} //fi jatm is_hydrogen ?
					} // loop for jatm
				} // fi iatm is_hydrogen ?
			} // loop for iatm

		} // loop for jaa
	}// loop for iaa

	// finalize
	Size tot( 0 );
	for ( Size i=1 ; i<=max_ssele; ++i ) {
		for ( Size j=1 ; j<=max_ssele; ++j ) {
			if ( calc_sselements[ i ][ j ] ) {
				tot += ncon_sselements[i][j];
			}
		}
	}
	ss_entrpy_ = 0.0;
	for ( Size i=1 ; i<=max_ssele; ++i ) {
		for ( Size j=1 ; j<=max_ssele; ++j ) {
			if ( calc_sselements[ i ][ j ] ) {
				Real prob = Real(ncon_sselements[i][j])/Real(tot);
				ss_entrpy_ += -prob*( std::log( prob )/std::log( 2.0 ) );
			}
		}
	}
	nc_allatm_ = Real(nc_allatm)/Real(nres);
	nc_hpatm_  = Real(nc_hpatm)/Real(nres);
	nc_hpres_  = Real(nc_hpres)/Real(nres);

} // recompute

} // fldsgn
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::fldsgn::NcontactsCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( condist_ ) ); // Real
	arc( CEREAL_NVP( isep_residue_ ) ); // Size
	arc( CEREAL_NVP( ignore_loops_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_same_sselement_ ) ); // _Bool
	arc( CEREAL_NVP( ignore_same_sheet_ ) ); // _Bool
	arc( CEREAL_NVP( use_only_calpha_ ) ); // _Bool
	arc( CEREAL_NVP( nc_allatm_ ) ); // Real
	arc( CEREAL_NVP( nc_hpatm_ ) ); // Real
	arc( CEREAL_NVP( nc_hpres_ ) ); // Real
	arc( CEREAL_NVP( ss_entrpy_ ) ); // Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::fldsgn::NcontactsCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( condist_ ); // Real
	arc( isep_residue_ ); // Size
	arc( ignore_loops_ ); // _Bool
	arc( ignore_same_sselement_ ); // _Bool
	arc( ignore_same_sheet_ ); // _Bool
	arc( use_only_calpha_ ); // _Bool
	arc( nc_allatm_ ); // Real
	arc( nc_hpatm_ ); // Real
	arc( nc_hpres_ ); // Real
	arc( ss_entrpy_ ); // Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::fldsgn::NcontactsCalculator );
CEREAL_REGISTER_TYPE( protocols::fldsgn::NcontactsCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_fldsgn_NcontactsCalculator )
#endif // SERIALIZATION
