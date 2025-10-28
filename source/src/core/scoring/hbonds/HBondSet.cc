// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBondSet.hh
/// @brief  Hydrogen bond set class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/HBondSet.hh>


// Package headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

// Utility headers
#include <utility>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

#include <numeric/xyz.functions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <map>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace hbonds {

/// @details Auto-generated virtual destructor
HBond::~HBond() = default;

using namespace ObjexxFCL::format;

static basic::Tracer t( "core.scoring.hbonds.HBondSet" );

HBond::HBond(
	Size const dhatm,
	bool const dhatm_is_protein_backbone,
	bool const dres_is_protein,
	bool const dres_is_dna,
	bool const dhatm_is_backbone,
	Size const dres,
	Size const aatm,
	bool const aatm_is_protein_backbone,
	bool const ares_is_protein,
	bool const ares_is_dna,
	bool const aatm_is_backbone,
	Size const ares,
	HBEvalTuple const & hbe_tuple,
	Real const energy_in, // unweighted
	Real const weight_in,
	HBondDerivs const & derivs_in
):
	don_hatm_( dhatm ),
	don_hatm_is_protein_backbone_( dhatm_is_protein_backbone ),
	don_res_is_protein_( dres_is_protein ),
	don_res_is_dna_( dres_is_dna ),
	don_hatm_is_backbone_( dhatm_is_backbone ),
	don_res_( dres ),
	acc_atm_( aatm ),
	acc_atm_is_protein_backbone_( aatm_is_protein_backbone ),
	acc_res_is_protein_( ares_is_protein ),
	acc_res_is_dna_( ares_is_dna ),
	acc_atm_is_backbone_( aatm_is_backbone ),
	acc_res_( ares ),
	eval_tuple_( hbe_tuple ),
	energy_( energy_in ),
	weight_( weight_in ),
	don_npd_weight_( 1.0 ),
	acc_npd_weight_( 1.0 ),
	derivs_( derivs_in ),
	index_( 0 ),
	don_index_( 0 ),
	acc_index_( 0 )
{}


Size
HBond::don_res() const
{
	return don_res_;
}


Size
HBond::don_hatm() const
{
	return don_hatm_;
}

/// needed for silly allow logic
bool
HBond::don_hatm_is_protein_backbone() const
{
	return don_hatm_is_protein_backbone_;
}

bool
HBond::don_res_is_protein() const
{
	return don_res_is_protein_;
}

bool
HBond::don_res_is_dna() const
{
	return don_res_is_dna_;
}

/// needed for silly allow logic
bool
HBond::don_hatm_is_backbone() const
{
	return don_hatm_is_backbone_;
}


Size
HBond::acc_res() const
{
	return acc_res_;
}


Size
HBond::acc_atm() const
{
	return acc_atm_;
}

/// needed for silly allow logic
bool
HBond::acc_atm_is_protein_backbone() const
{
	return acc_atm_is_protein_backbone_;
}

bool
HBond::acc_res_is_protein() const
{
	return acc_res_is_protein_;
}

bool
HBond::acc_res_is_dna() const
{
	return acc_res_is_dna_;
}

/// needed for silly allow logic
bool
HBond::acc_atm_is_backbone() const
{
	return acc_atm_is_backbone_;
}

/// NOTE: this is unweighted energy, see weight() for the weight
Real
HBond::energy() const
{
	return energy_;
}


Real
HBond::weight() const
{
	return weight_;
}

Real
HBond::don_npd_weight() const
{
	return don_npd_weight_;
}

Real
HBond::acc_npd_weight() const
{
	return acc_npd_weight_;

}

void
HBond::don_npd_weight( Real setting )
{
	don_npd_weight_ = setting;
}

void
HBond::acc_npd_weight( Real setting )
{
	acc_npd_weight_ = setting;
}

Size HBond::index() const
{
	return index_;
}

void HBond::index( Size setting )
{
	index_ = setting;
}

Size HBond::don_index() const
{
	return don_index_;
}

void HBond::don_index( Size setting )
{
	don_index_ = setting;
}

Size HBond::acc_index() const
{
	return acc_index_;
}

void HBond::acc_index( Size setting )
{
	acc_index_ = setting;
}



HBondDerivs const &
HBond::derivs() const
{
	return derivs_;
}


HBEvalType
HBond::eval_type() const
{
	return eval_tuple_.eval_type();
}


HBEvalTuple const &
HBond::eval_tuple() const
{
	return eval_tuple_;
}


bool
HBond::atom_is_donorH( id::AtomID const & atom ) const
{
	return ( atom.rsd() == don_res_ && atom.atomno() == don_hatm_ );
}


bool
HBond::atom_is_acceptor( id::AtomID const & atom ) const
{
	return ( atom.rsd() == acc_res_ && atom.atomno() == acc_atm_ );
}

///////////////////////////////////////////////////////////////////////////////
bool
HBond::hbond_energy_comparer( HBondCOP a, HBondCOP b )
{
	return ( a->energy() * a->weight() < b->energy() * b->weight() );
}

///////////////////////////////////////////////////////////////////////////////

Real
HBond::get_AHDangle(core::pose::Pose const & pose) const {

	PointPosition const Dxyz(pose.residue(don_res()).xyz(pose.residue(don_res()).atom_base(don_hatm())));
	PointPosition const Hxyz(pose.residue(don_res()).xyz(don_hatm()));
	PointPosition const Axyz(pose.residue(acc_res()).xyz(acc_atm()));
	return numeric::angle_degrees(Axyz, Hxyz, Dxyz);

}

Real
HBond::get_BAHangle(core::pose::Pose const & pose) const {

	PointPosition const Bxyz(pose.residue(acc_res()).xyz(pose.residue(acc_res()).atom_base(acc_atm())));
	PointPosition const Hxyz(pose.residue(don_res()).xyz(don_hatm()));
	PointPosition const Axyz(pose.residue(acc_res()).xyz(acc_atm()));
	return numeric::angle_degrees(Bxyz, Axyz, Hxyz);
}

Real
HBond::get_BAtorsion(core::pose::Pose const & pose) const {

	Size const base_atomno(pose.residue(acc_res()).atom_base(acc_atm()));
	Size const bbase_atomno(pose.residue(acc_res()).atom_base(base_atomno));
	PointPosition const Bxyz(pose.residue(acc_res()).xyz(base_atomno));
	PointPosition const BBxyz(pose.residue(acc_res()).xyz(bbase_atomno));
	PointPosition const Hxyz(pose.residue(don_res()).xyz(don_hatm()));
	PointPosition const Axyz(pose.residue(acc_res()).xyz(acc_atm()));

	return dihedral_degrees(BBxyz, Bxyz, Axyz, Hxyz);
}

Real
HBond::get_HAdist(core::pose::Pose const & pose) const {

	PointPosition const Hxyz(pose.residue(don_res()).xyz(don_hatm()));
	PointPosition const Axyz(pose.residue(acc_res()).xyz(acc_atm()));
	return Hxyz.distance(Axyz);

}


void
HBond::show( std::ostream & out ) const
{
	out << eval_tuple_.eval_type() << " ";
	out << "don: ";
	out << (don_res_is_protein_ ? "protein " : "");
	out << (don_hatm_is_protein_backbone_ ? "backbone " : "");
	out << (don_res_is_dna_ ? "dna " : "");
	out << don_res_ << " ";
	out << don_hatm_ << " ";
	out << "acc: ";
	out << (acc_res_is_protein_ ? "protein " : "");
	out << (acc_atm_is_backbone_ ? "backbone " : "");
	out << (acc_res_is_dna_ ? "dna " : "");
	out << acc_res_ << " ";
	out << acc_atm_ << " ";
	out << energy_ << " ";
	out << weight_ << std::endl;
}

void
HBond::show(
	pose::Pose const & pose,
	bool print_header,
	std::ostream & out
) const {

	if ( print_header ) {
		out
			<< "#Dch Dn Dres Da  Ach An Ares Aa  length AHDang BAHang  BAtor weight energy"
			<< std::endl;
	}

	std::string const don_pdb_chain(pose.pdb_info()->chain(don_res()));
	Size const don_pdb_number(pose.pdb_info()->number(don_res()));
	char const don_pdb_icode(pose.pdb_info()->icode(don_res()));
	Size const don_atomno(pose.residue(don_res()).atom_base(don_hatm()));
	std::string const don_res_name(pose.residue(don_res()).type().name3());
	std::string const don_atom_name(pose.residue(don_res()).atom_name(don_atomno));

	std::string const acc_pdb_chain(pose.pdb_info()->chain(acc_res()));
	Size const acc_pdb_number(pose.pdb_info()->number(acc_res()));
	char const acc_pdb_icode(pose.pdb_info()->icode(acc_res()));
	std::string const acc_res_name(pose.residue(acc_res()).type().name3());
	std::string const acc_atom_name(pose.residue(acc_res()).atom_name(acc_atm()));

	Real const HAdist = get_HAdist(pose);
	Real const AHDangle = get_AHDangle(pose);
	Real const BAHangle = get_BAHangle(pose);
	Real const BAtorsion = get_BAtorsion(pose);

	// In the pdb format, a residue is uniquely identified by
	// the chain, the residue number and the insertion code
	out
		<< "#"

		// donor site
		<< A(1, don_pdb_chain) << I(5, don_pdb_number) << A(1, don_pdb_icode)
		<< I(4, don_res_name) << I(5, don_atom_name)

		// acceptor site
		<< A(1, acc_pdb_chain) << I(5, acc_pdb_number)  << A(1, acc_pdb_icode)
		<< I(4, acc_res_name) << I(5, acc_atom_name)

		// bond geometry:
		<< F(6, 2, HAdist)

		// These angles are go from [0, 180] in degrees
		<< F(7, 1, AHDangle) << F(7, 1, BAHangle)

		// The axis of rotation is around the base-acceptor bond and the
		// zero angle in the plane of the base of the base of the
		// acceptor. [-180, 180] in degrees.
		<< F(7, 1, BAtorsion)

		// The hbond score function weight and energy
		<< F(7, 3, weight()) << F(7, 3, energy()) << std::endl;
}

void HBond::show(pose::Pose const & pose, bool const print_header) const
{
	show(pose, print_header, t);
}

std::ostream &
operator<< ( std::ostream & out, const HBond & hbond ){
	hbond.show( out );
	return out;
}


bool
operator==(HBond const & a, HBond const & b)
{
	return (
		a.don_hatm_is_protein_backbone_ == b.don_hatm_is_protein_backbone_&&
		a.don_res_is_protein_           == b.don_res_is_protein_          &&
		a.don_res_is_dna_               == b.don_res_is_dna_              &&
		a.don_hatm_is_backbone_         == b.don_hatm_is_backbone_        &&
		a.don_res_                      == b.don_res_                     &&
		a.acc_atm_                      == b.acc_atm_                     &&
		a.acc_atm_is_protein_backbone_  == b.acc_atm_is_protein_backbone_ &&
		a.acc_res_is_protein_           == b.acc_res_is_protein_          &&
		a.acc_res_is_dna_               == b.acc_res_is_dna_              &&
		a.acc_atm_is_backbone_          == b.acc_atm_is_backbone_         &&
		a.acc_res_                      == b.acc_res_                     &&
		a.eval_tuple_                   == b.eval_tuple_                  &&
		a.energy_                       == b.energy_                      &&
		a.weight_                       == b.weight_                      );
}


HBondSet::HBondSet():
	options_( utility::pointer::make_shared< HBondOptions >() ),
	atom_map_init_( false )
{}

HBondSet::~HBondSet() = default;

HBondSet::HBondSet( Size const nres ):
	options_( utility::pointer::make_shared< HBondOptions >() ),
	atom_map_init_( false )
{
	resize_bb_donor_acceptor_arrays( nres );
}


// Note that the following constructors are duplicated to prevent
// automatic casting.  Since owning-ptr is castible to Size and the
// above constructor takes a Size, automatic casting would
// ambiguous.
HBondSet::HBondSet( HBondOptions const & opts ):
	options_( utility::pointer::make_shared< HBondOptions >( opts ) ),
	atom_map_init_( false )
{}

HBondSet::HBondSet( HBondOptions const & opts, Size const nres):
	options_( utility::pointer::make_shared< HBondOptions >( opts ) ),
	atom_map_init_( false )
{
	resize_bb_donor_acceptor_arrays( nres );
}

/// @brief Convenience constructor. If you need more controlled
///construction please use one of the other constructors
///
/// The pose must be non-const because the neighbor graph may need to
/// be initialized.
HBondSet::HBondSet(
	pose::Pose & pose,
	bool const bb_only /*true*/) :
	options_(utility::pointer::make_shared< HBondOptions >()),
	atom_map_init_(false)
{
	pose.update_residue_neighbors();
	setup_for_residue_pair_energies(pose, false, bb_only);
}

/// @brief Convenience constructor. If you need more controlled
///construction please use one of the other constructors
///
/// The pose must be non-const because the neighbor graph may need to
/// be initialized.
HBondSet::HBondSet(
	HBondOptions const & opts,
	pose::Pose & pose,
	bool const bb_only /*true*/) :
	options_( utility::pointer::make_shared< HBondOptions >(opts)),
	atom_map_init_(false)
{
	pose.update_residue_neighbors();
	setup_for_residue_pair_energies(pose, false, bb_only);
}

/// @brief Convenience constructor: Find all the hbonds between two residues. BB only default.
/// @details A bit inefficient, internally, since this creates an intermediate two-residue pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
HBondSet::HBondSet(
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2,
	core::scoring::hbonds::HBondDatabase const &database
) :
	options_(utility::pointer::make_shared< HBondOptions >()),
	atom_map_init_(false)
{
	core::scoring::hbonds::identify_hbonds_1way( database, res1, res2, 1, 1, false, false, false, false, false, *this );
	core::scoring::hbonds::identify_hbonds_1way( database, res2, res1, 1, 1, false, false, false, false, false, *this );
}

/// @brief Convenience constructor: Find all the hbonds between two residues. BB only default.
/// @details A bit inefficient, internally, since this creates an intermediate two-residue pose.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
HBondSet::HBondSet(
	HBondOptions const & options,
	core::conformation::Residue const &res1,
	core::conformation::Residue const &res2,
	core::scoring::hbonds::HBondDatabase const &database
) :
	options_( utility::pointer::make_shared< HBondOptions >(options) ),
	atom_map_init_(false)
{
	core::scoring::hbonds::identify_hbonds_1way( database, res1, res2, 1, 1, false, false, false, false, false, *this );
	core::scoring::hbonds::identify_hbonds_1way( database, res2, res1, 1, 1, false, false, false, false, false, *this );
}


/// copy ctor
HBondSet::HBondSet( HBondSet const & src ):
	basic::datacache::CacheableData( src ),
	options_( utility::pointer::make_shared< HBondOptions >( *src.options_ ))
{
	for ( Size i=1; i<= src.hbonds_.size(); ++i ) {
		hbonds_.push_back( utility::pointer::make_shared< HBond >( *src.hbonds_[i] ) );
	}
	backbone_backbone_donor_ = src.backbone_backbone_donor_;
	backbone_backbone_acceptor_ = src.backbone_backbone_acceptor_;
	nbrs_ = src.nbrs_;
	// set this flag since we still have to setup the atommap
	atom_map_init_ = false;
}

/// copy ctor
HBondSet::HBondSet( HBondSet const & src, utility::vector1< core::Size > exclude_list ):
	basic::datacache::CacheableData(),
	options_( utility::pointer::make_shared< HBondOptions >( *src.options_ ))
{
	//bool exclude=false;
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		int const dres = hbond.don_res();
		int const ares = hbond.acc_res();
		bool exclude = false;
		for ( Size k = 1; k <= exclude_list.size(); k ++ ) {
			if ( exclude_list[k] == (core::Size)dres ) { exclude = true; break; };
			if ( exclude_list[k] == (core::Size)ares ) { exclude = true; break; };
		}
		if ( exclude ) continue;

		hbonds_.push_back( utility::pointer::make_shared< HBond >( *src.hbonds_[i] ) );
	}
	backbone_backbone_donor_ = src.backbone_backbone_donor_;
	backbone_backbone_acceptor_ = src.backbone_backbone_acceptor_;
	nbrs_ = src.nbrs_;
	// set this flag since we still have to setup the atommap
	atom_map_init_ = false;
}

/// copy ctor
HBondSet::HBondSet( HBondSet const & src, utility::vector1< bool > residue_mask ):
	basic::datacache::CacheableData(),
	options_( utility::pointer::make_shared< HBondOptions >( *src.options_ ))
{
	//bool exclude=false;
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		if ( residue_mask[hbond.don_res()] && residue_mask[hbond.acc_res()] ) {
			hbonds_.push_back( utility::pointer::make_shared< HBond >( *src.hbonds_[i] ) );
		}
	}
	backbone_backbone_donor_ = src.backbone_backbone_donor_;
	backbone_backbone_acceptor_ = src.backbone_backbone_acceptor_;
	nbrs_ = src.nbrs_;
	// set this flag since we still have to setup the atommap
	atom_map_init_ = false;
}


/// copy ctor
HBondSet::HBondSet(
	HBondSet const & src,
	Size seqpos
):
	basic::datacache::CacheableData(),
	options_( utility::pointer::make_shared< HBondOptions >( *src.options_ ))
{
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		if ( hbond.don_res() == seqpos || hbond.acc_res() == seqpos ) {
			hbonds_.push_back( utility::pointer::make_shared< HBond >( *src.hbonds_[i] ));
		}
	}
	backbone_backbone_donor_ = src.backbone_backbone_donor_;
	backbone_backbone_acceptor_ = src.backbone_backbone_acceptor_;
	nbrs_ = src.nbrs_;
	// set this flag since we still have to setup the atommap
	atom_map_init_ = false;
}


/// @brief clone this object
basic::datacache::CacheableDataOP
HBondSet::clone() const
{
	return utility::pointer::make_shared< HBondSet >( *this );
}


// for accessing the nbrs, allows hacky non-tenA_nbr count
int
HBondSet::nbrs( Size const seqpos ) const
{
	return nbrs_[ seqpos ];
}

void
HBondSet::set_nbrs( Size const seqpos, Size value )
{
	nbrs_[ seqpos ] = value;
}


/// \brief  Access hbond
HBond const &
HBondSet::hbond( Size const number ) const
{
	return *(hbonds_[ number ]);
}

HBondCOP
HBondSet::hbond_cop( Size const number ) const
{
	return hbonds_[ number ];
}

///////////////////////////////////////////////////////////////////////////////
/// \brief  Add a new hbond to the list
/// updates the "hbchk" array as necessary
void
HBondSet::append_hbond(
	Size const dhatm,
	conformation::Residue const & don_rsd,
	Size const aatm,
	conformation::Residue const & acc_rsd,
	HBEvalTuple const & hbe_tuple,
	Real const energy,
	Real const weight,
	HBondDerivs const & derivs
){
	atom_map_init_ = false;

	// note that the derivative may not have been evaluated
	bool const dhatm_is_protein_backbone
		( don_rsd.is_protein() && don_rsd.atom_is_backbone( dhatm ) );
	bool const aatm_is_protein_backbone
		( acc_rsd.is_protein() && acc_rsd.atom_is_backbone( aatm ) );

	bool const dhatm_is_backbone
		( don_rsd.atom_is_backbone( dhatm ) );
	bool const aatm_is_backbone
		( acc_rsd.atom_is_backbone( aatm ) );

	Size const don_pos( don_rsd.seqpos() );
	Size const acc_pos( acc_rsd.seqpos() );

	hbonds_.push_back( utility::pointer::make_shared< HBond >(
		dhatm, dhatm_is_protein_backbone, don_rsd.is_protein(), don_rsd.is_DNA(),
		dhatm_is_backbone, don_pos,
		aatm , aatm_is_protein_backbone, acc_rsd.is_protein(),
		acc_rsd.is_DNA(), aatm_is_backbone, acc_pos,
		hbe_tuple, energy, weight, derivs ) );
	hbonds_[ hbonds_.size() ]->index( hbonds_.size() );

	// update hbcheck
	// note: these dimension checks could be removed completely & replaced by pose-aware dimensioning
	if ( options_->bb_donor_acceptor_check() /* actually we could save a LOT of computation if we implement the removal of bb_donor_acceptor_check correctly... */ ) {
		if ( dhatm_is_protein_backbone || aatm_is_protein_backbone ) {
			Size const max_pos( std::max( acc_pos, don_pos ) );
			if ( max_pos > backbone_backbone_donor_.size() ||
					max_pos > backbone_backbone_acceptor_.size() ) {
				// this could be more efficient if we require specifying nres ahead of time
				// Specify nres in setup_for_residue_pair_energy
				resize_bb_donor_acceptor_arrays( max_pos );
			}
			if ( dhatm_is_protein_backbone && aatm_is_protein_backbone ) {
				backbone_backbone_donor_   [ don_pos ] = true;
				backbone_backbone_acceptor_[ acc_pos ] = true;
			}
		}
	}
}

/// \brief  Is this hbond allowed under the bb-bb exclusion scheme?
bool
HBondSet::allow_hbond( Size const index ) const
{
	return allow_hbond( *hbonds_[ index ] );
}

bool
HBondSet::allow_hbond( HBond const & hbond ) const
{
	// donor is backbone atom already making bb-bb hbond
	// acceptor is sidechain
	if ( hbond.don_hatm_is_protein_backbone() &&
			!hbond.acc_atm_is_protein_backbone() &&
			backbone_backbone_donor_[ hbond.don_res() ] ) return false;

	if ( hbond.acc_atm_is_protein_backbone() &&
			!hbond.don_hatm_is_protein_backbone() &&
			backbone_backbone_acceptor_[ hbond.acc_res() ] ) return false;

	return true;
}

/// @brief is the backbone bone acceptor group in a bb/bb hydrogen bond?
bool
HBondSet::acc_bbg_in_bb_bb_hbond( Size const residue ) const
{
	return backbone_backbone_acceptor_[ residue ];
}

/// @brief is the backbone bone donor group in a bb/bb hydrogen bond?
bool
HBondSet::don_bbg_in_bb_bb_hbond( Size const residue ) const
{
	return backbone_backbone_donor_[ residue ];
}

/// \brief  Number of hbonds
Size
HBondSet::nhbonds() const
{
	return hbonds_.size();
}

// Return the vector of HBond pointers -- super useful in PyRosetta
utility::vector1< HBondOP > const &
HBondSet::hbonds() const
{
	return hbonds_;
}

Size
HBondSet::nhbonds(Size const seqpos, bool include_only_allowed /* true */) const
{
	utility::vector1< HBondCOP > const bonds = residue_hbonds(seqpos, include_only_allowed);
	return bonds.size();
}

Size
HBondSet::nhbonds( AtomID const & atom, bool include_only_allowed /* true */) const
{
	utility::vector1< HBondCOP > const & bonds = atom_hbonds(atom, include_only_allowed);
	return bonds.size();
}

utility::vector1< HBondCOP >
HBondSet::atom_hbonds( AtomID const & atom, bool include_only_allowed /* true */ ) const
{
	setup_atom_map();
	utility::vector1< HBondCOP > result_bonds;
	HBondAtomMap::const_iterator iter( atom_map_.find( atom ) );
	if ( iter != atom_map_.end() ) {
		utility::vector1< HBondCOP > bonds = iter->second;
		for ( Size i = 1; i <= bonds.size(); ++i ) {
			if ( (include_only_allowed && allow_hbond(*bonds[i])) || !include_only_allowed ) {
				result_bonds.push_back(bonds[i]);
			} else continue;
		}
	}

	return result_bonds;
}

utility::vector1< HBondCOP > const &
HBondSet::atom_hbonds_all( AtomID const & atom ) const
{
	setup_atom_map();
	HBondAtomMap::const_iterator iter( atom_map_.find( atom ) );
	if ( iter != atom_map_.end() ) {
		return iter->second;
	} else {
		return empty_hbond_vector_;
	}
}

utility::vector1< HBondCOP >
HBondSet::residue_hbonds(const Size seqpos, bool include_only_allowed /* true */ ) const
{
	utility::vector1< HBondCOP > bonds;

	for ( Size i=1; i<= nhbonds(); ++i ) {
		HBondCOP bond( hbond_cop(i) );
		if ( (bond->don_res() == seqpos || bond->acc_res() == seqpos)
				&& ((include_only_allowed && allow_hbond(*bond)) || !include_only_allowed) ) {
			bonds.push_back( bond );
		}
	}
	return bonds;
}

/// @brief  Delete all the data
void
HBondSet::clear()
{
	atom_map_init_ = false;
	hbonds_.clear();
	backbone_backbone_donor_.clear();
	backbone_backbone_acceptor_.clear();
	nbrs_.clear();
}

/// \brief Resize bb info arrays
void
HBondSet::resize_bb_donor_acceptor_arrays( Size const new_dimension )
{
	backbone_backbone_donor_.resize( new_dimension, false );
	backbone_backbone_acceptor_.resize( new_dimension, false );
}

void
HBondSet::copy_bb_donor_acceptor_arrays( HBondSet const & src ) {
	backbone_backbone_donor_ = src.backbone_backbone_donor_;
	backbone_backbone_acceptor_ = src.backbone_backbone_acceptor_;
}

void
HBondSet::sort_by_weighted_energy()
{
	std::sort( hbonds_.begin(), hbonds_.end(), HBond::hbond_energy_comparer );
}

///////////////////////////////////////////////////////////////////////////////
/// \brief  Setup the mapping from atoms to hbonds
void
HBondSet::setup_atom_map() const // lazy updating
{
	if ( atom_map_init_ ) return;
	atom_map_init_ = true;
	for ( Size i=1; i<= hbonds_.size(); ++i ) {
		// You can get an OP without changing an object
		HBondOP hb( hbonds_[i] );
		AtomID const don_atom( hb->don_hatm(), hb->don_res() );
		AtomID const acc_atom( hb->acc_atm() , hb->acc_res() );
		for ( Size r=1; r<= 2; ++r ) {
			AtomID const & atom( r == 1 ? don_atom : acc_atom );
			auto iter( atom_map_.find( atom ) );
			if ( iter == atom_map_.end() ) {
				atom_map_.insert
					( std::make_pair( atom, utility::vector1< HBondCOP >() ));
				iter = atom_map_.find( atom );
			} // atom not already present
			iter->second.push_back( hb );

			if ( r == 1 ) hb->don_index( iter->second.size() );
			else          hb->acc_index( iter->second.size() );
		} // repeat for donor and acceptor
	} // i=1.nhbonds
}


void
HBondSet::setup_for_residue_pair_energies(
	pose::Pose const & pose,
	bool const calculate_derivative /*= false*/,
	bool const backbone_only /*= true*/
) {
	// now fill the hbond set with only bb-bb hbonds
	// in the process we fill the arrays that keep track of which protein bb groups are making a
	// bb-bb hbond
	//
	// sc-bb hbonds with these groups are disallowed if the
	// hbondoptions.bb_donor_acceptor_check_ boolean is true.

	SSWeightParameters ssdep;
	ssdep.ssdep_ = options_->length_dependent_srbb();
	ssdep.l_ = options_->length_dependent_srbb_lowscale();
	ssdep.h_ = options_->length_dependent_srbb_highscale();
	ssdep.len_l_ = options_->length_dependent_srbb_minlength();
	ssdep.len_h_ = options_->length_dependent_srbb_maxlength();

	HBondSet hbond_set(*this);
	fill_hbond_set(pose, calculate_derivative, *this, ssdep,
		false, backbone_only, backbone_only, backbone_only);

	// stash the nbr info
	TenANeighborGraph const & tenA_neighbor_graph
		( pose.energies().tenA_neighbor_graph() );

	nbrs_.resize( pose.size() );
	for ( Size i=1; i<= pose.size(); ++i ) {

		{ // hacky thing here until we update to current svn
			nbrs_[i] = 1;
			for ( utility::graph::Graph::EdgeListConstIter
					ir  = tenA_neighbor_graph.get_node( i )->const_edge_list_begin(),
					ire = tenA_neighbor_graph.get_node( i )->const_edge_list_end();
					ir != ire; ++ir ) {
				Size const neighbor_id( (*ir)->get_other_ind( i ) );
				// Don't include VRTs in neighbor count
				if ( pose.residue( neighbor_id ).aa() == core::chemical::aa_vrt ) continue;
				nbrs_[i] += 1;

			}
			//    std::cout << "nbrdiff: old= " << tenA_neighbor_graph.get_node(i)->num_neighbors_counting_self() << " new= " <<
			//     nbrs_[i] << std::endl;
		}
	}

	// ensure full size for the the HbondSet bb_acc/bb_don arrays
	// the bb arrays ought to be resized once at the beginning of this function
	// and not resized to 0 for later push_back operations.... such a waste.
	resize_bb_donor_acceptor_arrays( pose.size() );
}

HBondOptions const &
HBondSet::hbond_options() const {
	return *options_;
}

/// @brief set the hbond options for this hbond set; clears all hbonds already stored
void HBondSet::set_hbond_options( HBondOptions const & options )
{
	clear();
	options_ = utility::pointer::make_shared< HBondOptions >( options ); // make a copy of the input to guarantee its integrity
}

std::ostream &
operator<< ( std::ostream & out, const HBondSet & hbond_set ){
	hbond_set.show( out );
	return out;
}

void
HBondSet::show(
	std::ostream & out
) const {
	for ( core::Size i=1; i<=nhbonds(); ++i ) {
		out << hbond(i);
	}
}

// PyRosetta friendly version
void HBondSet::show() const { show(std::cout); }

void HBondSet::show(pose::Pose const & pose, bool const print_header) const
{
	show(pose, print_header, t);
}

void
HBondSet::show(
	pose::Pose const & pose,
	Size const residue,
	bool const print_header) const
{
	show(pose, residue, print_header, std::cout);
}



bool
operator==(HBondSet const & a, HBondSet const & b)
{
	if ( a.options_ != b.options_ ) return false;

	if ( a.atom_map_ != b.atom_map_ ) return false;
	if ( a.atom_map_init_ != b.atom_map_init_ ) return false;

	if ( a.backbone_backbone_donor_.size() != b.backbone_backbone_donor_.size() ) return false;
	for ( Size i=1; i<=a.backbone_backbone_donor_.size(); ++i ) {
		if ( a.backbone_backbone_donor_[i] != b.backbone_backbone_donor_[i] ) return false;
	}

	if ( a.backbone_backbone_acceptor_.size() != b.backbone_backbone_acceptor_.size() ) return false;
	for ( Size i=1; i<=a.backbone_backbone_donor_.size(); ++i ) {
		if ( a.backbone_backbone_acceptor_[i] != b.backbone_backbone_acceptor_[i] ) return false;
	}

	if ( a.nhbonds() != b.nhbonds() ) {
		t << "HBondSets not equal!" << std::endl;
		t << a << std::endl << b << std::endl;
		return false;
	}

	for ( Size i=1; i<=a.nhbonds(); ++i ) {
		if ( !(a.hbond(i) == b.hbond(i)) ) return false;
	}

	return true;
}

/// @detail Optionally print a header, and then a row for each hydrogen bond in the set using the iterprable version of the hbond show format formatted for easy parsing by R, Excel, etc.
void
HBondSet::show(
	pose::Pose const & pose,
	bool print_header/*=true*/,
	std::ostream & out/*=std::cout*/
) const {

	for ( Size i=1; i<= nhbonds(); ++i ) {
		hbond(i).show(pose, i==1 && print_header, out);
	}

}


/// @detail Optionally print a header, and then a row for each hydrogen
///bond in the set using the iterprable version of the hbond show
///format formatted for easy parsing by R, Excel, etc.
void
HBondSet::show(
	pose::Pose const & pose,
	Size const residue_number,
	bool const print_header /*=true*/,
	std::ostream & out /*=std::cout*/
) const {

	bool first_found(true);
	for ( Size i=1; i<=nhbonds(); ++i ) {
		if ( (hbond(i).don_res() == residue_number)
				|| (hbond(i).acc_res() == residue_number) ) {

			hbond(i).show(pose, print_header && first_found, out);
			first_found = false;
		}
	}
}

void HBondSet::hbond_don_npd_weight( Size index, Real setting )
{
	hbonds_[ index ]->don_npd_weight( setting );
}
void HBondSet::hbond_acc_npd_weight( Size index, Real setting )
{
	hbonds_[ index ]->acc_npd_weight( setting );
}

}
}
}


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::hbonds::HBond::HBond() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::HBond::save( Archive & arc ) const {
	arc( CEREAL_NVP( don_hatm_ ) ); // Size
	arc( CEREAL_NVP( don_hatm_is_protein_backbone_ ) ); // _Bool
	arc( CEREAL_NVP( don_res_is_protein_ ) ); // _Bool
	arc( CEREAL_NVP( don_res_is_dna_ ) ); // _Bool
	arc( CEREAL_NVP( don_hatm_is_backbone_ ) ); // _Bool
	arc( CEREAL_NVP( don_res_ ) ); // Size
	arc( CEREAL_NVP( acc_atm_ ) ); // Size
	arc( CEREAL_NVP( acc_atm_is_protein_backbone_ ) ); // _Bool
	arc( CEREAL_NVP( acc_res_is_protein_ ) ); // _Bool
	arc( CEREAL_NVP( acc_res_is_dna_ ) ); // _Bool
	arc( CEREAL_NVP( acc_atm_is_backbone_ ) ); // _Bool
	arc( CEREAL_NVP( acc_res_ ) ); // Size
	arc( CEREAL_NVP( eval_tuple_ ) ); // class core::scoring::hbonds::HBEvalTuple
	arc( CEREAL_NVP( energy_ ) ); // Real
	arc( CEREAL_NVP( weight_ ) ); // Real
	arc( CEREAL_NVP( don_npd_weight_ ) ); // Real
	arc( CEREAL_NVP( acc_npd_weight_ ) ); // Real
	arc( CEREAL_NVP( derivs_ ) ); // struct core::scoring::hbonds::HBondDerivs
	arc( CEREAL_NVP( index_ ) ); // Size
	arc( CEREAL_NVP( don_index_ ) ); // Size
	arc( CEREAL_NVP( acc_index_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::HBond::load( Archive & arc ) {
	arc( don_hatm_ ); // Size
	arc( don_hatm_is_protein_backbone_ ); // _Bool
	arc( don_res_is_protein_ ); // _Bool
	arc( don_res_is_dna_ ); // _Bool
	arc( don_hatm_is_backbone_ ); // _Bool
	arc( don_res_ ); // Size
	arc( acc_atm_ ); // Size
	arc( acc_atm_is_protein_backbone_ ); // _Bool
	arc( acc_res_is_protein_ ); // _Bool
	arc( acc_res_is_dna_ ); // _Bool
	arc( acc_atm_is_backbone_ ); // _Bool
	arc( acc_res_ ); // Size
	arc( eval_tuple_ ); // class core::scoring::hbonds::HBEvalTuple
	arc( energy_ ); // Real
	arc( weight_ ); // Real
	arc( don_npd_weight_ ); // Real
	arc( acc_npd_weight_ ); // Real
	arc( derivs_ ); // struct core::scoring::hbonds::HBondDerivs
	arc( index_ ); // Size
	arc( don_index_ ); // Size
	arc( acc_index_ ); // Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::HBond );
CEREAL_REGISTER_TYPE( core::scoring::hbonds::HBond )


/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::hbonds::HBondSet::save( Archive & arc ) const {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	arc( CEREAL_NVP( options_ ) ); // HBondOptionsCOP
	arc( CEREAL_NVP( hbonds_ ) ); // utility::vector1<HBondOP>
	arc( CEREAL_NVP( backbone_backbone_donor_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( backbone_backbone_acceptor_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( nbrs_ ) ); // utility::vector1<int>
	arc( CEREAL_NVP( atom_map_ ) ); // HBondAtomMap
	arc( CEREAL_NVP( atom_map_init_ ) ); // _Bool
	arc( CEREAL_NVP( empty_hbond_vector_ ) ); // utility::vector1<HBondCOP>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::hbonds::HBondSet::load( Archive & arc ) {
	arc( cereal::base_class< basic::datacache::CacheableData >( this ) );
	std::shared_ptr< core::scoring::hbonds::HBondOptions > local_options;
	arc( local_options ); // HBondOptionsCOP
	options_ = local_options; // copy the non-const pointer(s) into the const pointer(s)
	arc( hbonds_ ); // utility::vector1<HBondOP>
	arc( backbone_backbone_donor_ ); // utility::vector1<_Bool>
	arc( backbone_backbone_acceptor_ ); // utility::vector1<_Bool>
	arc( nbrs_ ); // utility::vector1<int>
	std::map< core::id::AtomID, utility::vector1< std::shared_ptr< core::scoring::hbonds::HBond > >, struct std::less< core::id::AtomID > > local_atom_map;
	arc( local_atom_map ); // HBondAtomMap
	// atom_map_ = local_atom_map; // copy the non-const pointer(s) into the const pointer(s)
	for ( auto iter : local_atom_map ) {
		atom_map_[ iter.first ] = iter.second;
	}

	arc( atom_map_init_ ); // _Bool
	utility::vector1< std::shared_ptr< core::scoring::hbonds::HBond > > local_empty_hbond_vector;
	arc( local_empty_hbond_vector ); // utility::vector1<HBondCOP>
	empty_hbond_vector_ = local_empty_hbond_vector; // copy the non-const pointer(s) into the const pointer(s)
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::hbonds::HBondSet );
CEREAL_REGISTER_TYPE( core::scoring::hbonds::HBondSet )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_hbonds_HBondSet )
#endif // SERIALIZATION
