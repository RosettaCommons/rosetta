// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/HBondSet.hh
/// @brief  Hydrogen bond set class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/scoring/hbonds/HBondSet.hh>


// Package headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TenANeighborGraph.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

// Utility headers
//#include <utility/exit.hh>
// #include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/conversions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++
#include <map> // what is the right header for std::pair ?

// option key includes
// AUTO-REMOVED #include <utility/options/keys/BooleanOptionKey.hh>

#include <utility/vector1.hh>




//#include <set>

namespace core {
namespace scoring {
namespace hbonds {

using namespace ObjexxFCL::fmt;

static basic::Tracer t("core.scoring.hbonds.HBondSet");

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
	HBEvalType const hbe_type,
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
	eval_type_( hbe_type ),
	energy_( energy_in ),
	weight_( weight_in ),
	derivs_( derivs_in )
{}

///
Size
HBond::don_res() const
{
	return don_res_;
}

///
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

///
Size
HBond::acc_res() const
{
	return acc_res_;
}

///
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

///
Real
HBond::weight() const
{
	return weight_;
}


HBondDerivs const &
HBond::derivs() const
{
	return derivs_;
}

///
HBEvalType const &
HBond::eval_type() const
{
	return eval_type_;
}

///
bool
HBond::atom_is_donorH( id::AtomID const & atom ) const
{
	return ( atom.rsd() == don_res_ && atom.atomno() == don_hatm_ );
}


///
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

///
void
HBond::show( std::ostream & out ) const
{
	out << eval_type_ << " ";
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
	//out << deriv_ << std::endl;
}

std::ostream &
operator<< ( std::ostream & out, const HBond & hbond ){
	hbond.show( out );
	return out;
}

///
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
		a.eval_type_                    == b.eval_type_                   &&
		a.energy_                       == b.energy_                      &&
		a.weight_                       == b.weight_                      ); /*&&
		a.deriv_.first[0]               == b.deriv_.first[0]              &&
		a.deriv_.first[1]               == b.deriv_.first[1]              &&
		a.deriv_.first[2]               == b.deriv_.first[2]              &&
		a.deriv_.second[0]              == b.deriv_.second[0]             &&
		a.deriv_.second[1]              == b.deriv_.second[1]             &&
		a.deriv_.second[2]              == b.deriv_.second[2]             );*/
}






HBondSet::HBondSet():
	options_( new HBondOptions() ),
	atom_map_init_( false )
{}

HBondSet::~HBondSet() {}

HBondSet::HBondSet( Size const nres ):
	options_( new HBondOptions() ),
	atom_map_init_( false )
{
	resize_bb_donor_acceptor_arrays( nres );
}


// Note that the following constructors are duplicated to prevent
// automatic casting.  Since owning-ptr is castible to Size and the
// above constructor takes a Size, automatic casting would
// ambiguous.
HBondSet::HBondSet( HBondOptions const & opts ):
	options_( new HBondOptions( opts ) ),
	atom_map_init_( false )
{}

HBondSet::HBondSet( HBondOptions const & opts, Size const nres):
	options_( new HBondOptions( opts ) ),
	atom_map_init_( false )
{
	resize_bb_donor_acceptor_arrays( nres );
}


/// copy ctor
HBondSet::HBondSet( HBondSet const & src ):
	basic::datacache::CacheableData(),
	options_( new HBondOptions( *src.options_ ))
{
	for ( Size i=1; i<= src.hbonds_.size(); ++i ) {
		hbonds_.push_back( new HBond( *src.hbonds_[i] ) );
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
	options_( new HBondOptions( *src.options_ ))
{

	bool exclude=false;
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		int const dres = hbond.don_res();
		int const ares = hbond.acc_res();
		exclude = false;
		for ( Size k = 1; k <= exclude_list.size(); k ++ ){
			if( exclude_list[k] == (core::Size)dres ){ exclude = true; break; };
			if( exclude_list[k] == (core::Size)ares ){ exclude = true; break; };
		}
		if( exclude ) continue;

		hbonds_.push_back( new HBond( *src.hbonds_[i] ) );
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
	options_( new HBondOptions( *src.options_ ))
{

	bool exclude=false;
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		if(residue_mask[hbond.don_res()] & residue_mask[hbond.acc_res()]){
			hbonds_.push_back( new HBond( *src.hbonds_[i] ) );
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
	options_( new HBondOptions( *src.options_ ))
{
	for ( Size i=1; i<= src.nhbonds(); ++i ) {
		HBond const & hbond(src.hbond(i));
		if (hbond.don_res() == seqpos || hbond.acc_res() == seqpos ){
			hbonds_.push_back( new HBond( *src.hbonds_[i] ));
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
	return new HBondSet( *this );
}

/// \brief  Number of hbonds
Size
HBondSet::nhbonds() const
{
	return hbonds_.size();
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

///////////////////////////////////////////////////////////////////////////////
/// \brief  Add a new hbond to the list
/// updates the "hbchk" array as necessary
void
HBondSet::append_hbond(
	Size const dhatm,
	conformation::Residue const & don_rsd,
	Size const aatm,
	conformation::Residue const & acc_rsd,
	HBEvalType const & hbe_type,
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

	hbonds_.push_back( new HBond(
		dhatm, dhatm_is_protein_backbone, don_rsd.is_protein(), don_rsd.is_DNA(),
		dhatm_is_backbone, don_pos,
		aatm , aatm_is_protein_backbone, acc_rsd.is_protein(),
		acc_rsd.is_DNA(), aatm_is_backbone, acc_pos,
		hbe_type, energy, weight, derivs ) );

	// update hbcheck
	// note: these dimension checks could be removed completely & replaced by pose-aware dimensioning
	if ( options_->bb_donor_acceptor_check() /* actually we could save a LOT of computation if we implement the removal of bb_donor_acceptor_check correctly... */ ){
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

/// \brief  Get a vector of all the hbonds involving this atom
utility::vector1< HBondCOP > const &
HBondSet::atom_hbonds( AtomID const & atom ) const
{
	setup_atom_map();
	HBondAtomMap::const_iterator iter( atom_map_.find( atom ) );
	if ( iter == atom_map_.end() ) {
		return empty_list_of_hbonds_;
	} else return iter->second;
}

/// \brief  Delete all the data
void
HBondSet::clear()
{
	atom_map_init_ = false;
	hbonds_.clear();
	backbone_backbone_donor_.clear();
	backbone_backbone_acceptor_.clear();
	nbrs_.clear();
}

/*
/// @details Leave the backbone_backbone_donor_ and backbone_backbone_acceptor_ arrays
/// intact
void
HBondSet::clear_hbonds()
{
	hbonds_.clear();
}
*/

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
		if ( !allow_hbond(i) ) continue;
		HBond const & hb( hbond(i) );
		AtomID const don_atom( hb.don_hatm(), hb.don_res() );
		AtomID const acc_atom( hb.acc_atm() , hb.acc_res() );
		for ( Size r=1; r<= 2; ++r ) {
			AtomID const & atom( r == 1 ? don_atom : acc_atom );
			HBondAtomMap::iterator iter( atom_map_.find( atom ) );
			if ( iter == atom_map_.end() ) {
				atom_map_.insert
					( std::make_pair( atom, utility::vector1< HBondCOP >() ));
				iter = atom_map_.find( atom );
			} // atom not already present
			iter->second.push_back( &hb );
		} // repeat for donor and acceptor
	} // i=1.nhbonds
}

///
void
HBondSet::setup_for_residue_pair_energies(
	pose::Pose const & pose,
	bool const calculate_derivative /*= false*/,
	bool const backbone_only /*= true*/
)
{

	// now fill the hbond set with only bb-bb hbonds
	// in the process we fill the arrays that keep track of which protein bb groups are making a
	// bb-bb hbond
	//
	// sc-bb hbonds with these groups are disallowed
	//

	HBondSet hbond_set(*this);
	fill_hbond_set(pose, calculate_derivative, *this,
		false, backbone_only, backbone_only, backbone_only);

	// stash the nbr info
	TenANeighborGraph const & tenA_neighbor_graph
		( pose.energies().tenA_neighbor_graph() );

	nbrs_.resize( pose.total_residue() );
	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		//nbrs_[i] = tenA_neighbor_graph.get_node(i)->num_neighbors_counting_self();

		{ // hacky thing here until we update to current svn
			nbrs_[i] = 1;
			for ( graph::Graph::EdgeListConstIter
				ir  = tenA_neighbor_graph.get_node( i )->const_edge_list_begin(),
				ire = tenA_neighbor_graph.get_node( i )->const_edge_list_end();
				ir != ire; ++ir ) {
				Size const neighbor_id( (*ir)->get_other_ind( i ) );
				// Don't include VRTs in neighbor count
				if ( pose.residue( neighbor_id ).aa() == core::chemical::aa_vrt ) continue;
				chemical::ResidueType const & nbr_rsd( pose.residue_type( neighbor_id ) );
				if ( nbr_rsd.is_protein() ) {
					nbrs_[i] += 1;
				} else if ( nbr_rsd.is_DNA() ) {
					nbrs_[i] += 1;
// 					nbrs_[i] += 3;
				} else {
					nbrs_[i] += 1;
				}
			}
// 			std::cout << "nbrdiff: old= " << tenA_neighbor_graph.get_node(i)->num_neighbors_counting_self() << " new= " <<
// 				nbrs_[i] << std::endl;
		}
	}

	// ensure full size for the the HbondSet bb_acc/bb_don arrays
	// the bb arrays ought to be resized once at the beginning of this function
	// and not resized to 0 for later push_back operations.... such a waste.
	resize_bb_donor_acceptor_arrays( pose.total_residue() );

	// this is a bug:
	//std::fill( backbone_backbone_donor_.begin(), backbone_backbone_donor_.end(), false );
	//std::fill( backbone_backbone_acceptor_.begin(), backbone_backbone_acceptor_.end(), false );

}

HBondOptions const &
HBondSet::hbond_options() const {
	return *options_;
}

/// @breif set the hbond options for this hbond set; clears all hbonds already stored
void HBondSet::set_hbond_options( HBondOptions const & options )
{
	clear();
	options_ = new HBondOptions( options ); // make a copy of the input to guarantee its integrity
}

std::ostream &
operator<< ( std::ostream & out, const HBondSet & hbond_set ){
	hbond_set.show( out );
	return out;
}

/// @brief various show functions
void
HBondSet::show(std::ostream & out) const
{
	for ( core::Size i=1; i<=(core::Size)nhbonds(); ++i )
		out << hbond((int)i) << std::endl;
}

///
bool
operator==(HBondSet const & a, HBondSet const & b)
{
	if(a.options_ != b.options_) return false;

	if(a.atom_map_ != b.atom_map_) return false;
	if(a.atom_map_init_ != b.atom_map_init_) return false;

	if(a.backbone_backbone_donor_.size() != b.backbone_backbone_donor_.size()) return false;
	for(Size i=1; i<=a.backbone_backbone_donor_.size(); ++i){
		if (a.backbone_backbone_donor_[i] != b.backbone_backbone_donor_[i]) return false;
	}

	if(a.backbone_backbone_acceptor_.size() != b.backbone_backbone_acceptor_.size()) return false;
	for(Size i=1; i<=a.backbone_backbone_donor_.size(); ++i){
		if (a.backbone_backbone_acceptor_[i] != b.backbone_backbone_acceptor_[i]) return false;
	}

	if(a.nhbonds() != b.nhbonds()) {
		t << "HBondSets not equal!" << std::endl;
		t << a << std::endl << b << std::endl;
		return false;
	}

	for(Size i=1; i<=a.nhbonds(); ++i){
		if (!(a.hbond(i) == b.hbond(i))) return false;
	}

	return true;
}

void
HBondSet::show(pose::Pose & pose, std::ostream & out) const
{
	core::Size don_pos, acc_pos, don_res, acc_res;
	std::string don_name, acc_name, don_atom_name, acc_atom_name;
	core::Real weight;
	core::Real hbE;

	// This version is formatted for easy parsing by R, Excel, etc.

	out << "# donor res atom acpt res atom  dist angleD angleA weight energy" << std::endl;
	for ( core::Size i=1; i<=(core::Size)nhbonds(); ++i ) {
		if ( allow_hbond((int)i) ) {
			don_pos = hbond((int)i).don_res();
			don_res = pose.pdb_info()->number(don_pos);
			don_name = pose.residue((int)don_pos).type().name3();
			core::id::AtomID const datm ( hbond((int)i).don_hatm(), don_pos );
			don_atom_name = pose.residue((int)don_pos).atom_name(datm.atomno());

			acc_pos = hbond((int)i).acc_res();
			acc_res = pose.pdb_info()->number(acc_pos);
			acc_name = pose.residue((int)acc_pos).type().name3();
			core::id::AtomID const aatm ( hbond((int)i).acc_atm(), acc_pos );
			acc_atom_name = pose.residue((int)acc_pos).atom_name(aatm.atomno());
			core::Size don_hv_atomno = pose.residue((int)don_pos).atom_base( hbond((int)i).don_hatm() );
			core::Size base_atomno   = pose.residue((int)acc_pos).atom_base( hbond((int)i).acc_atm() );

			numeric::xyzVector < core::Real> Hxyz = pose.residue((int)don_pos).xyz(hbond((int)i).don_hatm());
			numeric::xyzVector < core::Real> Axyz = pose.residue((int)acc_pos).xyz(hbond((int)i).acc_atm());
			numeric::xyzVector < core::Real> Dxyz = pose.residue((int)don_pos).xyz(don_hv_atomno);
			numeric::xyzVector < core::Real> Bxyz = pose.residue((int)acc_pos).xyz(base_atomno);

			core::Real HAdist = Hxyz.distance(Axyz);
			core::Real DHAangle = numeric::angle_radians(Dxyz, Hxyz, Axyz);
			numeric::conversions::to_degrees(DHAangle);
			core::Real HABangle = numeric::angle_radians(Hxyz, Axyz, Bxyz);
			numeric::conversions::to_degrees(HABangle);

			weight = hbond((int)i).weight();
			hbE = hbond((int)i).energy();

			out << "# " << I(5,don_res) << I(4,don_name) << I(5,don_atom_name)
				<< I(5,acc_res) << I(4,acc_name) << I(5,acc_atom_name)
				<< F(6,2,HAdist)
				<< F(7,1,DHAangle) << F(7,1,HABangle)
				<< F(7,3,weight) << F(7,3,hbE) << std::endl;
		}
	}
	out << "# end hbond info " << std::endl;
}

utility::vector1< HBondCOP > HBondSet::empty_list_of_hbonds_;

}
}
}
