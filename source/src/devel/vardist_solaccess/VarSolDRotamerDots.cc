// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/RotamerDots.cc
/// @brief  RotamerDots classes files - ported from rosetta++
/// @author Andrew Leaver-Fay
/// @author Ron Jacak

// Unit Headers
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <devel/vardist_solaccess/LoadVarSolDistSasaCalculatorMover.hh>

// Project headers
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

#include <core/scoring/sasa.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/pack/interaction_graph/RotamerDots.hh>

#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/ubyte.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1.hh>

// Numeric Headers
#include <numeric/constants.hh> // pi
#include <numeric/xyzVector.hh> // to get distance
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>

// Utility Headers
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/vector1.functions.hh>
#include <utility/string_util.hh>

// C++ Headers
// AUTO-REMOVED #include <cstring>
#include <vector>
// AUTO-REMOVED #include <fstream>
#include <iostream>

#include <core/chemical/AtomType.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
static basic::Tracer TR("devel.vardist_solaccess");
//static basic::Tracer TR_DS("core.pack.interaction_graph.RotamerDots.DotSphere");
//static basic::Tracer TR_RD("core.pack.interaction_graph.RotamerDots.RotamerDots");
//static basic::Tracer TR_RDC("core.pack.interaction_graph.RotamerDots.RotamerDotsCache");
//static basic::Tracer TR_RDRD("core.pack.interaction_graph.RotamerDots.RotamerDotsRadiusData");


using namespace ObjexxFCL::format;
using namespace core;

namespace devel {
namespace vardist_solaccess {



//----------------------------------------------------------------------------//
//---------------------------- Rotamer Dots Class ----------------------------//
//----------------------------------------------------------------------------//

Size const VarSolDRotamerDots::num_bytes_ = 21;
Real VarSolDRotamerDots::probe_radius_ = 1.4;

bool VarSolDRotamerDots::sasa_arrays_initialized_ = false;

utility::vector1< utility::vector1< core::Real > > VarSolDRotamerDots::radii_;
utility::vector1< core::Real > VarSolDRotamerDots::msas_radii_;
utility::vector1< core::Real > VarSolDRotamerDots::coll_radii_; // collision radii
utility::vector1< core::Real > VarSolDRotamerDots::int_radii_;  // interaction radii -- larger than coll radii for N and O
utility::vector1< utility::vector1< core::Real > > VarSolDRotamerDots::int_radii_sum_; // atom-type pair interaction radii sums
utility::vector1< utility::vector1< core::Real > > VarSolDRotamerDots::int_radii_sum2_; // atom-type pair interaction radii sums, squared


ObjexxFCL::FArray2D_int const *   VarSolDRotamerDots::lg_angles_( 0 );
ObjexxFCL::FArray2D_ubyte const * VarSolDRotamerDots::lg_masks_( 0 );


///
/// @begin RotamerDots::RotamerDots
///
VarSolDRotamerDots::VarSolDRotamerDots():
	rotamer_(0),
	num_atoms_(0)
{}

///
/// @begin RotamerDots::RotamerDots
///
/// @brief
/// Custom constructor for a RotamerDots object
///
/// @detailed
/// One RotamerDots object get allocated for every state of a first class IG Node, for all first class IG Nodes of a
/// protein being designed. That's potentially a very large number of states. This class should only hold the information
/// it needs to hold to do its job.
///
VarSolDRotamerDots::VarSolDRotamerDots(
	conformation::ResidueCOP rotamer,
	bool all_atom
) :
	rotamer_(rotamer)
{
	if ( ! sasa_arrays_initialized_ ) {
		initialize_sasa_arrays();
		//that function will set sasa_arrays_initialized_ to true;
	}

	if ( all_atom ) {
		num_atoms_ = rotamer->natoms();
	} else {
		num_atoms_ = rotamer_->nheavyatoms();
	}
	atom_coverage_.resize( num_atoms_ );
	for ( core::Size ii = 1; ii <= num_atoms_; ++ii ) {
		atom_coverage_[ ii ].resize( radii_[ rotamer_->atom( ii ).type() ].size() );
	}

}

///
/// @begin RotamerDots::~RotamerDots
///
VarSolDRotamerDots::~VarSolDRotamerDots() {
	//TR_RD << "called destructor" << std::endl;
}


///
/// @begin RotamerDots::RotamerDots
///
/// @brief
/// copy constructor
///
VarSolDRotamerDots::VarSolDRotamerDots( VarSolDRotamerDots const & rhs ) :
	utility::pointer::ReferenceCount(),
	rotamer_( rhs.rotamer_ ),
	num_atoms_( rhs.num_atoms_ ),
	atom_coverage_( rhs.atom_coverage_ )
{}


///
/// @begin RotamerDots::copy
///
/// @brief
/// Copy method for the RotamerDots class. Also used by the assignment operator.
///
void VarSolDRotamerDots::copy( VarSolDRotamerDots const & rhs ) {
	rotamer_ = rhs.rotamer_;
	num_atoms_ = rhs.num_atoms_;
	atom_coverage_ = rhs.atom_coverage_;
}

///
/// @begin RotamerDots::operator=
///
VarSolDRotamerDots const &
VarSolDRotamerDots::operator= ( VarSolDRotamerDots const & rhs ) {
	if ( this != & rhs ) {
		copy( rhs );
	}
	return *this;
}


///
/// @begin RotamerDots::operator!=
///
/// @brief
/// Used during debugging of the HPatchIG.  Some extra information is printed if current state dots is NOT EQUAL to
/// alternate state dots at a Node/BGNode.
///
//bool RotamerDots::operator!=( RotamerDots const & rhs ) {
//
//	if ( state_unassigned() || rhs.state_unassigned() ) {
//		return true; // I guess they could both be unassigned, but better to return not equal than equal
//	}
//
//	if ( rotamer_ != rhs.rotamer_ || num_atoms_ != rhs.num_atoms_ ) return true;
//	if ( sasa_ != rhs.sasa_ || sasa_is_current_ != rhs.sasa_is_current_ ) return true;
//	if ( atom_counts_.size() != rhs.atom_counts_.size() ) return true;
//	if ( atom_sasa_.size() != rhs.atom_sasa_.size() ) return true;
//	if ( radii_ != rhs.radii_ ) return true;
//
//	//TR_RD << "operator!=() score info same. checking if dot sphere objects match." << std::endl;
//	for ( Size ii=1; ii <= atom_counts_.size(); ++ii ) {
//		if ( atom_counts_[ ii ] != rhs.atom_counts_[ ii ] )
//			return true;
//	}
//
//	return false;
//}


///
/// @begin RotamerDots::zero
///
/// @brief
/// Zeros out all of the contained data except the rotamer pointer and the radii array.
///
/// @detailed
/// So far, this function only gets called by the BGNode::prep_for_simA() call so that multiple runs through an
/// interaction graph can be done. If the rotamer dots object on the BGNodes isn't "cleared" after a run, then the run
/// immediately following will have the incorrect counts.
///
//void RotamerDots::zero() {
//
//	// leave the rotamer and num_atoms_ variables untouched
//	//rotamer_ = 0;
//	//num_atoms_ = 0;
//
//	sasa_ = 0.0;
//	sasa_is_current_ = false;
//
//	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
//		atom_sasa_[ ii ] = 0.0;
//		atom_counts_[ ii ].zero(); // calls zero() on each DotSphere instance
//	}
//
//}


///
/// @begin RotamerDots::overlaps
///
/// @brief
/// Returns true if this RotamerDots object has any sphere overlap with the passed in RotamerDots object.
///
/// @detailed
/// This method only checks to see if two RotamerDots objects are within touching distance of each other. It is used
/// to determine whether Edges or BGEdges should be created in the IG. Calculate this using the expanded polar atom
/// radii. If we don't, there's a chance that a state substitution on a Node may cause SASA changes (when expanded polars
/// are used) on a BGNode, but if we didn't assume expanded radii in this method, there would be no edge between the two
/// nodes.
///
bool VarSolDRotamerDots::overlaps( VarSolDRotamerDots const & other ) const {

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= other.get_num_atoms(); ++jj ) {
			Real const distance_squared = get_atom_coords_xyz( ii ).distance_squared( other.get_atom_coords_xyz( jj ) );
			if ( distance_squared <= interaction_radii_squared( rotamer_->atom(ii).type(), other.rotamer_->atom(jj).type() )) return true;
		}
	}
	return false;
}


///
/// @begin RotamerDots::rotamer
///
core::conformation::ResidueCOP
VarSolDRotamerDots::rotamer() const {
	return rotamer_;
}

///
/// @begin RotamerDots::state_unassigned
///
/// @brief
/// Is the state of this RotamerDots object unassigned?
///
//bool RotamerDots::state_unassigned() const {
//	if ( rotamer_ == 0 )
//		return true;
//	return false;
//}

///
/// @begin RotamerDots::get_num_atoms
///
/// @brief
/// Returns the number of atoms this RotamerDots object is keeping SASA for.
///
Size VarSolDRotamerDots::get_num_atoms() const {
	return num_atoms_;
}

///
/// @begin RotamerDots::get_atom_coords_xyz
///
/// @brief
/// Return the xyz coordinates of an atom in this RotamerDots instance.
///
core::Vector
VarSolDRotamerDots::get_atom_coords_xyz( Size atom_index ) const {
	if ( rotamer_ == 0 )
		return numeric::xyzVector< Real >(0,0,0);

	return rotamer_->xyz( atom_index );
}


///
/// @begin RotamerDots::get_atom_radius
///
/// @brief
/// Returns the SASA radius for the passed in atom type. The DB file should have been read in at construct time.
///
/// @detailed
/// Many of the functions in this class iterate over 1 .. num_atoms_.
/// That's not the same thing as an atom type index which is what the radii vector is indexed with. So before we can return
/// the radius, we have to convert the passed in atom_index into the right atom in the residue and then use that to get the
/// right type.
///
//Real RotamerDots::get_atom_radius( Size atom_index ) const {
//	if ( rotamer_ == 0 )
//		return 0.0;
//
//	conformation::Atom const & atom( rotamer_->atom( atom_index ) );
//	return (*radii_)[ atom.type() ];
//
//}

core::Real
VarSolDRotamerDots::get_atom_collision_radius( Size atom_index ) const
{
	return radii_[ rotamer_->atom(atom_index).type() ][ 1 ];
}

core::Real
VarSolDRotamerDots::get_atom_interaction_radius( Size atom_index ) const
{
	return radii_[ rotamer_->atom(atom_index).type() ][ radii_[ rotamer_->atom(atom_index).type() ].size() ];
}


///
/// @begin RotamerDots::radius_for_attype
///
/// @brief
/// Same as the above, but skips the conversion from atom index to atom type index.
///
//core::Real
//RotamerDots::radius_for_attype( Size const attype_index ) {
//	return (*radii_)[ attype_index ];
//}


///
/// @begin RotamerDots::max_atom_radius
///
/// @brief
/// Returns the maximum atom radius. Used only by the SurfacePotential class.
///
//core::Real
//RotamerDots::max_atom_radius() {
//	return utility::max( *radii_ );
//}


///
/// @begin RotamerDots::get_radii
///
/// @brief
/// Returns a pointer to the radii vector. Used only by the InvRotamerDots class.
///
//utility::vector1< Real >*
//RotamerDots::get_radii() const {
//	return radii_;
//}

///
/// @begin RotamerDots::invert_to_boolmasks
///
/// @brief
/// Inverts the current dot counts and saves them to the passed in vector.
///
//void RotamerDots::invert_to_boolmasks( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots ) const {
//
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
//		atom_counts_[ ii ].invert_to_compact_array( inv_dots[ ii ] );
//	}
//}

/// @brief invert the current dot counts for a subset of the atoms in this rotamer.
//void
//RotamerDots::invert_to_boolmasks(
//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots,
//	utility::vector1< Size > const & ats_to_update
//) const
//{
//	for ( Size ii = 1, iiend = ats_to_update.size(); ii <= iiend; ++ii ) {
//		atom_counts_[ ats_to_update[ ii ] ].invert_to_compact_array( inv_dots[ ats_to_update[ ii ] ] );
//	}
//}


//core::Vector
//VarSolDRotamerDots::dot_coord( Size index ) {
//	if ( !sasa_arrays_initialized_ )
//		initialize_sasa_arrays();
//	return dot_coords_[ index ];
//}


///
/// @begin RotamerDots::initialize_sasa_arrays
///
/// @brief
/// Initializes the pointers to the angles and masks FArrays used by sasa.cc and inits the dot sphere coordinates.
///
/// @detailed
/// This call should only occur once (when the first RotamerDots object get constructed) and never again.
///
void VarSolDRotamerDots::initialize_sasa_arrays() {

	if ( sasa_arrays_initialized_ ) return;
	sasa_arrays_initialized_ = true;

	lg_angles_ = ( & core::scoring::get_angles() );
	lg_masks_  = ( & core::scoring::get_masks()  );

	using namespace core::chemical;

	AtomTypeSetCAP atset = ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	core::Size const SASA_RADIUS_INDEX = atset->extra_parameter_index( "SASA_RADIUS" );

	radii_.resize( atset->n_atomtypes() );
	msas_radii_.resize( atset->n_atomtypes() );
	for ( Size ii = 1; ii <= atset->n_atomtypes(); ++ii ) {
		AtomType const & iiattype = (*atset)[ ii ];
		msas_radii_[ii] = iiattype.extra_parameter( SASA_RADIUS_INDEX );

		if ( iiattype.element() == "O" ) {
			radii_[ ii ].resize( 5 );
			radii_[ ii ][ 1 ] = 2.6; // inner shell is collision radius
			radii_[ ii ][ 2 ] = 2.7;
			radii_[ ii ][ 3 ] = 2.8;
			radii_[ ii ][ 4 ] = 2.9;
			radii_[ ii ][ 5 ] = 3.0;
		} else if ( iiattype.element() == "N" ) {
			radii_[ ii ].resize( 5 );
			radii_[ ii ][ 1 ] = 2.7; // inner shell is collision radius
			radii_[ ii ][ 2 ] = 2.8;
			radii_[ ii ][ 3 ] = 2.9;
			radii_[ ii ][ 4 ] = 3.0;
			radii_[ ii ][ 5 ] = 3.1;
		} else if ( iiattype.element() == "C" ) {
			radii_[ ii ].resize( 1 );
			if ( iiattype.atom_type_name() == "COO" || iiattype.atom_type_name() == "CObb" ) {
				radii_[ ii ][ 1 ] = 1.65 + 1.4;
			} else {
				radii_[ ii ][ 1 ] = 1.75 + 1.4;
			}
		} else if ( iiattype.element() == "S" ) {
			radii_[ ii ].resize( 1 );
			radii_[ ii ][ 1 ] = 1.85 + 1.4;
		} else if ( iiattype.element() == "P" ) {
			radii_[ ii ].resize( 1 );
			radii_[ ii ][ 1 ] = 1.9 + 1.4;
		} else if ( iiattype.element() == "H" ) {
			if ( iiattype.atom_type_name() == "Hpol" || iiattype.atom_type_name() == "HNbb" ) {
				radii_[ ii ].resize( 5 );
				radii_[ ii ][ 1 ] = 1.7;
				radii_[ ii ][ 2 ] = 1.8;
				radii_[ ii ][ 3 ] = 1.9;
				radii_[ ii ][ 4 ] = 2.0;
				radii_[ ii ][ 5 ] = 2.1;
			} else if ( iiattype.atom_type_name() == "Haro" ) {
				radii_[ ii ].resize( 1 );
				radii_[ ii ][ 1 ] = 1.0 + 1.4;
			} else {
				radii_[ ii ].resize( 1 );
				radii_[ ii ][ 1 ] = 1.1 + 1.4;
			}
		}
		// else radii_[ ii ].size() == 0
	}
	coll_radii_.resize( atset->n_atomtypes() );
	int_radii_.resize( atset->n_atomtypes() );
	for ( Size ii = 1; ii <= atset->n_atomtypes(); ++ii ) {
		if ( radii_[ii].size() == 0 ) { coll_radii_[ ii ] = int_radii_[ ii ] = 0; continue; }
		coll_radii_[ii] = radii_[ ii ][ 1 ];
		int_radii_[ ii ] = radii_[ ii ][ radii_[ ii ].size() ]; // largest radius in the last position.
	}
	int_radii_sum_.resize(  atset->n_atomtypes() );
	int_radii_sum2_.resize( atset->n_atomtypes() );
 	for ( Size ii = 1; ii <= atset->n_atomtypes(); ++ii ) {
		int_radii_sum_[ ii ].resize( atset->n_atomtypes(), 0.0 );
		int_radii_sum2_[ ii ].resize( atset->n_atomtypes(), 0.0 );
		Real ii_col = coll_radii_[ ii ];
		if ( ii_col == 0 ) continue;
		Real ii_int = int_radii_[ ii ];
		for ( Size jj = 1; jj <= atset->n_atomtypes(); ++jj ) {
			Real jj_col = coll_radii_[ jj ];
			Real jj_int = int_radii_[ jj ];
			if ( ii_col == 0 ) continue;

			if ( ii_col + jj_int < ii_int + jj_col ) {
				int_radii_sum_[ ii ][ jj ] = ii_int + jj_col;
			} else {
				int_radii_sum_[ ii ][ jj ] = ii_col + jj_int;
			}
			int_radii_sum2_[ ii ][ jj ] = int_radii_sum_[ ii ][ jj ] * int_radii_sum_[ ii ][ jj ];
		}
	}

}

core::Real
VarSolDRotamerDots::interaction_radii_squared(
	Size attype1,
	Size attype2
)
{
	return int_radii_sum2_[ attype1 ][ attype2 ];
}

/// @brief
/// computes and stores self-induced dot coverage. uses a vector1 of vector1s of vector1s of
/// ubytes to store the calculated overlap information.
///
/// @detailed
/// uses overlap_atoms() which in turn uses get_atom_overlap_masks()
void VarSolDRotamerDots::increment_self_overlap() {

	using namespace utility; // for utility::vector1

	vector1< vector1< vector1< ObjexxFCL::ubyte > > > self_overlap( num_atoms_ );
	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		self_overlap[ ii ].resize( atom_coverage_[ ii ].size(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		// only have to iterate over the higher indexed atoms for INTRA res overlaps
		for ( Size jj = ii+1; jj <= num_atoms_; ++jj ) {
			overlap_atoms( *this, ii, jj, self_overlap[ii], self_overlap[jj] );
		}
	}

	//TR_RD << "increment_self_overlap(): incrementing with counts: " << std::endl;
	//for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
	//	RotamerDotsCache rdc;
	//	rdc.print_bit_string( self_overlap[ ii ] );
	//}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= self_overlap[ ii ].size(); ++jj ) {
			//TR_RD << "increment_self_overlap(): calling increment_count() on atom dotsphere " << ii << std::endl;
			atom_coverage_[ ii ][ jj ].increment_count( self_overlap[ ii ][ jj ] );
		}
	}
}

void VarSolDRotamerDots::intersect_residues( VarSolDRotamerDots & other )
{
	using namespace utility;
	vector1< vector1< vector1< ObjexxFCL::ubyte > > > this_overlap( num_atoms_ );
	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		this_overlap[ ii ].resize( atom_coverage_[ ii ].size(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	}

	vector1< vector1< vector1< ObjexxFCL::ubyte > > > other_overlap( other.num_atoms_ );
	for ( Size ii = 1; ii <= other.num_atoms_; ++ii ) {
		other_overlap[ ii ].resize( other.atom_coverage_[ ii ].size(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
	}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= other.num_atoms_; ++jj ) {
			overlap_atoms( other, ii, jj, this_overlap[ii], other_overlap[jj] );
		}
	}

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= this_overlap[ ii ].size(); ++jj ) {
			atom_coverage_[ii][jj].increment_count( this_overlap[ii][ jj ] );
		}
	}
	for ( Size ii = 1; ii <= other.num_atoms_; ++ii ) {
		for ( Size jj = 1; jj <= other_overlap[ ii ].size(); ++jj ) {
			other.atom_coverage_[ii][jj].increment_count( other_overlap[ii][ jj ] );
		}
	}

}

bool
VarSolDRotamerDots::any_exposed_dots( Size atom ) const
{
	for ( Size jj = 1; jj <= atom_coverage_[ atom ].size(); ++jj ) {
		if ( atom_coverage_[ atom ][ jj ].get_num_uncovered() != 0 ) return true;
	}
	return false;
}

Real
VarSolDRotamerDots::msas_for_atom( Size atom_index ) const
{
	using namespace core::pack::interaction_graph;
	Size const jj_end = atom_coverage_[ atom_index ].size();
	Size count_exposed = 0;
	for ( Size ii = 1; ii <= DotSphere::NUM_DOTS_TOTAL; ++ii ) {
		for ( Size jj = 1; jj <= jj_end; ++jj ) {
			if ( ! atom_coverage_[ atom_index ][ jj ].get_dot_covered( ii ) ) {
				++count_exposed;
				break;
			}
		}
	}
	return msas_radii_[ rotamer_->atom(atom_index).type() ] * msas_radii_[ rotamer_->atom(atom_index).type() ] *
		4 * numeric::constants::d::pi * ((double) count_exposed / DotSphere::NUM_DOTS_TOTAL );
}

///
/// @begin RotamerDots::increment_this_and_cache
///
/// @brief
/// Add rotamer coverage counts for dots on this object only, leaving rhs unaltered.
///
/// @detailed
/// In the context of the HPatchIG, this method is called by all FCNodes to increment the overlap a BG residue has on
/// the FCNode. It is called by all BG Edges, to make sure that all FCNodes that are connected to a BG residue get this
/// method called by them. 'other' in this case is a BG residue, and 'this_overlap_on_other' is the overlap that is
/// caused by all-states-possible-at-this-node on the BG node. Yes, that's keeping the same information in two places,
/// but it makes updating hpatch score calculations later faster. (ronj)
///
//void RotamerDots::increment_this_and_cache(
//	RotamerDots const & other,
//	RotamerDotsCache & this_overlap_on_other,
//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
//) {
//
//	RotamerDotsCache this_dots_covered_by_other( num_atoms_ );
//	// can't make the variable above be static because this function gets called by all kinds of FC nodes, which will
//	// have varying numbers of atoms.  Then the member vector inside RDC would not be sized correctly and errors would
//	// occur.
//
//	/*TR_RD << "atom_atom_overlaps_cache.size(): " << atom_atom_overlaps_cache.size() << std::endl;
//	for ( Size ii=1; ii <= atom_atom_overlaps_cache.size(); ++ii ) {
//		TR_RD << "atom_atom_overlaps_cache[ " << ii << " ].size(): " << atom_atom_overlaps_cache[ ii ].size() << std::endl;
//	}
//	TR_RD << std::endl;*/
//
//	get_overlap_cache( other, this_overlap_on_other, this_dots_covered_by_other, atom_atom_overlaps_cache );
//
//	// 'increment_count_for_some' is a class method in RotamerDotsCache objects
//	//TR_RD << "increment_this_and_cache(): this_dots_covered_by_other: " << std::endl;
//	//this_dots_covered_by_other.print( std::cout );
//
//	//TR_RD << "increment_this_and_cache(): others_dots_covered_by_this: " << std::endl;
//	//this_overlap_on_other.print( std::cout );
//
//	// increment_from_cached will update both regular and expanded polar SASA. the hard part is getting get_overlap_cache()
//	// to do it right.
//	increment_from_cached( this_dots_covered_by_other );
//
//}


///
/// @begin RotamerDots::get_overlap_cache
///
/// @brief
/// computes the overlap each rotamer (this & other) induce on each other and stores that information in the RotamerDotsCache objects
///
/// @param
/// other - [in] - the other RotamerDots object
/// other_dots_covered_by_this - [out] - the Cache for the dots on the surface of other that are covered by the atoms on this rotamer
/// this_dots_covered_by_other - [out] - the Cache for the dots on the surface of this that are covered by the atoms on the other rotamer
/// atom_atom_overlaps_cache   - [out] - holds a boolean indicating whether two atoms have overlapping, solvent-exposed surface area
///
//void RotamerDots::get_overlap_cache(
//	RotamerDots const & other,
//	RotamerDotsCache & others_dots_covered_by_this,
//	RotamerDotsCache & this_dots_covered_by_other,
//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
//) const {
//
//	using namespace utility;
//
//	//apl static so that they will be allocated exactly once - we don't want to use new and delete in the inner most loop.
//	static vector1< vector1< ObjexxFCL::ubyte > > vv_this_covered_by_other;
//	static vector1< vector1< ObjexxFCL::ubyte > > vv_other_covered_by_this;
//
//	// clear and resize takes alot of time; just make resize calls instead and zero out indices as necessary
//	//vv_this_covered_by_other.clear();
//	//vv_this_covered_by_other.resize( num_atoms_, vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
//	//vv_other_covered_by_this.clear();
//	//vv_other_covered_by_this.resize( other.get_num_atoms(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
//
//	if ( vv_this_covered_by_other.size() < num_atoms_ ) {
//		vv_this_covered_by_other.resize( num_atoms_, vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
//	}
//	if ( vv_other_covered_by_this.size() < other.get_num_atoms() ) {
//		vv_other_covered_by_this.resize( other.get_num_atoms(), vector1< ObjexxFCL::ubyte >( num_bytes_, ObjexxFCL::ubyte( 0 ) ) );
//	}
//
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
//		for ( Size jj = 1; jj <= num_bytes_; ++jj ) {
//			vv_this_covered_by_other[ ii ][ jj ] = 0;
//		}
//	}
//	for ( Size ii = 1; ii <= other.num_atoms_; ++ii ) {
//		for ( Size jj = 1; jj <= num_bytes_; ++jj ) {
//			vv_other_covered_by_this[ ii ][ jj ] = 0;
//		}
//	}
//
//	// calls the zero() class method on both of the RotamerDotsCache references passed in.  I guess we do that to invalidate the cache
//	// before storing the new cache values.
//	others_dots_covered_by_this.zero();
//	this_dots_covered_by_other.zero();
//
//	get_res_res_overlap( other, vv_this_covered_by_other, vv_other_covered_by_this, atom_atom_overlaps_cache );
//
//	this_dots_covered_by_other.increment_count( vv_this_covered_by_other );
//	others_dots_covered_by_this.increment_count( vv_other_covered_by_this );
//
//}


///
/// @begin RotamerDots::get_res_res_overlap()
///
/// @brief
/// Calls get_atom_atom_coverage for all atom pairs between res1 and res2. This method gets called by RotamerDots::get_overlap_cache().
///
//void RotamerDots::get_res_res_overlap(
//	RotamerDots const & other_res,
//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res1_covered_by_res2, // may include more entries than atoms
//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res2_covered_by_res1, // may include more entries than atoms
//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
//) const {
//
//	Real square_distance = 0.0f;
//	bool overlap;
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) { // 'this'/self residues atoms
//
//		Vector ii_atom_xyz = get_atom_coords_xyz( ii );
//		Real ii_atom_radius;
//		ii_atom_radius = get_atom_radius( ii );
//
//		for ( Size jj = 1; jj <= other_res.get_num_atoms(); ++jj ) {
//
//			Vector jj_atom_xyz = other_res.get_atom_coords_xyz( jj );
//			Real jj_atom_radius;
//			jj_atom_radius = other_res.get_atom_radius( jj );
//
//			overlap = get_atom_atom_coverage( ii_atom_xyz, ii_atom_radius, jj_atom_xyz, jj_atom_radius, res1_covered_by_res2[ ii ], res2_covered_by_res1[ jj ], square_distance );
//
//			// set the overlaps bool if overlap was found. the outer vector holds the changing node's atoms (or other_res, as it
//			// is referred to in this function) and the inner vector is for this RD object's atoms.
//			if ( overlap )
//				atom_atom_overlaps_cache[ jj ][ ii ] = true;
//
//			if ( square_distance > (ii_atom_radius + jj_atom_radius) * (ii_atom_radius + jj_atom_radius) ) {
//				break;
//			}
//		}
//	}
//}


///
/// @begin RotamerDots::get_atom_atom_coverage()
///
/// @brief
/// returns false if the two spheres do not overlap at all. otherwise, saves the overlap masks to the input vectors.
///
//bool RotamerDots::get_atom_atom_coverage( Vector const & at1_xyz, Real at1_base_radius,
//	Vector const & at2_xyz, Real at2_base_radius,
//	utility::vector1< ObjexxFCL::ubyte > & at1_sphere_covered,
//	utility::vector1< ObjexxFCL::ubyte > & at2_sphere_covered, Real dist_sq ) const {
//
//	int degree_of_overlap;
//	int aphi_1_2, aphi_2_1;
//	int theta_1_2, theta_2_1;
//	int masknum;
//
//	Real at1_radius = at1_base_radius + probe_radius_;
//	Real at2_radius = at2_base_radius + probe_radius_;
//
//	// exit if large probe radii do not touch
//	dist_sq = at1_xyz.distance_squared( at2_xyz );
//	if ( dist_sq > (at1_radius + at2_radius) * (at1_radius + at2_radius) ) {
//		return false;
//	}
//	Real const distance = std::sqrt( dist_sq );
//
//	//ronj this block represents the amount of surface area covered up on atom1 by atom2
//	core::scoring::get_overlap( at1_radius, at2_radius, distance, degree_of_overlap );
//	core::scoring::get_2way_orientation( at1_xyz, at2_xyz, aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );
//
//	Size closest_dot1 = (*lg_angles_)( aphi_1_2, theta_1_2 );
//	masknum = ( closest_dot1 * 100 ) + degree_of_overlap;
//	for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
//		at1_sphere_covered[ bb ] |= (*lg_masks_)[ bbli ];
//	}
//
//	//ronj the amount of surface area covered up on atom2 by atom1
//	core::scoring::get_overlap( at2_radius, at1_radius, distance, degree_of_overlap );
//
//	Size closest_dot2 = (*lg_angles_)( aphi_2_1, theta_2_1 );
//	masknum = ( closest_dot2 * 100 ) + degree_of_overlap;
//	for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
//		at2_sphere_covered[ bb ] |= (*lg_masks_)[ bbli ];
//	}
//
//	return true;
//}

bool
VarSolDRotamerDots::get_atom_overlap_masks(
	VarSolDRotamerDots const & other,
	Size at_this,
	Size at_other,
	Real & distance,
	Size & closest_dot1,
	Size & closest_dot2
) const
{
	core::Vector const & at1_xyz( rotamer_->xyz( at_this ));
	core::Vector const & at2_xyz( other.rotamer_->xyz( at_other ));
	Real dist_sq = at1_xyz.distance_squared( at2_xyz );
	if ( dist_sq > interaction_radii_squared( rotamer_->atom(at_this).type(), other.rotamer_->atom(at_other).type() ) ) return false;

	// int degree_of_overlap; // Unused variable causes warning.
	int aphi_1_2, aphi_2_1;
	int theta_1_2, theta_2_1;
	// int masknum; // Unused variable causes warning.

	distance = std::sqrt( dist_sq );

	//ronj this block represents the amount of surface area covered up on atom1 by atom2
	//core::scoring::get_overlap( at1_radius, at2_radius, distance, degree_of_overlap );
	core::scoring::get_2way_orientation( at1_xyz, at2_xyz, aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );
	closest_dot1 = (*lg_angles_)( aphi_1_2, theta_1_2 );
	closest_dot2 = (*lg_angles_)( aphi_2_1, theta_2_1 );
	return true;
}

bool
VarSolDRotamerDots::overlap_atoms(
	VarSolDRotamerDots const & other,
	Size at_this,
	Size at_other,
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & at_this_coverage,
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & at_other_coverage
) const
{
	assert( atom_coverage_[ at_this ].size() == at_this_coverage.size() );
	assert( other.atom_coverage_[ at_other ].size() == at_other_coverage.size() );
	Real distance;
	Size closest_dot1, closest_dot2;
	if ( ! get_atom_overlap_masks( other, at_this, at_other, distance, closest_dot1, closest_dot2 ) ) return false;
	Size const attype1( rotamer_->atom( at_this ).type() );
	Size const attype2( other.rotamer_->atom( at_other ).type() );
	Size const ii_nradii = radii_[ attype1 ].size(), jj_nradii = radii_[ attype2 ].size();
	for ( Size ii = 1; ii <= ii_nradii; ++ii ) {
		for ( Size jj = 1; jj <= jj_nradii; ++jj ) {
			if ( distance > radii_[ attype1 ][ ii ] + radii_[ attype2 ][ jj ] ) continue;
			if ( ii == 1 ) {
				int degree_of_overlap;
				core::scoring::get_overlap( radii_[ attype2 ][ jj ], radii_[ attype1 ][ ii ], distance, degree_of_overlap );
				Size masknum = ( closest_dot2 * 100 ) + degree_of_overlap;
				for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
					at_other_coverage[ jj ][ bb ] |= (*lg_masks_)[ bbli ];
				}
			}
			if ( jj == 1 ) {
				int degree_of_overlap;
				core::scoring::get_overlap( radii_[ attype1 ][ ii ], radii_[ attype2 ][ jj ], distance, degree_of_overlap );
				Size masknum = ( closest_dot1 * 100 ) + degree_of_overlap;
				for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
					at_this_coverage[ ii ][ bb ] |= (*lg_masks_)[ bbli ];
				}
			}
		}
	}
	return true;
}


///
/// @begin RotamerDots::increment_from_cached
///
/// @brief
/// Increments the dot coverage count for this rotamer from a coverage cache
///
//void RotamerDots::increment_from_cached( RotamerDotsCache const & cache ) {
//
//	if ( cache.atom_counts_.size() != atom_counts_.size() )
//		TR_RD << "increment_from_cached(): cache.size(): " << cache.atom_counts_.size() << ", atom_counts_.size(): " << atom_counts_.size() << std::endl;
//	assert( cache.atom_counts_.size() == atom_counts_.size() );
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
//		atom_counts_[ ii ] += cache.atom_counts_[ ii ];
//	}
//
//	sasa_is_current_ = false;
//}
//
//
///
/// @begin RotamerDots::decrement_from_cached
///
/// @brief
/// decrements the dot coverage count for this by the coverage stored in the input RotamerDotsCache object
///
//void RotamerDots::decrement_from_cached( RotamerDotsCache const & cache ) {
//
//	assert( cache.atom_counts_.size() == atom_counts_.size() );
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
//		atom_counts_[ ii ] -= cache.atom_counts_[ ii ];
//	}
//
//	sasa_is_current_ = false;
//}
//
//
///
/// @begin RotamerDots::increment_both
///
/// @brief
/// Add rotamer coverage counts for dots on both this and other. sets sasa_is_current_ to false on both this and rhs
///
/// @detailed
/// One use case involves BGNodes initializing overlap with other BGNodes. This is a brute force all BGNode v all BGNode
/// pairwise increment on the RotamerDots objects each Node holds. That's why we call increment *both*.
/// We don't care about the cache values in this case, but to use increment_both_and_cache, we have to create Cache variables
/// to use as references.
///
//void RotamerDots::increment_both( RotamerDots & other ) {
//
//	RotamerDotsCache others_dots_covered_by_this( other.get_num_atoms() );
//	others_dots_covered_by_this.zero();
//
//	RotamerDotsCache this_dots_covered_by_other( num_atoms_ );
//	this_dots_covered_by_other.zero();
//
//	// can't make the variables above be static because this function gets called by all kinds of FC nodes, which will
//	// have varying numbers of atoms.  Then the member vector inside RDC would not be sized correctly and errors would
//	// occur. But RDC object are pretty lightweight, so it shouldn't result in a big performance hit.
//	utility::vector1< utility::vector1< bool > > dummy( other.get_num_atoms(), utility::vector1< bool >( num_atoms_, false ) );
//	increment_both_and_cache( other, others_dots_covered_by_this, this_dots_covered_by_other, dummy );
//}
//
///
/// @begin RotamerDots::increment_both_and_cache
///
/// @detailed
/// Add rotamer coverage counts for dots on both this and other. Cache's the overlap this and other have with each other for greater efficiency.
/// The second parameter are the dots on the surface of other_rotamer that this overlaps. The third parameter are the dots on the surface of this
/// that other_rotamer overlaps. The fourth parameter lives on the Edges of the IG and stores a boolean indicating whether two atoms have
/// overlapping, solvent-exposed surface area. Instead of recalculating this for all residue pairs every substitution, keep it in this cache.
/// The vectors are already sized. The outer vector has the "other_rotamer"'s atoms and the inner vector is for the atoms in this class's residue.
/// The structure just gets passed down to get_overlap_cache() and it eventually gets filled in get_atom_atom_coverage().
///
/// If the class is keeping the expanded polar atom SASA also, then this function makes an extra call to get_overlap_cache()
/// to get the extra overlap that happens when polar atom radii are expanded.
///
//void RotamerDots::increment_both_and_cache( RotamerDots & other_rotamer, RotamerDotsCache & others_dots_covered_by_this,
//	RotamerDotsCache & this_dots_covered_by_other, utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache ) {
//
//	get_overlap_cache( other_rotamer, others_dots_covered_by_this, this_dots_covered_by_other, atom_atom_overlaps_cache );
//
//	//TR_RD << "increment_both_and_cache(): overlap found:" << std::endl;
//	//this_dots_covered_by_other.print( std::cout );
//	//others_dots_covered_by_this.print( std::cout );
//
//	increment_from_cached( this_dots_covered_by_other );
//	other_rotamer.increment_from_cached( others_dots_covered_by_this );
//
//}
//
//
///
/// @begin RotamerDots::get_sasa
///
/// @brief
/// Given the current dot coverage counts, returns the total SASA of this residue.
///
/// @detailed
/// This method does not do any work figuring out who is overlapping with this rotamer. It assumes that work has been
/// done. Instead, it returns the SASA of the dot counts currently held.  If the dot counts have not changed since the
/// last time get_sasa() got called then sasa_is_current_ will be true.  In that case, the method will just return the
/// the variable sasa_. If the counts have changed, it iterates over all the atoms and recalculates the total SASA.
/// That value is stored in sasa_ and sasa_is_current_ is set to true.
///
/// The reason both sasa_ and sasa_is_current_ are "mutable" is so that they can be modified inside this const function.
///
//Real RotamerDots::get_sasa() const {
//
//	if ( state_unassigned() )
//		return 0.0f;
//
//	if ( sasa_is_current_ )
//		return sasa_;
//
//	Real fraction_uncovered = 0.0;
//	Real atom_radius = 0.0;
//	Real const four_pi = 4.0 * Real( numeric::constants::d::pi );
//	Real atom_area_exposed = 0.0;
//
//	sasa_ = 0.0;
//	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
//		fraction_uncovered = static_cast< Real >( get_num_uncovered( ii ) ) / atom_counts_[ii].get_total_dots();
//		atom_radius = get_atom_radius( ii ) + probe_radius_;
//		atom_area_exposed = four_pi * ( atom_radius * atom_radius ) * fraction_uncovered;
//		//std::cout << "get_sasa(): atom: " << ii << ", rad: " << atom_radius << ", %-uncovered: " << fraction_uncovered << ", exposed: " << atom_area_exposed << std::endl;
//
//		atom_sasa_[ ii ] = atom_area_exposed;
//		sasa_ += atom_area_exposed;
//	}
//
//	sasa_is_current_ = true;
//
//	return sasa_;
//}
//
//
///
/// @begin RotamerDots::get_atom_sasa
///
/// @brief
/// Given the current dot coverage counts, returns the total SASA for a particular atom index.
/// Assumes that get_atom_sasa() will never be called when the object is in the unassigned state.
///
//Real RotamerDots::get_atom_sasa( Size atom_index ) const {
//	if ( ! sasa_is_current_ )
//		get_sasa();
//
//	return atom_sasa_[ atom_index ];
//}
//
//
///
/// @begin RotamerDots::get_num_uncovered
///
/// @brief
/// Returns the number of uncovered dots on the given atom, when using standard SASA radii.
/// Note: no expanded polars version of this method.
///
//Size RotamerDots::get_num_uncovered( Size atom ) const {
//	return atom_counts_[ atom ].get_num_uncovered();
//}
//
//
///
/// @begin RotamerDots::get_num_covered
///
/// @brief
/// Note: no expanded polars version of this method.
///
//Size RotamerDots::get_num_covered_total() const {
//
//	Size total_num_covered = 0;
//	for ( Size ii=1; ii <= atom_counts_.size(); ++ii ) {
//		total_num_covered += atom_counts_[ ii ].get_num_covered();
//	}
//
//	return total_num_covered;
//}
//
//
///
/// @begin RotamerDots::write_dot_kinemage
///
///*void RotamerDots::write_dot_kinemage( std::ofstream & kinfile ) {
//
//	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
//		write_dotlist_header( kinfile, "1.4 exposed dots", "red");
//		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
//			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
//				write_dot( kinfile, ii, jj, probe_radius_ );
//			}
//		}
//		write_dotlist_header( kinfile, "SA", "red");
//		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
//			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
//				write_dot( kinfile, ii, jj, 0 );
//			}
//		}
//		write_dotlist_header( kinfile, "1.0 A probe Accessible", "green");
//		for ( Size jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj ) {
//			if ( ! atom_counts_[ii].get_dot_covered( jj ) ) {
//				write_dot( kinfile, ii, jj, 0 );
//			}
//		}
//		//write_dotlist_header( kinfile, "void surface", "blue");
//		//for (int jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj) {
//		//	if ( atom_counts_[ii].get_dot_covered( jj ) && ! atom_counts_small_probe_[ii].get_dot_covered( jj ) ) {
//		//		write_dot( kinfile, ii, jj, 0 );
//		//	}
//		//}
//		//write_dotlist_header( kinfile, "all dots", "white");
//		//for (int jj = 1; jj <= DotSphere::NUM_DOTS_TOTAL; ++jj)
//		//{
//		//	write_dot( kinfile, ii, jj, 0 );
//		//}
//	}
//}
//
//
///
/// @begin RotamerDots::write_dotlist_header
///
//void RotamerDots::write_dotlist_header( std::ofstream & kinfile, std::string master_name, std::string color ) {
//	kinfile << "@dotlist master= {" << master_name << "} color= " << color << "\n";
//}
//
///
/// @begin RotamerDots::write_dot
///
//void RotamerDots::write_dot( std::ofstream & kinfile, Size atom, Size dot, Real radius ) {
//	static numeric::xyzVector< Real > coord;
//	coord = get_atom_coords_xyz( atom );
//	coord += (radius + get_atom_radius( atom )) * get_dot_coord( dot );
//
//	write_dot( kinfile, coord, "dot" );
//}
//
///
/// @begin RotamerDots::write_dot
///
//void RotamerDots::write_dot( std::ofstream & kinfile, numeric::xyzVector< Real > const & coord, std::string atname ) {
//	static std::string last_atname = "";
//	if ( last_atname == atname ) {
//		kinfile << "{\"} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
//	} else {
//		kinfile << "{" << atname << "} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
//		last_atname = atname;
//	}
//}
//
//
///
/// @begin RotamerDots::get_dot_coord
///
//numeric::xyzVector<Real> const & RotamerDots::get_dot_coord( Size dot_id ) {
//	if ( ! dot_sphere_coordinates_initialized_ ) {
//		initialize_dot_sphere_coordinates_from_file();
//	}
//	return dot_sphere_coordinates_[ dot_id ];
//}
//
//
///
/// @begin RotamerDots::initialize_dot_sphere_coordinates_from_file
///
//void RotamerDots::initialize_dot_sphere_coordinates_from_file() {
//
//	dot_sphere_coordinates_.resize( DotSphere::NUM_DOTS_TOTAL );
//
//	std::ifstream dotfile("sphere.txt");
//	for ( Size ii = 1; ii <= DotSphere::NUM_DOTS_TOTAL; ++ii ) {
//		Real x, y, z;
//		dotfile >> x >> y >> z;
//		dot_sphere_coordinates_[ ii ] = -1 * numeric::xyzVector< Real >( x,y,z );
//	}
//	dotfile.close();
//
//	dot_sphere_coordinates_initialized_ = true;
//}
//*/
//
///
/// @begin RotamerDots::print
///
//void RotamerDots::print( std::ostream & os ) const {
//
//	if ( state_unassigned() ) {
//		os << "dots: unassigned" << std::endl;
//		return;
//	}
//
//	os << "dots: " << rotamer_->name3() << rotamer_->seqpos() << std::endl;
//	for ( Size ii=1; ii <= num_atoms_; ++ii ) {
//		os << "atom " << rotamer_->atom_name( ii ) << ", rad: " << get_atom_radius( ii )
//			<< ", covered: " << ObjexxFCL::format::I(3,atom_counts_[ii].get_num_covered()) << ", counts: ";
//		atom_counts_[ ii ].print( os );
//	}
//	os << "num_covered: " << get_num_covered_total() << ", sasa_is_current_: " << sasa_is_current_ << ", sasa: " << get_sasa() << std::endl;
//
//}
//
//
//std::string RotamerDots::name3() const { return rotamer_->name3(); } // for operator<<
//core::Size  RotamerDots::seqpos() const { return rotamer_->seqpos(); } // for operator<<
//
///
/// @begin operator<< ( ostream, RotamerDots )
///
/// @brief
/// invokes print on the input RotamerDots object
///
//std::ostream & operator<< ( std::ostream & os, RotamerDots const & rd ) {
//	if ( rd.state_unassigned() ) {
//		os << "unassigned";
//		return os;
//	}
//	os << "rotamer: " << rd.name3() << rd.seqpos() << ", sasa: " << rd.get_sasa();
//	return os;
//}

void VarSolDRotamerDots::write_dot_kinemage( std::ostream & kinfile ) const
{
	using namespace core::pack::interaction_graph;

	for ( Size ii = 1; ii <= num_atoms_; ++ii ) {
		Size const ii_attype = rotamer_->atom( ii ).type();

		// Color hydrogen dots by their base atom
		Size iirepatom = ii;
		if ( rotamer_->atom_is_hydrogen( ii ) ) iirepatom = rotamer_->atom_base( ii );

		if ( rotamer_->type().atom_type( iirepatom ).element() == "O" ) {
			write_dotlist_header( kinfile, "exposed polar dots", "red");
		} else if ( rotamer_->type().atom_type( iirepatom ).element() == "N" ) {
			write_dotlist_header( kinfile, "exposed polar dots", "blue");
		} else {
			write_dotlist_header( kinfile, "exposed hydrophobic dots", "gray");
		}

		for ( Size jj =1; jj <= atom_coverage_[ ii ].size(); ++jj ) {
			for ( Size kk = 1; kk <= DotSphere::NUM_DOTS_TOTAL; ++kk ) {
				if ( ! atom_coverage_[ii][jj].get_dot_covered( kk ) ) {
					write_dot( kinfile, ii, kk, radii_[ ii_attype ][ jj ] );
				}
			}
		}
	}
}

void
VarSolDRotamerDots::write_dotlist_header(
	std::ostream & kinfile,
	std::string const & master_name,
	std::string const & color
) const
{
	kinfile << "@dotlist master= {" << master_name << "} color= " << color << "\n";
}

void
VarSolDRotamerDots::write_dot(
	std::ostream & kinfile,
	core::Size atom,
	core::Size dot,
	core::Real radius
) const
{
	numeric::xyzVector< Real > coord;
	coord = rotamer_->xyz( atom );
	coord += radius * core::pack::interaction_graph::RotamerDots::dot_coord( dot );

	write_dot( kinfile, coord, "dot" );
}

void VarSolDRotamerDots::write_dot(
	std::ostream & kinfile,
	core::Vector const & coord,
	std::string const & atname
) const
{
	kinfile << "{" << atname << "} P " << coord.x() << " " << coord.y() << " " << coord.z() << "\n";
}



VarSolDistSasaCalculator::VarSolDistSasaCalculator()
{
}

core::pose::metrics::PoseMetricCalculatorOP
VarSolDistSasaCalculator::clone() const{
	return new VarSolDistSasaCalculator;
}

void
VarSolDistSasaCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const
{
	// STOLEN from SasaCalculator.cc
	TR << "VarSolDistSasaCalculator::lookup" << std::endl;

	if ( key == "total_sasa" ) {
		basic::check_cast( valptr, &total_sasa_, "total_sasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_sasa_ );

	} else if ( key == "atom_sasa" ) {
		basic::check_cast( valptr, &atom_sasa_, "atom_sasa expects to return a id::AtomID_Map< Real >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< Real > > *>(valptr))->set( atom_sasa_ );

	} else if ( key == "residue_sasa" ) {
		basic::check_cast( valptr, &residue_sasa_, "residue_sasa expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_sasa_ );

	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}

std::string
VarSolDistSasaCalculator::print( std::string const & key ) const
{
	// STOLEN from SasaCalculator.cc

	if ( key == "total_sasa" ) {
		return utility::to_string( total_sasa_ );
	} else if ( key == "atom_sasa" ) {
		basic::Error() << "id::AtomID_Map< Real > has no output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "residue_sasa" ) {
		return utility::to_string( residue_sasa_ );
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";
}

void
VarSolDistSasaCalculator::recompute( core::pose::Pose const & this_pose )
{
	TR << "VarSolDistSasaCalculator::recompute" << std::endl;
	core::pose::initialize_atomid_map( atom_sasa_, this_pose, 0.0 );

	rotamer_dots_.resize( this_pose.total_residue() );
	residue_sasa_.resize( this_pose.total_residue() );
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		rotamer_dots_[ ii ] = new VarSolDRotamerDots( new core::conformation::Residue( this_pose.residue( ii ) ), true );
		rotamer_dots_[ ii ]->increment_self_overlap();
	}
	core::scoring::EnergyGraph const & energy_graph( this_pose.energies().energy_graph() );
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		for ( core::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			rotamer_dots_[ ii ]->intersect_residues( *rotamer_dots_[ (*iru)->get_second_node_ind() ] );
		}
	}
	total_sasa_ = 0.0;
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		residue_sasa_[ ii ] = 0.0;
		for ( Size jj = 1; jj <= this_pose.residue(ii).natoms(); ++jj ){
			core::id::AtomID at( jj, ii );
			residue_sasa_[ ii ] += atom_sasa_[ at ] = rotamer_dots_[ ii ]->msas_for_atom( jj );
		}
		total_sasa_ += residue_sasa_[ ii ];
	}

}



//----------------------------------------------------------------------------//
//------------------------- RotamerDotsRadiusData  ---------------------------//
//----------------------------------------------------------------------------//

/// @brief set initial value as no instance
//RotamerDotsRadiusData* RotamerDotsRadiusData::instance_( 0 );
//
/// @brief static function to get the instance of (pointer to) this singleton class
//RotamerDotsRadiusData* RotamerDotsRadiusData::get_instance() {
//	if ( instance_ == 0 ) {
//		instance_ = new RotamerDotsRadiusData();
//	}
//	return instance_;
//}
//
/// @brief private constructor to guarantee the singleton
//RotamerDotsRadiusData::RotamerDotsRadiusData() {}
//
//utility::vector1< Real >*
//RotamerDotsRadiusData::get_ROSETTA_SASA_radii() {
//
//	using namespace core::chemical;
//
//	if ( ROSETTA_SASA_radii_.size() == 0 ) {
//		// need to size and set the values in the vector
//
//		TR_RDRD << "get_ROSETTA_SASA_radii(): reading in sasa radii database file" << std::endl;
//
//		//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
//		//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
//		//ronj for each atom type into the radii array. each index of the radii array corresponds to some atom type.
//		AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
//		ROSETTA_SASA_radii_.resize( atom_type_set.n_atomtypes(), 0.0 );
//
//		core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "SASA_RADIUS" );
//
//		TR_RDRD << "ROSETTA_SASA_radii_: [ ";
//		for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
//			ROSETTA_SASA_radii_[ ii ] = atom_type_set[ ii ].extra_parameter( SASA_RADIUS_INDEX );
//			TR_RDRD << ROSETTA_SASA_radii_[ ii ] << ", ";
//		}
//		TR_RDRD << "]" << std::endl;
//	}
//
//	return &ROSETTA_SASA_radii_;
//
//}
//
//utility::vector1< Real >*
//RotamerDotsRadiusData::get_NACCESS_SASA_radii() {
//
//	using namespace core::chemical;
//
//	if ( NACCESS_SASA_radii_.size() == 0 ) {
//		// need to size and set the values in the vector
//
//		TR_RDRD << "get_NACCESS_SASA_radii(): reading in sasa radii database file" << std::endl;
//
//		AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
//		NACCESS_SASA_radii_.resize( atom_type_set.n_atomtypes(), 0.0 );
//
//		core::Size NACCESS_SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "NACCESS_SASA_RADIUS" );
//
//		TR_RDRD << "NACCESS_SASA_radii_: [ ";
//		for ( core::Size ii=1; ii <= atom_type_set.n_atomtypes(); ++ii ) {
//			NACCESS_SASA_radii_[ ii ] = atom_type_set[ ii ].extra_parameter( NACCESS_SASA_RADIUS_INDEX );
//			TR_RDRD << NACCESS_SASA_radii_[ ii ] << ", ";
//		}
//		TR_RDRD << "]" << std::endl;
//	}
//
//	return &NACCESS_SASA_radii_;
//
//}
//
//utility::vector1< Real >*
//RotamerDotsRadiusData::get_NACCESS_SASA_radii_with_expanded_polars( Real polar_expansion_radius ) {
//
//	using namespace core::chemical;
//
//	if ( NACCESS_SASA_radii_with_expanded_polars.size() == 0 ) {
//
//		TR_RDRD << "get_NACCESS_SASA_radii_with_expanded_polars(): reading in sasa radii database file" << std::endl;
//
//		utility::vector1< Real >* radii = get_NACCESS_SASA_radii();
//		NACCESS_SASA_radii_with_expanded_polars.resize( (*radii).size() );
//
//		// used to figure out what atom_type a given index is
//		AtomTypeSet const & atom_type_set = * ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
//
//		TR_RDRD << "NACCESS_SASA_radii_with_expanded_polars_: [ ";
//		for ( Size ii=1; ii <= radii->size(); ++ii ) {
//			NACCESS_SASA_radii_with_expanded_polars[ ii ] = (*radii)[ ii ];
//
//			core::chemical::AtomType const & at( atom_type_set[ ii ] );
//			if ( at.element() == "N" || at.element() == "O" ) {
//				NACCESS_SASA_radii_with_expanded_polars[ ii ] += polar_expansion_radius;
//			}
//			TR_RDRD << NACCESS_SASA_radii_with_expanded_polars[ ii ] << ", ";
//		}
//		TR_RDRD << "]" << std::endl;
//
//	}
//
//	return &NACCESS_SASA_radii_with_expanded_polars;
//
//}
//
//
////----------------------------------------------------------------------------//
////------------------------- Rotamer Dots Cache Class -------------------------//
////----------------------------------------------------------------------------//
//
//
///
/// @begin RotamerDotsCache::RotamerDotsCache
///
//RotamerDotsCache::RotamerDotsCache() {}
//
///
/// @begin RotamerDotsCache::RotamerDotsCache
///
//RotamerDotsCache::RotamerDotsCache( Size num_atoms ) {
//	atom_counts_.resize( num_atoms );
//}
//
//
///
/// @begin RotamerDotsCache::RotamerDotsCache
///
/// @brief
/// copy constructor
///
//RotamerDotsCache::RotamerDotsCache( RotamerDotsCache const & rhs ) :
//	atom_counts_( rhs.atom_counts_ )
//{}
//
//
///
/// @begin RotamerDotsCache::~RotamerDotsCache
///
//RotamerDotsCache::~RotamerDotsCache() {}
//
//
///
/// @begin RotamerDotsCache::operator=
///
/// @brief
/// assignment operator
///
//RotamerDotsCache const & RotamerDotsCache::operator=( RotamerDotsCache const & rhs ) {
//	atom_counts_ = rhs.atom_counts_;
//	return *this;
//}
//
//
///
/// @begin RotamerDotsCache::resize
///
//void RotamerDotsCache::resize( Size num_atoms ) {
//	//TR_RDC << "resize() called with num_atoms: " << num_atoms << std::endl;
//	atom_counts_.clear();
//	atom_counts_.resize( num_atoms );
//}
//
//
///
/// @begin RotamerDotsCache::zero
///
/// @brief
/// sets the dot counts to zero for all atoms
///
/// @detailed
/// if the cache already knows the atom's counts are uniformly 0, it skips it
///
//void
//RotamerDotsCache::zero() {
//	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
//		//if ( non_zero_overlap_[ ii ] ) {
//			atom_counts_[ ii ].zero(); // calls zero() on each DotSphere instance
//		//}
//	}
//	//atom_counts_.clear(); // leaving this in causes the vector to be sized down to 0, causing problems
//}
//
//
///
/// @begin RotamerDotsCache::increment_count
///
/// @brief
/// increments the dot coverage counts for all atoms in the cache
///
/// @param
/// covered - [in] - compact ubyte array of dots; '1' for any covered dot for the vdw + 1.4 A sphere
///
//void RotamerDotsCache::increment_count( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & covered ) {
//
//	//TR_RDC << "increment_count(): atom_counts_.size(): " << atom_counts_.size() << ", covered.size(): " << covered.size() << std::endl;
//	// do not assert this -- let there be more entries, possibly, in covered array;
//	// assert( atom_counts_.size() == covered.size() );
//	assert( atom_counts_.size() <= covered.size() );
//	for ( Size ii = 1, ii_end = atom_counts_.size(); ii <= ii_end; ++ii ) {
//		//utility::vector1< ObjexxFCL::ubyte > const & atom_mask = covered[ ii ];
//		atom_counts_[ ii ].increment_count( covered[ ii ] );
//	}
//
//}
//
//
///
/// @begin RotamerDotsCache::write_to_compact_array
///
/// @detailed
/// Called by BGEdges with a reference to a vector1 of vector1 of ubytes representing where to put the compact (count)
/// based representation of all the atom overlap masks.
///
/// Note: This method only writes the standard SASA dot counts to the passed in array. No version of this function
/// exists for expanded polar atom SASA dot counts.
///
//void RotamerDotsCache::write_to_compact_array( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & compact ) const {
//	assert( compact.size() != 0 );
//	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
//		atom_counts_[ ii ].write_to_compact_array( compact[ ii ] );
//	}
//}
//
//
///
/// @begin RotamerDotsCache::print()
///
//void RotamerDotsCache::print_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const {
//	for ( Size bb = 1; bb <= values.size(); ++bb ) {
//		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
//		for ( int index=7; index >= 0; index-- ) {
//			std::cout << ( ( (int)values[ bb ] >> index ) & 1 );
//		}
//		std::cout << " ";
//	}
//	std::cout << std::endl;
//}
//
//
///
/// @begin RotamerDotsCache::print()
///
//void RotamerDotsCache::print( std::ostream & os ) const {
//	for ( Size ii = 1; ii <= atom_counts_.size(); ++ii ) {
//		os << "atom " << I(2,ii) << ": ";
//		atom_counts_[ ii ].print( os );
//	}
//}
//
////----------------------------------------------------------------------------//
////------------------------------ Inverted Dots Class -------------------------//
////----------------------------------------------------------------------------//
//
//Real const InvRotamerDots::max_dist_from_dot_to_intersection = 0.8;
//
//InvRotamerDots::InvRotamerDots() :
//	rotamer_( 0 )
//{}
//
//
//InvRotamerDots::InvRotamerDots( InvRotamerDots const & src ) :
//	utility::pointer::ReferenceCount(),
//	rotamer_( src.rotamer_ ),
//	inv_dots_( src.inv_dots_ ),
//	radii_( src.radii_ )
//{}
//
//InvRotamerDots::~InvRotamerDots() {}
//
//InvRotamerDots const &
//InvRotamerDots::operator= ( InvRotamerDots const & rhs ) {
//	if ( this != & rhs ) {
//		rotamer_  = rhs.rotamer_;
//		inv_dots_ = rhs.inv_dots_;
//		radii_ = rhs.radii_;
//	}
//	return *this;
//}
//
//void
//InvRotamerDots::setup_from_rotamer_dots( RotamerDots const & rdots ) {
//	rotamer_ = rdots.rotamer();
//	if ( ! rotamer_ ) return;
//	if ( inv_dots_.size() < rdots.get_num_atoms() ) {
//		inv_dots_.resize( rdots.get_num_atoms(), utility::vector1< ObjexxFCL::ubyte >( RotamerDots::num_bytes_, ObjexxFCL::ubyte(0) )  );
//	}
//	rdots.invert_to_boolmasks( inv_dots_ );
//	radii_ = rdots.get_radii();
//}
//
//void
//InvRotamerDots::setup_from_rotamer_dots(
//	RotamerDots const & rdots,
//	utility::vector1< Size > const & ats_to_update
//)
//{
//	rotamer_ = rdots.rotamer();
//	if ( ! rotamer_ ) return;
//	if ( inv_dots_.size() < rdots.get_num_atoms() ) {
//		inv_dots_.resize( rdots.get_num_atoms(), utility::vector1< ObjexxFCL::ubyte >( RotamerDots::num_bytes_, ObjexxFCL::ubyte(0) )  );
//	}
//	rdots.invert_to_boolmasks( inv_dots_, ats_to_update );
//	radii_ = rdots.get_radii();
//}
//
//
//core::conformation::ResidueCOP
//InvRotamerDots::rotamer() const
//{
//	return rotamer_;
//}
//
/// @brief Is the intersection between two atoms on this inv-rotamer-dots object exposed?
//bool
//InvRotamerDots::atom_overlap_is_exposed( Size at1, Size at2 ) const {
//	assert( rotamer_ );
//	return overlap_exposed( rotamer_->atom( at1 ), inv_dots_[ at1 ], rotamer_->atom( at2 ), inv_dots_[ at2 ] );
//}
//
/// @brief Is the intersection between one atom on this inv-rotamer-dots object,
/// and one atom on another inv-rotamer-dots object exposed?
//bool
//InvRotamerDots::atom_overlap_is_exposed(
//	Size at_this,
//	InvRotamerDots const & other,
//	Size at_other
//) const
//{
//	return overlap_exposed( rotamer_->atom( at_this ), inv_dots_[ at_this ], other.rotamer_->atom( at_other ), other.inv_dots_[ at_other ] );
//}
//
//
//bool InvRotamerDots::dot_exposed( Size atomid, Size dot_index ) const {
//	assert( dot_index > 0 && dot_index <= 162 );
//	dot_index -= 1;
//	Size const which_byte = dot_index / 8;
//	Size const which_bit  = dot_index - which_byte * 8 ;
//	return unpack_ubyte( inv_dots_[ atomid ][ which_byte + 1 ], which_bit );
//}
//
//void InvRotamerDots::write_exposed_dots_to_kinemage(
//	std::ostream & ostr,
//	bool group
//) const
//{
//	if ( ! rotamer_ ) return;
//
//	if ( group ) {
//		ostr << "@group { invdots " << rotamer_->seqpos() << "} dominant\n";
//	}
//
//	for ( Size ii = 1; ii <= rotamer_->nheavyatoms(); ++ii ) {
//		Real const ii_rad = (*radii_)[ rotamer_->atom(ii).type() ] + RotamerDots::probe_radius_;
//		write_sphere_list_uv1( ostr, rotamer_->atom_name( ii ), "gray", rotamer_->xyz( ii ), ii_rad, inv_dots_[ ii ] );
//	}
//}
//
//void
//InvRotamerDots::write_circle_intersection_mask_to_kinemage(
//	std::ostream & ostr,
//	Size const atom_this,
//	InvRotamerDots const & invdots_other,
//	Size const atom_other,
//	bool group
//) const
//{
//	core::conformation::Atom const & at1( rotamer_->atom( atom_this ) );
//	core::conformation::Atom const & at2( invdots_other.rotamer_->atom( atom_other ) );
//
//	Real const rad1 = (*radii_)[ at1.type() ] + RotamerDots::probe_radius_;
//	Real const rad2 = (*radii_)[ at2.type() ] + RotamerDots::probe_radius_;
//
//	Real const dist_sq = at1.xyz().distance_squared( at2.xyz() );
//	assert( dist_sq < (rad1 + rad2) * (rad1 + rad2) );
//	Real const distance = std::sqrt( dist_sq );
//
//	Real const step_size1 = rad1 * 0.02;
//	Real const step_size2 = rad2 * 0.02;
//
//	Size const nsteps1 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size1 ));
//	Size const nsteps2 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size2 ));
//
//	int degree_of_overlap1, degree_of_overlap1_stepped, degree_of_overlap2, degree_of_overlap2_stepped;
//	int aphi_1_2, aphi_2_1;
//	int theta_1_2, theta_2_1;
//
//	//ronj this block represents the amount of surface area covered up on atom1 by atom2
//	core::scoring::get_overlap( rad1, rad2, distance, degree_of_overlap1 );
//	core::scoring::get_2way_orientation( at1.xyz(), at2.xyz(), aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );
//
//	utility::vector1< ObjexxFCL::ubyte > ring1( 21, ObjexxFCL::ubyte(0) );
//	utility::vector1< ObjexxFCL::ubyte > ring2( 21, ObjexxFCL::ubyte(0) );
//
//	utility::vector1< ObjexxFCL::ubyte > hit_ring1( 21, ObjexxFCL::ubyte(0) );
//	utility::vector1< ObjexxFCL::ubyte > hit_ring2( 21, ObjexxFCL::ubyte(0) );
//
//	Size closest_dot1 = (*RotamerDots::lg_angles_)( aphi_1_2, theta_1_2 );
//	if ( degree_of_overlap1 + nsteps1 > 100 ) {
//		degree_of_overlap1_stepped = 100;
//	} else {
//		degree_of_overlap1_stepped = degree_of_overlap1 + nsteps1;
//	}
//
//	int const masknum1a = ( closest_dot1 * 100 ) + degree_of_overlap1;
//	int const masknum1b = ( closest_dot1 * 100 ) + degree_of_overlap1_stepped;
//
//	// so we take two "offsets" into the "masks" table: 1) the one that normally gets used to figure out which dots are covered
//	// by the neighboring atom, and 2) one that's one "step" in, which represents the "ring" of dots that's just past the ones
//	// that are covered by the neighboring atom. if we then take the inverse (negate) the normally used mask, we'll get 0's
//	// whereever there are dots covered by the other atom (instead of 1's) and 1's everywhere else. if we logical AND that
//	// result with the dots that one "step" in, we'll get 1's at just the ring of dots that's next to the ones that are covered.
//
//	// if we then logical AND the ring of dots with all of the exposed dots on this atom, we can determine if there are any
//	// exposed dots adjacent to the intersection circle.
//
//	for ( Size bb = 1, bblia = (*RotamerDots::lg_masks_).index( bb, masknum1a ), bblib = (*RotamerDots::lg_masks_).index( bb, masknum1b );
//			bb <= RotamerDots::num_bytes_; ++bb, ++bblia, ++bblib ) {
//		ring1[ bb ] = (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ];
//		hit_ring1[ bb ] = inv_dots_[ atom_this ][ bb ] & ( (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ] );
//	}
//	std::cout << "exposed dots: ";
//	utility::vector1< ObjexxFCL::ubyte > exposed_copy1 = inv_dots_[ atom_this ];
//	print_dot_bit_string( exposed_copy1 );
//	std::cout << "ring1: ";
//	print_dot_bit_string( ring1 );
//	std::cout << "hit_ring1: ";
//	print_dot_bit_string( hit_ring1 );
//
//	//ronj the amount of surface area covered up on atom2 by atom1
//	core::scoring::get_overlap( rad2, rad1, distance, degree_of_overlap2 );
//
//	Size closest_dot2 = (*RotamerDots::lg_angles_)( aphi_2_1, theta_2_1 );
//	if ( degree_of_overlap2 + nsteps2 > 100 ) {
//		degree_of_overlap2_stepped = 100;
//	} else {
//		degree_of_overlap2_stepped = degree_of_overlap2 + nsteps2;
//	}
//
//	int const masknum2a = ( closest_dot2 * 100 ) + degree_of_overlap2;
//	int const masknum2b = ( closest_dot2 * 100 ) + degree_of_overlap2_stepped;
//
//	for ( Size bb = 1, bblia = (*RotamerDots::lg_masks_).index( bb, masknum2a ), bblib = (*RotamerDots::lg_masks_).index( bb, masknum2b );
//			bb <= RotamerDots::num_bytes_; ++bb, ++bblia, ++bblib ) {
//		ring2[ bb ] = (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ];
//		hit_ring2[ bb ] = invdots_other.inv_dots_[ atom_other ][ bb ] & ( (*RotamerDots::lg_masks_)[ bblib ] & ~ (*RotamerDots::lg_masks_)[ bblia ] );
//	}
//	std::cout << "exposed dots: ";
//	utility::vector1< ObjexxFCL::ubyte > exposed_copy2 = invdots_other.inv_dots_[ atom_other ];
//	print_dot_bit_string( exposed_copy2 );
//	std::cout << "ring2: ";
//	print_dot_bit_string( ring2 );
//	std::cout << "hit_ring2: ";
//	print_dot_bit_string( hit_ring2 );
//
//	if ( group ) {
//		ostr << "@group { invdots " << rotamer_->seqpos() << " " << rotamer_->atom_name( atom_this ) << "} dominant\n";
//	}
//	write_sphere_list_uv1( ostr, rotamer_->atom_name( atom_this ), "gray", at1.xyz(), rad1, ring1 );
//
//	if ( group ) {
//		ostr << "@group { invdots " << invdots_other.rotamer_->seqpos() << " " << invdots_other.rotamer_->atom_name( atom_other ) << "} dominant\n";
//	}
//	write_sphere_list_uv1( ostr,  invdots_other.rotamer_->atom_name( atom_other ), "gray", at2.xyz(), rad2, ring2 );
//
//	if ( group ) {
//		ostr << "@group { olap " << rotamer_->seqpos() << " " << rotamer_->atom_name( atom_this ) << "} dominant\n";
//	}
//	write_sphere_list_uv1( ostr, rotamer_->atom_name( atom_this ), "blue", at1.xyz(), rad1, hit_ring1 );
//
//	if ( group ) {
//		ostr << "@group { olap " << invdots_other.rotamer_->seqpos() << " " << invdots_other.rotamer_->atom_name( atom_other ) << "} dominant\n";
//	}
//	write_sphere_list_uv1( ostr,  invdots_other.rotamer_->atom_name( atom_other ), "blue", at2.xyz(), rad2, hit_ring2 );
//
//}
//
//
//bool
//InvRotamerDots::overlap_exposed(
//	core::conformation::Atom const & at1,
//	utility::vector1< ObjexxFCL::ubyte > const & at1exposed_dots,
//	core::conformation::Atom const & at2,
//	utility::vector1< ObjexxFCL::ubyte > const & at2exposed_dots
//) const
//{
//
//	Real const rad1 = (*radii_)[ at1.type() ] + RotamerDots::probe_radius_;
//	Real const rad2 = (*radii_)[ at2.type() ] + RotamerDots::probe_radius_;
//
//	Real const dist_sq = at1.xyz().distance_squared( at2.xyz() );
//	assert( dist_sq <= (rad1 + rad2) * (rad1 + rad2) );
//	Real const distance = std::sqrt( dist_sq );
//
//	Real const step_size1 = rad1 * 0.02;
//	Real const step_size2 = rad2 * 0.02;
//
//	Size const nsteps1 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size1 ));
//	Size const nsteps2 = (Size)(ceil( max_dist_from_dot_to_intersection / step_size2 ));
//
//	int degree_of_overlap1, degree_of_overlap2;
//	int aphi_1_2, aphi_2_1;
//	int theta_1_2, theta_2_1;
//
//	//ronj this block represents the amount of surface area covered up on atom1 by atom2
//	core::scoring::get_overlap( rad1, rad2, distance, degree_of_overlap1 );
//	core::scoring::get_2way_orientation( at1.xyz(), at2.xyz(), aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );
//
//	bool at1_intersection_exposed = false;
//	Size closest_dot1 = (*RotamerDots::lg_angles_)( aphi_1_2, theta_1_2 );
//	if ( degree_of_overlap1 + nsteps1 > 100 ) {
//		for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ){
//			if ( at1exposed_dots[ bb ] ) {
//				at1_intersection_exposed = true;
//				break;
//			}
//		}
//	} else {
//		int masknum = ( closest_dot1 * 100 ) + degree_of_overlap1 + nsteps1;
//		for ( Size bb = 1, bbli = (*RotamerDots::lg_masks_).index( bb, masknum ); bb <= RotamerDots::num_bytes_; ++bb, ++bbli ){
//			if ( at1exposed_dots[ bb ] & (*RotamerDots::lg_masks_)[ bbli ] ) {
//				at1_intersection_exposed = true;
//				break;
//			}
//		}
//	}
//
//	if ( ! at1_intersection_exposed ) return false;
//
//	//ronj the amount of surface area covered up on atom2 by atom1
//	core::scoring::get_overlap( rad2, rad1, distance, degree_of_overlap2 );
//
//	bool at2_intersection_exposed = false;
//
//	Size closest_dot2 = (*RotamerDots::lg_angles_)( aphi_2_1, theta_2_1 );
//	if ( degree_of_overlap2 + nsteps2 > 100 ) {
//		for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ){
//			if ( at2exposed_dots[ bb ] ) {
//				at2_intersection_exposed = true;
//				break;
//			}
//		}
//	} else {
//		int masknum = ( closest_dot2 * 100 ) + degree_of_overlap2 + nsteps2;
//		for ( Size bb = 1, bbli = (*RotamerDots::lg_masks_).index( bb, masknum ); bb <= RotamerDots::num_bytes_; ++bb, ++bbli ){
//			if ( at2exposed_dots[ bb ] & (*RotamerDots::lg_masks_)[ bbli ] ) {
//				at2_intersection_exposed = true;
//				break;
//			}
//		}
//	}
//
//	return at2_intersection_exposed;
//}
//
///
/// @begin InvRotamerDots::print_dot_bit_string
///
/// @brief
/// Helper method I am using to confirm that the dots are being overlapped and bits are being set correctly.
///
//void InvRotamerDots::print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const {
//	for ( Size bb = 1; bb <= RotamerDots::num_bytes_; ++bb ) {
//		int bit;
//		if ( (bb-1)*8 % 16 == 0 ) std::cout << (bb-1) * 8 << ":";
//		for ( int index=7; index >= 0; index-- ) {
//			bit = ( ( (int)values[ bb ] >> index ) & 1 );
//			std::cout << bit;
//		}
//		std::cout << " ";
//	}
//	std::cout << std::endl;
//}

protocols::moves::MoverOP
LoadVarSolDistSasaCalculatorMoverCreator::create_mover() const
{
	return new LoadVarSolDistSasaCalculatorMover;
}
std::string LoadVarSolDistSasaCalculatorMoverCreator::keyname() const
{
	return "LoadVarSolDistSasaCalculatorMover";
}

LoadVarSolDistSasaCalculatorMover::LoadVarSolDistSasaCalculatorMover() :
	protocols::moves::Mover( "LoadVarSolDistSasaCalculatorMover" )
{}

LoadVarSolDistSasaCalculatorMover::~LoadVarSolDistSasaCalculatorMover()
{}

protocols::moves::MoverOP
LoadVarSolDistSasaCalculatorMover::clone() const {
	return new LoadVarSolDistSasaCalculatorMover;
}

std::string
LoadVarSolDistSasaCalculatorMover::get_name() const { return "LoadVarSolDistSasaCalculatorMover"; }

void
LoadVarSolDistSasaCalculatorMover::apply( core::pose::Pose & )
{
	using core::pose::metrics::CalculatorFactory;
	if ( CalculatorFactory::Instance().check_calculator_exists( "bur_unsat_calc_default_sasa_calc" ) ) {
		CalculatorFactory::Instance().remove_calculator( "bur_unsat_calc_default_sasa_calc" );
	}
	CalculatorFactory::Instance().register_calculator( "bur_unsat_calc_default_sasa_calc", new VarSolDistSasaCalculator );
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void LoadVarSolDistSasaCalculatorMover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & )
{}



} // vardist_solaccess
} // devel

