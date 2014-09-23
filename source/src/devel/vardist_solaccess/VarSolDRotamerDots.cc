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

//#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/util.hh>
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
#include <algorithm>

#include <core/chemical/AtomType.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
static thread_local basic::Tracer TR( "devel.vardist_solaccess" );
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

//Real VarSolDRotamerDots::probe_radius_ = 1.4;

// No longer static -KH
//bool VarSolDRotamerDots::sasa_arrays_initialized_ = false;

/// @detailed
/// One RotamerDots object get allocated for every state of a first class IG Node, for all first class IG Nodes of a
/// protein being designed. That's potentially a very large number of states. This class should only hold the information
/// it needs to hold to do its job.
///
/*
VarSolDRotamerDots::VarSolDRotamerDots(
	conformation::ResidueCOP rotamer,
	bool all_atom,
	Real probe_radius,
	Real wobble
) :
	rotamer_(rotamer),
	probe_radius_(probe_radius),
	wobble_(wobble),
	lg_masks_( 0 )
{
	// this work now done in constructor, no need to check
	//if ( ! sasa_arrays_initialized_ ) {
	//	initialize_sasa_arrays();
	//	//that function will set sasa_arrays_initialized_ to true;
	//}

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
 */

VarSolDRotamerDots::VarSolDRotamerDots(
	conformation::ResidueCOP rotamer,
	VarSolDistSasaCalculator const & vsasa_calc
) :
	rotamer_(rotamer),
	radii_(vsasa_calc.radii_),
	msas_radii_(vsasa_calc.msas_radii_),
	coll_radii_(vsasa_calc.coll_radii_),
	int_radii_(vsasa_calc.int_radii_),
	int_radii_sum_(vsasa_calc.int_radii_sum_),
	int_radii_sum2_(vsasa_calc.int_radii_sum2_),
	lg_angles_(vsasa_calc.lg_angles_),
	lg_masks_(vsasa_calc.lg_masks_),
	num_bytes_(vsasa_calc.num_bytes_),
	polar_expansion_radius_(vsasa_calc.polar_expansion_radius_)
{
	num_atoms_ = rotamer->natoms();
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
	rotamer_(rhs.rotamer_),
	radii_(rhs.radii_),
	msas_radii_(rhs.msas_radii_),
	coll_radii_(rhs.coll_radii_),
	int_radii_(rhs.int_radii_),
	int_radii_sum_(rhs.int_radii_sum_),
	int_radii_sum2_(rhs.int_radii_sum2_),
	lg_angles_(rhs.lg_angles_),
	lg_masks_(rhs.lg_masks_),
	num_bytes_(rhs.num_bytes_),
	polar_expansion_radius_(rhs.polar_expansion_radius_),
	num_atoms_(rhs.num_atoms_),
	atom_coverage_(rhs.atom_coverage_)
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
/// @begin RotamerDots::get_atom_collision_radius
///
/// @brief
/// Returns the collision radius for the passed in atom type.
///
/// @detailed
/// Many of the functions in this class iterate over 1 .. num_atoms_.
/// That's not the same thing as an atom type index which is what the radii vector is indexed with. So before we can return
/// the radius, we have to convert the passed in atom_index into the right atom in the residue and then use that to get the
/// right type.
///
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

core::Size VarSolDRotamerDots::nshells_for_atom( core::Size atom_index ) const
{
	return radii_[ rotamer_->atom(atom_index).type() ].size();
}

core::Real VarSolDRotamerDots::shell_radius_for_atom( core::Size atom_index, core::Size shell_index ) const
{
	return radii_[ rotamer_->atom( atom_index ).type() ][ shell_index ];
}

core::Size VarSolDRotamerDots::ndots() const
{
	return atom_coverage_[ num_atoms_ ][ 1 ].get_total_dots();
}

bool
VarSolDRotamerDots::get_dot_covered(
	core::Size atom_index,
	core::Size shell_index,
	core::Size dot_index
) const
{
	return atom_coverage_[ atom_index ][ shell_index ].get_dot_covered( dot_index );
}

/// @brief
/// Initializes the pointers to the angles and masks FArrays used by sasa.cc and inits the dot sphere coordinates.
///
/// @detailed
/// This call should only occur once (when the first RotamerDots object get constructed) and never again.
///
void VarSolDistSasaCalculator::initialize_sasa_arrays() {

	lg_angles_ = ( & core::scoring::sasa::get_legrand_sasa_angles() );
	lg_masks_  = ( & core::scoring::sasa::get_legrand_sasa_masks()  );

	using namespace core::chemical;

	AtomTypeSetCOP atset = ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	core::Size const SASA_RADIUS_INDEX = atset->extra_parameter_index( "REDUCE_SASA_RADIUS" );

	radii_.resize( atset->n_atomtypes() );
	msas_radii_.resize( atset->n_atomtypes() );
	for ( Size ii = 1; ii <= atset->n_atomtypes(); ++ii ) {
		AtomType const & iiattype = (*atset)[ ii ];
		msas_radii_[ii] = iiattype.extra_parameter( SASA_RADIUS_INDEX );

		Size steps=0;
		Real coll_radius = 0;
		Real int_radius = 0;
		if ( iiattype.element() == "O" ) {
			steps = 5;
			//coll_radius = std::max(0., 1.2 + probe_radius_ - wobble_);
			//int_radius = std::max(0., 3.0 + wobble_);
			coll_radius = 2.6;
			int_radius = 3.0;
		} else if ( iiattype.element() == "N" ) {
			steps = 5;
			//coll_radius = std::max(0., 1.3 + probe_radius_ - wobble_);
			//int_radius = std::max(0., 3.1 + wobble_);
			coll_radius = 2.7;
			int_radius = 3.1;
		} else if ( iiattype.element() == "C" ) {
			steps = 1;
			if ( iiattype.atom_type_name() == "COO" || iiattype.atom_type_name() == "CObb" ) {
				//coll_radius = std::max(0., 1.65 + probe_radius_ - wobble_);
				//int_radius = std::max(0., 3.05 + wobble_);
				coll_radius = 3.05;
			} else {
				//coll_radius = std::max(0., 1.75 + probe_radius_ - wobble_);
				coll_radius = 3.15;
			}
		} else if ( iiattype.element() == "S" ) {
			steps = 1;
			//coll_radius = std::max(0., 1.85 + probe_radius_ - wobble_);
			//int_radius = std::max(0., 3.25 + wobble_);
			coll_radius = 3.25;
		} else if ( iiattype.element() == "P" ) {
			steps = 1;
			//coll_radius = std::max(0., 1.9 + probe_radius_ - wobble_);
			//int_radius = std::max(0., 3.3 + wobble_);
			coll_radius = 3.3;
		} else if ( iiattype.element() == "H" ) {
			if ( iiattype.atom_type_name() == "Hpol" || iiattype.atom_type_name() == "HNbb" ) {
				steps = 5;
				//coll_radius = std::max(0., 0.3 + probe_radius_ - wobble_);
				//int_radius = std::max(0., 2.1 + wobble_);
				coll_radius = 1.7;
				int_radius = 2.1;
			} else if ( iiattype.atom_type_name() == "Haro" ) {
				steps = 1;
				//coll_radius = std::max(0., 1.0 + probe_radius_ - wobble_);
				//int_radius = std::max(0., 2.4 + wobble_);
				coll_radius = 2.4;
			} else {
				steps = 1;
				//coll_radius = std::max(0., 1.1 + probe_radius_ - wobble_);
				//int_radius = std::max(0., 2.5 + wobble_);
				coll_radius = 2.5;
			}
		}

		/*
        if ( iiattype.element() == "O" ) {
            steps = 5;
            coll_radius = 2.413468;
            int_radius = 3.078613;
        } else if ( iiattype.element() == "N" ) {
            steps = 5;
            coll_radius = 2.572567;
            int_radius = 3.265344;
        } else if ( iiattype.element() == "C" ) {
            steps = 1;
            coll_radius = 2.983538;
        } else if ( iiattype.element() == "S" ) {
            steps = 1;
            coll_radius = 2.767480;
        } else if ( iiattype.element() == "P" ) {
            steps = 1;
            // not statistically derived
            coll_radius = std::max(0., 1.9 + probe_radius_ - wobble_);
        } else if ( iiattype.element() == "H" ) {
            if ( iiattype.atom_type_name() == "Hpol" || iiattype.atom_type_name() == "HNbb" ) {
                steps = 5;
                coll_radius = 1.44593;
                int_radius = 2.45833;
            } else {
                steps = 1;
                coll_radius = 2.072345;
            }
        }
		 */

		// coll_radius is collision radius, int_radius is interaction radius
		// create evenly spaced shells inbetween, number of shells total = steps
		radii_[ ii ].resize( steps );
		if (steps > 0) {
			radii_[ii][1] = std::max(0., coll_radius);
		}
		for (Size j=2; j<steps; j++) {
			radii_[ii][j] = coll_radius + ( (int_radius - coll_radius) * ((j-1) / static_cast<Real>(steps-1)) );
		}
		if (steps > 1) {
			radii_[ii][steps] = int_radius;
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
) const
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
	//core::scoring::sasa::get_legrand_atomic_overlap( at1_radius, at2_radius, distance, degree_of_overlap );
	core::scoring::sasa::get_legrand_2way_orientation( at1_xyz, at2_xyz, aphi_1_2, theta_1_2, aphi_2_1, theta_2_1, distance );
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
				core::scoring::sasa::get_legrand_atomic_overlap( radii_[ attype2 ][ jj ], radii_[ attype1 ][ ii ], distance, degree_of_overlap );
				Size masknum = ( closest_dot2 * 100 ) + degree_of_overlap;
				for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
					at_other_coverage[ jj ][ bb ] |= (*lg_masks_)[ bbli ];
				}
			}
			if ( jj == 1 ) {
				int degree_of_overlap;
				core::scoring::sasa::get_legrand_atomic_overlap( radii_[ attype1 ][ ii ], radii_[ attype2 ][ jj ], distance, degree_of_overlap );
				Size masknum = ( closest_dot1 * 100 ) + degree_of_overlap;
				for ( Size bb = 1, bbli = (*lg_masks_).index( bb, masknum ); bb <= num_bytes_; ++bb, ++bbli ) {
					at_this_coverage[ ii ][ bb ] |= (*lg_masks_)[ bbli ];
				}
			}
		}
	}
	return true;
}


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


VarSolDistSasaCalculator::VarSolDistSasaCalculator():
	// probe_radius_(0),
	// wobble_(0),
	num_bytes_(21),
	lg_masks_(0),
	lg_angles_(0)
{
	initialize_sasa_arrays();
}

core::pose::metrics::PoseMetricCalculatorOP
VarSolDistSasaCalculator::clone() const{
	return new VarSolDistSasaCalculator();
}

void VarSolDistSasaCalculator::set_atom_type_radii(std::string atype_name, Real coll_radius, Real int_radius, Size nshells) {

	using namespace core::chemical;
	AtomTypeSetCOP atset = ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	Size i = atset->atom_type_index(atype_name);
	radii_[i].resize( nshells );
	if (nshells > 0) {
		radii_[i][1] = std::max(0., coll_radius);
	}
	for (Size j=2; j<nshells; j++) {
		radii_[i][j] = coll_radius + ( (int_radius - coll_radius) * ((j-1) / static_cast<Real>(nshells-1)) );
	}
	if (nshells > 0) {
		radii_[i][nshells] = int_radius;
	}

	coll_radii_[i] = coll_radius;
	int_radii_[i] = int_radius;
	for (Size j = 1; j <= atset->n_atomtypes(); j++) {
		Real other_coll_radius = coll_radii_[j];
		Real other_int_radius = int_radii_[j];
		int_radii_sum_[i][j] = std::max(coll_radius + other_int_radius, int_radius + other_coll_radius);
		int_radii_sum_[j][i] = int_radii_sum_[i][j];
		int_radii_sum2_[i][j] = std::pow(int_radii_sum_[i][j], 2);
		int_radii_sum2_[j][i] = int_radii_sum2_[i][j];
	}
}


void VarSolDistSasaCalculator::set_element_radii(std::string elem, Real coll_radius, Real int_radius, Size nshells) {

	using namespace core::chemical;
	AtomTypeSetCOP atset = ChemicalManager::get_instance()->atom_type_set( FA_STANDARD );
	// find index of atomtype
	for (Size i=1; i <= atset->n_atomtypes(); i++) {
		if ( (*atset)[i].element() == elem ) {
			set_atom_type_radii( (*atset)[i].atom_type_name(), coll_radius, int_radius, nshells);
		}
	}
}

id::AtomID_Map< core::Real >
VarSolDistSasaCalculator::calculate(const pose::Pose& pose) {
	recompute(pose);
	return atom_sasa_;
}

void
VarSolDistSasaCalculator::lookup( std::string const & key, basic::MetricValueBase * valptr ) const
{
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

	rotamer_dots_vec_.resize( this_pose.total_residue() );
	residue_sasa_.resize( this_pose.total_residue() );
	// TR << "Initializing vSASA arrays with probe radius = " << probe_radius_ << " and wobble = " << wobble_ << std::endl;
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		rotamer_dots_vec_[ ii ] = new VarSolDRotamerDots( core::conformation::ResidueOP( new core::conformation::Residue(this_pose.residue( ii ) ) ), *this );
		rotamer_dots_vec_[ ii ]->increment_self_overlap();
	}
	core::scoring::EnergyGraph const & energy_graph( this_pose.energies().energy_graph() );
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		for ( core::graph::Graph::EdgeListConstIter
				iru  = energy_graph.get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph.get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			rotamer_dots_vec_[ ii ]->intersect_residues( *rotamer_dots_vec_[ (*iru)->get_second_node_ind() ] );
		}
	}
	total_sasa_ = 0.0;
	for ( Size ii = 1; ii <= this_pose.total_residue(); ++ii ) {
		residue_sasa_[ ii ] = 0.0;
		for ( Size jj = 1; jj <= this_pose.residue(ii).natoms(); ++jj ){
			core::id::AtomID at( jj, ii );
			residue_sasa_[ ii ] += atom_sasa_[ at ] = rotamer_dots_vec_[ ii ]->msas_for_atom( jj );
		}
		total_sasa_ += residue_sasa_[ ii ];
	}
}

protocols::moves::MoverOP
LoadVarSolDistSasaCalculatorMoverCreator::create_mover() const
{
	return new LoadVarSolDistSasaCalculatorMover;
}
std::string LoadVarSolDistSasaCalculatorMoverCreator::keyname() const
{
	return "LoadVarSolDistSasaCalculatorMover";
}

LoadVarSolDistSasaCalculatorMover::LoadVarSolDistSasaCalculatorMover(Real /*probe_radius*/, Real /*wobble*/) :
	protocols::moves::Mover( "LoadVarSolDistSasaCalculatorMover" )
	// probe_radius_(probe_radius),
	// wobble_(wobble)
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
	CalculatorFactory::Instance().register_calculator( "bur_unsat_calc_default_sasa_calc", core::pose::metrics::PoseMetricCalculatorOP( new VarSolDistSasaCalculator() ) );
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
