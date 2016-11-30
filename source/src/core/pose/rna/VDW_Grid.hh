// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/VDW_Grid.hh
/// @brief
/// @details
/// @author Caleb Geniesse, geniesse@stanford.edu


#ifndef INCLUDED_core_pose_rna_VDW_Grid_HH
#define INCLUDED_core_pose_rna_VDW_Grid_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/pose/rna/VDW_Grid.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace rna {


struct Atom_Bin {
	int x;
	int y;
	int z;
};


class VDW_Grid : public utility::pointer::ReferenceCount {

public:

	VDW_Grid();

	VDW_Grid( VDW_Grid const & src );

	virtual ~VDW_Grid();

	void
	setup( int const bin_max ) const;

	void
	reset() const;

	void
	reset( utility::vector1< Atom_Bin > const & occupied_xyz_bins_ ) const;

	core::Size
	size() const;

	bool
	get_bin( int const x, int const y, int const z ) const;

	void
	set_bin( int const x, int const y, int const z, bool const value ) const;

	bool
	get_xyz_bin( Atom_Bin const & xyz_bin ) const;

	void
	set_xyz_bin( Atom_Bin const & xyz_bin, bool const value ) const;

	bool
	is_occupied() const;

	void
	set_bin_max( int const value ) const;

	int
	get_bin_max() const;

	void
	set_atom_bin_size( core::Real const value ) const;

	core::Real
	get_atom_bin_size() const;

	void
	set_bin_offset( int const value ) const;

	int
	get_bin_offset() const;

	void
	set_ref_xyz( numeric::xyzVector< core::Real > const & value ) const;

	numeric::xyzVector< core::Real >
	get_ref_xyz() const;

private:

	mutable utility::vector1< utility::vector1< utility::vector1< bool > > > bins_;
	mutable bool is_occupied_;
	mutable int bin_max_;
	mutable core::Real atom_bin_size_;
	mutable int bin_offset_;
	mutable numeric::xyzVector< core::Real > ref_xyz_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //rna
} //pose
} //core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pose_rna_VDW_Grid )
#endif // SERIALIZATION


#endif
