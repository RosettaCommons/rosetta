// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/VDW_Grid.hh
/// @brief
/// @detailed
/// @author Caleb Geniesse, geniesse@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_rna_checker_VDW_Grid_HH
#define INCLUDED_protocols_stepwise_modeler_rna_checker_VDW_Grid_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/rna/checker/VDW_Grid.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


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
	setup( int const & bin_max ) const;

	void
	reset() const;

	void
	reset( utility::vector1< Atom_Bin > const & occupied_xyz_bins_ ) const;

	core::Size
	size() const;

	bool
	get_bin( int const & x, int const & y, int const & z ) const;

	void
	set_bin( int const & x, int const & y, int const & z, bool const & value ) const;

	bool
	get_xyz_bin( Atom_Bin const & xyz_bin ) const;

	void
	set_xyz_bin( Atom_Bin const & xyz_bin, bool const & value ) const;

	bool
	is_occupied() const;

private:

	mutable utility::vector1< utility::vector1< utility::vector1< bool > > > bins_;
	mutable bool is_occupied_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_stepwise_modeler_rna_checker_VDW_Grid )
#endif // SERIALIZATION


#endif
