// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/Basis.hh
/// @brief Class for storing residue information needed to generate alignments
/// @author Minnie Langlois (minnie@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_Basis_hh
#define INCLUDED_protocols_sewing_data_storage_Basis_hh

#include <protocols/sewing/data_storage/Basis.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>

namespace protocols {
namespace sewing {
namespace data_storage {


///@brief Class for storing residue information needed to generate alignments
class Basis : public utility::pointer::ReferenceCount {

public:

	Basis();
	Basis( core::Size segment_id, core::Size resnum );
	Basis(Basis const & src);

	virtual ~Basis();

	BasisOP
	clone() const;

	friend
	bool
	operator< (Basis const & a, Basis const & b);

	void
	segment_id( core::Size segment_id);

	void
	resnum( core::Size resnum );

	core::Size
	segment_id() const { return segment_id_; }

	core::Size
	resnum() const { return resnum_; }

	Basis &
	operator=( Basis const & other );

private:
	core::Size segment_id_;
	core::Size resnum_;

};

typedef std::pair< Basis, Basis > BasisPair;

} //data_storage
} //sewing
} //protocols



#endif //INCLUDED_protocols_sewing_data_storage_Basis_hh





