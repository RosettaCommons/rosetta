// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaInterfaceFinder_hh
#define INCLUDED_protocols_dna_DnaInterfaceFinder_hh

#include <protocols/dna/DnaInterfaceFinder.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class DnaNeighbor {
public:
	DnaNeighbor() : close_(false), contact_(false) {}
	void close( bool val ) { close_ |= val; }
	void contact( bool val ) { contact_ |= val; }
	bool close() const { return close_; }
	bool contact() const { return contact_; }
private:
	bool close_, contact_;
};

typedef std::map< core::Size, DnaNeighbor > DnaNeighbors;

class DnaInterfaceFinder : public utility::pointer::ReferenceCount {
public:

	DnaInterfaceFinder(
		core::Real close_threshold = 10 * 10,
		core::Real contact_threshold = 3.7 * 3.7,
		core::Real z_cutoff = 10.,
		bool base_only = false
	) : utility::pointer::ReferenceCount(),
		close_threshold_( close_threshold ),
		contact_threshold_( contact_threshold ),
		z_cutoff_( z_cutoff ),
		base_only_( base_only ),
		initialized_( false )
	{}

	virtual ~DnaInterfaceFinder(){}

	DnaInterfaceFinderOP clone() const { return DnaInterfaceFinderOP( new DnaInterfaceFinder(*this) ); }

	void determine_protein_interface( core::pose::Pose const & pose );

	void
	determine_interface(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & protein_positions,
		utility::vector1< core::Size > const & dna_positions
	);

	void
	determine_interface( core::pose::Pose const & pose );

	void
	determine_protein_interface(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & protein_positions,
		utility::vector1< core::Size > const & dna_positions
	);

	void
	determine_dna_interface(
		core::pose::Pose const & pose,
		utility::vector1< core::Size > const & protein_positions,
		utility::vector1< core::Size > const & dna_positions
	);

	bool initialized() const { return initialized_; }

	DnaNeighbors const & protein_neighbors() const;
	DnaNeighbors const & dna_neighbors() const;

private:
	core::Real close_threshold_, contact_threshold_, z_cutoff_;
	bool base_only_;
	bool initialized_;
	DnaNeighbors protein_neighbors_;
	DnaNeighbors dna_neighbors_;

};

} // namespace dna
} // namespace protocols

#endif
