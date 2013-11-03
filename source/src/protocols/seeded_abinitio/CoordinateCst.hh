// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/seeded_abinitio/CoordinateCst.hh
/// @brief
/// @author Eva-Maria Strauch (evas01@u.washington.edu)

#ifndef INCLUDED_protocols_seeded_abinitio_CoordinateCst_hh
#define INCLUDED_protocols_seeded_abinitio_CoordinateCst_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace seeded_abinitio {

class CoordinateCst : public protocols::moves::Mover
{
public:
	typedef core::pose::Pose Pose;

public:
	CoordinateCst();

	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	void parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	virtual ~CoordinateCst();

private:
	/// stedev for the constraint function
	core::Real stddev_;
	/// residue number to which the coordinate constraints are anchored to
	/// stored as string, parsed at runtime
	std::string anchor_res_;
	/// whether to use the jump to determine the anchor residue
	bool use_jumps_;
	/// vector with list of residues numbers that will getcoordinate constraints.
	/// they are parsed during runtime to be compatible with pose length changes
	std::string unparsed_residue_;
	///container for spans of residues that will get coordinate constraints
	utility::vector1 < std::pair < std::string, std::string > >  span_vector_;
	///which jump atom pair to choose when the anchor residue should be part of the jump atoms
	core::Size jump_;

	/// which atom on the anchor residue to place the constraint onto
	std::string anchor_atom_id_;

	///which atom id on the moving residue to place the constraint onto
	std::string atom_id_;

};
} //seeded_abinitio
} // protocols

#endif
