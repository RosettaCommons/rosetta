// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/membrane/benchmark/MakeCanonicalHelix.hh
/// @brief Creates an ideal a-helix from sequence given a range of residues
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelix_hh
#define INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelix_hh

// Unit headers
#include <protocols/membrane/benchmark/MakeCanonicalHelix.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>

#include <protocols/filters/Filter.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace membrane {
namespace benchmark {

///@brief Creates an ideal a-helix from sequence given a range of residues
class MakeCanonicalHelix : public protocols::moves::Mover {

public:

	MakeCanonicalHelix();

	MakeCanonicalHelix( core::Size helix_start, core::Size helix_end ); 

	// copy constructor
	MakeCanonicalHelix( MakeCanonicalHelix const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~MakeCanonicalHelix();

	virtual void
	apply( core::pose::Pose & pose );

public:

	virtual void
	show( std::ostream & output=std::cout ) const;

	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

private: 

	/// @brief Check that the start & end positons are valid given the pose
	void is_valid( core::pose::Pose & pose ); 

private:

	// Ideal Helix Dihedral Angles
	core::Real phi_; 
	core::Real psi_; 
	core::Real omega_; 

	// Position where Helix Begins and Ends
	core::Real helix_start_; 
	core::Real helix_end_; 
};

} // benchmark
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_benchmark_MakeCanonicalHelix_hh







