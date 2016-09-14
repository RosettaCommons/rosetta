// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/farna/RNAIdealizeMover.hh
/// @brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_RNAIdealizeMover_hh
#define INCLUDED_protocols_farna_RNAIdealizeMover_hh

// Unit headers
#include <protocols/farna/RNAIdealizeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/AtomID.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace farna {

///@brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
class RNAIdealizeMover : public protocols::moves::Mover {

public:

	RNAIdealizeMover();
	
	RNAIdealizeMover( Size const iterations, bool const noise, bool const final_minimize ):
		protocols::moves::Mover( RNAIdealizeMover::class_name() ),
		iterations_( iterations ),
		noise_( noise ),
		final_minimize_( final_minimize )
	{}

	// copy constructor (not needed unless you need deep copies)
	//RNAIdealizeMover( RNAIdealizeMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~RNAIdealizeMover();

	static std::string
	class_name();

public:
	// mover virtual API
	virtual void
	apply( core::pose::Pose & pose );

	virtual void
	show( std::ostream & output = std::cout ) const;

	virtual std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	//RNAIdealizeMover & operator=( RNAIdealizeMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;
	
	Size get_iterations() const { return iterations_; }
	
	void set_iterations( Size const iterations ) { iterations_ = iterations; }
	
private:
	
	void perturb_pose( core::pose::Pose & pose ) const;
	void constrain_to_ideal( core::pose::Pose & pose ) const;
	
	void
	add_bond_angle_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::id::AtomID const & atom_id3,
		core::pose::Pose & pose
	) const;
	
	void 
	add_bond_constraint(
		core::id::AtomID const & atom_id1,
		core::id::AtomID const & atom_id2,
		core::pose::Pose & pose
	) const;
		
	core::pose::Pose ref_pose_;

	/// @brief Number of iterations made (1/n idealization per iteration).
	Size iterations_ = 100;
	bool noise_ = false;
	bool final_minimize_ = false;
};

std::ostream &
operator<<( std::ostream & os, RNAIdealizeMover const & mover );

} //protocols
} //farna

#endif //protocols/farna_RNAIdealizeMover_hh
