// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/rna/movers/RNAIdealizeMover.hh
/// @brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
/// @author Andy Watkins (amw579@nyu.edu)

#ifndef INCLUDED_protocols_farna_RNAIdealizeMover_hh
#define INCLUDED_protocols_farna_RNAIdealizeMover_hh

// Unit headers
#include <protocols/rna/movers/RNAIdealizeMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/AtomID.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace rna {
namespace movers {

///@brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
class RNAIdealizeMover : public protocols::moves::Mover {

public:

	RNAIdealizeMover();

	RNAIdealizeMover( Size const iterations, bool const noise, bool const final_minimize, core::Real const ang_significance_threshold, bool const handle_suites ):
		protocols::moves::Mover( RNAIdealizeMover::mover_name() ),
		iterations_( iterations ),
		noise_( noise ),
		final_minimize_( final_minimize ),
		ang_significance_threshold_( ang_significance_threshold ),
		handle_suites_( handle_suites )
	{}

	// copy constructor (not needed unless you need deep copies)
	//RNAIdealizeMover( RNAIdealizeMover const & src );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~RNAIdealizeMover();

	// XRW TEMP  static std::string
	// XRW TEMP  class_name();

public:
	// mover virtual API
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	// XRW TEMP  virtual std::string
	// XRW TEMP  get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//RNAIdealizeMover & operator=( RNAIdealizeMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	Size get_iterations() const { return iterations_; }

	void set_iterations( Size const iterations ) { iterations_ = iterations; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	void perturb_pose( core::pose::Pose & pose ) const;
	void constrain_to_ideal( core::pose::Pose & pose, utility::vector1< core::Size > const & bad_suite_res ) const;

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
	core::Real ang_significance_threshold_ = 5;
	bool handle_suites_ = false;
};

std::ostream &
operator<<( std::ostream & os, RNAIdealizeMover const & mover );

} //movers
} //rna
} //protocols

#endif //protocols/farna_RNAIdealizeMover_hh
