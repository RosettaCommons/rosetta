// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/scientific_tests/PDBDiagnosticMover.hh
/// @brief does some very lightweight modeling on a PDB.  Meant to be run against the whole PDB as a scientific test
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_devel_scientific_tests_PDBDiagnosticMover_HH
#define INCLUDED_devel_scientific_tests_PDBDiagnosticMover_HH

// Unit headers
#include <devel/scientific_tests/PDBDiagnosticMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

#include <protocols/jd2/Job.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace devel {
namespace scientific_tests {

///@brief does some very lightweight modeling on a PDB.  Meant to be run against the whole PDB as a scientific test.  It is meant to only ever be run on one PDB and then shut down; it pairs with JD0 in the tools repo.  It uses JD2 but it can't be properly job distributed by MPI, etc, because there will be segfaults due to wonky weird stuff in the PDB.
class PDBDiagnosticMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	PDBDiagnosticMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	PDBDiagnosticMover( PDBDiagnosticMover const & src );

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~PDBDiagnosticMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	//PDBDiagnosticMover & operator=( PDBDiagnosticMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	bool
	reinitialize_for_each_job() const override { return true; }

	bool
	reinitialize_for_new_input() const override { return true; }


private: // methods

	/// @brief loops over all residues and all ResidueProperties - or the is_DNA, etc, in the ResidueType, anyway - and dumps that into the scorefile output.  This implementation is clearly tied to JD2; if someone wants to change that go ahead, but since this tool is un-job-distributable in the broad sense due to segfaults in the broad PDB test it's a waste of time.
	void residue_type_statistics( core::pose::Pose const & pose, protocols::jd2::JobOP job_me, core::Size const nres );


private: // data

};

std::ostream &
operator<<( std::ostream & os, PDBDiagnosticMover const & mover );

} //devel
} //scientific_tests

#endif //devel_scientific_tests_PDBDiagnosticMover_HH
