// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelConstraintGenerator.hh
///
/// @brief  abstract base class that generates constraints during forge loop remodelling
/// @author Florian Richter, floric@u.washington.edu, april 2009
/// @modified Tom Linsky, tlinsky@uw.edu


#ifndef INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_hh
#define INCLUDED_protocols_forge_remodel_RemodelConstraintGenerator_hh

//unit headers
#include <protocols/moves/ConstraintGenerator.hh>
#include <protocols/forge/components/VarLengthBuild.fwd.hh>
#include <protocols/forge/remodel/RemodelConstraintGenerator.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh> // WIN32 INCLUDE
#endif
#include <core/id/SequenceMapping.fwd.hh>

//utility headers

namespace protocols {
namespace forge {
namespace remodel {

/// @brief pure virtual base class
class RemodelConstraintGenerator : public protocols::moves::ConstraintGenerator
{
public: // typedefs
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~RemodelConstraintGenerator();

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & );

	typedef protocols::forge::components::VarLengthBuildAP VarLengthBuildAP;


public: //constructors

	RemodelConstraintGenerator();

public:

	VarLengthBuildAP
	vlb() const;

	void
	set_vlb( VarLengthBuildAP vlb );

	void
	set_seqmap( core::id::SequenceMappingCOP seqmap );

	core::id::SequenceMappingCOP
	seqmap() const;

private:
	core::id::SequenceMappingCOP seqmap_;
	VarLengthBuildAP vlb_;
}; //class RemodelConstraintGenerator


} //namespace remodel
} //namespace forge
} //namespace protocols


#endif
