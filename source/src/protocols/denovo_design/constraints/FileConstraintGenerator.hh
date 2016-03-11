// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/denovo_design/constraints/FileConstraintGenerator.hh
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky ( tlinsky@uw.edu ), Nov 2012

#ifndef INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_hh
#define INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_hh

// Unit Header
#include <protocols/denovo_design/constraints/FileConstraintGenerator.fwd.hh>

// Package Header
#include <protocols/forge/remodel/RemodelConstraintGenerator.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

// Proeject Header
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace denovo_design {
namespace constraints {

class FileConstraintGenerator : public protocols::forge::remodel::RemodelConstraintGenerator{
public:

	typedef std::string String;
	typedef core::pose::Pose Pose;

public:
	FileConstraintGenerator();

	FileConstraintGenerator( String const & filename );

	virtual ~FileConstraintGenerator();

	void set_cstfile( String const & filename );

	virtual void
	parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual std::string
	get_name() const;

	virtual protocols::moves::MoverOP
	fresh_instance() const;

	virtual protocols::moves::MoverOP
	clone() const;

	virtual core::scoring::constraints::ConstraintCOPs
	generate_constraints( Pose const & pose );

private:

	String filename_;

}; //class FileConstraintGenerator


} //namespace constraints
} //namespace denovo_design
} //namespace protocols


#endif // INCLUDED_protocols_denovo_design_constraints_FileConstraintGenerator_HH
