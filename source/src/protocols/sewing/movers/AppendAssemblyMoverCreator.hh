// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/movers/AppendAssemblyMoverCreator.hh
/// @brief an AssemblyMover for adding to existing poses
/// @author frankdt (frankdt@email.unc.edu)

#ifndef INCLUDED_protocols_sewing_movers_AppendAssemblyMoverCreator_hh
#define INCLUDED_protocols_sewing_movers_AppendAssemblyMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>
#include <protocols/sewing/movers/AssemblyMoverCreator.hh>
namespace protocols {
namespace sewing {
namespace movers {

class AppendAssemblyMoverCreator : public protocols::moves::MoverCreator {

public:

	virtual protocols::moves::MoverOP
	create_mover() const;

	virtual std::string
	keyname() const;

	static std::string class_name();

	static std::string mover_name();

	virtual void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const;
};

} //protocols
} //sewing
} //movers

#endif //INCLUDED_protocols/sewing/movers_AppendAssemblyMover_fwd_hh
