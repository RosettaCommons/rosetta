// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/WriteSSEMover.hh
/// @brief Writes SSE assignation from DSSP or prediction from PSIPRED as REMARK.
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#ifndef INCLUDED_protocols_simple_filters_ReportFilter_hh
#define INCLUDED_protocols_simple_filters_ReportFilter_hh

//unit headers
#include <protocols/simple_moves/WriteSSEMover.fwd.hh>

// Project Headers
#include <core/io/external/PsiPredInterface.fwd.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

namespace protocols {
namespace simple_moves {

/// @brief Writes SSE assignation from DSSP or prediction from PSIPRED as REMARK.
class WriteSSEMover : public moves::Mover
{
public:
	WriteSSEMover();
	~WriteSSEMover() override;

	void apply( core::pose::Pose & ) override;
	void dssp( bool pick ) { dssp_ = pick; };
	void cmd( std::string const & cmd );

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new WriteSSEMover( *this ));
	}

	protocols::moves::MoverOP fresh_instance() const override{
		return protocols::moves::MoverOP( new WriteSSEMover() );
	}

	void parse_my_tag(utility::tag::TagCOP,basic::datacache::DataMap &,protocols::filters::Filters_map const &f,protocols::moves::Movers_map const &,core::pose::Pose const &) override;
	std::string get_name() const override;

	static
	std::string
	mover_name();

	static
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	/// @brief the object which communicates with psipred and interprets its output
	core::io::external::PsiPredInterfaceOP psipred_interface_;
	/// @brief path to PSIPRED
	std::string cmd_;
	/// @brief use dssp
	bool dssp_;

};
}
}

#endif
