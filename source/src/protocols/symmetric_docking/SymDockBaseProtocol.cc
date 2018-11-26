// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/symmetric_docking/SymDataBaseProtocol.cc
///
/// @brief
/// @author Ingemar Andre


#include <protocols/symmetric_docking/SymDockBaseProtocol.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <basic/Tracer.hh>


#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>


namespace protocols {
namespace symmetric_docking {

static basic::Tracer TR( "protocols.symmetric_docking.SymDockBaseProtocol" );

SymDockBaseProtocol::SymDockBaseProtocol() :
	Mover()
{
	Mover::type( "SymDockBaseProtocol" );

	using namespace basic::options;
	using namespace core::scoring;
	scorefxn_lowres_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
	scorefxn_hires_ =  ScoreFunctionFactory::create_score_function( core::scoring::PRE_TALARIS_2013_STANDARD_WTS, core::scoring::DOCK_PATCH );

}

SymDockBaseProtocol::~SymDockBaseProtocol() = default;

std::string
SymDockBaseProtocol::get_name() const {
	return "SymDockBaseProtocol";
}

} // symmetric_docking
} // protocols
