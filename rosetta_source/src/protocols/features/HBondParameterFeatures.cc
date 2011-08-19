// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/HBondParameterFeatures.cc
/// @brief  report HBond Parameters to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/HBondParameterFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

namespace protocols{
namespace features{

using std::string;
using std::endl;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::chemical::BondName;
using core::chemical::AtomIndices;
using core::chemical::ResidueType;
using core::scoring::ScoreFunctionOP;
using core::scoring::getScoreFunction;
using core::scoring::hbonds::HBondDatabase;
using core::scoring::hbonds::HBondDatabaseCOP;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;
using basic::Tracer;

static Tracer TR("protocols.features.HBondParameterFeatures");

HBondParameterFeatures::HBondParameterFeatures() :
	FeaturesReporter(),
	scfxn_(getScoreFunction())
{}

HBondParameterFeatures::HBondParameterFeatures(
	ScoreFunctionOP scfxn) :
	FeaturesReporter(),
	scfxn_(getScoreFunction())
{}

HBondParameterFeatures::HBondParameterFeatures(
	HBondParameterFeatures const & src) :
	FeaturesReporter(),
	scfxn_(src.scfxn_)
{}

HBondParameterFeatures::~HBondParameterFeatures() {}

string
HBondParameterFeatures::type_name() const { return "HBondParameterFeatures"; }

string
HBondParameterFeatures::schema() const {
	HBondDatabaseCOP hb_database(HBondDatabase::get_database());
	return hb_database->report_parameter_features_schema();
}

Size
HBondParameterFeatures::report_features(
	Pose const &,
	vector1< bool > const &,
	Size const,
	sessionOP db_session
){

	string const & database_tag(
		scfxn_->energy_method_options().hbond_options().params_database_tag());
	HBondDatabaseCOP hb_database(HBondDatabase::get_database(database_tag));
	return hb_database->report_parameter_features(db_session);
}

} // namesapce
} // namespace
