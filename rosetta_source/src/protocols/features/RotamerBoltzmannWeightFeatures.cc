// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RotamerBoltzmannWeightFeatures.cc
/// @brief  report residue burial to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/RotamerBoltzmannWeightFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
#include <protocols/protein_interface_design/filters/RotamerBoltzmannWeight.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

namespace protocols{
namespace features{

using std::string;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using core::scoring::ScoreFunctionOP;
using protocols::protein_interface_design::filters::RotamerBoltzmannWeight;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;

RotamerBoltzmannWeightFeatures::RotamerBoltzmannWeightFeatures() :
	rotamer_boltzmann_weight_(new RotamerBoltzmannWeight())
{}

RotamerBoltzmannWeightFeatures::RotamerBoltzmannWeightFeatures(
	ScoreFunctionOP scfxn) :
	rotamer_boltzmann_weight_(new RotamerBoltzmannWeight())
{
	rotamer_boltzmann_weight_->scorefxn(scfxn);
	rotamer_boltzmann_weight_->type("monomer");
}

RotamerBoltzmannWeightFeatures::RotamerBoltzmannWeightFeatures(RotamerBoltzmannWeightFeatures const & src) :
	FeaturesReporter(),
	rotamer_boltzmann_weight_(src.rotamer_boltzmann_weight_)
{}

RotamerBoltzmannWeightFeatures::~RotamerBoltzmannWeightFeatures(){}

string
RotamerBoltzmannWeightFeatures::type_name() const { return "RotamerBoltzmannWeightFeatures"; }

string
RotamerBoltzmannWeightFeatures::schema() const {
	return
		"CREATE TABLE IF NOT EXISTS rotamer_boltzmann_weight (\n"
		"	struct_id INTEGER,\n"
		"	resNum INTEGER,\n"
		"	boltzmann_weight REAL,\n"
		"	FOREIGN KEY (struct_id, resNum)\n"
		"		REFERENCES residues (struct_id, resNum)\n"
		"		DEFERRABLE INITIALLY DEFERRED,\n"
		"	PRIMARY KEY (struct_id, resNum));\n";
}
Size
RotamerBoltzmannWeightFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size const struct_id,
	sessionOP db_session
){

	for(Size resNum=1; resNum <= pose.total_residue(); ++resNum){
		if(!relevant_residues[resNum]) continue;
		Real const boltzmann_weight(
			rotamer_boltzmann_weight_->compute_Boltzmann_weight(pose, resNum));

		statement stmt = (*db_session)
			<< "INSERT INTO rotamer_boltzmann_weight VALUES (?,?,?);"
			<< struct_id
			<< resNum
			<< boltzmann_weight;
		stmt.exec();
	}
	return 0;
}

} // namespace
} // namespace
