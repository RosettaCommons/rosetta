// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/LoopAnchorFeatures.cc
/// @brief  report comments stored with each pose
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/LoopAnchorFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Platform Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <protocols/moves/DataMap.hh>

// Numeric Headers
#include <numeric/HomogeneousTransform.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

// External Headers
#include <cppdb/frontend.h>

// Boost Headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers
#include <string>
#include <map>
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using std::stringstream;
using std::endl;
using std::map;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using core::Size;
using core::SSize;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using protocols::moves::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;
using numeric::HomogeneousTransform;
using numeric::xyzVector;
using utility::sql_database::sessionOP;
using utility::vector1;
using utility::tag::TagPtr;
using cppdb::statement;
using cppdb::result;

static basic::Tracer TR( "protocols.features.LoopAnchorFeatures" );

LoopAnchorFeatures::LoopAnchorFeatures() :
	FeaturesReporter(),
	min_loop_length_(5),
	max_loop_length_(30)
{}

LoopAnchorFeatures::LoopAnchorFeatures( LoopAnchorFeatures const & src) :
	FeaturesReporter(),
	min_loop_length_(src.min_loop_length_),
	max_loop_length_(src.max_loop_length_)
{}

LoopAnchorFeatures::~LoopAnchorFeatures() {}



string
LoopAnchorFeatures::type_name() const { return "LoopAnchorFeatures"; }

string
LoopAnchorFeatures::schema() const {
	std::string db_mode(basic::options::option[basic::options::OptionKeys::inout::database_mode]);

	if(db_mode == "sqlite3")
	{
		return
			"CREATE TABLE IF NOT EXISTS loop_anchors (\n"
			"	struct_id INTEGER,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	FOREIGN KEY (struct_id, residue_begin)\n"
			"		REFERENCES residues (struct_id, resNum)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	FOREIGN KEY (struct_id, residue_end)\n"
			"		REFERENCES residues (struct_id, resNum)\n"
			"		DEFERRABLE INITIALLY DEFERRED,\n"
			"	PRIMARY KEY(struct_id, residue_begin, residue_end));\n"
			"\n"
			"CREATE TABLE IF NOT EXISTS loop_anchor_transforms (\n"
			"	struct_id INTEGER,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
			"       REFERENCES loop_anchors (struct_id, residue_begin, residue_end)\n"
			"       DEFERRABLE INITIALLY DEFERRED,\n"
			"   PRIMARY KEY(struct_id, residue_begin, residue_end));";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS loop_anchors (\n"
			"	struct_id INTEGER,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
            "   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
            "       REFERENCES residues (struct_id, resNum),\n"
			"	PRIMARY KEY(struct_id, residue_begin, residue_end));"
			"\n"
			"CREATE TABLE IF NOT EXISTS loop_anchor_transforms (\n"
			"	struct_id INTEGER,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
			"       REFERENCES loop_anchors (struct_id, residue_begin, residue_end),\n"
			"   PRIMARY KEY(struct_id, residue_begin, residue_end));";
	} else {
		return "";
	}

}

utility::vector1<std::string>
LoopAnchorFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
LoopAnchorFeatures::parse_my_tag(
	TagPtr const tag,
	DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	min_loop_length_ = tag->getOption<Size>("min_loop_length", 5);
	max_loop_length_ = tag->getOption<Size>("max_loop_length", 30);

	if(max_loop_length_ < min_loop_length_){
		stringstream error_msg;
		error_msg
			<< "The max_loop_length, '" << max_loop_length_ << "',"
			<< " must be longer than the min_loop_length, '" << min_loop_length_ << "'." << endl;
		utility_exit_with_message(error_msg.str());
	}
}



/// @details
/// An anchor is a take off and landing for a loop.
/// Every residue in the loop must be relevant in order for the loop to be stored.
Size
LoopAnchorFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	Size struct_id,
	sessionOP db_session
){
	string loop_anchors_stmt_string = "INSERT INTO loop_anchors VALUES (?,?,?);";
	statement loop_anchors_stmt(
		safely_prepare_statement(loop_anchors_stmt_string, db_session));

	string loop_anchor_transforms_stmt_string =
		"INSERT INTO loop_anchor_transforms VALUES (?,?,?,?,?,?,?,?,?);";
	statement loop_anchor_transforms_stmt(
		safely_prepare_statement(loop_anchor_transforms_stmt_string, db_session));

	vector1<Size>::const_iterator chain_ending(pose.conformation().chain_endings().begin());
	vector1<Size>::const_iterator chain_ending_end(pose.conformation().chain_endings().end());
    
	for(SSize begin=1; begin < SSize( pose.total_residue() - min_loop_length_ ); ++begin){

		for(Size end=begin;
				(end <= begin + max_loop_length_ + 1 && end <= pose.total_residue());
				++end){

			// Note chain_endings does not have the last residue in the pose
			// in the chain endings vector. So handle this separately
			if(((chain_ending != chain_ending_end) && end == *chain_ending) ||
				!relevant_residues[end]){
				++chain_ending;
				if(end - begin < min_loop_length_){
					begin = end;
				}
				break;
			}

			if(end - begin > min_loop_length_){
				loop_anchors_stmt.bind(1,struct_id);
				loop_anchors_stmt.bind(2,begin);
				loop_anchors_stmt.bind(3,end);
				basic::database::safely_write_to_database(loop_anchors_stmt);


				HomogeneousTransform<Real> anchor_transform(
					compute_anchor_transform(pose, begin, end));

				xyzVector<Real> t(anchor_transform.point());
				xyzVector<Real> r(anchor_transform.euler_angles_rad());

				loop_anchor_transforms_stmt.bind(1,struct_id);
				loop_anchor_transforms_stmt.bind(2,begin);
				loop_anchor_transforms_stmt.bind(3,end);
				loop_anchor_transforms_stmt.bind(4,t.x());
				loop_anchor_transforms_stmt.bind(5,t.y());
				loop_anchor_transforms_stmt.bind(6,t.z());
				loop_anchor_transforms_stmt.bind(7,r.x());
				loop_anchor_transforms_stmt.bind(8,r.y());
				loop_anchor_transforms_stmt.bind(9,r.z());
				basic::database::safely_write_to_database(loop_anchor_transforms_stmt);
			}
		}
	}
	return 0;
}

HomogeneousTransform<Real>
LoopAnchorFeatures::compute_anchor_transform(
	Pose const & pose,
	Size residue_begin,
	Size residue_end){

  HomogeneousTransform<Real> take_off_frame(
		pose.residue(residue_begin).xyz(1),
		pose.residue(residue_begin).xyz(2),
		pose.residue(residue_begin).xyz(3));

	HomogeneousTransform<Real> landing_frame(
		pose.residue(residue_end).xyz(1),
		pose.residue(residue_end).xyz(2),
		pose.residue(residue_end).xyz(3));

	HomogeneousTransform<Real> anchor_transform(
		take_off_frame.inverse() * landing_frame);

	return anchor_transform;
}

} //namesapce
} //namespace
