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
#include <core/conformation/Conformation.hh>

//External
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_io.hpp>

// Utility Headers
#include <numeric/HomogeneousTransform.hh>
 #include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>

namespace protocols{
namespace features{

using std::string;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using core::Size;
using core::SSize;
using core::Real;
using core::pose::Pose;
using numeric::HomogeneousTransform;
using numeric::xyzVector;
using utility::sql_database::sessionOP;
using utility::vector1;
using cppdb::statement;
using cppdb::result;

static basic::Tracer TR( "protocols.features.LoopAnchorFeatures" );

LoopAnchorFeatures::LoopAnchorFeatures() :
	FeaturesReporter(),
	use_single_residue_to_define_anchor_transfrom_(true),
	min_loop_length_(5),
	max_loop_length_(30)
{}

LoopAnchorFeatures::LoopAnchorFeatures( LoopAnchorFeatures const & src) :
	FeaturesReporter(),
	use_single_residue_to_define_anchor_transfrom_(src.use_single_residue_to_define_anchor_transfrom_),
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
			"	struct_id BLOB,\n"
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
			"	struct_id BLOB,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"	alpha REAL,\n"
			"	omega REAL,\n"
			"   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
			"       REFERENCES loop_anchors (struct_id, residue_begin, residue_end)\n"
			"       DEFERRABLE INITIALLY DEFERRED,\n"
			"   PRIMARY KEY(struct_id, residue_begin, residue_end));"
			"\n"
			"CREATE TABLE IF NOT EXISTS loop_anchor_transforms_three_res (\n"
			"	struct_id BLOB,\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"	alpha REAL,\n"
			"	omega REAL,\n"
			"   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
			"       REFERENCES loop_anchors (struct_id, residue_begin, residue_end)\n"
			"       DEFERRABLE INITIALLY DEFERRED,\n"
			"   PRIMARY KEY(struct_id, residue_begin, residue_end));";
	}else if(db_mode == "mysql")
	{
		return
			"CREATE TABLE IF NOT EXISTS loop_anchors (\n"
			"	struct_id BINARY(36),\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
            "   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
            "       REFERENCES residues (struct_id, resNum),\n"
			"	PRIMARY KEY(struct_id, residue_begin, residue_end));"
			"\n"
			"CREATE TABLE IF NOT EXISTS loop_anchor_transforms (\n"
			"	struct_id BINARY(36),\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"	alpha REAL,\n"
			"	omega REAL,\n"
			"   FOREIGN KEY (struct_id, residue_begin, residue_end)\n"
			"       REFERENCES loop_anchors (struct_id, residue_begin, residue_end),\n"
			"   PRIMARY KEY(struct_id, residue_begin, residue_end));"
			"\n"
			"CREATE TABLE IF NOT EXISTS loop_anchor_transforms_three_res (\n"
			"	struct_id BINARY(36),\n"
			"	residue_begin INTEGER,\n"
			"	residue_end INTEGER,\n"
			"	x REAL,\n"
			"	y REAL,\n"
			"	z REAL,\n"
			"	phi REAL,\n"
			"	psi REAL,\n"
			"	theta REAL,\n"
			"	alpha REAL,\n"
			"	omega REAL,\n"
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
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	min_loop_length_ = tag->getOption<Size>("min_loop_length", 5);
	max_loop_length_ = tag->getOption<Size>("max_loop_length", 30);

	set_use_relevant_residues_as_loop_length( tag->getOption<bool>("use_relevant_residues_as_loop_length", 0) );

	if(max_loop_length_ < min_loop_length_){
		std::stringstream error_msg;
		error_msg
			<< "The max_loop_length, '" << max_loop_length_ << "',"
		<< " must be longer than the min_loop_length, '" << min_loop_length_ << "'." << std::endl;
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
	boost::uuids::uuid struct_id,
	sessionOP db_session
){
	string loop_anchors_stmt_string = "INSERT INTO loop_anchors VALUES (?,?,?);";
	statement loop_anchors_stmt(
		safely_prepare_statement(loop_anchors_stmt_string, db_session));

	string loop_anchor_transforms_stmt_string =
		"INSERT INTO loop_anchor_transforms VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	statement loop_anchor_transforms_stmt(
		safely_prepare_statement(loop_anchor_transforms_stmt_string, db_session));

	string loop_anchor_transforms_three_res_stmt_string =
		"INSERT INTO loop_anchor_transforms_three_res VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	statement loop_anchor_transforms_three_res_stmt(
		safely_prepare_statement(loop_anchor_transforms_three_res_stmt_string, db_session));
	
	vector1<Size>::const_iterator chain_ending(pose.conformation().chain_endings().begin());
	vector1<Size>::const_iterator chain_ending_end(pose.conformation().chain_endings().end());
	
	Size local_min_loop_length = min_loop_length( relevant_residues );
	Size local_max_loop_length = max_loop_length( relevant_residues );
	
	for(SSize begin=1; begin <= SSize( pose.total_residue() - local_min_loop_length ); ++begin){

		for(Size end=begin;
				(end <= begin + local_max_loop_length && end <= pose.total_residue());
				++end){

			bool bail_out = !relevant_residues[end];

			// Note chain_endings does not have the last residue in the pose
			// in the chain endings vector. So handle this separately
			if((chain_ending != chain_ending_end) && (end == *chain_ending))
			{
				++chain_ending;
				bail_out = true;
			}
			if ( bail_out )
			{
				if( (end - begin + 1) < local_min_loop_length){
					begin = end;
				}
				break;				
			}

			if( (end - begin + 1) >= local_min_loop_length){
				loop_anchors_stmt.bind(1,struct_id);
				loop_anchors_stmt.bind(2,begin);
				loop_anchors_stmt.bind(3,end);
				basic::database::safely_write_to_database(loop_anchors_stmt);

				set_use_single_residue_to_define_anchor_transfrom(true);
				compute_transform_and_write_to_db(struct_id, begin, end, pose, loop_anchor_transforms_stmt);

				set_use_single_residue_to_define_anchor_transfrom(false);
				compute_transform_and_write_to_db(struct_id, begin, end, pose, loop_anchor_transforms_three_res_stmt);
				
			}
		}
	}
	return 0;
}

void LoopAnchorFeatures::set_use_relevant_residues_as_loop_length( bool const use_relevant_residues_as_loop_length )
{
    use_relevant_residues_as_loop_length_ = use_relevant_residues_as_loop_length;
}

void LoopAnchorFeatures::set_use_single_residue_to_define_anchor_transfrom( bool const use_single_residue_to_define_anchor_transfrom )
{
    use_single_residue_to_define_anchor_transfrom_ = use_single_residue_to_define_anchor_transfrom;
}

Size LoopAnchorFeatures::min_loop_length( vector1< bool > const & relevant_residues )
{
	return determine_correct_length( relevant_residues, min_loop_length_ );
}

Size LoopAnchorFeatures::max_loop_length( vector1< bool > const & relevant_residues )
{
	return determine_correct_length( relevant_residues, max_loop_length_ );
}

Size LoopAnchorFeatures::determine_correct_length( vector1< bool > const & relevant_residue, Size default_length )
{
	if ( use_relevant_residues_as_loop_length_ )
	{
		Size number_of_residues = 0;
		for ( vector1< bool >::const_iterator it = relevant_residue.begin(); it != relevant_residue.end(); ++it )
		{
			if ( *it ) ++number_of_residues;
		}
		return number_of_residues;
	}
	return default_length;
}

vector1<Size>
LoopAnchorFeatures::start_residue(Size resNo){
	vector1<Size> residue_vector;
	if (use_single_residue_to_define_anchor_transfrom_){
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo);
	} else{
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo + 1);
		residue_vector.push_back(resNo + 2);
	}
	return residue_vector;
}

vector1<Size>
LoopAnchorFeatures::end_residue(Size resNo){
	vector1<Size> residue_vector;
	if (use_single_residue_to_define_anchor_transfrom_){
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo);
	} else{
		residue_vector.push_back(resNo);
		residue_vector.push_back(resNo - 1);
		residue_vector.push_back(resNo - 2);
	}
	return residue_vector;
}

vector1<Size>
LoopAnchorFeatures::atoms(){
	vector1<Size> atom_vector;
	// 1, 2, 3 corresponds to N, CA, C, respectively
	if (use_single_residue_to_define_anchor_transfrom_) {
		atom_vector.push_back(1);
		atom_vector.push_back(2);
		atom_vector.push_back(3);
	} else{
		// Multiple residues being used to define transform - only use CA coordinates
		atom_vector.push_back(2);
		atom_vector.push_back(2);
		atom_vector.push_back(2);
	}
	return atom_vector;
}


HomogeneousTransform<Real>
LoopAnchorFeatures::compute_anchor_transform(
	Pose const & pose,
	vector1<Size> const & residue_begin,
	vector1<Size> const & residue_end,
	vector1<Size> const & atoms){
	
	runtime_assert_string_msg( (residue_begin.size() == atoms.size()) && (residue_begin.size() == residue_end.size()) && 
		(residue_begin.size() == 3), "The length of the residues and atoms vectors must be exactly 3.");
	
	HomogeneousTransform<Real> take_off_frame(
		pose.residue(residue_begin[1]).xyz(atoms[1]),
		pose.residue(residue_begin[2]).xyz(atoms[2]),
		pose.residue(residue_begin[3]).xyz(atoms[3]));

	HomogeneousTransform<Real> landing_frame(
		pose.residue(residue_end[1]).xyz(atoms[1]),
		pose.residue(residue_end[2]).xyz(atoms[2]),
		pose.residue(residue_end[3]).xyz(atoms[3]));

	HomogeneousTransform<Real> anchor_transform(
		take_off_frame.inverse() * landing_frame);

	return anchor_transform;
}

void
LoopAnchorFeatures::compute_transform_and_write_to_db(
	boost::uuids::uuid struct_id,
	Size begin,
	Size end,
	Pose const & pose,
	statement & stmt){

	vector1<Size> start_res = start_residue(begin);
	vector1<Size> end_res = end_residue(end);
	vector1<Size> atom_set = atoms();

	HomogeneousTransform<Real> anchor_transform(
		compute_anchor_transform(pose, start_res, end_res, atom_set));

	xyzVector<Real> t(anchor_transform.point());
	xyzVector<Real> r(anchor_transform.euler_angles_rad());

	Real alpha = compute_atom_angles(pose, start_res, atom_set);
	Real omega = compute_atom_angles(pose, end_res, atom_set);

	stmt.bind(1, struct_id);
	stmt.bind(2, begin);
	stmt.bind(3, end);
	stmt.bind(4, t.x());
	stmt.bind(5, t.y());
	stmt.bind(6, t.z());
	stmt.bind(7, r.x());
	stmt.bind(8, r.y());
	stmt.bind(9, r.z());
	stmt.bind(10, alpha);
	stmt.bind(11, omega);
	basic::database::safely_write_to_database(stmt);
}

Real
LoopAnchorFeatures::compute_atom_angles(
	Pose const & pose,
	vector1<Size> const & residues,
	vector1<Size> const & atoms){
	
	runtime_assert_string_msg( residues.size() == atoms.size() && residues.size() == 3, 
		"The length of the residues and atoms vectors must be exactly 3.");
	
	xyzVector<Real> first_vector(
		pose.residue(residues[1]).xyz(atoms[1]) - pose.residue(residues[2]).xyz(atoms[2]));
	xyzVector<Real> second_vector(
		pose.residue(residues[3]).xyz(atoms[3])  - pose.residue(residues[2]).xyz(atoms[2]));
	
	return numeric::arccos(first_vector.dot(second_vector)/(first_vector.norm() * second_vector.norm()));
}

} //namesapce
} //namespace
