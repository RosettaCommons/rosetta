// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/LoopAnchorFeatures.cc
/// @brief  report loop anchor features to a features database
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// @author Brian Weitzner

// Unit Headers
#include <protocols/features/LoopAnchorFeatures.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

//External

// Utility Headers
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyz.functions.hh>
#include <utility/exit.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <utility/excn/Exceptions.hh>
#include <sstream>

namespace protocols {
namespace features {

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

static THREAD_LOCAL basic::Tracer TR( "protocols.features.LoopAnchorFeatures" );

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

void
LoopAnchorFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_loop_anchors_table_schema(db_session);
	write_loop_anchor_transforms_table_schema(db_session);
}

void
LoopAnchorFeatures::write_loop_anchors_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column residue_begin("residue_begin", DbDataTypeOP( new DbInteger() ));
	Column residue_end("residue_end", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_begin);
	primary_key_columns.push_back(residue_end);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	foreign_key_columns1.push_back(residue_begin);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	reference_columns1.push_back("resNum");
	ForeignKey foreign_key1(foreign_key_columns1, "residues", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(residue_end);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("resNum");
	ForeignKey foreign_key2(foreign_key_columns2, "residues", reference_columns2, true);

	Schema table("loop_anchors", primary_key);
	table.add_foreign_key(foreign_key1);
	table.add_foreign_key(foreign_key2);

	table.write(db_session);
}

void
LoopAnchorFeatures::write_loop_anchor_transforms_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column residue_begin("residue_begin", DbDataTypeOP( new DbInteger() ));
	Column residue_end("residue_end", DbDataTypeOP( new DbInteger() ));
	Column x("x", DbDataTypeOP( new DbReal() ));
	Column y("y", DbDataTypeOP( new DbReal() ));
	Column z("z", DbDataTypeOP( new DbReal() ));
	Column phi("phi", DbDataTypeOP( new DbReal() ));
	Column psi("psi", DbDataTypeOP( new DbReal() ));
	Column theta("theta", DbDataTypeOP( new DbReal() ));
	Column alpha101("alpha101", DbDataTypeOP( new DbReal() ));
	Column tau101("tau101", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(residue_begin);
	primary_key_columns.push_back(residue_end);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	foreign_key_columns.push_back(residue_begin);
	foreign_key_columns.push_back(residue_end);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	reference_columns.push_back("residue_begin");
	reference_columns.push_back("residue_end");
	ForeignKey foreign_key(foreign_key_columns, "loop_anchors", reference_columns, true);

	Schema table("loop_anchor_transforms", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(x);
	table.add_column(y);
	table.add_column(z);
	table.add_column(phi);
	table.add_column(psi);
	table.add_column(theta);
	table.add_column(alpha101);
	table.add_column(tau101);

	table.write(db_session);
}

utility::vector1<std::string>
LoopAnchorFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	return dependencies;
}

void
LoopAnchorFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	min_loop_length_ = tag->getOption<Size>("min_loop_length", 5);
	max_loop_length_ = tag->getOption<Size>("max_loop_length", 30);

	set_use_relevant_residues_as_loop_length( tag->getOption<bool>("use_relevant_residues_as_loop_length", 0) );

	if ( max_loop_length_ < min_loop_length_ ) {
		std::stringstream error_msg;
		error_msg
			<< "The max_loop_length, '" << max_loop_length_ << "',"
			<< " must be longer than the min_loop_length, '" << min_loop_length_ << "'." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption(error_msg.str());
	}
}

/// @details
/// An anchor is a take off and landing for a loop.
/// Every residue in the loop must be relevant in order for the loop to be stored.
Size
LoopAnchorFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID struct_id,
	sessionOP db_session
){
	string loop_anchors_stmt_string = "INSERT INTO loop_anchors (struct_id, residue_begin, residue_end) VALUES (?,?,?);";
	statement loop_anchors_stmt(
		safely_prepare_statement(loop_anchors_stmt_string, db_session));

	string loop_anchor_transforms_stmt_string =
		"INSERT INTO loop_anchor_transforms (struct_id, residue_begin, residue_end, x, y, z, phi, psi, theta, alpha101, tau101) VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	statement loop_anchor_transforms_stmt(
		safely_prepare_statement(loop_anchor_transforms_stmt_string, db_session));

	vector1<Size>::const_iterator chain_ending(pose.conformation().chain_endings().begin());
	vector1<Size>::const_iterator chain_ending_end(pose.conformation().chain_endings().end());

	Size local_min_loop_length = min_loop_length( relevant_residues );
	Size local_max_loop_length = max_loop_length( relevant_residues );

	for ( SSize begin=1; begin <= SSize( pose.total_residue() - local_min_loop_length ); ++begin ) {
		// check if end+1 is a chain ending or the last residue because we will use end+1 to calculate a pseudodihedral
		for ( Size end=begin;
				(end <= begin + local_max_loop_length && end + 1 <= pose.total_residue());
				++end ) {

			bool bail_out = !relevant_residues[end];

			// Note: chain_endings does not have the last residue in the pose
			// in the chain endings vector, so we handle this separately.
			if ( (chain_ending != chain_ending_end) && (end + 1 == *chain_ending) ) {
				++chain_ending;
				bail_out = true;
			}
			if ( bail_out ) {
				if ( (end - begin + 1) < local_min_loop_length ) {
					begin = end;
				}
				break;
			}

			if ( (end - begin + 1) >= local_min_loop_length ) {
				loop_anchors_stmt.bind(1, struct_id);
				loop_anchors_stmt.bind(2, begin);
				loop_anchors_stmt.bind(3, end);
				basic::database::safely_write_to_database(loop_anchors_stmt);

				set_use_single_residue_to_define_anchor_transfrom(true);
				compute_transform_and_write_to_db(struct_id, begin, end, pose, loop_anchor_transforms_stmt);
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

Size LoopAnchorFeatures::min_loop_length( vector1< bool > const & relevant_residues ) const
{
	return determine_correct_length( relevant_residues, min_loop_length_ );
}

Size LoopAnchorFeatures::max_loop_length( vector1< bool > const & relevant_residues ) const
{
	return determine_correct_length( relevant_residues, max_loop_length_ );
}

Size LoopAnchorFeatures::determine_correct_length( vector1< bool > const & relevant_residue, Size default_length ) const
{
	if ( use_relevant_residues_as_loop_length_ ) {
		Size number_of_residues = 0;
		for ( vector1< bool >::const_iterator it = relevant_residue.begin(); it != relevant_residue.end(); ++it ) {
			if ( *it ) ++number_of_residues;
		}
		return number_of_residues;
	}
	return default_length;
}

HomogeneousTransform<Real> const
LoopAnchorFeatures::frame_for_residue(core::conformation::Residue const & residue) const
{
	return HomogeneousTransform<Real>(residue.xyz("N"), residue.xyz("CA"), residue.xyz("C"));
}

HomogeneousTransform<Real> const
LoopAnchorFeatures::compute_anchor_transform(
	Pose const & pose,
	Size const residue_begin,
	Size const residue_end) const
{
	HomogeneousTransform<Real> take_off_frame = frame_for_residue(pose.residue(residue_begin));
	HomogeneousTransform<Real> landing_frame = frame_for_residue(pose.residue(residue_end));

	HomogeneousTransform<Real> anchor_transform(
		take_off_frame.inverse() * landing_frame);

	return anchor_transform;
}

void
LoopAnchorFeatures::compute_transform_and_write_to_db(
	StructureID struct_id,
	Size begin,
	Size end,
	Pose const & pose,
	statement & stmt) const {

	HomogeneousTransform<Real> anchor_transform(
		compute_anchor_transform(pose, begin, end));

	xyzVector<Real> t(anchor_transform.point());
	xyzVector<Real> r(anchor_transform.euler_angles_rad());

	// tau and alpha definitions from phipsi program (http://bioserv.rpbs.jussieu.fr/Yakusa/download/README_phipsi)

	// compute tau101 -- tauN defined as CA(N-1)--CA(N)--CA(N+1) pseudo-bond angle
	Real tau101 = numeric::angle_degrees(
		pose.residue(end - 2).xyz("CA"),
		pose.residue(end - 1).xyz("CA"),
		pose.residue(end).xyz("CA"));

	// compute alpha101 -- alphaN defined as CA(N-1)--CA(N)--CA(N+1)--CA(N+2) pseudo-dihedral angle
	Real alpha101 = numeric::dihedral_degrees(
		pose.residue(end - 2).xyz("CA"),
		pose.residue(end - 1).xyz("CA"),
		pose.residue(end).xyz("CA"),
		pose.residue(end + 1).xyz("CA"));

	stmt.bind(1, struct_id);
	stmt.bind(2, begin);
	stmt.bind(3, end);
	stmt.bind(4, t.x());
	stmt.bind(5, t.y());
	stmt.bind(6, t.z());
	stmt.bind(7, r.x());
	stmt.bind(8, r.y());
	stmt.bind(9, r.z());
	stmt.bind(10, alpha101);
	stmt.bind(11, tau101);
	basic::database::safely_write_to_database(stmt);
}

} //namesapce
} //namespace
