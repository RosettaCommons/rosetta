// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SmotifFeatures.cc
///
/// @brief Calculate 4 geometric values used to classify smotifs. Reference (http://www.ncbi.nlm.nih.gov/pubmed/20421995)
/// @author Tim Jacobs

// Unit Headers
#include <protocols/features/SmotifFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>

#ifdef WIN32
#include <core/scoring/ScoreFunction.hh>
#endif

//Basic Headers
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/PCA.hh>
#include <numeric/xyzVector.io.hh>


// Utility Headers
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/SmotifFeaturesCreator.hh>

namespace protocols {
namespace features {

static basic::Tracer TR( "protocols.features.SmotifFeatures" );


SmotifFeatures::SmotifFeatures() = default;

SmotifFeatures::SmotifFeatures(core::scoring::ScoreFunctionOP /*scfxn*/)
{
}

SmotifFeatures::SmotifFeatures( SmotifFeatures const & /*src*/ ) : FeaturesReporter(/*src*/)
{
}

SmotifFeatures::~SmotifFeatures()= default;

/// @brief return string with class name
// XRW TEMP std::string
// XRW TEMP SmotifFeatures::type_name() const
// XRW TEMP {
// XRW TEMP  return "SmotifFeatures";
// XRW TEMP }

/// @brief generate the table schemas and write them to the database
void
SmotifFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;
	using namespace utility;

	//******smotifs******//
	Column smotif_id("smotif_id", DbDataTypeOP(new DbInteger()), true);
	Column struct_id("struct_id", DbDataTypeOP(new DbBigInt()), false);
	Column ss1("secondary_struct_segment_id_1", DbDataTypeOP(new DbInteger()), false);
	Column ss2("secondary_struct_segment_id_2", DbDataTypeOP(new DbInteger()), false);
	Column loop("loop_segment_id", DbDataTypeOP(new DbInteger()), false);

	//Distance, in angstroms, between C-term of SS1 (P1) and N term of SS2 (P2)
	Column distance("distance", DbDataTypeOP(new DbReal()), false);

	//Angle between L (vector between P1 and P2) and M1 (principal moment of inertia for SS1)
	Column hoist("hoist", DbDataTypeOP(new DbReal()), false);

	//Angle between M1 and M2 (principal moment of inertia for SS1 and SS2)
	Column packing("packing", DbDataTypeOP(new DbReal()), false);

	//Angle between M2 and Gamma (defined by M1 and normal to plane Pie, which is defined
	Column meridian("meridian", DbDataTypeOP(new DbReal()), false);

	Columns foreign_key_columns1;
	foreign_key_columns1.push_back(struct_id);
	vector1< std::string > reference_columns1;
	reference_columns1.push_back("struct_id");
	ForeignKey foreign_key1(foreign_key_columns1, "structures", reference_columns1, true);

	Columns foreign_key_columns2;
	foreign_key_columns2.push_back(struct_id);
	foreign_key_columns2.push_back(ss1);
	vector1< std::string > reference_columns2;
	reference_columns2.push_back("struct_id");
	reference_columns2.push_back("segment_id");
	ForeignKey foreign_key2(foreign_key_columns2, "secondary_structure_segments", reference_columns2, true);

	Columns foreign_key_columns3;
	foreign_key_columns3.push_back(struct_id);
	foreign_key_columns3.push_back(ss2);
	vector1< std::string > reference_columns3;
	reference_columns3.push_back("struct_id");
	reference_columns3.push_back("segment_id");
	ForeignKey foreign_key3(foreign_key_columns3, "secondary_structure_segments", reference_columns3, true);

	Columns foreign_key_columns4;
	foreign_key_columns4.push_back(struct_id);
	foreign_key_columns4.push_back(loop);
	vector1< std::string > reference_columns4;
	reference_columns4.push_back("struct_id");
	reference_columns4.push_back("segment_id");
	ForeignKey foreign_key4(foreign_key_columns4, "secondary_structure_segments", reference_columns4, true);

	Schema smotifs("smotifs", smotif_id);
	smotifs.add_foreign_key(foreign_key1);
	smotifs.add_foreign_key(foreign_key2);
	smotifs.add_foreign_key(foreign_key3);
	smotifs.add_foreign_key(foreign_key4);
	smotifs.add_column(distance);
	smotifs.add_column(hoist);
	smotifs.add_column(packing);
	smotifs.add_column(meridian);

	smotifs.write(db_session);
}

/// @brief return the set of features reporters that are required to
///also already be extracted by the time this one is used.
utility::vector1<std::string>
SmotifFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("SecondaryStructureSegmentFeatures");
	return dependencies;
}

void
SmotifFeatures::parse_my_tag(
	utility::tag::TagCOP const /*tag*/,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/)
{
}

void
SmotifFeatures::calculate_angles(
	utility::vector1< numeric::xyzVector<core::Real> > const & ss1_coords,
	utility::vector1< numeric::xyzVector<core::Real> > const & ss2_coords,
	core::Real & distance, /*output*/
	core::Real & hoist_angle_degrees, /*output*/
	core::Real & packing_angle_degrees, /*output*/
	core::Real & meridian_angle_degrees /*output*/
) {
	using namespace numeric;

	xyzVector<core::Real> ss1_com = center_of_mass(ss1_coords);
	xyzVector<core::Real> ss1_first_principal_component =
		first_principal_component(ss1_coords);

	//Draw principal component vector from COM instead of origin
	xyzVector<core::Real> ss1_com_principal_component=
		ss1_first_principal_component+=ss1_com;

	//p0 is point on principal component vector closest to N-term of SS1
	xyzVector<core::Real> p0 =
		closest_point_on_line(ss1_com, ss1_com_principal_component, ss1_coords[1]);

	//p1 is point on principal component vector closest to C-term of SS1
	xyzVector<core::Real> p1 =
		closest_point_on_line(ss1_com, ss1_com_principal_component, ss1_coords.back());

	//SS2
	xyzVector<core::Real> ss2_com = center_of_mass(ss2_coords);
	xyzVector<core::Real> ss2_first_principal_component =
		first_principal_component(ss2_coords);

	xyzVector<core::Real> ss2_com_principal_component=
		ss2_first_principal_component+=ss2_com;

	//p2 is point on principal component vector closest to N-term of SS2
	xyzVector<core::Real> p2 =
		closest_point_on_line(ss2_com, ss2_com_principal_component, ss2_coords[1]);

	//p3 is point on principal component vector closest to C-term of SS2
	xyzVector<core::Real> p3 =
		closest_point_on_line(ss2_com, ss2_com_principal_component, ss2_coords.back());

	distance = p1.distance(p2);
	hoist_angle_degrees = angle_degrees(p0,p1,p2);
	packing_angle_degrees = angle_degrees(p0,p1,p2,p3);

	meridian_angle_degrees = dihedral_degrees(p0,p1,p2,p3);
}

/// @brief collect all the feature data for the pose
core::Size
SmotifFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1< bool > const & /*relevant_residues*/,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{
	using core::Size;
	using core::Real;

	utility::vector1<SecondaryStructureSegment> ss_segments =
		get_ss_segments(struct_id, db_session);

	std::string smotif_insert_string =
		"INSERT INTO smotifs (struct_id, secondary_struct_segment_id_1,\n"
		"\tsecondary_struct_segment_id_2, loop_segment_id, distance, hoist, packing,\n"
		"\tmeridian) VALUES(?,?,?,?,?,?,?,?)";

	cppdb::statement smotif_insert_stmt =
		basic::database::safely_prepare_statement(smotif_insert_string, db_session);

	if ( ss_segments.size()>2 ) {
		for ( Size i=1; i<=ss_segments.size()-2; ++i ) {
			SecondaryStructureSegment ss1=ss_segments[i];
			SecondaryStructureSegment loop=ss_segments[i+1];
			SecondaryStructureSegment ss2=ss_segments[i+2];

			if ( (ss1.dssp=="H" || ss1.dssp=="E") &&
					(ss2.dssp=="H" || ss2.dssp=="E") &&
					//       (ss1.dssp=="H" && ss2.dssp=="H") &&
					loop.dssp == "L" ) { //make sure loop is loop
				utility::vector1< numeric::xyzVector<core::Real> > ss1_coords;
				for ( core::Size i=ss1.residue_begin; i<=ss1.residue_end; ++i ) {
					ss1_coords.push_back(pose.residue(i).atom("CA").xyz());
				}

				utility::vector1< numeric::xyzVector<core::Real> > ss2_coords;
				for ( core::Size i=ss2.residue_begin; i<=ss2.residue_end; ++i ) {
					ss2_coords.push_back(pose.residue(i).atom("CA").xyz());
				}

				Real distance, hoist_angle_degrees, packing_angle_degrees, meridian_angle_degrees;
				calculate_angles(ss1_coords, ss2_coords, distance, hoist_angle_degrees, packing_angle_degrees, meridian_angle_degrees);

				smotif_insert_stmt.bind(1, struct_id);
				smotif_insert_stmt.bind(2, ss1.segment_id);
				smotif_insert_stmt.bind(3, ss2.segment_id);
				smotif_insert_stmt.bind(4, loop.segment_id);
				smotif_insert_stmt.bind(5, distance);
				smotif_insert_stmt.bind(6, hoist_angle_degrees);
				smotif_insert_stmt.bind(7, packing_angle_degrees);
				smotif_insert_stmt.bind(8, meridian_angle_degrees);

				basic::database::safely_write_to_database(smotif_insert_stmt);
			}
		}
	}
	return 0;  // Added to remove warning; should this be returning an ID of some sort? ~Labonte
}

utility::vector1<SecondaryStructureSegment>
SmotifFeatures::get_ss_segments(
	StructureID struct_id,
	utility::sql_database::sessionOP db_session)
{

	using cppdb::statement;
	using cppdb::result;
	using core::Size;

	std::string ss_segment_select_string=
		"SELECT\n"
		"\tsegment_id, residue_begin, residue_end, dssp\n"
		"FROM\n"
		"\tsecondary_structure_segments\n"
		"WHERE\n"
		"\tstruct_id = ?\n"
		"ORDER BY segment_id";

	statement ss_segment_select_statement(
		basic::database::safely_prepare_statement(ss_segment_select_string,db_session));
	ss_segment_select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(ss_segment_select_statement));

	utility::vector1<SecondaryStructureSegment> ss_segments;
	while ( res.next() ) {

		Size segment_id, residue_begin, residue_end;
		std::string dssp;
		res >> segment_id >> residue_begin >> residue_end >> dssp;

		SecondaryStructureSegment ss_segment;
		ss_segment.segment_id = segment_id;
		ss_segment.residue_begin = residue_begin;
		ss_segment.residue_end = residue_end;
		ss_segment.dssp = dssp;

		ss_segments.push_back(ss_segment);
	}
	return ss_segments;

}

std::string SmotifFeatures::type_name() const {
	return class_name();
}

std::string SmotifFeatures::class_name() {
	return "SmotifFeatures";
}

void SmotifFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::features::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Calculate 4 geometric values used to classify smotifs. "
		"Reference (http://www.ncbi.nlm.nih.gov/pubmed/20421995)",
		attlist );
}

std::string SmotifFeaturesCreator::type_name() const {
	return SmotifFeatures::class_name();
}

protocols::features::FeaturesReporterOP
SmotifFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new SmotifFeatures );
}

void SmotifFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SmotifFeatures::provide_xml_schema( xsd );
}


} // namespace
} // namespace
