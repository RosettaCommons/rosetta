// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/AtomAtomPairFeatures.cc
/// @brief  report atom-atom pair geometry and scores to features statistics scientific benchmark
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/features/AtomAtomPairFeatures.hh>

// Project Headers
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <utility/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/types.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>

#include <basic/datacache/DataMap.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray3D.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// External Headers
#include <cppdb/frontend.h>

// C++ Headers
#include <algorithm>
#include <utility/excn/Exceptions.hh>
#include <map>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/AtomAtomPairFeaturesCreator.hh>

namespace protocols {
namespace features {

using std::map;
using std::string;
using std::endl;
using std::upper_bound;
using core::chemical::num_canonical_aas;
using core::chemical::AtomTypeSetCOP;
using core::chemical::AtomIndices;
using core::chemical::ChemicalManager;
using core::pose::Pose;
using core::Size;
using core::Real;
using core::Distance;
using core::Vector;
using utility::graph::Graph;
using core::conformation::Residue;
using core::scoring::TenANeighborGraph;
using protocols::filters::Filters_map;
using basic::datacache::DataMap;
using protocols::moves::Movers_map;
using utility::tag::TagCOP;
using utility::sql_database::sessionOP;
using utility::vector1;
using basic::Tracer;
using basic::database::safely_write_to_database;
using basic::database::safely_prepare_statement;
using ObjexxFCL::FArray3D;
using cppdb::statement;

static Tracer TR("protocols.features.AtomAtomPairFeatures");

AtomAtomPairFeatures::AtomAtomPairFeatures() :
	min_dist_(0.0),
	max_dist_(10.0),
	nbins_(15)
{
	relevant_atom_names_.push_back("CAbb");
	relevant_atom_names_.push_back("CObb");
	relevant_atom_names_.push_back("OCbb");
	relevant_atom_names_.push_back("CNH2");
	relevant_atom_names_.push_back("COO" );
	relevant_atom_names_.push_back("CH1" );
	relevant_atom_names_.push_back("CH2" );
	relevant_atom_names_.push_back("CH3" );
	relevant_atom_names_.push_back("aroC");
	relevant_atom_names_.push_back("Nbb" );
	relevant_atom_names_.push_back("Ntrp");
	relevant_atom_names_.push_back("Nhis");
	relevant_atom_names_.push_back("NH2O");
	relevant_atom_names_.push_back("Nlys");
	relevant_atom_names_.push_back("Narg");
	relevant_atom_names_.push_back("Npro");
	relevant_atom_names_.push_back("OH"  );
	relevant_atom_names_.push_back("ONH2");
	relevant_atom_names_.push_back("OOC" );
	relevant_atom_names_.push_back("Oaro");
	relevant_atom_names_.push_back("Hpol");
	relevant_atom_names_.push_back("Hapo");
	relevant_atom_names_.push_back("Haro");
	relevant_atom_names_.push_back("HNbb");
	relevant_atom_names_.push_back("HOH" );
	relevant_atom_names_.push_back("S"   );

	AtomTypeSetCOP atom_type_set(ChemicalManager::get_instance()->atom_type_set("fa_standard"));

	for ( Size i=1; i <= relevant_atom_names_.size(); ++i ) {
		atom_index_to_relevant_atom_index_[
			atom_type_set->atom_type_index(relevant_atom_names_[i])] = i;
	}

	relevant_elements_["C"] = 1;
	relevant_elements_["N"] = 2;
	relevant_elements_["O"] = 3;
	relevant_elements_["H"] = 4;

}

AtomAtomPairFeatures::AtomAtomPairFeatures(AtomAtomPairFeatures const & ) = default;

AtomAtomPairFeatures::~AtomAtomPairFeatures()= default;

// XRW TEMP string
// XRW TEMP AtomAtomPairFeatures::type_name() const { return "AtomAtomPairFeatures"; }

void
AtomAtomPairFeatures::write_schema_to_db(
	sessionOP db_session
) const {
	write_atom_pairs_table_schema(db_session);
}

void
AtomAtomPairFeatures::write_atom_pairs_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP( new DbBigInt() ));
	Column atom_type("atom_type", DbDataTypeOP( new DbText() ));
	Column element("element", DbDataTypeOP( new DbText() ));
	Column lower_break("lower_break", DbDataTypeOP( new DbReal() ));
	Column upper_break("upper_break", DbDataTypeOP( new DbReal() ));
	Column count("count", DbDataTypeOP( new DbInteger() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(struct_id);
	primary_key_columns.push_back(atom_type);
	primary_key_columns.push_back(element);
	primary_key_columns.push_back(lower_break);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns;
	foreign_key_columns.push_back(struct_id);
	vector1< std::string > reference_columns;
	reference_columns.push_back("struct_id");
	ForeignKey foreign_key(foreign_key_columns, "structures", reference_columns, true);

	Schema table("atom_pairs", primary_key);
	table.add_foreign_key(foreign_key);
	table.add_column(upper_break);
	table.add_column(count);

	table.write(db_session);
}

utility::vector1<std::string>
AtomAtomPairFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

void
AtomAtomPairFeatures::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	min_dist_ = tag->getOption<Real>("min_dist", 0.0);
	max_dist_ = tag->getOption<Real>("max_dist", 10.0);
	nbins_ = tag->getOption<Size>("nbins", 15);
	if ( nbins_ < 1 ) {
		throw utility::excn::EXCN_RosettaScriptsOption("The parameter 'nbins' must be an integer greater than 0.");
	}
}


Size
AtomAtomPairFeatures::report_features(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){
	report_atom_pairs(pose, relevant_residues, struct_id, db_session);
	return 0;
}


void
AtomAtomPairFeatures::report_atom_pairs(
	Pose const & pose,
	vector1< bool > const & relevant_residues,
	StructureID const struct_id,
	sessionOP db_session
){

	// assert pose.update_residue_neighbors() has been called:
	runtime_assert(
		!pose.conformation().structure_moved() &&
		pose.energies().residue_neighbors_updated());

	if ( pose.size() ==0 ) {
		return;
	}

	if ( pose.residue(1).type().atom_type_set().name() != "fa_standard" ) {
		TR.Warning
			<< "Currently AtomAtomPairFeatures only works "
			<< "for the 'fa_standard' AtomTypeSet. This pose has AtomTypeSet '"
			<< pose.residue(1).type().atom_type_set().name() << "'.";
		utility_exit();
	}


	vector1<Distance> bin_breaks;
	Distance const bin_width((max_dist_-min_dist_)/nbins_);
	for ( Size i=0; i <= nbins_; ++i ) {
		bin_breaks.push_back(i*bin_width + min_dist_);
	}

	int const dim1(relevant_atom_names_.size());
	int const dim2(relevant_elements_.size());
	int const dim3(nbins_);
	Size const initial_value(0);
	FArray3D< Size > counts;
	counts.dimension(dim1, dim2, dim3, initial_value);

	for ( Size res_num1=1; res_num1 <= pose.size(); ++res_num1 ) {
		Residue res1(pose.residue(res_num1));

		for ( Size atom_num1=1; atom_num1 <= res1.natoms(); ++atom_num1 ) {
			Vector const & atom1_xyz( res1.xyz(atom_num1) );

			Size const atom_index1(res1.type().atom(atom_num1).atom_type_index());
			map<Size, Size>::const_iterator const i_relevant_atom_index1(
				atom_index_to_relevant_atom_index_.find(atom_index1));
			if ( i_relevant_atom_index1 == atom_index_to_relevant_atom_index_.end() ) {
				continue;
			}

			for ( Size res_num2=1; res_num2 <= pose.size(); ++res_num2 ) {
				if ( !check_relevant_residues(
						relevant_residues, res_num1, res_num2) ) continue;
				Residue res2( pose.residue(res_num2) );

				for ( Size atom_num2=1; atom_num2 <= res2.natoms(); ++atom_num2 ) {
					string const elem_name2(res2.type().atom_type(atom_num2).element());
					map< string, Size>::const_iterator i_elem2(
						relevant_elements_.find(elem_name2));
					if ( i_elem2 == relevant_elements_.end() ) continue;

					Vector const & atom2_xyz( res2.xyz(atom_num2) );
					Distance dist(atom1_xyz.distance(atom2_xyz));
					if ( dist <= min_dist_ || dist > max_dist_ ) continue;


					Size const dist_bin(
						static_cast<Size>(ceil(
						(dist-min_dist_)*nbins_/(max_dist_-min_dist_))));
					counts(i_relevant_atom_index1->second, i_elem2->second, dist_bin) += 1;
				}
			}
		}
	}

	string stmt_string = "INSERT INTO atom_pairs (struct_id, atom_type, element, lower_break, upper_break, count) VALUES (?,?,?,?,?,?);";
	statement stmt(safely_prepare_statement(stmt_string,db_session));


	for ( Size i_atom1=1; i_atom1 <= relevant_atom_names_.size(); ++i_atom1 ) {
		for ( map<string, Size>::const_iterator
				i_elem2=relevant_elements_.begin(),
				ie_elem2=relevant_elements_.end(); i_elem2 != ie_elem2; ++i_elem2 ) {
			for ( Size dist_bin=1; dist_bin <= nbins_; ++dist_bin ) {
				Real const lower_break(min_dist_ + (dist_bin - 1)*bin_width);
				Real const upper_break(min_dist_ + dist_bin * bin_width);
				Size const count(counts(i_atom1,i_elem2->second, dist_bin));
				stmt.bind(1,struct_id);
				stmt.bind(2, relevant_atom_names_[i_atom1]);
				stmt.bind(3, i_elem2->first);
				stmt.bind(4, lower_break);
				stmt.bind(5, upper_break);
				stmt.bind(6, count);
				safely_write_to_database(stmt);
			}
		}
	}
}

std::string AtomAtomPairFeatures::type_name() const {
	return class_name();
}

std::string AtomAtomPairFeatures::class_name() {
	return "AtomAtomPairFeatures";
}

void AtomAtomPairFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "min_dist", xsct_real, "Minimum distance of interest", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default( "max_dist", xsct_real, "Maximum distance of interest", "10.0" )
		+ XMLSchemaAttribute::attribute_w_default( "nbins", xsct_positive_integer, "Number of bins to subdivide the above interval", "15" );

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Atom pair features for any two atoms in a pose", attlist );
}

std::string AtomAtomPairFeaturesCreator::type_name() const {
	return AtomAtomPairFeatures::class_name();
}

protocols::features::FeaturesReporterOP
AtomAtomPairFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new AtomAtomPairFeatures );
}

void AtomAtomPairFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AtomAtomPairFeatures::provide_xml_schema( xsd );
}


} // namesapce
} // namespace
