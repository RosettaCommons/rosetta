// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/features/helixAssembly/HelixBundleFeatures.cc
/// @brief Search through a pose for sets of 3 helices and generate RMSD calculations between all pairs of them
/// @author Tim Jacobs

//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

//Devel
#include <protocols/features/helixAssembly/HelixBundleFeatures.hh>
#include <protocols/features/helixAssembly/HelicalFragment.hh>

//Utility
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/tag/Tag.hh>
#include <utility/LexicographicalIterator.hh>

//Core
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <utility/graph/Graph.hh>

//C++
#include <string>
#include <cmath>

//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>
#include <basic/database/schema_generator/DbDataType.hh>

//Numeric
#include <numeric/PCA.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/features/feature_schemas.hh>
#include <protocols/features/helixAssembly/HelixBundleFeaturesCreator.hh>

static basic::Tracer TR( "protocols.features.helixAssembly.HelixBundleFeatures" );

namespace protocols {
namespace features {
namespace helixAssembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

HelixBundleFeatures::HelixBundleFeatures() :
	helix_cap_dist_cutoff_(18.0),
	num_helices_per_bundle_(3),
	min_helix_size_(14)
{
	scorefxn_ = core::scoring::get_score_function();
}

utility::vector1<std::string>
HelixBundleFeatures::features_reporter_dependencies() const
{
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	dependencies.push_back("SecondaryStructureSegmentFeatures");
	return dependencies;
}

void
HelixBundleFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const
{
	using namespace basic::database::schema_generator;

	Column struct_id("struct_id", DbDataTypeOP(new DbBigInt()), false /*not null*/, false /*don't autoincrement*/);

	/******helix_bundles******/
	Column bundle_id("bundle_id", DbDataTypeOP(new DbInteger()), false /*not null*/, true /*autoincrement*/);
	Column num_helices_per_bundle("num_helices_per_bundle", DbDataTypeOP(new DbInteger()), false, false);

	Schema helix_bundles("helix_bundles", bundle_id);
	helix_bundles.add_column(num_helices_per_bundle);
	helix_bundles.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", true /*defer*/));

	helix_bundles.write(db_session);

	/******bundle_helices******/
	Column helix_id("helix_id", DbDataTypeOP(new DbInteger()), false /*not null*/, true /*autoincrement*/);
	Column bundle_id_fkey("bundle_id", DbDataTypeOP(new DbInteger()), false, false);

	Column residue_begin("residue_begin", DbDataTypeOP(new DbInteger()));
	Column residue_end("residue_end", DbDataTypeOP(new DbInteger()));
	Column reversed("reversed", DbDataTypeOP(new DbInteger()));

	utility::vector1<std::string> fkey_res_reference_cols;
	fkey_res_reference_cols.push_back("struct_id");
	fkey_res_reference_cols.push_back("resNum");

	utility::vector1<Column> fkey_cols_begin;
	fkey_cols_begin.push_back(struct_id);
	fkey_cols_begin.push_back(residue_begin);

	utility::vector1<Column> fkey_cols_end;
	fkey_cols_end.push_back(struct_id);
	fkey_cols_end.push_back(residue_end);

	Schema bundle_helices("bundle_helices", helix_id);
	bundle_helices.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", false /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(bundle_id_fkey, "helix_bundles", "bundle_id", false /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_begin, "residues", fkey_res_reference_cols, false /*defer*/));
	bundle_helices.add_foreign_key(ForeignKey(fkey_cols_end, "residues", fkey_res_reference_cols, false /*defer*/));
	bundle_helices.add_column(reversed);

	bundle_helices.write(db_session);

	/******helix_pairs******/

	Column pair_id("pair_id", DbDataTypeOP(new DbInteger()), false /*not null*/, true /*autoincrement*/);
	Column helix_id_1("helix_id_1", DbDataTypeOP(new DbInteger()), false, false);
	Column helix_id_2("helix_id_2", DbDataTypeOP(new DbInteger()), false, false);

	Schema helix_pairs("helix_pairs", pair_id);
	helix_pairs.add_foreign_key(ForeignKey(struct_id, "structures", "struct_id", false /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(bundle_id_fkey, "helix_bundles", "bundle_id", false /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(helix_id_1, "bundle_helices", "helix_id", false /*defer*/));
	helix_pairs.add_foreign_key(ForeignKey(helix_id_2, "bundle_helices", "helix_id", false /*defer*/));

	helix_pairs.write(db_session);
}

//Select all helical segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<HelicalFragmentOP> HelixBundleFeatures::get_helices(StructureID struct_id, sessionOP db_session) {
	std::string select_string =
		"SELECT\n"
		"\tresidue_begin,\n"
		"\tresidue_end\n"
		"FROM\n"
		"\tsecondary_structure_segments\n"
		"WHERE\n"
		"\tstruct_id = ?\n AND "
		"\tdssp = 'H'\n"
		"ORDER BY residue_begin;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<HelicalFragmentOP> all_helices;

	core::Size prev_residue_begin = 0;
	core::Size prev_residue_end = 0;
	while ( res.next() )
			{
		Size residue_begin, residue_end;
		res >> residue_begin >> residue_end;

		if ( residue_begin - prev_residue_end == 2 && all_helices.size() > 0 ) {
			all_helices.pop_back();
			all_helices.push_back(HelicalFragmentOP(new HelicalFragment(prev_residue_begin, residue_end)));
			TR  << "combining helix segments: " << prev_residue_begin << "-" << prev_residue_end
				<< " and " << residue_begin << "-" << residue_end << std::endl;
		} else {
			all_helices.push_back(HelicalFragmentOP(new HelicalFragment(residue_begin, residue_end)));
		}
		prev_residue_begin=residue_begin;
		prev_residue_end=residue_end;
	}

	utility::vector1<HelicalFragmentOP> all_valid_helices;
	for ( core::Size i=1; i<=all_helices.size(); ++i ) {
		if ( all_helices[i]->size() >= min_helix_size_ ) {
			all_valid_helices.push_back(all_helices[i]);
		}
	}

	TR << "number of uninterrupted helix_segments with at least " << min_helix_size_ << " residues: " << all_valid_helices.size() << std::endl;

	return all_valid_helices;
}

/// @brief calculate the shared fa_attr for each pair of helices
/// in the bundle
//void
//HelixBundleFeatures::calc_fa_energy(
// core::pose::Pose const & pose,
// FragmentPair & helix_pair
//){
// HelicalFragmentOP frag_1=helix_pair.fragment_1;
// HelicalFragmentOP frag_2=helix_pair.fragment_2;
//
// core::Real helix_pair_fa_attr=0;
// std::set<core::Size> interacting_residues; //residues that make favorable inter-fragment contacts (negative fa_attr)
// for(core::Size ii=frag_1->seq_start(); ii<=frag_1->seq_end(); ++ii)
// {
//  //iterate through interaction graph of current residue and add any fa_attr
//  //between this residue and any of the residues from fragment j
//  core::scoring::EnergyGraph const & eg =
//  pose.energies().energy_graph();
//  for (utility::graph::Node::EdgeListConstIter
//    iter = eg.get_node(ii)->const_upper_edge_list_begin(),
//    iter_end = eg.get_node(ii)->const_upper_edge_list_end();
//    iter != iter_end; ++iter )
//  {
//   core::Size jj = (*iter)->get_other_ind(ii);
//   if(jj >= frag_2->seq_start() && jj <= frag_2->seq_end())
//   {
//    const core::scoring::EnergyEdge* ee =
//    dynamic_cast< const core::scoring::EnergyEdge* > (*iter);
//    helix_pair_fa_attr += (*ee)[core::scoring::fa_atr];
//    if((*ee)[core::scoring::fa_atr] < 0)
//    {
//     interacting_residues.insert(ii);
//     interacting_residues.insert(jj);
//    }
//   }
//  }
// }
//
// //normalize fragment pair fa attr by number of residues and check for per-residue threshold
// helix_pair.fa_fraction = interacting_residues.size()/(min_helix_size_*2.0);
// helix_pair.fa_attr = helix_pair_fa_attr/(min_helix_size_*2.0);
//
//// TR << "residue-normalized FA attractive between: " << (*frag_1)
////  << " " << (*frag_2) << ": " << per_residue_fa_attr << std::endl;
////
//// TR << "interacting set fraction: " << (*frag_1)
////  << " " << (*frag_2) << ": " << interacting_set_fraction << std::endl;
//}

void
HelixBundleFeatures::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	protocols::filters::Filters_map const & /*filters*/,
	protocols::moves::Movers_map const & /*movers*/,
	core::pose::Pose const & /*pose*/
){
	runtime_assert(tag->getName() == type_name());

	if ( tag->hasOption("num_helices_per_bundle") ) {
		num_helices_per_bundle_ = tag->getOption<core::Size>("num_helices_per_bundle");
	}
	if ( tag->hasOption("min_helix_size") ) {
		min_helix_size_ = tag->getOption<core::Size>("min_helix_size");
	}

	if ( tag->hasOption("helix_cap_dist_cutoff") ) {
		helix_cap_dist_cutoff_ = tag->getOption<core::Real>("helix_cap_dist_cutoff");
	}

	min_per_residue_fa_attr_ = tag->getOption<core::Real>("min_per_residue_fa_attr", 0.2);
	min_interacting_set_fraction_ = tag->getOption<core::Real>("min_interacting_set_fraction", 0.4);
	// max_degrees_off_parallel_ = tag->getOption<core::Real>("max_degrees_off_parallel", 40);

}

core::Size
HelixBundleFeatures::report_features(
	core::pose::Pose const & pose,
	utility::vector1<bool> const &,
	StructureID struct_id,
	utility::sql_database::sessionOP db_session
){

	//Get all the helices from the database
	utility::vector1<HelicalFragmentOP> helices = get_helices(struct_id, db_session);
	if ( helices.size()<num_helices_per_bundle_ ) { return 0; }


	//Iterate through all sets of n-helices (where n is number of helices per bundle)
	core::Size bundle_counter=0;
	utility::vector1<core::Size> dim_sizes(num_helices_per_bundle_, helices.size());
	utility::LexicographicalIterator lex( dim_sizes );
	while ( ! lex.at_end() ) {
		//Only evaluate the upper diagonal
		bool broke=false;
		for ( core::Size i=1; i<num_helices_per_bundle_; ++i ) {
			for ( core::Size j=i+1; j<=num_helices_per_bundle_; ++j ) {
				if ( lex[i] >= lex[j] ) {
					lex.continue_at_dimension(j);
					broke=true;
					break;
				}
			}
			if ( broke ) {
				break;
			}
		}
		if ( broke ) {
			continue;
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Working with helix group " << std::endl;
			for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
				TR.Debug << *helices[i] << std::endl;;
			}
		}

		utility::vector1<HelicalFragmentOP> bundle_helices(num_helices_per_bundle_);
		utility::vector1<core::Size> offset_dims;
		for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
			bundle_helices[i] = helices[lex[i]];
			offset_dims.push_back(helices[lex[i]]->size()-1);
			offset_dims.push_back(helices[lex[i]]->size()-1);
		}

		//Now, make another lex iterator  to iterate through the possible start and end positions
		core::Size max_size=0;
		utility::vector1<HelicalFragment> largest_bundle;
		utility::LexicographicalIterator offset_lex( offset_dims );
		while ( ! offset_lex.at_end() ) {

			utility::vector1< std::pair<core::Size, core::Size> > helix_offsets(num_helices_per_bundle_);
			core::Size counter=0;
			for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
				++counter;
				helix_offsets[i].first = offset_lex[counter]-1;

				++counter;
				helix_offsets[i].second = offset_lex[counter]-1;
			}

			if ( TR.Debug.visible() ) {
				TR.Debug << "Offsets are " << std::endl;
				for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
					TR.Debug << helix_offsets[i].first << " " << helix_offsets[i].second << std::endl;;
				}
			}

			utility::vector1< HelicalFragment > trimmed_helices(num_helices_per_bundle_);

			core::Size invalid_index=0;
			for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
				core::Size trimmed_start = bundle_helices[i]->start()+helix_offsets[i].first;
				core::Size trimmed_stop = bundle_helices[i]->end()-helix_offsets[i].second;
				HelicalFragment trimmed_helix(trimmed_start, trimmed_stop);
				if ( trimmed_helix.size() < min_helix_size_ ) {
					//index is i+2 because we always want to continue at the second offset for a given helix
					invalid_index=i*2;
					break;
				}
				trimmed_helices[i] = trimmed_helix;
			}
			if ( invalid_index!=0 ) {
				if ( TR.Debug.visible() ) { TR << "Lex is continuing at dimension " << invalid_index << std::endl; }
				offset_lex.continue_at_dimension(invalid_index);
				continue;
			}

			//If we found a bundle where each helix meets the size cutoff, validate the bundle as a whole. If
			//it passed, find out how big it is and possible save it
			if ( validate_bundle(pose, trimmed_helices) ) {
				core::Size bundle_size=0;
				for ( core::Size i=1; i<=trimmed_helices.size(); ++i ) {
					bundle_size+=trimmed_helices[i].size();
				}
				if ( bundle_size > max_size ) {
					max_size=bundle_size;
					largest_bundle = trimmed_helices;
				}
			}

			++offset_lex;
		}
		if ( max_size != 0 ) {
			TR << "Found bundle:" << std::endl;
			for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
				TR << largest_bundle[i] << std::endl;
			}
			++bundle_counter;
			write_bundle_to_db(struct_id, db_session, largest_bundle);
		}
		++lex;
	}

	TR << "Found " << bundle_counter << " total bundles" << std::endl;
	return 0;
}



/// @details validate that the bundle fits all the requirements. Currently, the only
///requirement is that each helix in the bundle has a c-alpha, c-alpha distance less
///than the specified cutoff with at least 2 other helices in the same bundle
bool
HelixBundleFeatures::validate_bundle(
	core::pose::Pose const & pose,
	utility::vector1<HelicalFragment> const & helices
){
	core::Size dist_sq = helix_cap_dist_cutoff_*helix_cap_dist_cutoff_;

	if ( TR.Debug.visible() ) {
		TR.Debug << "Validating bundle:" << std::endl;
		for ( core::Size i=1; i<=num_helices_per_bundle_; ++i ) {
			TR.Debug << helices[i] << std::endl;;
		}

		TR.Debug << "Requirements: " << std::endl;
		TR.Debug << "\tCap distance cutoff: " << helix_cap_dist_cutoff_ << std::endl;
	}

	//Now check that each helix cap is within a certain distance of at least two other helices
	utility::vector1<core::Size> cap_contacts(helices.size(), 0);
	for ( core::Size i=1; i<=helices.size(); ++i ) {
		for ( core::Size j=i+1; j<=helices.size(); ++j ) {
			core::Size start_distance = pose.residue(helices[i].start()).atom("CA").xyz().distance_squared(pose.residue(helices[j].start()).atom("CA").xyz());
			core::Size end_distance = pose.residue(helices[i].end()).atom("CA").xyz().distance_squared(pose.residue(helices[j].end()).atom("CA").xyz());

			if ( start_distance < dist_sq && end_distance < dist_sq ) {
				cap_contacts[i]++;
				cap_contacts[j]++;
			}
		}
	}

	for ( core::Size i=1; i<=cap_contacts.size(); ++i ) {
		if ( cap_contacts[i] < 2 ) {
			return false;
		}
	}

	return true;
}


/// @details Write the given bundle to the database. This involves writing
///to three tables:
///helix_bundles - the top-level table. Has a unique bundle_id. FK to the structure
///bundle_helices - the individual helical segments that make up the bundle. FK to the bundle table
///bundle_pairs - the pairs of helices that make up the bundle. FKs to the bundle table and the bundle_helices table
void
HelixBundleFeatures::write_bundle_to_db(
	protocols::features::StructureID const & struct_id,
	utility::sql_database::sessionOP db_session,
	utility::vector1<HelicalFragment> const & bundle
){
	string bundle_insert =
		"INSERT INTO helix_bundles (struct_id, num_helices_per_bundle)  VALUES (?,?);";
	statement bundle_insert_stmt(basic::database::safely_prepare_statement(bundle_insert,db_session));

	string bundle_helix_insert =
		"INSERT INTO bundle_helices (struct_id, bundle_id, residue_begin, residue_end, reversed) VALUES (?,?,?,?,?);";
	statement bundle_helix_insert_stmt(basic::database::safely_prepare_statement(bundle_helix_insert,db_session));

	string pair_insert =
		"INSERT INTO helix_pairs (struct_id, bundle_id, helix_id_1, helix_id_2)\n"
		"VALUES (?,?,?,?)";
	statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert, db_session));

	bundle_insert_stmt.bind(1, struct_id);
	bundle_insert_stmt.bind(2, num_helices_per_bundle_);
	basic::database::safely_write_to_database( bundle_insert_stmt );
	core::Size bundle_id = bundle_insert_stmt.sequence_last("helix_bundles_bundle_id_seq");

	utility::vector1<core::Size> helix_ids;
	for ( core::Size i=1; i<=bundle.size(); ++i ) {
		bundle_helix_insert_stmt.bind(1, struct_id);
		bundle_helix_insert_stmt.bind(2, bundle_id);
		bundle_helix_insert_stmt.bind(3, bundle[i].seq_start());
		bundle_helix_insert_stmt.bind(4, bundle[i].seq_end());
		bundle_helix_insert_stmt.bind(5, bundle[i].reversed());
		basic::database::safely_write_to_database( bundle_helix_insert_stmt );
		helix_ids.push_back(bundle_helix_insert_stmt.sequence_last("bundle_helices_helix_id_seq"));
	}

	for ( core::Size i=1; i<=bundle.size(); ++i ) {
		for ( core::Size j=i+1; j<=bundle.size(); ++j ) {
			pair_insert_stmt.bind(1, struct_id);
			pair_insert_stmt.bind(2, bundle_id);
			pair_insert_stmt.bind(3, helix_ids[i]);
			pair_insert_stmt.bind(4, helix_ids[j]);
			basic::database::safely_write_to_database( pair_insert_stmt );
		}
	}

}

std::string HelixBundleFeatures::type_name() const {
	return class_name();
}

std::string HelixBundleFeatures::class_name() {
	return "HelixBundleFeatures";
}

void HelixBundleFeatures::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute("num_helices_per_bundle", xsct_non_negative_integer, "Number of helices to put in each bundle")
		+ XMLSchemaAttribute("min_helix_size", xsct_non_negative_integer, "Minimum size of helix to include in a bundle")
		+ XMLSchemaAttribute("helix_cap_dist_cutoff", xsct_real, "Maximum distance between helix caps in the bundle")
		+ XMLSchemaAttribute::attribute_w_default("min_per_residue_fa_attr", xsct_real, "Minimum fa_attr score per residue for a bundle", "0.2")
		+ XMLSchemaAttribute::attribute_w_default("min_interacting_set_fraction", xsct_real, "Minimum fraction of residues in the bundle that interact with another helix", "0.4");

	protocols::features::xsd_type_definition_w_attributes( xsd, class_name(), "Forms helical bundles from helices in a database and outputs bundles to a database", attlist );
}

std::string HelixBundleFeaturesCreator::type_name() const {
	return HelixBundleFeatures::class_name();
}

protocols::features::FeaturesReporterOP
HelixBundleFeaturesCreator::create_features_reporter() const {
	return protocols::features::FeaturesReporterOP( new HelixBundleFeatures );
}

void HelixBundleFeaturesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HelixBundleFeatures::provide_xml_schema( xsd );
}


} //namespace helixAssembly
} //namespace features
} //namespace protocols
