// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/StrandBundleFeatures.cc
/// @brief Search through a pose for sets of 2 strands 
/// @author Doo Nam Kim (based on Tim Jacobs' helixAssembly)


//Core
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

//External
#include <boost/uuid/uuid.hpp>

//Devel
#include <protocols/features/strand_assembly/StrandBundleFeatures.hh>
#include <protocols/features/strand_assembly/StrandFragment.hh>

//Utility and basic
#include <basic/database/sql_utils.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <numeric/xyz.functions.hh> // for torsion calculations

//C library
#include <string>
#include <math.h>


//External Headers
#include <cppdb/frontend.h>

//Basic
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/strand_assembly.OptionKeys.gen.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/Constraint.hh>

#include <protocols/analysis/InterfaceAnalyzerMover.hh> // for SASA
#include <core/scoring/ScoreFunction.hh> // ScoreFunction.hh seems required for compilation of InterfaceAnalyzerMover.hh

static basic::Tracer TR("protocols.features.strand_assembly.StrandBundleFeatures");

namespace protocols {
namespace features {
namespace strand_assembly {

using namespace std;
using namespace core;
using core::pose::Pose;
using utility::vector1;
using utility::sql_database::sessionOP;
using cppdb::statement;
using cppdb::result;

StrandBundleFeatures::StrandBundleFeatures() :
min_strand_size_(4),
max_strand_size_(15),
min_O_N_dis_(2.6),
max_O_N_dis_(3.1),
min_sheet_dis_(7.0),
max_sheet_dis_(14.0), // 15 Angstrom may seem OK, but I set it to 14, because 15 cases are not easy to make assemblies as of now
min_sheet_torsion_(-40.0), // according to swissmodel.expasy.org/course/text/chapter4.htm, -20 < torsion < -50, but I set it to -40, because -50 cases are not easy to make assemblies  as of now
max_sheet_torsion_(-20.0)
{
	init_from_options();
}

void StrandBundleFeatures::init_from_options(){
	using namespace basic::options;

	if(option[OptionKeys::strand_assembly::min_strand_size].user()){
		min_strand_size_ = option[OptionKeys::strand_assembly::min_strand_size];
	}
	if(option[OptionKeys::strand_assembly::max_strand_size].user()){
		max_strand_size_ = option[OptionKeys::strand_assembly::max_strand_size];
	}
	if(option[OptionKeys::strand_assembly::min_O_N_dis].user()){
		min_O_N_dis_ = option[OptionKeys::strand_assembly::min_O_N_dis];
	}
	if(option[OptionKeys::strand_assembly::max_O_N_dis].user()){
		max_O_N_dis_ = option[OptionKeys::strand_assembly::max_O_N_dis];
	}
	if(option[OptionKeys::strand_assembly::min_sheet_dis].user()){
		min_sheet_dis_ = option[OptionKeys::strand_assembly::min_sheet_dis];
	}
	if(option[OptionKeys::strand_assembly::max_sheet_dis].user()){
		max_sheet_dis_ = option[OptionKeys::strand_assembly::max_sheet_dis];
	}
	if(option[OptionKeys::strand_assembly::min_sheet_torsion].user()){
		min_sheet_torsion_ = option[OptionKeys::strand_assembly::min_sheet_torsion];
	}
	if(option[OptionKeys::strand_assembly::max_sheet_torsion].user()){
		max_sheet_torsion_ = option[OptionKeys::strand_assembly::max_sheet_torsion];
	}
}

utility::vector1<std::string>
StrandBundleFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("ResidueFeatures");
	dependencies.push_back("ProteinResidueConformationFeatures");
	dependencies.push_back("ResidueSecondaryStructureFeatures");
	return dependencies;
}

void
StrandBundleFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const{
	using namespace basic::database::schema_generator;

	/****** write strand_pairs ******/

		// PrimaryKey
			// id of strand_pairs
			Column pairs_id              ("pairs_id",               DbInteger(), false /*not null*/, true /*autoincrement*/);

			// bool_parallel is PrimaryKey just because it doesn't point to any foreign values
			Column bool_parallel         ("bool_parallel",          DbInteger(), false /*not null*/, false /*don't autoincrement*/);

		// unique key of original PDB file
			Column struct_id             ("struct_id",              DbUUID(),    false /*not null*/, false /*don't autoincrement*/);

		// ForeignKey
			Column i_strand_residue_start("i_strand_residue_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column i_strand_residue_end  ("i_strand_residue_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column j_strand_residue_start("j_strand_residue_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column j_strand_residue_end  ("j_strand_residue_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);

		// Schema
		// PrimaryKey
			Schema strand_pairs("strand_pairs",  PrimaryKey(pairs_id));
//			Schema strand_pairs("bool_parallel", PrimaryKey(bool_parallel));

		// ForeignKey

			strand_pairs.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));
			// (reference) wiki.rosettacommons.org/index.php/MultiBodyFeaturesReporters#StructureFeatures

			// fkey_reference_cols -> one definition of fkey_reference_cols is enough since it will be used repeatedly
				utility::vector1<std::string> fkey_reference_cols;
				fkey_reference_cols.push_back("struct_id");
				fkey_reference_cols.push_back("resNum");

			// i_start_fkey_cols
				utility::vector1<Column> i_start_fkey_cols;
				i_start_fkey_cols.push_back(struct_id);
				i_start_fkey_cols.push_back(i_strand_residue_start);

				strand_pairs.add_foreign_key(ForeignKey(i_start_fkey_cols, "residues",	fkey_reference_cols,	true /*defer*/));
				// (reference) wiki.rosettacommons.org/index.php/OneBodyFeaturesReporters#ResidueFeatures

			// i_end_fkey_cols
				utility::vector1<Column> i_end_fkey_cols;
				i_end_fkey_cols.push_back(struct_id);
				i_end_fkey_cols.push_back(i_strand_residue_end);

				strand_pairs.add_foreign_key(ForeignKey(i_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

			// j_start_fkey_cols
				utility::vector1<Column> j_start_fkey_cols;
				j_start_fkey_cols.push_back(struct_id);
				j_start_fkey_cols.push_back(j_strand_residue_start);

				strand_pairs.add_foreign_key(ForeignKey(j_start_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

			// j_end_fkey_cols
				utility::vector1<Column> j_end_fkey_cols;
				j_end_fkey_cols.push_back(struct_id);
				j_end_fkey_cols.push_back(j_strand_residue_end);

				strand_pairs.add_foreign_key(ForeignKey(j_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		strand_pairs.write(db_session);


	/****** write sheet_pairs ******/

		// PrimaryKey
			Column sheet_pairs_id              ("sheet_pairs_id", DbInteger(), false /*not null*/, true /*autoincrement*/);

		// group1 (g1) should be close pair with each other, group2 should be close pair with each other
		// ForeignKey
			// group 1
				// strand 1
			Column g1_strand_1_res_start  ("g1_strand_1_res_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column g1_strand_1_res_end    ("g1_strand_1_res_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);
				// strand 2			
			Column g1_strand_2_res_start  ("g1_strand_2_res_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column g1_strand_2_res_end    ("g1_strand_2_res_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			
			// group 2
			Column g2_strand_1_res_start  ("g2_strand_1_res_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column g2_strand_1_res_end    ("g2_strand_1_res_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column g2_strand_2_res_start  ("g2_strand_2_res_start", DbInteger(), false /*not null*/, false /*don't autoincrement*/);
			Column g2_strand_2_res_end    ("g2_strand_2_res_end",   DbInteger(), false /*not null*/, false /*don't autoincrement*/);

			Column shortest_sc_dis	("shortest_sc_dis",	DbDouble(), false /*not null*/, false /*don't autoincrement*/);

		// Schema
		// PrimaryKey
			Schema sheet_pairs("sheet_pairs", PrimaryKey(sheet_pairs_id));

		// ForeignKey
			sheet_pairs.add_foreign_key(ForeignKey(struct_id,	"structures",	"struct_id",	true /*defer*/));

		// g1_s1_start_fkey_cols
			utility::vector1<Column> g1_s1_start_fkey_cols;
			g1_s1_start_fkey_cols.push_back(struct_id);
			g1_s1_start_fkey_cols.push_back(g1_strand_1_res_start);

			sheet_pairs.add_foreign_key(ForeignKey(g1_s1_start_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		// g1_s1_end_fkey_cols
			utility::vector1<Column> g1_s1_end_fkey_cols;
			g1_s1_end_fkey_cols.push_back(struct_id);
			g1_s1_end_fkey_cols.push_back(g1_strand_1_res_end);

			sheet_pairs.add_foreign_key(ForeignKey(g1_s1_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

	
		// g1_s2_start_fkey_cols
			utility::vector1<Column> g1_s2_start_fkey_cols;
			g1_s2_start_fkey_cols.push_back(struct_id);
			g1_s2_start_fkey_cols.push_back(g1_strand_2_res_start);

			sheet_pairs.add_foreign_key(ForeignKey(g1_s2_start_fkey_cols,	"residues",	fkey_reference_cols, true /*defer*/));

		// g1_s2_end_fkey_cols
			utility::vector1<Column> g1_s2_end_fkey_cols;
			g1_s2_end_fkey_cols.push_back(struct_id);
			g1_s2_end_fkey_cols.push_back(g1_strand_2_res_end);

			sheet_pairs.add_foreign_key(ForeignKey(g1_s2_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		// g2_s1_start_fkey_cols
			utility::vector1<Column> g2_s1_start_fkey_cols;
			g2_s1_start_fkey_cols.push_back(struct_id);
			g2_s1_start_fkey_cols.push_back(g2_strand_1_res_start);

			sheet_pairs.add_foreign_key(ForeignKey(g2_s1_start_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		// g2_s1_end_fkey_cols
			utility::vector1<Column> g2_s1_end_fkey_cols;
			g2_s1_end_fkey_cols.push_back(struct_id);
			g2_s1_end_fkey_cols.push_back(g2_strand_1_res_end);

			sheet_pairs.add_foreign_key(ForeignKey(g2_s1_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		// g2_s2_start_fkey_cols
			utility::vector1<Column> g2_s2_start_fkey_cols;
			g2_s2_start_fkey_cols.push_back(struct_id);
			g2_s2_start_fkey_cols.push_back(g2_strand_2_res_start);

			sheet_pairs.add_foreign_key(ForeignKey(g2_s2_start_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		// g2_s2_end_fkey_cols
			utility::vector1<Column> g2_s2_end_fkey_cols;
			g2_s2_end_fkey_cols.push_back(struct_id);
			g2_s2_end_fkey_cols.push_back(g2_strand_2_res_end);

			sheet_pairs.add_foreign_key(ForeignKey(g2_s2_end_fkey_cols,	"residues",	fkey_reference_cols,	true /*defer*/));

		sheet_pairs.add_column(shortest_sc_dis);

		sheet_pairs.write(db_session);

}

//Select all strand segments reported by the ResidueSecondaryStructureFeatures and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_full_strands(boost::uuids::uuid struct_id, sessionOP db_session){
	std::string select_string =
	"SELECT\n"
	"	strands.beta_id,\n"
	"	strands.residue_begin,\n"
	"	strands.residue_end\n"
	"FROM\n"
	"	beta_segments as strands\n"
	"WHERE\n"
	"	strands.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> all_strands;
	while(res.next()){
		Size strand_id,     residue_begin,   residue_end;
		res >> strand_id >> residue_begin >> residue_end;
		all_strands.push_back(StrandFragment(residue_begin, residue_end));
	}

	return all_strands;
}

//Select all i strand reported by the StrandBundleFeatures::get_full_strands and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_i_strand_from_full_strand_pairs(boost::uuids::uuid struct_id, sessionOP db_session){
	std::string select_string =
	"SELECT\n"
	"	strand_pairs_table.pairs_id,\n"
	"	strand_pairs_table.i_strand_residue_start,\n"
	"	strand_pairs_table.i_strand_residue_end\n"
	"FROM\n"
	"	strand_pairs as strand_pairs_table\n"
	"WHERE\n"
	"	strand_pairs_table.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> i_strand;
	while(res.next()){
		Size   strand_pairs_id,   i_strand_residue_start,   i_strand_residue_end;
		res >> strand_pairs_id >> i_strand_residue_start >> i_strand_residue_end;
		i_strand.push_back(StrandFragment(i_strand_residue_start, i_strand_residue_end));
	}
	return i_strand;
}

//Select all j strand reported by the StrandBundleFeatures::get_full_strands and save them in a vector
utility::vector1<StrandFragment> StrandBundleFeatures::get_j_strand_from_full_strand_pairs(boost::uuids::uuid struct_id, sessionOP db_session){
	std::string select_string =
	"SELECT\n"
	"	strand_pairs_table.pairs_id,\n"
	"	strand_pairs_table.j_strand_residue_start,\n"
	"	strand_pairs_table.j_strand_residue_end\n"
	"FROM\n"
	"	strand_pairs as strand_pairs_table\n"
	"WHERE\n"
	"	strand_pairs_table.struct_id = ?;";

	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,struct_id);
	result res(basic::database::safely_read_from_database(select_statement));

	utility::vector1<StrandFragment> j_strand;
	while(res.next()){
		Size   strand_pairs_id,   j_strand_residue_start,   j_strand_residue_end;
		res >> strand_pairs_id >> j_strand_residue_start >> j_strand_residue_end;
		j_strand.push_back(StrandFragment(j_strand_residue_start, j_strand_residue_end));
	}
	return j_strand;
}



bool
StrandBundleFeatures::find_antiparallel(
									Pose const & pose,
									StrandFragment strand_i,
									StrandFragment strand_j
									)
{
//	TR.Info << "strand_i.get_start() : " << strand_i.get_start() << endl;
//	TR.Info << "strand_i.get_end() : " << strand_i.get_end() << endl;
//	TR.Info << "strand_j.get_start() : " << strand_j.get_start() << endl;
//	TR.Info << "strand_j.get_end() : " << strand_j.get_end() << endl;

	// seeing distances between 'O' of strand "i" and 'N' of strand "j"
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{

			Size j_resnum = strand_j.get_start()+strand_j_res;
			//			 TR.Info << "residue number (strand_i.get_start()+strand_i_res): " << i_resnum << endl;
			//			 TR.Info << "residue number (strand_j.get_start()+strand_j_res): " << j_resnum << endl;
		//				 TR.Info << "pose.residue(i_resnum): " << pose.residue(i_resnum) << endl;
		//				 TR.Info << "pose.residue(j_resnum): " << pose.residue(j_resnum) << endl;

			Real dis_N_O = pose.residue(i_resnum).atom("N").xyz().distance(pose.residue(j_resnum).atom("O").xyz());
//			TR.Info << "distance between resnum("<< i_resnum << ")'s N and resnum(" << j_resnum << ")'s O = " << dis_N_O << endl;
			if (dis_N_O > min_O_N_dis_ && dis_N_O < max_O_N_dis_) {
//			TR.Info << "<first check of anti-parallel> distance between N and O is within " <<min_O_N_dis_ << " and " << max_O_N_dis_ << endl;
				Real dis_O_N_confirm_antiparallel = pose.residue(i_resnum).atom("O").xyz().distance(pose.residue(j_resnum).atom("N").xyz());
//				TR.Info << "distance between resnum("<< i_resnum << ")'s O and resnum(" << j_resnum << ")'s N = " << dis_O_N_confirm_antiparallel << endl;
				if (dis_O_N_confirm_antiparallel > min_O_N_dis_ && dis_O_N_confirm_antiparallel < max_O_N_dis_) {
//					TR.Info << "<second check of antiparallel> distance between O and N is within " <<min_O_N_dis_ << " and " << max_O_N_dis_ << " too, anti-parallel found!" << endl;
					return true;
				}
			}
		}
	}
	return false;
} //StrandBundleFeatures::find_antiparallel


bool
StrandBundleFeatures::find_parallel(
									Pose const & pose,
									StrandFragment strand_i,
									StrandFragment strand_j
									)
{
	 // seeing distances between 'O' of strand "i" and 'N' of strand "j"
	 for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	 {
		 Size i_resnum = strand_i.get_start()+strand_i_res;
		 for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		 {

			 Size j_resnum = strand_j.get_start()+strand_j_res;
//			 TR.Info << "residue number (strand_i.get_start()+strand_i_res): " << i_resnum << endl;
//			 TR.Info << "residue number (strand_j.get_start()+strand_j_res): " << j_resnum << endl;
//			 TR.Info << "pose.residue(i_resnum): " << pose.residue(i_resnum) << endl;
//			 TR.Info << "pose.residue(j_resnum): " << pose.residue(j_resnum) << endl;

			 Real dis_O_N = pose.residue(i_resnum).atom("O").xyz().distance(pose.residue(j_resnum).atom("N").xyz());
		//	 TR.Info << "distance between resnum("<< i_resnum << ")'s O and resnum(" << j_resnum << ")'s N = " << dis_O_N << endl;
			 if (dis_O_N > min_O_N_dis_ && dis_O_N < max_O_N_dis_) {
		//		 TR.Info << "<first check of parallel> distance between O and N within min_O_N_dis_ and max_O_N_dis_" << endl;
				 Real dis_N_O_confirm_parallel = pose.residue(i_resnum+2).atom("N").xyz().distance(pose.residue(j_resnum).atom("O").xyz());
	 	//		 TR.Info << "distance between resnum("<< i_resnum+2 << ")'s N and resnum(" << j_resnum << ")'s O = " << dis_N_O_confirm_parallel << endl;
				 if (dis_N_O_confirm_parallel > min_O_N_dis_ && dis_N_O_confirm_parallel < max_O_N_dis_) {
		//			 TR.Info << "<second check of parallel> distance between N and O within min_O_N_dis_ and max_O_N_dis_ too, parallel found!" << endl;
					 return true;
				 }
			 }
		  }
	 }
	return false;
} //StrandBundleFeatures::find_parallel




Real
StrandBundleFeatures::check_sheet_dis_antiparallel(
											   Pose const & pose,
											   StrandFragment strand_i,
											   StrandFragment strand_j
											   )
{
	// check anti-parallel sheet distance

//	TR.Info << "strand_i size                  : " << strand_i.get_size() << endl;
//	TR.Info << "strand_i first  residue number : " << strand_i.get_start() << endl;
//	TR.Info << "strand_i ending residue number : " << strand_i.get_end() << endl;

//	TR.Info << "strand_j size                  : " << strand_j.get_size() << endl;
//	TR.Info << "strand_j first  residue number : " << strand_j.get_start() << endl;
//	TR.Info << "strand_j ending residue number : " << strand_j.get_end() << endl;

	// first, check the shortest distance between the two strand_pairs
	// seeing distances between 'CA' of strand "i" and 'CA' of strand "j"
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;

			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
		//	TR.Info << "distance between resnum("<< i_resnum << ")'s CA and resnum(" << j_resnum << ")'s CA = " << dis_CA_CA << endl;

			if (dis_CA_CA < 6.0)
			{
			//	TR.Info << "these two pair of strands are within 6 Angstrom <check_sheet_dis_antiparallel>" << endl;
				return -99;
			}
		}
	}

//	TR.Info << "OK, these two strand_pairs are farther than 6 Angstrom <check_sheet_dis_antiparallel>" << endl;

//	TR.Info << "let me see distances between strands " << endl;
//	TR.Info << "pose: " << pose << endl;
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			//TR.Info << "residue number (strand_i.get_start()+strand_i_res): " << i_resnum << endl;
			//TR.Info << "residue number (strand_j.get_start()+strand_j_res): " << j_resnum << endl;
			//TR.Info << "pose.residue(i_resnum): " << pose.residue(i_resnum) << endl;
			//TR.Info << "pose.residue(j_resnum): " << pose.residue(j_resnum) << endl;


			Real dis_CA_CA_1 = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			//TR.Info << "distance between resnum("<< i_resnum << ")'s CA and resnum(" << j_resnum << ")'s CA = " << dis_CA_CA << endl;
			if (dis_CA_CA_1 > min_sheet_dis_ && dis_CA_CA_1 < max_sheet_dis_)
			{
				//TR.Info << "<1st check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
				Size current_i_resnum = i_resnum+1;
				Size current_j_resnum = j_resnum-1;
				if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
				{
					Real dis_CA_CA_2 = pose.residue(i_resnum+1).atom("CA").xyz().distance(pose.residue(j_resnum-1).atom("CA").xyz());
					//TR.Info << "distance between resnum("<< i_resnum+1 << ")'s CA and resnum(" << j_resnum-1 << ")'s CA = " << dis_CA_CA << endl;
					if (dis_CA_CA_2 > min_sheet_dis_ && dis_CA_CA_2 < max_sheet_dis_)
					{
					//	TR.Info << "<2nd check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
						Size current_i_resnum = i_resnum+2;
						Size current_j_resnum = j_resnum-2;
						if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
						{
							Real dis_CA_CA_3 = pose.residue(i_resnum+2).atom("CA").xyz().distance(pose.residue(j_resnum-2).atom("CA").xyz());
					//		TR.Info << "distance between resnum("<< i_resnum+2 << ")'s CA and resnum(" << j_resnum-2 << ")'s CA = " << dis_CA_CA << endl;
							if (dis_CA_CA_3 > min_sheet_dis_ && dis_CA_CA_2 < max_sheet_dis_)
							{
					//			TR.Info << "<3rd check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
								Size current_i_resnum = i_resnum+3;
								Size current_j_resnum = j_resnum-3;
								if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
								{
									Real dis_CA_CA_4 = pose.residue(i_resnum+3).atom("CA").xyz().distance(pose.residue(j_resnum-3).atom("CA").xyz());
							//		TR.Info << "distance between resnum("<< i_resnum+3 << ")'s CA and resnum(" << j_resnum-3 << ")'s CA = " << dis_CA_CA << endl;
									if (dis_CA_CA_4 > min_sheet_dis_ && dis_CA_CA_4 < max_sheet_dis_)
									{
						//				TR.Info << "<4th check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
										Real avg_dis_CA_CA = (dis_CA_CA_1 + dis_CA_CA_2 + dis_CA_CA_3 + dis_CA_CA_4)/4;
										return avg_dis_CA_CA;
									}
								}
								else
								{
									TR.Info << "maybe no residue here" << endl;
									return -99;
								}
							}
						}
						else
						{
							TR.Info << "maybe no residue here" << endl;
							return -99;
						}
					}
				}
				else
				{
					TR.Info << "maybe no residue here" << endl;
					return -99;
				}
			}
		}
	}
	return -99;
} //StrandBundleFeatures::check_sheet_dis_antiparallel



Real
StrandBundleFeatures::check_sheet_dis_parallel(
												   Pose const & pose,
												   StrandFragment strand_i,
												   StrandFragment strand_j
												   ) // check parallel sheet distance
{


//	TR.Info << "strand_i size                  : " << strand_i.get_size() << endl;
//	TR.Info << "strand_i first  residue number : " << strand_i.get_start() << endl;
//	TR.Info << "strand_i ending residue number : " << strand_i.get_end() << endl;

//	TR.Info << "strand_j size                  : " << strand_j.get_size() << endl;
//	TR.Info << "strand_j first  residue number : " << strand_j.get_start() << endl;
//	TR.Info << "strand_j ending residue number : " << strand_j.get_end() << endl;


	// first, check the shortest distance between the two strand_pairs
	// seeing distances between 'CA' of strand "i" and 'CA' of strand "j"
	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;

			Real dis_CA_CA = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
		//	TR.Info << "distance between resnum("<< i_resnum << ")'s CA and resnum(" << j_resnum << ")'s CA = " << dis_CA_CA << endl;

			if (dis_CA_CA < 6.0)
			{
		//		TR.Info << "these two pair of strands are within 6 Angstrom <check_sheet_dis_parallel>" << endl;
				return -99;
			}
		}
	}

//	TR.Info << "OK, these two strand_pairs are farther than 6 Angstrom <check_sheet_dis_parallel>" << endl;

	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;
		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			Size j_resnum = strand_j.get_start()+strand_j_res;
			//	TR.Info << "residue number (strand_i.get_start()+strand_i_res): " << i_resnum << endl;
			//	TR.Info << "residue number (strand_j.get_start()+strand_j_res): " << j_resnum << endl;
			//	TR.Info << "pose.residue(i_resnum): " << pose.residue(i_resnum) << endl;
			//	TR.Info << "pose.residue(j_resnum): " << pose.residue(j_resnum) << endl;



			Real dis_CA_CA_1 = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
		//	TR.Info << "distance between resnum("<< i_resnum << ")'s CA and resnum(" << j_resnum << ")'s CA = " << dis_CA_CA << endl;
			if (dis_CA_CA_1 > min_sheet_dis_ && dis_CA_CA_1 < max_sheet_dis_)
			{
		//		TR.Info << "<1st check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
				Size current_i_resnum = i_resnum+1;
				Size current_j_resnum = j_resnum+1;
				if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
				{
					Real dis_CA_CA_2 = pose.residue(i_resnum+1).atom("CA").xyz().distance(pose.residue(j_resnum+1).atom("CA").xyz());
		//			TR.Info << "distance between resnum("<< i_resnum+1 << ")'s CA and resnum(" << j_resnum+1 << ")'s CA = " << dis_CA_CA << endl;
					if (dis_CA_CA_2 > min_sheet_dis_ && dis_CA_CA_2 < max_sheet_dis_)
					{
		//				TR.Info << "<2nd check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
						Size current_i_resnum = i_resnum+2;
						Size current_j_resnum = j_resnum+2;
						if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
						{
							Real dis_CA_CA_3 = pose.residue(i_resnum+2).atom("CA").xyz().distance(pose.residue(j_resnum+2).atom("CA").xyz());
					//		TR.Info << "distance between resnum("<< i_resnum+2 << ")'s CA and resnum(" << j_resnum+2 << ")'s CA = " << dis_CA_CA << endl;
							if (dis_CA_CA_3 > min_sheet_dis_ && dis_CA_CA_3 < max_sheet_dis_)
							{
					//			TR.Info << "<3rd check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;
								Size current_i_resnum = i_resnum+3;
								Size current_j_resnum = j_resnum+3;
								if (current_i_resnum > 0 && current_i_resnum <= pose.total_residue() && current_j_resnum > 0 && current_j_resnum <= pose.total_residue())
								{
									Real dis_CA_CA_4 = pose.residue(i_resnum+3).atom("CA").xyz().distance(pose.residue(j_resnum+3).atom("CA").xyz());
						//			TR.Info << "distance between resnum("<< i_resnum+3 << ")'s CA and resnum(" << j_resnum+3 << ")'s CA = " << dis_CA_CA << endl;
									if (dis_CA_CA_4 > min_sheet_dis_ && dis_CA_CA_4 < max_sheet_dis_)
									{
						//				TR.Info << "<4th check of sheet> distance between CA and CA is within " << min_sheet_dis_ << " and " << max_sheet_dis_ << endl;

										Real avg_dis_CA_CA = (dis_CA_CA_1 + dis_CA_CA_2 + dis_CA_CA_3 + dis_CA_CA_4)/4;
										return avg_dis_CA_CA;
									}
								}
								else
								{
									TR.Info << "maybe no residue here" << endl;
									return -99;
								}
							}
						}
						else
						{
							TR.Info << "maybe no residue here" << endl;
							return -99;
						}
					}
				}
				else
				{
					TR.Info << "maybe no residue here" << endl;
					return -99;
				}
			}
		}
	}
	return -99;
} //StrandBundleFeatures::check_sheet_dis_parallel


Real
StrandBundleFeatures::sheet_torsion(
									Pose const & pose,
									StrandFragment strand_i,
									StrandFragment strand_j
									)
{
	// check anti-parallel sheet
	//TR.Info << "strand_i size                  : " << strand_i.get_size() << endl;
	//	TR.Info << "strand_i first residue number  : " << strand_i.get_start() << endl;
	//TR.Info << "pose.residue(strand_i.get_start()): " << pose.residue(strand_i.get_start()) << endl;
	//	TR.Info << "strand_i ending residue number : " << strand_i.get_end() << endl;
	//TR.Info << "pose.residue(strand_i.get_end()): " << pose.residue(strand_i.get_end()) << endl;

	//	TR.Info << "strand_j size                  : " << strand_j.get_size() << endl;
	//	TR.Info << "strand_j first residue number  : " << strand_j.get_start() << endl;
	//TR.Info << "pose.residue(strand_j.get_start()): " << pose.residue(strand_j.get_start()) << endl;
	//	TR.Info << "strand_j ending residue number : " << strand_j.get_end() << endl;
	//TR.Info << "pose.residue(strand_j.get_end()): " << pose.residue(strand_j.get_end()) << endl;

	//I need to define 4 terminal residues correctly for torsion calculation
	Real dis_i_end_and_j_start = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_start()).atom("CA").xyz());
	Real dis_i_end_and_j_end   = pose.residue(strand_i.get_end()).atom("CA").xyz().distance(pose.residue(strand_j.get_end()).atom("CA").xyz());
	if (dis_i_end_and_j_start < dis_i_end_and_j_end)
	{
		//TR.Info << "dis_i_end_and_j_start < dis_i_end_and_j_end" << endl;

		Vector const& first_xyz    ( pose.residue(strand_i.get_start()).xyz("CA") );
		//		TR.Info << "first_xyz : " << first_xyz << endl;
		//		TR.Info << "pose.residue(strand_i.get_start()).xyz(CA) : " << pose.residue(strand_i.get_start()).xyz("CA") << endl;

		Vector const& second_xyz   ( pose.residue(strand_i.get_end()).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(strand_j.get_start()).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(strand_j.get_end()).xyz("CA") );

		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		Real torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
		//		TR.Info << "torsion_i_j : " << torsion_i_j << endl;

		return torsion_i_j;

	}

	else
	{
		//TR.Info << "dis_i_end_and_j_start >= dis_i_end_and_j_end" << endl;
		Vector const& first_xyz    ( pose.residue(strand_i.get_start()).xyz("CA") );
		Vector const& second_xyz   ( pose.residue(strand_i.get_end()).xyz("CA") );
		Vector const& third_xyz    ( pose.residue(strand_j.get_end()).xyz("CA") );
		Vector const& fourth_xyz   ( pose.residue(strand_j.get_start()).xyz("CA") );


		// calculates a torsion angles between four atoms of 'CA' of strand "i" and 'CA' of strand "j"
		Real torsion_i_j = numeric::dihedral_degrees(first_xyz, second_xyz, third_xyz, fourth_xyz);
		//TR.Info << "torsion_i_j : " << torsion_i_j << endl;

		return torsion_i_j;
	}

} //StrandBundleFeatures::sheet_torsion

Real
StrandBundleFeatures::shortest_dis_sidechain(
									Pose const & pose,
									StrandFragment strand_i,
									StrandFragment strand_j
									) // calculate shortest distances for each pair
{
	Real temp_shortest_dis = 9999;

	//TR.Info << "strand_i size                  : " << strand_i.get_size() << endl;
	//TR.Info << "strand_i first  residue number : " << strand_i.get_start() << endl;
	//TR.Info << "strand_i ending residue number : " << strand_i.get_end() << endl;

	//TR.Info << "strand_j size                  : " << strand_j.get_size() << endl;
	//TR.Info << "strand_j first  residue number : " << strand_j.get_start() << endl;
	//TR.Info << "strand_j ending residue number : " << strand_j.get_end() << endl;


	for(Size strand_i_res=0; strand_i_res < strand_i.get_size(); strand_i_res++)
	{
		Size i_resnum = strand_i.get_start()+strand_i_res;

		for(Size strand_j_res=0; strand_j_res < strand_j.get_size(); strand_j_res++)
		{
			//TR.Info << "pose.residue_type(i_resnum).all_sc_atoms(): " << pose.residue_type(i_resnum).all_sc_atoms() << endl;

			Size j_resnum = strand_j.get_start()+strand_j_res;

			//TR.Info << "i_resnum : " << i_resnum << endl;
			//TR.Info << "pose.residue_type(i_resnum).name3(): " << pose.residue_type(i_resnum).name3() << endl;

			//TR.Info << "j_resnum : " << j_resnum << endl;
			//TR.Info << "pose.residue_type(j_resnum).name3(): " << pose.residue_type(j_resnum).name3() << endl;

			//TR.Info << "pose.residue(i_resnum).xyz(1): " << pose.residue(i_resnum).xyz(1) << endl;
			//TR.Info << "pose.residue(i_resnum).xyz(2): " << pose.residue(i_resnum).xyz(2) << endl;

			//TR.Info << " pose.residue(i_resnum).natoms(): " << pose.residue(i_resnum).natoms() << endl;
			//TR.Info << " pose.residue(i_resnum).nheavyatoms(): " << pose.residue(i_resnum).nheavyatoms() << endl;

			//TR.Info << " pose.residue(j_resnum).natoms(): " << pose.residue(j_resnum).natoms() << endl;
			//TR.Info << " pose.residue(j_resnum).nheavyatoms(): " << pose.residue(j_resnum).nheavyatoms() << endl;

			for (Size i_AtomNum=1; i_AtomNum < pose.residue(i_resnum).natoms(); i_AtomNum++)
			{

				for (Size j_AtomNum=1; j_AtomNum < pose.residue(j_resnum).natoms(); j_AtomNum++)
				{

					//TR.Info << "i_AtomNum : " << i_AtomNum << endl;
					//TR.Info << "j_AtomNum : " << j_AtomNum << endl;

					//TR.Info << "pose.residue(i_resnum).atom_type_index(i_AtomNum) : " << pose.residue(i_resnum).atom_type_index(i_AtomNum) << endl;
					//TR.Info << "pose.residue(j_resnum).atom_type_index(j_AtomNum) : " << pose.residue(j_resnum).atom_type_index(j_AtomNum) << endl;
					Real dis_sc_sc = pose.residue(i_resnum).xyz(i_AtomNum).distance(pose.residue(j_resnum).xyz(j_AtomNum));
					//TR.Info << "dis_sc_sc: " << dis_sc_sc << endl;
					if (temp_shortest_dis > dis_sc_sc)
					{
						temp_shortest_dis = dis_sc_sc;
					}

				}

			}
			//TR.Info << "temp_shortest_dis between these sheet pairs until " << pose.residue_type(i_resnum).name3() << i_resnum << " and " << pose.residue_type(j_resnum).name3() << j_resnum << ": " << temp_shortest_dis << endl;


//			for (Size i_AtomNum=0; pose.residue_type(i_resnum).nheavyatoms(); i_AtomNum++)
//			{
			//	Real dis_sc_sc = pose.residue_type(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
//			}


			/*Real dis_sc_sc = pose.residue(i_resnum).atom("CA").xyz().distance(pose.residue(j_resnum).atom("CA").xyz());
			TR.Info << "dis_sc_sc (CA): " << dis_sc_sc << endl;
			if (temp_shortest_dis > dis_sc_sc)
				temp_shortest_dis = dis_sc_sc;

			Real dis_sc_sc = pose.residue(i_resnum).atom("C").xyz().distance(pose.residue(j_resnum).atom("C").xyz());
			TR.Info << "dis_sc_sc (C): " << dis_sc_sc << endl;
			if (temp_shortest_dis > dis_sc_sc)
				temp_shortest_dis = dis_sc_sc;*/
		}
	}

	return temp_shortest_dis;

} //StrandBundleFeatures::shortest_dis_sidechain


Real
StrandBundleFeatures::shortest_dis_sidechain(
											 Real val_shortest_dis_sidechain_1,
											 Real val_shortest_dis_sidechain_2,
											 Real val_shortest_dis_sidechain_3,
											 Real val_shortest_dis_sidechain_4
											 ) // shortest of shortest distances
{
	Real temp_shortest_dis = val_shortest_dis_sidechain_1;

	if (temp_shortest_dis > val_shortest_dis_sidechain_2)
	{
		{
			temp_shortest_dis = val_shortest_dis_sidechain_2;
			if (temp_shortest_dis > val_shortest_dis_sidechain_3)
			{
				temp_shortest_dis = val_shortest_dis_sidechain_3;
				if (temp_shortest_dis > val_shortest_dis_sidechain_4)
				{
					temp_shortest_dis = val_shortest_dis_sidechain_4;
				}

			}
			else
			{
				if (temp_shortest_dis > val_shortest_dis_sidechain_4)
				{
					temp_shortest_dis = val_shortest_dis_sidechain_4;
				}
			}

		}
	}

	else
	{
		if (temp_shortest_dis > val_shortest_dis_sidechain_3)
		{
			temp_shortest_dis = val_shortest_dis_sidechain_3;
			if (temp_shortest_dis > val_shortest_dis_sidechain_4)
			{
				temp_shortest_dis = val_shortest_dis_sidechain_4;
			}
		}
		else
		{
			if (temp_shortest_dis > val_shortest_dis_sidechain_4)
			{
				temp_shortest_dis = val_shortest_dis_sidechain_4;
			}
		}
	}

	return temp_shortest_dis;

} //StrandBundleFeatures::shortest_dis_sidechain (with four parameters)



Real
StrandBundleFeatures::shortest_dis_pairs(
											 Real return_of_check_sheet_dis_antiparallel_1,
											 Real return_of_check_sheet_dis_antiparallel_2,
											 Real return_of_check_sheet_dis_antiparallel_3,
											 Real return_of_check_sheet_dis_antiparallel_4
											 ) // find the shortest distances between strands for Tim's way of finding nearest strand pairs
{
	Real temp_shortest_dis = return_of_check_sheet_dis_antiparallel_1;

	if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_2)
	{
		{
			temp_shortest_dis = return_of_check_sheet_dis_antiparallel_2;
			if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_3)
			{
				temp_shortest_dis = return_of_check_sheet_dis_antiparallel_3;
				if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_4)
				{
					temp_shortest_dis = return_of_check_sheet_dis_antiparallel_4;
				}

			}
			else
			{
				if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_4)
				{
					temp_shortest_dis = return_of_check_sheet_dis_antiparallel_4;
				}
			}

		}
	}

	else
	{
		if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_3)
		{
			temp_shortest_dis = return_of_check_sheet_dis_antiparallel_3;
			if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_4)
			{
				temp_shortest_dis = return_of_check_sheet_dis_antiparallel_4;
			}
		}
		else
		{
			if (temp_shortest_dis > return_of_check_sheet_dis_antiparallel_4)
			{
				temp_shortest_dis = return_of_check_sheet_dis_antiparallel_4;
			}
		}
	}

	return temp_shortest_dis;

} //StrandBundleFeatures::shortest_dis_pairs (to find closer group of strand pairs)


///@brief collect all the feature data for the pose
core::Size
StrandBundleFeatures::report_features(
									  core::pose::Pose const & pose,
									  utility::vector1<bool> const &,
									  boost::uuids::uuid struct_id,
									  utility::sql_database::sessionOP db_session)
{
	TR.Info << "======================= <report_features begin> =========================" << endl;
//	TR.Info << "pose: " << pose << endl;



	utility::vector1<StrandFragment> all_strands = get_full_strands(struct_id, db_session);
	//TR.Info << "--------------- <sql> all strands are read --------------- " << endl;

	//TR.Info << "all_strands.size(): " << all_strands.size() << endl;

	// strand pairing begins
	for(Size i=1; i<all_strands.size(); ++i) // I don't need the last strand since this double for loops are exhaustive search for all pairs of strands
	{
//		TR.Info << "i: " << i << endl;
//		TR.Info << "all_strands[i].get_size(): " << all_strands[i].get_size() << endl;
//		TR.Info << "all_strands[i].get_start(): " << all_strands[i].get_start() << endl;
//		TR.Info << "all_strands[i].get_end(): " << all_strands[i].get_end() << endl;
		if (all_strands[i].get_size() >= min_strand_size_ && all_strands[i].get_size() <= max_strand_size_) // the legnth of this beta strand is between min_strand_size_ and max_strand_size_
		{
			for(Size j=i+1; j<=all_strands.size(); ++j) // I need the last strand for this second for loop
			{
			//	TR.Info << "j: " << j << endl;
			//	TR.Info << "all_strands[j].get_size(): " << all_strands[j].get_size() << endl;
				if (all_strands[j].get_size() >= min_strand_size_ && all_strands[j].get_size() <= max_strand_size_) // the legnth of this beta strand is between min_strand_size_ and max_strand_size_ too
				{
	//				TR.Info << "-- both strands are long enough to be considered for strands pairing --" << endl;

					StrandFragment temp_strand_i(all_strands[i].get_start(), all_strands[i].get_end());
					StrandFragment temp_strand_j(all_strands[j].get_start(), all_strands[j].get_end());

					//TR.Info << "-- anti-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand begins --- " << endl;
					bool return_of_find_antiparallel = find_antiparallel (pose, temp_strand_i, temp_strand_j);
					if (return_of_find_antiparallel)
					{

					//	TR.Info << "!!!!! i(" << i << ") strand and j(" << j << ") strand are antiparallel to each other !!!!! " << endl;
						//					TR.Info << "=========== saving antiparallel strand pairs =========== " << endl;
						string pair_insert =  "INSERT INTO strand_pairs (struct_id, i_strand_residue_start, i_strand_residue_end, j_strand_residue_start, j_strand_residue_end)  VALUES (?,?,?,?,?);";
						statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,db_session));
						pair_insert_stmt.bind(1,struct_id);
						pair_insert_stmt.bind(2,all_strands[i].get_start());
						pair_insert_stmt.bind(3,all_strands[i].get_end());
						pair_insert_stmt.bind(4,all_strands[j].get_start());
						pair_insert_stmt.bind(5,all_strands[j].get_end());
						basic::database::safely_write_to_database(pair_insert_stmt);
					}
					else
					{
					//	TR.Info << "<parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand begins>" << endl;
						bool return_of_find_parallel = find_parallel (pose, temp_strand_i, temp_strand_j);
						if (return_of_find_parallel)
						{

					//		TR.Info << "!!!!! i(" << i << ") strand and j(" << j << ") strand are parallel to each other !!!!! " << endl;
							//					TR.Info << "=========== saving parallel strand pairs =========== " << endl;
							string pair_insert =  "INSERT INTO strand_pairs (struct_id, i_strand_residue_start, i_strand_residue_end, j_strand_residue_start, j_strand_residue_end)  VALUES (?,?,?,?,?);";
							statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,db_session));
							pair_insert_stmt.bind(1,struct_id);
							pair_insert_stmt.bind(2,all_strands[i].get_start());
							pair_insert_stmt.bind(3,all_strands[i].get_end());
							pair_insert_stmt.bind(4,all_strands[j].get_start());
							pair_insert_stmt.bind(5,all_strands[j].get_end());
							basic::database::safely_write_to_database(pair_insert_stmt);
						}

					}
/* archive to be kept for bool_parallel
						TR.Info << "!!!!! i(" << i << ") strand and j(" << j << ") strand are antiparallel to each other !!!!! " << endl;
	//					TR.Info << "=========== saving antiparallel strand pairs =========== " << endl;
						string pair_insert =  "INSERT INTO strand_pairs (struct_id, bool_parallel, i_strand_residue_start, i_strand_residue_end, j_strand_residue_start, j_strand_residue_end)  VALUES (?,?,?,?,?,?);";
						statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,db_session));
						pair_insert_stmt.bind(1,struct_id);
						pair_insert_stmt.bind(2,"FALSE");
						pair_insert_stmt.bind(3,all_strands[i].get_start());
						pair_insert_stmt.bind(4,all_strands[i].get_end());
						pair_insert_stmt.bind(5,all_strands[j].get_start());
						pair_insert_stmt.bind(6,all_strands[j].get_end());
						basic::database::safely_write_to_database(pair_insert_stmt);
					}
					else{
		//				TR.Info << "-parallel check between i(" << i << ")'s strand and j(" << j << ")'s strand begins- " << endl;
						bool return_of_find_parallel = find_parallel (pose, temp_strand_i, temp_strand_j);
						if (return_of_find_parallel) {
							TR.Info << "!!!!! i(" << i << ") strand and j(" << j << ") strand are parallel to each other !!!!!! " << endl;
		//					TR.Info << "=========== saving parallel strand pairs =========== " << endl;
							string pair_insert =  "INSERT INTO strand_pairs (struct_id, bool_parallel, i_strand_residue_start, i_strand_residue_end, j_strand_residue_start, j_strand_residue_end)  VALUES (?,?,?,?,?,?);";
							statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,db_session));
							pair_insert_stmt.bind(1,struct_id);
							pair_insert_stmt.bind(2,"TRUE");
							pair_insert_stmt.bind(3,all_strands[i].get_start());
							pair_insert_stmt.bind(4,all_strands[i].get_end());
							pair_insert_stmt.bind(5,all_strands[j].get_start());
							pair_insert_stmt.bind(6,all_strands[j].get_end());
							basic::database::safely_write_to_database(pair_insert_stmt);
						}
					}
 */
				}
			}
		}
	}
	TR << "============== (Done) saving pairs of strands ==========" << endl;
	// strand pairing ends

	// sheet pairing begins
	utility::vector1<StrandFragment> i_strand_from_full_strand_pairs = get_i_strand_from_full_strand_pairs(struct_id, db_session);
	//TR.Info << "--------------- <sql> i_strand_from_full_strand_pairs are read --------------- " << endl;

	for(Size ii=1; ii<i_strand_from_full_strand_pairs.size(); ++ii) // I don't need the last pair of strands in this first 'for' loop
	{
		//TR.Info << "ii: " << ii << endl;
		//TR.Info << "i_strand_from_full_strand_pairs[ii].get_size() : " << i_strand_from_full_strand_pairs[ii].get_size()  << endl;
		//TR.Info << "i_strand_from_full_strand_pairs[ii].get_start(): " << i_strand_from_full_strand_pairs[ii].get_start() << endl;
		//TR.Info << "i_strand_from_full_strand_pairs[ii].get_end()  : " << i_strand_from_full_strand_pairs[ii].get_end()   << endl;

		utility::vector1<StrandFragment> j_strand_from_full_strand_pairs = get_j_strand_from_full_strand_pairs(struct_id, db_session);

		for(Size jj=ii+1; jj<=i_strand_from_full_strand_pairs.size(); ++jj) // I need the last pair of strands in this second 'for' loop
		{

			//TR.Info << "jj: " << jj << endl;
			//TR.Info << "j_strand_from_full_strand_pairs[jj].get_size() : " << j_strand_from_full_strand_pairs[jj].get_size()  << endl;
			//TR.Info << "j_strand_from_full_strand_pairs[jj].get_start(): " << j_strand_from_full_strand_pairs[jj].get_start() << endl;
			//TR.Info << "j_strand_from_full_strand_pairs[jj].get_end()  : " << j_strand_from_full_strand_pairs[jj].get_end()   << endl;

			//TR.Info << "<anti-parallel sheet check by distance> begins between " << ii << " th pair of strand_pair (i) and " << jj << " th strand_pair (j)" << endl;

			// all four possible combinations search
			StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
			StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
			Real return_of_check_sheet_dis_antiparallel_1 = check_sheet_dis_antiparallel (pose, temp_strand_i, temp_strand_j); // first check of distance between strand i and strand j
			if (return_of_check_sheet_dis_antiparallel_1 != -99)
			{
				StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
				StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
				Real return_of_check_sheet_dis_antiparallel_2 = check_sheet_dis_antiparallel (pose, temp_strand_i, temp_strand_j); // since the first distance between strands is within range, calculate the second distance between strands interchangeably
				if (return_of_check_sheet_dis_antiparallel_2 != -99)
				{
					StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
					StrandFragment temp_strand_j(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
					Real return_of_check_sheet_dis_antiparallel_3 = check_sheet_dis_antiparallel (pose, temp_strand_i, temp_strand_j); // since the 1st and 2nd distances between strands are within range, calculate the 3rd distance between i strands
					if (return_of_check_sheet_dis_antiparallel_3 != -99)
					{
						StrandFragment temp_strand_i(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
						StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
						Real return_of_check_sheet_dis_antiparallel_4 = check_sheet_dis_antiparallel (pose, temp_strand_i, temp_strand_j); // since 1st,2nd,3rd distances between strands are within range, calculate the 4th distance between j strands
						if (return_of_check_sheet_dis_antiparallel_4 != -99)
						{
							//TR.Info << "return_of_check_sheet_dis_antiparallel_1~4: " << return_of_check_sheet_dis_antiparallel_1 << ", " << return_of_check_sheet_dis_antiparallel_2 <<  ", " << return_of_check_sheet_dis_antiparallel_3 <<  ", " <<  return_of_check_sheet_dis_antiparallel_4 << endl;
							//TR.Info << "<sheet by distance found (by anti-parallel way)> " << ii << " th strand_pair (i) and " << jj << " th strand_pair (j) are within " << min_sheet_dis_ << " Angstrom and " << max_sheet_dis_ << " Angstrom" << endl;

							//TR.Info << "<check sheet by torsion begins (by anti-parallel way)> between " << ii << " th pair of strand_pair (i) and " << jj << " th strand_pair (j)" << endl;

							StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
							StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
							Real sheet_torsion_1 = sheet_torsion(pose, temp_strand_i, temp_strand_j);

							StrandFragment temp_strand_k(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
							StrandFragment temp_strand_l(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
							Real sheet_torsion_2 = sheet_torsion(pose, temp_strand_k, temp_strand_l);

							StrandFragment temp_strand_m(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
							StrandFragment temp_strand_n(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
							Real sheet_torsion_3 = sheet_torsion(pose, temp_strand_m, temp_strand_n);

							StrandFragment temp_strand_o(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
							StrandFragment temp_strand_p(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
							Real sheet_torsion_4 = sheet_torsion(pose, temp_strand_o, temp_strand_p);

							if (sheet_torsion_1 < 0 && sheet_torsion_2 < 0 && sheet_torsion_3 < 0 && sheet_torsion_4 < 0)
							{
							//	TR.Info << "torsions are : " << sheet_torsion_1 << ", " << sheet_torsion_2 << ", " << sheet_torsion_3 << ", " << sheet_torsion_4 << endl;
								Real sheet_torsion_avg = (sheet_torsion_1 + sheet_torsion_2 + sheet_torsion_3 + sheet_torsion_4)/4;
							//	TR.Info << "sheet_torsion_avg (by anti-parallel way) : " << sheet_torsion_avg << endl;

								if (sheet_torsion_avg > min_sheet_torsion_ && sheet_torsion_avg < max_sheet_torsion_)
								{

									TR.Info << "<sheet by torsion found (by anti-parallel way)> the average torsion between " << ii << " th strand_pair (i) and " << jj << " th strand_pair (j) is within " << min_sheet_torsion_ << " and " << max_sheet_torsion_ << endl;

									StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
									StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());

									Real val_shortest_dis_sidechain_1 = shortest_dis_sidechain (pose, temp_strand_i, temp_strand_j);
									//TR.Info << "one shortest_distance_sidechain_1 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_1 << endl;


									StrandFragment temp_strand_k(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
									StrandFragment temp_strand_l(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());

									Real val_shortest_dis_sidechain_2 = shortest_dis_sidechain (pose, temp_strand_k, temp_strand_l);
									//TR.Info << "one shortest_distance_sidechain_2 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_2 << endl;


									StrandFragment temp_strand_m(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
									StrandFragment temp_strand_n(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());

									Real val_shortest_dis_sidechain_3 = shortest_dis_sidechain (pose, temp_strand_m, temp_strand_n);
									//TR.Info << "one shortest_distance_sidechain_3 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_3 << endl;


									StrandFragment temp_strand_o(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
									StrandFragment temp_strand_p(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());

									Real val_shortest_dis_sidechain_4 = shortest_dis_sidechain (pose, temp_strand_o, temp_strand_p);
									//TR.Info << "one shortest_distance_sidechain_4 between " << ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_4 << endl;



									Real val_shortest_dis_sidechain = shortest_dis_sidechain(val_shortest_dis_sidechain_1, val_shortest_dis_sidechain_2, val_shortest_dis_sidechain_3, val_shortest_dis_sidechain_4);

									//TR.Info << "the shortest_distance_sidechain between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain << endl;

									//TR.Info << "pose: " << pose << endl;

									//protocols::simple_filters::InterfaceSasaFilter ob_SASA_filter;
									//core::pose::Pose non_const_pose = pose;
									//ob_SASA_filter.apply(non_const_pose);


									/* InterfaceAnalyzerMover
									//protocols::analysis::InterfaceAnalyzerMover ob_SASA_Mover;
									//core::pose::Pose non_const_pose = pose;
									//ob_SASA_Mover.apply(non_const_pose);
								   	//TR.Info << "ob_SASA_Mover.get_interface_delta_sasa(): " << ob_SASA_Mover.get_interface_delta_sasa() << endl;
									 */

									// to find i-i and j-i by calculating distances between pairs
									Real val_shortest_dis_pairs = shortest_dis_pairs(return_of_check_sheet_dis_antiparallel_1, return_of_check_sheet_dis_antiparallel_2, return_of_check_sheet_dis_antiparallel_3, return_of_check_sheet_dis_antiparallel_4);



									string pair_insert =
									"INSERT INTO sheet_pairs (struct_id, g1_strand_1_res_start, g1_strand_1_res_end, g1_strand_2_res_start, g1_strand_2_res_end, g2_strand_1_res_start, g2_strand_1_res_end, g2_strand_2_res_start, g2_strand_2_res_end, shortest_sc_dis)  VALUES (?, ?,?,?,?, ?,?,?,?, ?);";
									statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,	db_session));
									pair_insert_stmt.bind(1,	struct_id);

									// I need to store closer strands together as a group for Tim's graph based assembly application
									
									if (val_shortest_dis_pairs == return_of_check_sheet_dis_antiparallel_1) // first check of distance between strand i and strand j
									{
										{ // these two strands should be closer
											pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
											pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
											pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
											pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand	
										}

										{ // these two strands should be closer
											pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
											pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
											pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
											pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
										}
									}


									else if (val_shortest_dis_pairs == return_of_check_sheet_dis_antiparallel_2) // 2nd check of distance between strand i and strand j
									{
										{ // these two strands should be closer
											pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
											pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
											pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
											pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand	
										}
										
										{ // these two strands should be closer
											pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
											pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
											pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
											pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
										}
									}
									
									else if (val_shortest_dis_pairs == return_of_check_sheet_dis_antiparallel_3) // 3rd check of distance between strand i and strand j
									{
										{ // these two strands should be closer
											pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
											pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
											pair_insert_stmt.bind(4,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
											pair_insert_stmt.bind(5,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
										}
										
										{ // these two strands should be closer
											pair_insert_stmt.bind(6,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
											pair_insert_stmt.bind(7,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
											pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
											pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
										}
									}
									
									else
									{
										{ // these two strands should be closer
											pair_insert_stmt.bind(2,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
											pair_insert_stmt.bind(3,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
											pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
											pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
										}

										{ // these two strands should be closer
											pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
											pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
											pair_insert_stmt.bind(8,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
											pair_insert_stmt.bind(9,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
										}
									}


									pair_insert_stmt.bind(10,	val_shortest_dis_sidechain);
									basic::database::safely_write_to_database(pair_insert_stmt);
								}
								
								/*else
								{
									TR.Info << "avg torsion angle is not within range " << endl;
								}*/

							}
							/*else
							{
								TR.Info << "<check sheet (by anti-parallel way)> not right-handed beta sheet" << endl;
							}*/

						}
					}
				}
			}



			else // check by parallel way
			{
			//TR.Info << "<parallel sheet check by distance> begins between " << i << " th pair of strand_pair (i) and " << j << " th strand_pair (j)" << endl;

				StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
				StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
				Real return_of_check_sheet_dis_parallel_1 = check_sheet_dis_parallel (pose, temp_strand_i, temp_strand_j);
				if (return_of_check_sheet_dis_parallel_1 != -99)
				{
					StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
					StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
					Real return_of_check_sheet_dis_parallel_2 = check_sheet_dis_parallel (pose, temp_strand_i, temp_strand_j);
					if (return_of_check_sheet_dis_parallel_2 != -99)
					{
						StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
						StrandFragment temp_strand_j(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
						Real return_of_check_sheet_dis_parallel_3 = check_sheet_dis_parallel (pose, temp_strand_i, temp_strand_j);
						if (return_of_check_sheet_dis_parallel_3 != -99)
						{
							StrandFragment temp_strand_i(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
							StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
							Real return_of_check_sheet_dis_parallel_4 = check_sheet_dis_parallel (pose, temp_strand_i, temp_strand_j);
							if (return_of_check_sheet_dis_parallel_4 != -99)
							{
							//	TR.Info << "<sheet by distance found (by parallel way)> " << ii << " th strand_pair (i) and " << jj << " th strand_pair (j) are within " << min_sheet_dis_ << " Angstrom and " << max_sheet_dis_ << " Angstrom" << endl;

							//	TR.Info << "<check sheet by torsion begins (by parallel way)> between " << ii << " th pair of strand_pair (i) and " << jj << " th strand_pair (j)" << endl;

								StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
								StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
								Real sheet_torsion_1 = sheet_torsion(pose, temp_strand_i, temp_strand_j);

								StrandFragment temp_strand_k(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
								StrandFragment temp_strand_l(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
								Real sheet_torsion_2 = sheet_torsion(pose, temp_strand_k, temp_strand_l);

								StrandFragment temp_strand_m(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
								StrandFragment temp_strand_n(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
								Real sheet_torsion_3 = sheet_torsion(pose, temp_strand_m, temp_strand_n);

								StrandFragment temp_strand_o(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
								StrandFragment temp_strand_p(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());
								Real sheet_torsion_4 = sheet_torsion(pose, temp_strand_o, temp_strand_p);

								if (sheet_torsion_1 < 0 && sheet_torsion_2 < 0 && sheet_torsion_3 < 0 && sheet_torsion_4 < 0)
								{

//									TR.Info << "torsions are : " << sheet_torsion_1 << ", " << sheet_torsion_2 << ", " << sheet_torsion_3 << ", " << sheet_torsion_4 << endl;
									Real sheet_torsion_avg = (sheet_torsion_1 + sheet_torsion_2 + sheet_torsion_3 + sheet_torsion_4)/4;
//									TR.Info << "sheet_torsion_avg (by parallel way) : " << sheet_torsion_avg << endl;

									if (sheet_torsion_avg > min_sheet_torsion_ && sheet_torsion_avg < max_sheet_torsion_)
									{

										//TR.Info << "<sheet by torsion found (by parallel way)> the average torsion between " << ii << " th strand_pair (i) and " << jj << " th strand_pair (j) is within " << min_sheet_torsion_ << " and " << max_sheet_torsion_ << endl;
										StrandFragment temp_strand_i(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
										StrandFragment temp_strand_j(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());

										Real val_shortest_dis_sidechain_1 = shortest_dis_sidechain (pose, temp_strand_i, temp_strand_j);
										//TR.Info << "one shortest_distance_sidechain_1 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_1 << endl;


										StrandFragment temp_strand_k(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());
										StrandFragment temp_strand_l(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());

										Real val_shortest_dis_sidechain_2 = shortest_dis_sidechain (pose, temp_strand_k, temp_strand_l);
										//TR.Info << "one shortest_distance_sidechain_2 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_2 << endl;


										StrandFragment temp_strand_m(i_strand_from_full_strand_pairs[ii].get_start(), i_strand_from_full_strand_pairs[ii].get_end());
										StrandFragment temp_strand_n(i_strand_from_full_strand_pairs[jj].get_start(), i_strand_from_full_strand_pairs[jj].get_end());

										Real val_shortest_dis_sidechain_3 = shortest_dis_sidechain (pose, temp_strand_m, temp_strand_n);
										//TR.Info << "one shortest_distance_sidechain_3 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " <<  val_shortest_dis_sidechain_3 << endl;


										StrandFragment temp_strand_o(j_strand_from_full_strand_pairs[jj].get_start(), j_strand_from_full_strand_pairs[jj].get_end());
										StrandFragment temp_strand_p(j_strand_from_full_strand_pairs[ii].get_start(), j_strand_from_full_strand_pairs[ii].get_end());

										Real val_shortest_dis_sidechain_4 = shortest_dis_sidechain (pose, temp_strand_o, temp_strand_p);
										//TR.Info << "one shortest_distance_sidechain_4 between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain_4 << endl;


										Real val_shortest_dis_sidechain = shortest_dis_sidechain(val_shortest_dis_sidechain_1, val_shortest_dis_sidechain_2, val_shortest_dis_sidechain_3, val_shortest_dis_sidechain_4);
										//TR.Info << "the shortest_distance_sidechain between " <<  ii << " th strand_pair and " << jj << " th strand_pair : " << val_shortest_dis_sidechain << endl;


										// to find i-i and j-i by calculating distances between pairs
										Real val_shortest_dis_pairs = shortest_dis_pairs(return_of_check_sheet_dis_parallel_1, return_of_check_sheet_dis_parallel_2, return_of_check_sheet_dis_parallel_3, return_of_check_sheet_dis_parallel_4);



										TR.Info << "=========== saving (parallel) sheet pairs =========== " << endl;
										string pair_insert =
										"INSERT INTO sheet_pairs (struct_id, g1_strand_1_res_start, g1_strand_1_res_end, g1_strand_2_res_start, g1_strand_2_res_end, g2_strand_1_res_start, g2_strand_1_res_end, g2_strand_2_res_start, g2_strand_2_res_end, shortest_sc_dis)  VALUES (?, ?,?,?,?, ?,?,?,?, ?);";
										statement pair_insert_stmt(basic::database::safely_prepare_statement(pair_insert,	db_session));
										pair_insert_stmt.bind(1,	struct_id);
										
										// I need to store closer strands together as a group for Tim's graph based assembly application
										
										if (val_shortest_dis_pairs == return_of_check_sheet_dis_parallel_1) // first check of distance between strand i and strand j
										{
											{ // these two strands should be closer
												pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
												pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
												pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
												pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand	
											}
											
											{ // these two strands should be closer
												pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
												pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
												pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
												pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
											}
										}
										
										
										else if (val_shortest_dis_pairs == return_of_check_sheet_dis_parallel_2) // 2nd check of distance between strand i and strand j
										{
											{ // these two strands should be closer
												pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
												pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
												pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
												pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand	
											}
											
											{ // these two strands should be closer
												pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
												pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
												pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
												pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
											}
										}
										
										else if (val_shortest_dis_pairs == return_of_check_sheet_dis_parallel_3) // 3rd check of distance between strand i and strand j
										{
											{ // these two strands should be closer
												pair_insert_stmt.bind(2,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
												pair_insert_stmt.bind(3,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
												pair_insert_stmt.bind(4,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
												pair_insert_stmt.bind(5,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
											}
											
											{ // these two strands should be closer
												pair_insert_stmt.bind(6,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
												pair_insert_stmt.bind(7,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
												pair_insert_stmt.bind(8,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
												pair_insert_stmt.bind(9,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
											}
										}
										
										else
										{
											{ // these two strands should be closer
												pair_insert_stmt.bind(2,	j_strand_from_full_strand_pairs[jj].get_start()); // start residue of j-j strand
												pair_insert_stmt.bind(3,	j_strand_from_full_strand_pairs[jj].get_end());   // end   residue of j-j strand
												pair_insert_stmt.bind(4,	j_strand_from_full_strand_pairs[ii].get_start()); // start residue of j-i strand
												pair_insert_stmt.bind(5,	j_strand_from_full_strand_pairs[ii].get_end());   // end   residue of j-i strand
											}
											
											{ // these two strands should be closer
												pair_insert_stmt.bind(6,	i_strand_from_full_strand_pairs[ii].get_start()); // start residue of i-i strand
												pair_insert_stmt.bind(7,	i_strand_from_full_strand_pairs[ii].get_end());   // end   residue of i-i strand
												pair_insert_stmt.bind(8,	i_strand_from_full_strand_pairs[jj].get_start()); // start residue of i-j strand
												pair_insert_stmt.bind(9,	i_strand_from_full_strand_pairs[jj].get_end());   // end   residue of i-j strand
											}
										}
										
										pair_insert_stmt.bind(10,	val_shortest_dis_sidechain);
										basic::database::safely_write_to_database(pair_insert_stmt);
									}
									/*else
									{
									//	TR.Info << "<check sheet (by parallel way)> : avg torsion angle is not within range " << endl;
									}*/
								}
								/*else
								{
								//	TR.Info << "not right-handed beta sheet" << endl;
								}*/
							}
						}
					}
				}
			}
		}
	}

	TR << "============== (Done) saving \"ideal\" pairs of sheets ==========" << endl;
	// sheet pairing ends

	return 0;
} //StrandBundleFeatures::report_features


} //namespace strand_assembly
} //namespace features
} //namespace protocols
