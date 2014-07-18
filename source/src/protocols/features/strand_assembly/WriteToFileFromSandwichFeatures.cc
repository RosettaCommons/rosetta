// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file protocols/features/strand_assembly/WriteToFileFromSW.cc
/// @brief Write to a file after SandwichFeatures
/// @author Doo Nam Kim
/// @overview

//Devel
#include <protocols/features/strand_assembly/WriteToFileFromSandwichFeatures.hh>

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

static basic::Tracer TR("protocols.features.strand_assembly.WriteToFileFromSandwichFeatures");


utility::vector1<Size>
get_vector_loop_AA_distribution (
	StructureID struct_id,
	sessionOP db_session,
	string loop_kind // like 'hairpin' or 'inter-sheet-loop'
	)
{
	string	sum_string =
		"SELECT\n"
		"	sum(A), sum(C), sum(D), sum(E), sum(F), \n"
		"	sum(G), sum(H), sum(I), sum(K), sum(L), \n"
		"	sum(M), sum(N), sum(P), sum(Q), sum(R), \n"
		"	sum(S), sum(T), sum(V), sum(W), sum(Y) \n"
		"FROM\n"
		"	sw_by_components \n"
		"WHERE\n"
		"	loop_kind = ? \n"
		"	AND struct_id = ? ;";

	statement sum_statement(basic::database::safely_prepare_statement(sum_string, db_session));
	sum_statement.bind(1,	loop_kind);
	sum_statement.bind(2,	struct_id);
	result res(basic::database::safely_read_from_database(sum_statement));

	Size num_A, num_C,	num_D,	num_E, num_F, num_G, num_H, num_I, num_K, num_L, num_M, num_N, num_P, num_Q, num_R, num_S, num_T, num_V, num_W, num_Y;

	utility::vector1<Size> vector_loop_AA_distribution;

	while(res.next())
	{
		res >> num_A >>  num_C >> 	num_D >> 	num_E >>  num_F >>  num_G >>  num_H >>  num_I >>  num_K >>  num_L >>  num_M >>  num_N >>  num_P >>  num_Q >>  num_R >>  num_S >>  num_T >>  num_V >>  num_W >>  num_Y;
		vector_loop_AA_distribution.push_back(num_A);
		vector_loop_AA_distribution.push_back(num_C);
		vector_loop_AA_distribution.push_back(num_D);
		vector_loop_AA_distribution.push_back(num_E);
		vector_loop_AA_distribution.push_back(num_F);

		vector_loop_AA_distribution.push_back(num_G);
		vector_loop_AA_distribution.push_back(num_H);
		vector_loop_AA_distribution.push_back(num_I);
		vector_loop_AA_distribution.push_back(num_K);
		vector_loop_AA_distribution.push_back(num_L);

		vector_loop_AA_distribution.push_back(num_M);
		vector_loop_AA_distribution.push_back(num_N);
		vector_loop_AA_distribution.push_back(num_P);
		vector_loop_AA_distribution.push_back(num_Q);
		vector_loop_AA_distribution.push_back(num_R);

		vector_loop_AA_distribution.push_back(num_S);
		vector_loop_AA_distribution.push_back(num_T);
		vector_loop_AA_distribution.push_back(num_V);
		vector_loop_AA_distribution.push_back(num_W);
		vector_loop_AA_distribution.push_back(num_Y);
	}

	return vector_loop_AA_distribution;
} // get_vector_loop_AA_distribution


core::Size
write_AA_distribution_without_direction_to_a_file(
	string	tag,
	StructureID	struct_id,
	utility::sql_database::sessionOP	db_session)
{
	Size tag_len = tag.length();
	string pdb_file_name = tag.substr(0, tag_len-5);
	string AA_dis_file_name = pdb_file_name + "_AA_distribution_of_loops_sorted_alphabetically.txt";
	ofstream AA_dis_file;

	AA_dis_file.open(AA_dis_file_name.c_str());
	utility::vector1<Size> vector_of_hairpin_AA = get_vector_loop_AA_distribution (struct_id,	db_session, "hairpin");
	utility::vector1<Size> vector_of_inter_sheet_loop_AA = get_vector_loop_AA_distribution (struct_id,	db_session, "loop_connecting_two_sheets");

	AA_dis_file << "hairpin_AA	inter_sheet_loop_AA" << endl;
	for (Size i =1; i<=(vector_of_hairpin_AA.size()); i++)
	{
		AA_dis_file << vector_of_hairpin_AA[i] << "	" << vector_of_inter_sheet_loop_AA[i] << endl;
	}
	AA_dis_file.close();

	return 0;
}


} //namespace strand_assembly
} //namespace features
} //namespace protocols
