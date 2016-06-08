// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/antibody/grafting/grafter.cc
/// @brief Grafter implementation: take SCS results and create antibody Pose with CDR loops from results
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/grafter.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>

#include <basic/Tracer.hh>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

using core::Size;

using namespace core::pose;


static THREAD_LOCAL basic::Tracer TR("protocols.antibody.grafting");



int find_chain(Pose const &pose, char pdb_chain_letter, string const &template_name)
{
	int chain = -1;
	// Finding right chain in SCS template and delete eveything else
	for(uint i=1; i<=pose.total_residue(); ++i) {
		if( pose.pdb_info()->chain(i) == pdb_chain_letter ) { chain=pose.chain(i);  break; }
	}

	if( chain >= 0 ) return chain;
	else throw _AE_grafting_failed_( string("Could not find chain: ") + pdb_chain_letter + " in template " + template_name );
}


struct PDB_N
{
	PDB_N(string s) {
		n = std::stoi(s);

		if( s.size()  and std::isalpha( s[s.size()-1])) icode=s[s.size()-1];
		else icode=' ';
	}

	int n;
	char icode;
};



/// @brief Construct Pose SCS_ResultSet as templates and superimpose 'orientation' template on results. Write results into specified output_prefix.
core::pose::PoseOP construct_antibody(AntibodySequence const &A, SCS_ResultSet const &scs, string const & prefix, string const & suffix, string const & database)
{
	string frh_pdb_name = database + "/antibody_database/pdb" + scs.frh->pdb + "_chothia.pdb";
	string frl_pdb_name = database + "/antibody_database/pdb" + scs.frl->pdb + "_chothia.pdb";
	string orientation_pdb_name = database + "/antibody_database/pdb" + scs.orientation->pdb + "_chothia.pdb";

	PoseOP frh = core::import_pose::pose_from_file( frh_pdb_name, core::import_pose::PDB_file);
	PoseOP frl = core::import_pose::pose_from_file( frl_pdb_name, core::import_pose::PDB_file);

	frh->dump_pdb(prefix + "frh" + suffix + ".pdb");
	frl->dump_pdb(prefix + "frl" + suffix + ".pdb");

	core::pose::PoseOP orientation = core::import_pose::pose_from_file( orientation_pdb_name , core::import_pose::PDB_file);
	orientation->dump_pdb(prefix + "orientation" + suffix + ".pdb");

	AntibodyFramework trimmed_heavy_fr = A.heavy_framework();
	AntibodyFramework trimmed_light_fr = A.light_framework();

	trim_framework(A, trimmed_heavy_fr, trimmed_light_fr);

	AntibodyNumbering an( Chothia_Numberer().number(A, trimmed_heavy_fr, trimmed_light_fr) );

	struct {
		char chain;
		core::pose::PoseOP &pose;
		AntibodyChainNumbering numbering;
		string trimmed_sequence;
		string color;
	} J[] {
		{'H', frh, an.heavy, trimmed_heavy_fr.fr1 + A.h1_sequence() + trimmed_heavy_fr.fr2 + A.h2_sequence() + trimmed_heavy_fr.fr3 + A.h3_sequence() + trimmed_heavy_fr.fr4, TR.Blue},
		{'L', frl, an.light, trimmed_light_fr.fr1 + A.l1_sequence() + trimmed_light_fr.fr2 + A.l2_sequence() + trimmed_light_fr.fr3 + A.l3_sequence() + trimmed_light_fr.fr4, TR.Green},
	};

	for(auto &j : J) {
		char chain_lower = char(std::tolower(j.chain));
		TR << "Adjusting fr" << chain_lower << " template sequence [" << j.pose->pdb_info()->name() << "]..." << std::endl;

		AntibodyChainNumbering::NumberingVector numbering = j.numbering.all();

		TR << "By using numbering: " << utility::join( numbering, " ") << std::endl;

		TR.Debug << "Sequence before (all chains): " << j.pose->sequence() << std::endl;

		j.pose = j.pose->split_by_chain( find_chain(*j.pose, j.chain, j.pose->pdb_info()->name() ) );

		TR << "Sequence before: " << j.color << j.pose->sequence() << TR.Reset << std::endl;

		for(uint i=j.pose->total_residue()-1; i>=1; --i) {
			string pdb_res_n = std::to_string(j.pose->pdb_info()->number(i)) + ( j.pose->pdb_info()->icode(i) == ' ' ?  "" : string(1, j.pose->pdb_info()->icode(i) ) );
			auto np = std::find(numbering.begin(), numbering.end(), pdb_res_n);

			if( np < numbering.end() ) {
				char aa = j.trimmed_sequence[ np-numbering.begin() ];
				TR.Trace << "Replacing pose residue " << i << " with " << aa << std::endl;
				protocols::simple_moves::MutateResidue(i, aa).apply( *j.pose.get() );

				j.trimmed_sequence.erase( j.trimmed_sequence.begin() + ( np - numbering.begin() ) );
				numbering.erase( np );
			}
		}
		if( numbering.size() ) {
			TR << TR.Red << "WARNING: Was not able to adjust all residue in chain: " << j.chain << "!!!" << TR.Reset << std::endl;
			TR << TR.Red << "Leftovers numbering: " << numbering  << "  AA: " << j.trimmed_sequence << TR.Reset << std::endl;
			TR << TR.Red << TR.Underline << "Original numbering was:" << TR.Reset << ' ' << j.numbering << std::endl;
		}

		TR << "Sequence after:  " << TR.Bold << j.color << j.pose->sequence() << TR.Reset << std::endl;

		j.pose->dump_pdb( string(prefix + "fr") + chain_lower + "_after_seqeunce_adjustment" + suffix + ".pdb" );

		PoseOP O = orientation->split_by_chain( find_chain(*orientation, j.chain, "orientation" ) );
		
		protocols::moves::MoverOP imposer( new protocols::simple_moves::SuperimposeMover( *O,
																																										  1 /*ref_start*/,
																																										  std::min( O->n_residue(), j.pose->n_residue() ) /*ref_end*/,
																																										  1 /*target_start*/,
																																										  std::min(O->n_residue(), j.pose->n_residue() ) /*target_end*/,
																																										  true /*CA_only*/) );
		imposer->apply(*j.pose);
	}

	frh->append_pose_by_jump(*frl, 1);
	frh->pdb_info()->obsolete(false);

	frh->dump_pdb(prefix + "frh_frl_oriented" + suffix + ".pdb");

	// for(uint i=1; i<=frh->total_residue(); ++i) {
	// 	TR << "New pose pdb info for res " << i << ":" << frh->pdb_info()->pose2pdb(i) << " i:" << frh->pdb_info()->icode(i) << std::endl;
	// }

	return frh;
}


/// @brief graft cdr-loops using best scs-results and write results into specified output_prefix
core::pose::PoseOP graft_cdr_loops(AntibodySequence const &A, SCS_ResultSet const &scs, string const & prefix, string const & suffix, string const & database)
{
	if ( !(scs.h1 and scs.h2 and scs.h3 and scs.l1 and scs.l2 and scs.l2 and scs.l3 and scs.frh and scs.frl and scs.orientation) ) throw _AE_grafting_failed_("SimpleGrafter::graft: not all nessesary SCS results is specified!");

	PoseOP result = construct_antibody(A, scs, prefix, suffix, database);

	AntibodyFramework trimmed_heavy_fr = A.heavy_framework();
	AntibodyFramework trimmed_light_fr = A.light_framework();

	trim_framework(A, trimmed_heavy_fr, trimmed_light_fr);

	AntibodyNumbering an( Chothia_Numberer().number(A, trimmed_heavy_fr, trimmed_light_fr) );

	struct {
		string name; char chain; string pdb;
		AntibodyChainNumbering::NumberingVector cdr_numbering;  // Numbering for framework regions just befor and after CDR
	} G[] {
		{"h1", 'H', scs.h1->pdb,  an.heavy.cdr1}, {"h2", 'H', scs.h2->pdb,  an.heavy.cdr2}, {"h3", 'H', scs.h3->pdb,  an.heavy.cdr3},
		{"l1", 'L', scs.l1->pdb,  an.light.cdr1}, {"l2", 'L', scs.l2->pdb,  an.light.cdr2}, {"l3", 'L', scs.l3->pdb,  an.light.cdr3},
	};

	for(auto &g : G) {
		
		TR << "Attaching CDR loop: " << TR.Bold << g.name << ", from pdb: " << g.pdb << std::endl;

		string pdb_name = database + "/antibody_database/pdb" + g.pdb + "_chothia.pdb";
		core::pose::PoseOP cdr = core::import_pose::pose_from_file(pdb_name, core::import_pose::PDB_file);

		const int overlap = 2;

		if( !g.cdr_numbering.size() ) throw _AE_grafting_failed_( string("Empty template:") + g.pdb +" supplied as cdr for region:" + g.name );

		PDB_N pdb_n_cdr_first( g.cdr_numbering.front() );
		PDB_N pdb_n_cdr_last(  g.cdr_numbering.back()  );

		Size pose_n_cdr_first = cdr->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_first.n, pdb_n_cdr_first.icode);
		Size pose_n_cdr_last  = cdr->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_last .n, pdb_n_cdr_last .icode);

		if( !pose_n_cdr_first ) throw _AE_grafting_failed_( string("Could not find residue:") + g.cdr_numbering.front() + " in template:" + g.pdb + " region:"+ g.name);
		if( !pose_n_cdr_last )  throw _AE_grafting_failed_( string("Could not find residue:") + g.cdr_numbering.back() + " in template:" + g.pdb  + " region:"+ g.name);

		pose_n_cdr_first -= overlap;
		pose_n_cdr_last += overlap;

		if( pose_n_cdr_first < 1  or  pose_n_cdr_last > cdr->n_residue() ) throw _AE_grafting_failed_( string("There is not enough overlap residue at:")+ g.pdb + " in template:" + g.pdb + " region:" + g.name);

		//TR << "Deleting residues: " << pose_n_cdr_last+1 << ":" << cdr->n_residue() << " from cdr template... [template size: " << cdr->n_residue() << "]" << std::endl;
		if( cdr->n_residue() > pose_n_cdr_last ) cdr->delete_residue_range_slow(pose_n_cdr_last+1, cdr->n_residue());
		//TR << "Deleting residues: " << "1:" << pose_n_cdr_first-1 << " from cdr template... [template size: " << cdr->n_residue() << "]" << std::endl;
		if( pose_n_cdr_first > 1 ) cdr->delete_residue_range_slow(1, pose_n_cdr_first-1);

		Size result_pose_cdr_first = result->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_first.n, pdb_n_cdr_first.icode);
		Size result_pose_cdr_last  = result->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_last .n, pdb_n_cdr_last .icode);

		if( !result_pose_cdr_first ) throw _AE_grafting_failed_( string("Could not find residue:") + g.cdr_numbering.front() + " in superimposed pdb. Region:" + g.name);
		if( !result_pose_cdr_last )  throw _AE_grafting_failed_( string("Could not find residue:") + g.cdr_numbering.back() +  " in superimposed pdb. Region:" + g.name);

		result_pose_cdr_first -= overlap;
		result_pose_cdr_last += overlap;

		if( result_pose_cdr_first < 1  or result_pose_cdr_last > result->n_residue() ) throw _AE_grafting_failed_( string("There is not enough overlap residue in superimposed template! region:") + g.name);

		TR << "Grafting..." << std::endl;

		// Reasons for +1/-1: Jared's code uses an "insert pose into pose" mover written by Steven Lewis where the "start"
		// residue is the residue before the insertion region, and the "end" residue is the residue after the insertion region;
		// he continues this convention through his code. Since we're inserting a whole CDR, this puts the required start/end
		// values one residue outside of the loop on either end. Since we define our start/end points as the first/last residues
		// of the loop, we have to subtract/add one to match up with the terminology.
		protocols::grafting::CCDEndsGraftMoverOP grafter(new protocols::grafting::CCDEndsGraftMover(result_pose_cdr_first+1, result_pose_cdr_last-1, *cdr, overlap, overlap, true) );
		grafter->stop_at_closure(true);  grafter->set_cycles(128);

		grafter->apply(*result);
	}

	// prior to dumping, restore proper sequence to CDRs as grafter copies over both structure and sequence from the template
	AntibodyInfoOP ab_info = AntibodyInfoOP( new AntibodyInfo(*result) );
	
	struct{
		string cdr_name; Size cdr_start; Size cdr_end; string cdr_seq;
	} H[] {
		{ "h1", ab_info->get_CDR_start(h1, *result), ab_info->get_CDR_end(h1, *result),  A.h1_sequence() },
		{ "h2", ab_info->get_CDR_start(h2, *result), ab_info->get_CDR_end(h2, *result),  A.h2_sequence() },
		{ "h3", ab_info->get_CDR_start(h3, *result), ab_info->get_CDR_end(h3, *result),  A.h3_sequence() },
		{ "l1", ab_info->get_CDR_start(l1, *result), ab_info->get_CDR_end(l1, *result),  A.l1_sequence() },
		{ "l2", ab_info->get_CDR_start(l2, *result), ab_info->get_CDR_end(l2, *result),  A.l2_sequence() },
		{ "l3", ab_info->get_CDR_start(l3, *result), ab_info->get_CDR_end(l3, *result),  A.l3_sequence() },
	};
	
	for (auto &h : H) {
		// check for matching lengths of cdr and cdr sequence
		// everything should be kosher since we're using the chothia definition throughout (AFAIK)
		if ( h.cdr_seq.size() != h.cdr_end - h.cdr_start + 1 ) throw _AE_grafting_failed_( string("Could revert sequence to query after grafting cdr ") + h.cdr_name + ". Length mismatch between CDR in grafted model (" + utility::to_string(h.cdr_end - h.cdr_start + 1) +  ") and query sequence (" + utility::to_string(h.cdr_seq.size()) + ").");

		
		for(Size i = h.cdr_start; i < h.cdr_end + 1; ++i) {
			protocols::simple_moves::MutateResidue( i, h.cdr_seq[i-h.cdr_start] ).apply( *result );
		}
		
	}

	result->dump_pdb(prefix + "model" + suffix + ".pdb");

	return result;
}



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
