// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/grafter.cc
/// @brief Grafter implementation: take SCS results and create antibody Pose with CDR loops from results
/// @author Sergey Lyskov
/// @author Jared Adolf-Bryfogle
/// @author Jeliazko Jeliazkov

#include <protocols/antibody/grafting/util.hh>

#include <protocols/antibody/grafting/chothia_numberer.hh> // AUTO IWYU For AntibodyChainNumbering::NumberingVector
#include <protocols/antibody/grafting/scs_blast.hh> // AUTO IWYU For SCS_ResultSet, SCS_Result, trim_fram...

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/grafter.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>

#include <core/import_pose/import_pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/init_id_map.hh>
#include <core/scoring/rms_util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
#include <protocols/antibody/AntibodyCDRGrafter.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>

namespace protocols {
namespace antibody {
namespace grafting {

using std::string;

using core::Size;

using namespace core::pose;

static basic::Tracer TR("protocols.antibody.grafting");



int find_chain(Pose const &pose, std::string const & pdb_chain_letter, string const &template_name)
{
	int chain = -1;
	// Finding right chain in SCS template and delete eveything else
	for(uint i=1; i<=pose.size(); ++i) {
		if( pose.pdb_info()->chain(i) == pdb_chain_letter ) { chain=pose.chain(i);  break; }
	}

	if( chain >= 0 ) return chain;
	else throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could not find chain: ") + pdb_chain_letter + " in template " + template_name );
}


struct PDB_N
{
	PDB_N(string s) {
		n = std::stoi(s);
		std::locale loc;
		if( s.size()  and std::isalpha( s[s.size()-1], loc)) icode=s[s.size()-1];
		else icode=' ';
	}

	int n;
	char icode;
};



/// @brief Construct Pose SCS_ResultSet as templates and superimpose 'orientation' template on results. Write results into specified output_prefix.
core::pose::PoseOP construct_antibody(AntibodySequence const &A, SCS_ResultSet const &scs, string const & prefix, string const & suffix, string const & database)
{

	string frh_pdb_name = database + "/antibody_database/" + scs.frh->pdb + "_trunc.pdb";
	string frl_pdb_name = database + "/antibody_database/" + scs.frl->pdb + "_trunc.pdb";
	string orientation_pdb_name = database + "/antibody_database/" + scs.orientation->pdb + "_trunc.pdb";

	PoseOP frh = core::import_pose::pose_from_file( frh_pdb_name, core::import_pose::PDB_file);
	PoseOP frl = core::import_pose::pose_from_file( frl_pdb_name, core::import_pose::PDB_file);

	frh->dump_pdb(prefix + "frh" + suffix + ".pdb");
	frl->dump_pdb(prefix + "frl" + suffix + ".pdb");

	core::pose::PoseOP orientation = core::import_pose::pose_from_file( orientation_pdb_name , core::import_pose::PDB_file);
	orientation->dump_pdb(prefix + "orientation" + suffix + ".pdb");

	AntibodyFramework trimmed_heavy_fr = A.heavy_framework();
	AntibodyFramework trimmed_light_fr = A.light_framework();

	trim_framework(A, trimmed_heavy_fr, trimmed_light_fr); // what happens if no light?

	AntibodyNumbering an( Chothia_Numberer().number(A, trimmed_heavy_fr, trimmed_light_fr) );

    utility::vector1< core::Size > conserved_frh_residues = get_conserved_residue_list('H');
    utility::vector1< core::Size > conserved_frl_residues = get_conserved_residue_list('L');

	TR.Debug << "Conserved FRH regions: " << conserved_frh_residues << std::endl;
	TR.Debug << "Conserved FRL regions: " << conserved_frl_residues << std::endl;

	struct {
		std::string chain;
		core::pose::PoseOP &pose;
		AntibodyChainNumbering numbering;
		string trimmed_sequence;
		utility::CSI::CSI_Enum color;
		utility::vector1< core::Size > conserved_fr_residues; // chothia numbering is assumed (and true for members of the database)
	} J[] {
		{"H", frh, an.heavy, trimmed_heavy_fr.fr1 + A.h1_sequence() + trimmed_heavy_fr.fr2 + A.h2_sequence() + trimmed_heavy_fr.fr3 + A.h3_sequence() + trimmed_heavy_fr.fr4, TR.Blue, conserved_frh_residues},
		{"L", frl, an.light, trimmed_light_fr.fr1 + A.l1_sequence() + trimmed_light_fr.fr2 + A.l2_sequence() + trimmed_light_fr.fr3 + A.l3_sequence() + trimmed_light_fr.fr4, TR.Green, conserved_frl_residues},
	};

	std::locale loc;
	for(auto &j : J) {
		auto chain_lower = char(std::tolower(j.chain[0],loc));

		TR << "Adjusting fr" << chain_lower << " template sequence [" << j.pose->pdb_info()->name() << "]..." << std::endl;

		AntibodyChainNumbering::NumberingVector numbering = j.numbering.all();

		TR << "By using numbering: " << utility::join( numbering, " ") << std::endl;

		TR.Debug << "Sequence before (all chains): " << j.pose->sequence() << std::endl;

		j.pose = j.pose->split_by_chain( find_chain(*j.pose, j.chain, j.pose->pdb_info()->name() ) );

		TR << "Sequence before: " << j.color << j.pose->sequence() << TR.Reset << std::endl;

		for(uint i=j.pose->size()-1; i>=1; --i) {
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
			TR.Warning << TR.Red << "Was not able to adjust all residue in chain: " << j.chain << "!!!" << TR.Reset << std::endl;
			TR.Warning << TR.Red << "Leftovers numbering: " << numbering  << "  AA: " << j.trimmed_sequence << TR.Reset << std::endl;
			TR.Warning << TR.Red << TR.Underline << "Original numbering was:" << TR.Reset << ' ' << j.numbering << std::endl;
		}

		TR << "Sequence after:  " << TR.Bold << j.color << j.pose->sequence() << TR.Reset << std::endl;

		j.pose->dump_pdb( string(prefix + "fr") + chain_lower + "_after_seqeunce_adjustment" + suffix + ".pdb" );

		PoseOP O = orientation->split_by_chain( find_chain(*orientation, j.chain, "orientation" ) );

		// super impose only the FR regions
		core::id::AtomID_Map< core::id::AtomID> atom_map; // map j.pose CAs of FR to orientation CAs of FR
		core::pose::initialize_atomid_map( atom_map, *j.pose, core::id::AtomID::BOGUS_ATOM_ID() );


		for (auto it = j.conserved_fr_residues.begin(); it != j.conserved_fr_residues.end(); ++it) {

			// convert to pose numbering, luckily no insertion codes to worry about!
			core::Size it1 = j.pose->pdb_info()->pdb2pose(j.chain ,*it);
			core::Size it2 = O->pdb_info()->pdb2pose(j.chain ,*it);

			core::id::AtomID const id1( j.pose->residue(it1).atom_index("CA"), it1);
			core::id::AtomID const id2( O->residue(it2).atom_index("CA"), it2);
			atom_map[ id1 ] = id2;
		}

		core::scoring::superimpose_pose( *j.pose, *O, atom_map );

		//j.pose->dump_pdb( prefix + j.chain + "_super_imposed" + suffix + ".pdb");
	}

	frh->append_pose_by_jump(*frl, 1);
	frh->pdb_info()->obsolete(false);

	frh->dump_pdb(prefix + "frh_frl_oriented" + suffix + ".pdb");

	// for(uint i=1; i<=frh->size(); ++i) {
	// 	TR << "New pose pdb info for res " << i << ":" << frh->pdb_info()->pose2pdb(i) << " i:" << frh->pdb_info()->icode(i) << std::endl;
	// }

	return frh;
}


/// @brief graft cdr-loops using best scs-results and write results into specified output_prefix
core::pose::PoseOP graft_cdr_loops(AntibodySequence const &A, SCS_ResultSet const &scs, string const & prefix, string const & suffix, string const & database, bool optimal_graft /* false */, bool optimize_cdrs /* false */)
{
	if ( !(scs.h1 and scs.h2 and scs.h3 and scs.l1 and scs.l2 and scs.l2 and scs.l3 and scs.frh and scs.frl and scs.orientation) ) throw CREATE_EXCEPTION(_AE_grafting_failed_, "SimpleGrafter::graft: not all nessesary SCS results is specified!");

	PoseOP result = construct_antibody(A, scs, prefix, suffix, database);

	AntibodyFramework trimmed_heavy_fr = A.heavy_framework();
	AntibodyFramework trimmed_light_fr = A.light_framework();

	trim_framework(A, trimmed_heavy_fr, trimmed_light_fr);

	AntibodyNumbering an( Chothia_Numberer().number(A, trimmed_heavy_fr, trimmed_light_fr) );

	struct {
		string name; std::string chain; string pdb;
		AntibodyChainNumbering::NumberingVector cdr_numbering;  // Numbering for framework regions just befor and after CDR
		core::Size length; // length of CDR to be inserted
	} G[] {
		{"h1", "H", scs.h1->pdb,  an.heavy.cdr1, A.h1_sequence().size()},
		{"h2", "H", scs.h2->pdb,  an.heavy.cdr2, A.h2_sequence().size()},
		{"h3", "H", scs.h3->pdb,  an.heavy.cdr3, A.h3_sequence().size()},
		{"l1", "L", scs.l1->pdb,  an.light.cdr1, A.l1_sequence().size()},
		{"l2", "L", scs.l2->pdb,  an.light.cdr2, A.l2_sequence().size()},
		{"l3", "L", scs.l3->pdb,  an.light.cdr3, A.l3_sequence().size()},
	};

	AntibodyEnumManager enum_manager = AntibodyEnumManager();

	//result->dump_pdb("result_pre_graft.pdb");

	for(auto &g : G) {

		TR << "Attaching CDR loop: " << TR.Bold << g.name << ", from pdb: " << g.pdb << std::endl;

		string pdb_name = database + "/antibody_database/" + g.pdb + "_trunc.pdb";
		core::pose::PoseOP cdr = core::import_pose::pose_from_file(pdb_name, core::import_pose::PDB_file);

		if( !g.cdr_numbering.size() ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Empty template:") + g.pdb +" supplied as cdr for region:" + g.name );

		PDB_N pdb_n_cdr_first( g.cdr_numbering.front() );
		PDB_N pdb_n_cdr_last(  g.cdr_numbering.back()  );

		core::Size pose_n_cdr_first = cdr->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_first.n, pdb_n_cdr_first.icode);
		core::Size pose_n_cdr_last  = cdr->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_last .n, pdb_n_cdr_last .icode);

		if( !pose_n_cdr_first ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could not find residue:") + g.cdr_numbering.front() + " in template:" + g.pdb + " region:"+ g.name);
		if( !pose_n_cdr_last )  throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could not find residue:") + g.cdr_numbering.back() + " in template:" + g.pdb  + " region:"+ g.name);

		if (optimal_graft){

			//AntibodyInfoOP scaffold_ab_info = utility::pointer::make_shared< AntibodyInfo >( *result ); //Only needs AbInfo for on-the-fly numbering conversion.
			AntibodyCDRGrafter grafter = AntibodyCDRGrafter();
			CDRNameEnum cdr_to_graft = enum_manager.cdr_name_string_to_enum(g.name);
			grafter.set_cdr_only(cdr_to_graft);
			grafter.set_donor_structure( *cdr );
			grafter.set_use_secondary_graft_mover_if_needed( true );
			grafter.set_optimize_cdrs( false ); //We will do this later.
			grafter.set_dihedral_constraint_weight( .3 ); //Value Used originally for ab design
			grafter.set_idealize_insert( false ); //Take bond angles and lengths from original grafted structure.  If you want to optimize these, use cart-min at the end of RosettaAntibody.

			std::cout << "Grafting "<< g.name << " Using mover" << std::endl;
			grafter.set_stop_after_closure(false);
			grafter.apply(*result);


		}
		else {

			const core::Size overlap = 2; // used for alignment of scaffold/loop
			core::Size insert_flexibility   = 2; // flexible regions in loop
			core::Size scaffold_flexibility = 2; // flexible regions in scaffold

			// test cdr loop length before inserting and update flexibility accordingly
			if (g.length == 4) {
				// decrease flexibility in the loop, so there is no overlap
				insert_flexibility   = 1;
				// increase flexibility in the scaffold to compensate
				scaffold_flexibility = 3;
			}

			pose_n_cdr_first -= overlap;
			pose_n_cdr_last += overlap;

			if( pose_n_cdr_first < 1  or  pose_n_cdr_last > cdr->total_residue() ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("There is not enough overlap residue at:")+ g.pdb + " in template:" + g.pdb + " region:" + g.name);

			//TR << "Deleting residues: " << pose_n_cdr_last+1 << ":" << cdr->total_residue() << " from cdr template... [template size: " << cdr->total_residue() << "]" << std::endl;
			if( cdr->total_residue() > pose_n_cdr_last ) cdr->delete_residue_range_slow(pose_n_cdr_last+1, cdr->total_residue());
			//TR << "Deleting residues: " << "1:" << pose_n_cdr_first-1 << " from cdr template... [template size: " << cdr->total_residue() << "]" << std::endl;
			if( pose_n_cdr_first > 1 ) cdr->delete_residue_range_slow(1, pose_n_cdr_first-1);

            // we must erase disulfide memory of the input as this break CDR grafting later on...
            cdr->conformation().detect_disulfides();

			core::Size result_pose_cdr_first = result->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_first.n, pdb_n_cdr_first.icode);
			core::Size result_pose_cdr_last  = result->pdb_info()->pdb2pose(g.chain, pdb_n_cdr_last .n, pdb_n_cdr_last .icode);

			if( !result_pose_cdr_first ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could not find residue:") + g.cdr_numbering.front() + " in superimposed pdb. Region:" + g.name);
			if( !result_pose_cdr_last )  throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could not find residue:") + g.cdr_numbering.back() +  " in superimposed pdb. Region:" + g.name);

			result_pose_cdr_first -= overlap;

			result_pose_cdr_last += overlap;

			if( result_pose_cdr_first < 1  or result_pose_cdr_last > result->total_residue() ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("There is not enough overlap residue in superimposed template! region:") + g.name);

			TR << "Grafting..." << std::endl;

			// Reasons for +1/-1: Jared's code uses an "insert pose into pose" mover written by Steven Lewis where the "start"
			// residue is the residue before the insertion region, and the "end" residue is the residue after the insertion region;
			// he continues this convention through his code. Since we're inserting a whole CDR, this puts the required start/end
			// values one residue outside of the loop on either end. Since we define our start/end points as the first/last residues
			// of the loop, we have to subtract/add one to match up with the terminology.
			protocols::grafting::CCDEndsGraftMoverOP grafter(new protocols::grafting::CCDEndsGraftMover(result_pose_cdr_first+1, result_pose_cdr_last-1, *cdr, overlap, overlap, true) );

			//JAB - increasing flexibility of stem of insert from 0 to 2 to fix disulfide issues in L1.
			//  This needs to be refactored to use both CCDEndsGraftMover and AnchoredGraftMover if if the loop is not closed,
			//   then a dihedral-constrained minimization of the CDR loops with at least a few residues into the stem.
			//   This is a drastic change, that I don't think is ready to be done now...

			grafter->set_cycles(128);
			grafter->set_scaffold_flexibility(scaffold_flexibility, scaffold_flexibility);
			grafter->set_insert_flexibility(insert_flexibility, insert_flexibility);
			grafter->stop_at_closure( false );
			grafter->copy_pdbinfo( true );
			grafter->apply(*result);

		}

		result->dump_pdb(prefix + "debug-model-A" + suffix + ".pdb");
	}


	// prior to dumping, restore proper sequence to CDRs as grafter copies over both structure and sequence from the template
	AntibodyInfoOP ab_info = utility::pointer::make_shared< AntibodyInfo >(*result);

	struct{
		string cdr_name; core::Size cdr_start; core::Size cdr_end; string cdr_seq;
	} H[] {
		{ "h1", ab_info->get_CDR_start(h1, *result), ab_info->get_CDR_end(h1, *result),  A.h1_sequence() },
		{ "h2", ab_info->get_CDR_start(h2, *result), ab_info->get_CDR_end(h2, *result),  A.h2_sequence() },
		{ "h3", ab_info->get_CDR_start(h3, *result), ab_info->get_CDR_end(h3, *result),  A.h3_sequence() },
		{ "l1", ab_info->get_CDR_start(l1, *result), ab_info->get_CDR_end(l1, *result),  A.l1_sequence() },
		{ "l2", ab_info->get_CDR_start(l2, *result), ab_info->get_CDR_end(l2, *result),  A.l2_sequence() },
		{ "l3", ab_info->get_CDR_start(l3, *result), ab_info->get_CDR_end(l3, *result),  A.l3_sequence() },
	};

	result->dump_pdb(prefix + "debug-model-B" + suffix + ".pdb");

	for (auto &h : H) {
		// check for matching lengths of cdr and cdr sequence
		// everything should be kosher since we're using the chothia definition throughout (AFAIK)
		if ( h.cdr_seq.size() != h.cdr_end - h.cdr_start + 1 ) throw CREATE_EXCEPTION(_AE_grafting_failed_,  string("Could revert sequence to query after grafting cdr ") + h.cdr_name + ". Length mismatch between CDR in grafted model (" + utility::to_string(h.cdr_end - h.cdr_start + 1) +  ") and query sequence (" + utility::to_string(h.cdr_seq.size()) + ").");


		for(core::Size i = h.cdr_start; i < h.cdr_end + 1; ++i) {
			protocols::simple_moves::MutateResidue( i, h.cdr_seq[i-h.cdr_start] ).apply( *result );
		}

	}

	result->dump_pdb(prefix + "debug-model-B" + suffix + ".pdb");

	if (optimal_graft || optimize_cdrs){
		//Optmize all CDRs and DE loop to incorporate the new CDRs using dihedral constraints (including stem residues into the framework).

		// Default ab_design scorefunction has dihedral constraints on at a decent weight.
		// Default relax is to relax CDRs and optimize neighbor side-chains.  The neighbors are updated every min step,
		//  This is exactly what we use for antibody design.  By default, no design is done.

		result->dump_pdb(prefix + "pre-model" + suffix + ".pdb");
		TR << "Optimizing post-grafted CDRs including CDR4" << std::endl;
		//AntibodyInfoOP ab_info = utility::pointer::make_shared< AntibodyInfo >(*result);

		constraints::CDRDihedralConstraintMover cst_mover = constraints::CDRDihedralConstraintMover(ab_info);
		cst_mover.set_use_cluster_csts( false );

		for (core::Size i = 1; i <= CDRNameEnum_proto_total; ++i){
			auto cdr = static_cast<CDRNameEnum>(i);
			cst_mover.set_cdr(cdr);
			cst_mover.apply(*result);
		}

		//Relax the CDRs and DE loop.
		design::GeneralAntibodyModeler modeler = design::GeneralAntibodyModeler(ab_info);
		modeler.set_cdr_range( CDRNameEnum_start, CDRNameEnum_proto_total, true );
		modeler.set_overhang( 2 ); //3 Residue overhang to account for 3 residue flexibility in CCDEndsGraftMover.
		modeler.relax_cdrs( * result );

		//Remove our added constraints.
		result->remove_constraints(); //Since we just created the antibody, this should be OK here.

	}

	result->dump_pdb(prefix + "model" + suffix + ".pdb");

	return result;
}



} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // __ANTIBODY_GRAFTING__
