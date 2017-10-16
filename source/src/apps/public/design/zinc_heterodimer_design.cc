// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/design/zinc_heterodimer_design.cc
/// @brief This is the application for a metal interface design project.  In its initial conception, the idea was to design an interface between ankyrin and ubc12.  First, RosettaMatch is used to design a histidine/cystine zinc binding site on ankyrin (3 residues).  The fourth residue to coordinate a tetrahedral zinc then comes from ubc12 (present natively).  This metal binding gets the interface going; the protocol searches rigid-body space to try to find a shape-complementary interaction and then designs the interface.
/// @author Steven Lewis, Ramesh Jha, Bryan Der

// Unit Headers
#include <devel/init.hh>
#include <protocols/metal_interface/FindClosestAtom.hh>
#include <protocols/metal_interface/ZincHeterodimerMover.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh> //for changing the HIS tautomer
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD

#include <protocols/jd2/JobDistributor.hh>
#include <utility/file/FileName.hh>

// Numeric Headers
#include <numeric/xyzVector.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/variant_util.hh>


//tracers
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.rjha" );

//local options
namespace local {
basic::options::FileOptionKey const partner1("partner1");
basic::options::FileOptionKey const partner2("partner2");
basic::options::FileOptionKey const match_pdb("match_pdb");
basic::options::IntegerOptionKey const partner2_residue("partner2_residue");
basic::options::BooleanOptionKey const skip_sitegraft_repack("skip_sitegraft_repack");
}//local


//local enum - this names the (hopefully 5) residues in the input from RosettaMatch
enum MatchPosition {
	p1_1 = 1, //partner 1, residue 1
	p1_2 = 2, //partner 1, residue 2
	p1_3 = 3, //partner 1, residue 3
	p2 = 4,   //residue on partner2 - replaces partner2_residue in that pose
	metal = 5
};

int
main( int argc, char* argv[] )
{
	try {

		//add extra options to option system
		///.def call sets the defaults
		using basic::options::option;
		option.add( local::partner1, "partner 1 for metal interface design (ankyrin)" ).def("2D4X.pdb");
		option.add( local::partner2, "partner 2 for metal interface design (ubc12)" ).def("1UBQ.pdb");
		option.add( local::match_pdb, "RosettaMatch for metal interface design" ).def("match.pdb");
		option.add( local::partner2_residue, "partner 2 residue to replace with the 4th residue in the match pdb").def(68);
		option.add( local::skip_sitegraft_repack, "skip repacking partner1 after grafting metal site (faster for debugging purposes)").def(false);

		//initialize options
		devel::init(argc, argv);
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].user() ) {
			utility_exit_with_message("do not use -s with zinc_heterodimer_design (program uses internally)");
		}

		//basic::options::option[ basic::options::OptionKeys::run::skip_set_reasonable_fold_tree ].value(true);
		//read in our starting structures
		core::pose::Pose match;
		core::import_pose::pose_from_file( match, basic::options::option[local::match_pdb].value() , core::import_pose::PDB_file);
		core::Size const matchlength = match.size();
		debug_assert(matchlength == /*MatchPosition*/metal); //5
		//basic::options::option[ basic::options::OptionKeys::run::skip_set_reasonable_fold_tree ].value(false);

		//debugging
		//match.dump_pdb("match_hydrogened.pdb");

		core::pose::Pose partner2;
		core::import_pose::pose_from_file( partner2, basic::options::option[local::partner2].value() , core::import_pose::PDB_file);
		core::Size const partner2length = partner2.size();

		core::pose::Pose partner1;
		core::import_pose::pose_from_file( partner1, basic::options::option[local::partner1].value() , core::import_pose::PDB_file);
		core::Size const partner1length = partner1.size();


		//establish these ahead of time for convenience
		core::chemical::ResidueTypeSetCOP typeset(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
		core::chemical::ResidueType const & CYZ(typeset->name_map("CYZ"));


		////////////////////////////////////////////////////////////////////////
		//replace ankyrin (partner1) residues with appropriate liganding residues from RosettaMatch
		//we will fill this vector with the positions in the combined PDB that are the metal site
		utility::vector1< core::Size > metal_site;
		//we are assuming the last two residues of the match are the second partner's single ligand residue and the metal
		for ( core::Size i = 1; i <= matchlength-2; ++i ) {
			int pdbnum = match.pdb_info()->number(i);
			char chain = partner1.pdb_info()->chain(1);

			core::Size partner1_resid(partner1.pdb_info()->pdb2pose(chain, pdbnum));
			//  TR << "pdbnum and chain and partner1_resid" << pdbnum << chain << partner1_resid << std::endl;
			partner1.replace_residue(partner1_resid, match.residue(i), true);
			//strip termini and clean up residue connections
			core::pose::remove_variant_type_from_pose_residue(partner1, core::chemical::LOWER_TERMINUS_VARIANT, partner1_resid);
			core::pose::remove_variant_type_from_pose_residue(partner1, core::chemical::UPPER_TERMINUS_VARIANT, partner1_resid);
			partner1.conformation().update_polymeric_connection(partner1_resid-1, true);
			partner1.conformation().update_polymeric_connection(partner1_resid, true);
			partner1.conformation().update_polymeric_connection(partner1_resid+1, true);
			if ( partner1.residue_type(partner1_resid).aa() == core::chemical::aa_cys ) {
				core::pose::replace_pose_residue_copying_existing_coordinates(partner1, partner1_resid, CYZ);
			}
			metal_site.push_back(partner1_resid);
		}

		//After this graft, repacking is needed to accomodate the new site - actually this should probably be upgraded to limited design
		if ( !basic::options::option[ local::skip_sitegraft_repack ].value() ) {
			core::pack::task::PackerTaskOP task(core::pack::task::TaskFactory::create_packer_task(partner1));
			utility::vector1_bool packable(partner1length, true); //false = nobody is packable
			for ( core::Size i(1); i <= metal_site.size(); ++i ) packable[metal_site[i]] = false;
			task->restrict_to_residues(packable);
			task->restrict_to_repacking();
			for ( core::Size i(1); i <= partner1length; ++i ) {
				task->nonconst_residue_task(i).or_ex1(true);
				task->nonconst_residue_task(i).or_ex2(true);
				task->nonconst_residue_task(i).or_include_current(true);
				task->nonconst_residue_task(i).and_extrachi_cutoff(1); //histidines are near surface
			}

			protocols::simple_moves::PackRotamersMover pack(score_fxn, task);
			pack.apply(partner1);
			//partner1.dump_scored_pdb("partner1_repack.pdb", *score_fxn);
		}

		/////////////////////////////////////////////////////////////////////////////////////
		//set up partner2 fold tree to allow residue replacement to align partner2 onto match
		//this tree isolates partner2_residue by jumps and makes it root
		//this means moving that residue will move the whole pose (our intention)
		core::kinematics::FoldTree tree(partner2.size());
		int const p2_res(basic::options::option[local::partner2_residue].value());
		//For ubc12 as initially designed, where partner2_residue = 62, we want:
		//FOLD_TREE  EDGE 62 61 1  EDGE 62 63 2  EDGE 66 160 -1  EDGE 64 1 -1
		using core::kinematics::Edge;

		tree.add_edge(Edge(p2_res-1, 1, Edge::PEPTIDE));//center to start peptide edge
		tree.add_edge(Edge(p2_res+1, partner2length, Edge::PEPTIDE));//center to end peptide edge
		tree.add_edge(Edge(p2_res, p2_res-1, 1, "N", "C", false));//jump mimics peptide bond (C to N)
		tree.add_edge(Edge(p2_res, p2_res+1, 2, "C", "N", false));//jump mimics peptide bond (N to C)
		tree.delete_unordered_edge(1, partner2length, Edge::PEPTIDE);//kill the default edge

		tree.reorder(p2_res);
		partner2.fold_tree(tree);
		TR << "partner2 fold tree isolating " << p2_res << ": " << partner2.fold_tree() << std::flush;

		///////////////////////////////////////////////////////////////////////////////////////////////
		//align partner2 onto match
		//store current partner2 jumps
		using core::kinematics::Jump;
		Jump const j1(partner2.jump(1));
		Jump const j2(partner2.jump(2));

		//replace residue partner2_residue
		partner2.replace_residue(p2_res, match.residue(/*MatchPosition*/p2), false);
		core::pose::remove_variant_type_from_pose_residue(partner2, core::chemical::LOWER_TERMINUS_VARIANT, p2_res);
		core::pose::remove_variant_type_from_pose_residue(partner2, core::chemical::UPPER_TERMINUS_VARIANT, p2_res);
		partner2.conformation().update_polymeric_connection(p2_res-1, true);
		partner2.conformation().update_polymeric_connection(p2_res, true);
		partner2.conformation().update_polymeric_connection(p2_res+1, true);
		if ( partner2.residue_type(p2_res).aa() == core::chemical::aa_cys ) {
			core::pose::replace_pose_residue_copying_existing_coordinates(partner2, p2_res, CYZ);
		}

		//restore old jump transforms, which causes remainder of partner2 to align itself properly
		partner2.set_jump(1, j1);
		partner2.set_jump(2, j2);

		//we're lying about one particular atom (backbone hydrogen, number 11, from histidine, number 62) being missing to force it to be rebuilt
		core::id::AtomID_Mask missing( false );
		core::pose::initialize_atomid_map( missing, partner2 ); // dimension the missing-atom mask
		missing[ core::id::AtomID(partner2.residue_type(p2_res).atom_index("H"), p2_res) ] = true;
		partner2.conformation().fill_missing_atoms( missing );
		//partner2.dump_pdb("partner2_moved.pdb");

		//////////////////////////////////////////////////////////////////////////////
		//build the combined pose
		core::pose::Pose combined(partner1);
		combined.append_residue_by_jump(match.residue(/*MatchPosition*/metal), combined.size());
		combined.append_residue_by_jump(partner2.residue(1), 1);

		for ( core::Size i=2; i<=partner2length; ++i ) {
			combined.append_residue_by_bond(partner2.residue(i));
		}

		//calculate where special residues are
		//core::Size const combinedlength(combined.size());
		core::Size const metal_res(partner1length + 1);
		core::Size const p2_res_comb(partner1length + 1 + p2_res);

		//fix the pose into three chains
		combined.conformation().insert_chain_ending(partner1length);
		combined.conformation().insert_chain_ending(metal_res);

		//finish filling vector of metal site residues - maintain order from match.pdb
		metal_site.push_back(p2_res_comb);
		metal_site.push_back(metal_res);

		/////////////////////////////////ensure proper tautomer for all histidines////////////////////////
		core::Size const metal_site_size(metal_site.size());
		debug_assert( metal_site_size == matchlength);
		for ( core::Size i(1); i <= metal_site_size; ++i ) {
			core::conformation::Residue const old_rsd(combined.residue(metal_site[i]));
			if ( combined.residue_type(metal_site[i]).aa() != core::chemical::aa_his ) continue; //matters for HIS/HIS_D
			core::Vector const & metal_atom(combined.residue(metal_res).atom(1).xyz());
			core::Real dis_ND1 = (old_rsd.atom("ND1").xyz()).distance(metal_atom); //ND1
			core::Real dis_NE2 = (old_rsd.atom("NE2").xyz()).distance(metal_atom); //NE2

			//four cases exist:
			//  (is HIS_D) and  (dis_ND1 < dis_NE2) = change to NE2 HIS
			// !(is HIS_D) and !(dis_ND1 < dis_NE2) = change to ND1 HIS_D
			// !(is HIS_D) and  (dis_ND1 < dis_NE2) = okay
			//  (is HIS_D) and !(dis_ND1 < dis_NE2) = okay

			bool const HIS_D(combined.residue_type(metal_site[i]).name() == "HIS_D"); //true if HIS_D
			bool const ND1(dis_ND1 < dis_NE2); //true if should be HIS

			//This commented-out method ensures perfect placement of histidine ring (at expense of backbone)
			//   utility::vector1< std::pair< std::string, std::string > > pairs;
			//   pairs.push_back(std::make_pair("CG", "CG"));
			//   pairs.push_back(std::make_pair("ND1", "ND1"));
			//   pairs.push_back(std::make_pair("NE2", "NE2"));

			//if the two bools match, we have case 1 or 2 above
			if ( HIS_D == ND1 ) {
				//using namespace core::chemical;   using namespace core::conformation;
				std::string const new_name(HIS_D ? "HIS" : "HIS_D");

				//get type set from original residue; query it for a residue type of the other name
				core::chemical::ResidueTypeCOP new_type( core::pose::get_restype_for_pose(combined, new_name, old_rsd.type().mode() ) );
				TR << "mutating from " << old_rsd.name() << " to " << new_type->name() << " at combined position "
					<< metal_site[i] << std::endl;
				core::conformation::ResidueOP new_rsd(core::conformation::ResidueFactory::create_residue(*new_type, old_rsd, combined.conformation()));
				new_rsd->set_chi( 1, old_rsd.chi(1) );
				new_rsd->set_chi( 2, old_rsd.chi(2) );
				combined.replace_residue( metal_site[i], *new_rsd, true );
				//Rebuild connection-dependent atoms:
				combined.conformation().update_polymeric_connection(metal_site[i]-1, true);
				combined.conformation().update_polymeric_connection(metal_site[i], true);
				combined.conformation().update_polymeric_connection(metal_site[i]+1, true);
			}
		}

		//determine a good name for the starting pdb
		//combined.dump_scored_pdb("partner1_metal_partner2.pdb", *score_fxn);
		utility::file::FileName match_name( basic::options::option[local::match_pdb].value() );
		utility::file::FileName partner2_name( basic::options::option[local::partner2].value() );
		utility::file::FileName partner1_name( basic::options::option[local::partner1].value() );
		utility::file::FileName whole_name( match_name.base() + "_" + partner2_name.base() + "_" + partner1_name.base() + ".pdb" );
		TR << "saving starting pdb to " << whole_name.name() << std::endl;
		combined.dump_scored_pdb(whole_name.name(), *score_fxn);

		///////////////////////////////build special bond-mimic jumps for combined pose////////////////////////////
		//we are assuming the first atom in the metal residue is the metal...
		std::string const & metal_atom_name( combined.residue_type(metal_res).atom_name(1) );
		core::Vector const & metal_atom_xyz( combined.residue(metal_res).atom(1).xyz() );

		std::string const & partner1_atom_name( protocols::metal_interface::find_closest_atom(combined.residue(metal_site[1]), metal_atom_xyz) );
		std::string const & partner2_atom_name( protocols::metal_interface::find_closest_atom(combined.residue(metal_site[p2]), metal_atom_xyz) );

		//using core::kinematics::Edge;
		Edge const partner1_to_metal(metal_site[p1_1], metal_res, 1, partner1_atom_name, metal_atom_name, false);
		Edge const metal_to_partner2(metal_res, p2_res_comb, 2, metal_atom_name, partner2_atom_name, false);

		//////////////////////////////////////////////////////////////////////////////////
		///run movers

		//HACK HACK HACK to get pose out so it can be read into job distributor
		core::import_pose::pose_from_file(combined, whole_name.name(), core::import_pose::PDB_file);
		basic::options::option[ basic::options::OptionKeys::in::file::s ].value(whole_name.name());
		//combined.dump_scored_pdb("reread_partner1_metal_partner2.pdb", *score_fxn);

		protocols::moves::MoverOP mover( new protocols::metal_interface::ZincHeterodimerMover( metal_site, partner1_to_metal, metal_to_partner2 ) );

		protocols::jd2::JobDistributor::get_instance()->go(mover);

		TR << "************************d**o**n**e**********************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}


	return 0;
}
