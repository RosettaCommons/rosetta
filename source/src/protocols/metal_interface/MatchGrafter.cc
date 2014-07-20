// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/metal_interface/MatchGrafter.cc
/// @brief  Takes a scaffold protein and a match pdb from RosettaMatch, grafts the match onto the protein.  For zinc homodimer design, it can then combine two grafted poses by overlaying the zinc atoms.
/// @author Bryan Der


#include <protocols/metal_interface/MatchGrafter.hh>
#include <protocols/metal_interface/FindClosestAtom.hh> // find closest atom to zinc

#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueTypeSet.hh> //for changing the HIS tautomer
#include <core/chemical/ResidueType.hh> //for CYZ residue
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/ChemicalManager.hh> //CENTROID, FA_STANDARD
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh> //create HIS of proper tautomer
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/PDB_Info.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.io.hh> //print 'metalsite_atoms_'
#include <utility/vector1.hh>

#include <basic/Tracer.hh>


// Forward declarations
using basic::T;
using basic::Error;
using basic::Warning;
static basic::Tracer TR( "protocols.metal_interface.MatchGrafter" );

typedef core::pose::Pose Pose;
typedef numeric::xyzVector<core::Real> point;

namespace protocols {
namespace metal_interface {


MatchGrafter::MatchGrafter() // default constructor
//	: metalsite_atoms_ ( 5, 0 ), metalsite_residues_ ( 5, 0 )
{
}

MatchGrafter::~MatchGrafter()
{
}


///@brief Takes match pose (2 residues + zinc) and partner pose, grafts match onto partner, returns grafted partner
Pose
MatchGrafter::graft( Pose & match,
						Pose & partner_ungrafted ) {

	core::chemical::ResidueTypeSetCAP typeset(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
	core::chemical::ResidueType const & CYZ(typeset->name_map("CYZ"));

	//graft match, excluding zinc, onto target
	for ( core::Size i = 1; i <= match.pdb_info()->nres() - 1; ++i ){

		//int pdbnum = match.pdb_info()->number(i);
		//core::Size chain_num = match.residue(i).chain();
		//char chain = match.pdb_info()->chain(chain_num);
		//core::Size partner_resid(partner_ungrafted.pdb_info()->pdb2pose(chain, pdbnum));
		// replace residue requires residue number and residue object
		//TR << "Replacing now: nres " << match.pdb_info()->nres() << " chain " << chain << " chain_num " << chain_num << " pdbnum " << pdbnum << " partner_resid " << partner_resid << " match.residue(i) " << i << std::endl;

		core::Size partner_resid( (core::Size) match.pdb_info()->number(i) );
		TR << "partner_resid " << partner_resid << std::endl;

		partner_ungrafted.replace_residue(partner_resid, match.residue(i), true);

		core::pose::remove_variant_type_from_pose_residue(partner_ungrafted, core::chemical::LOWER_TERMINUS_VARIANT, partner_resid);
		core::pose::remove_variant_type_from_pose_residue(partner_ungrafted, core::chemical::UPPER_TERMINUS_VARIANT, partner_resid);
		partner_ungrafted.conformation().update_polymeric_connection(partner_resid-1, true);
		partner_ungrafted.conformation().update_polymeric_connection(partner_resid, true);
		partner_ungrafted.conformation().update_polymeric_connection(partner_resid+1, true);
		if(partner_ungrafted.residue_type(partner_resid).aa() == core::chemical::aa_cys){
			core::pose::replace_pose_residue_copying_existing_coordinates(partner_ungrafted, partner_resid, CYZ);
		}
	}
	Pose partner(partner_ungrafted);

	partner.append_residue_by_jump(match.residue(match.pdb_info()->nres()/*metal*/), partner.total_residue());
	//	partner.dump_pdb("3UBQ_graft.pdb");

	return partner;
}//graft

///@brief Takes two matches (for residue number info) and two grafted poses, overlays the zinc atoms, and uses a partner2 fold tree to combine the two poses
Pose
MatchGrafter::build_combined_pose_with_zinc_overlay( Pose & partner1, Pose & partner2 ) { //pass match1, match2???

	// keep in mind that zinc as been appended to both partners, thus adding a residue
	// zinc of partner 1 will be retained in combined pose
	// zinc of partner 2 (p2_zinc) is only needed to overlay the two zinc atoms
	core::Size metal_res_num = partner1.total_residue();
	//core::Size partner1length = partner1.total_residue() - 1;
	core::Size partner2length = partner2.total_residue() - 1;

	/////////////////////////////////////////////////////////////////////////////////////
	//set up partner1 fold tree to allow zinc atoms to superimpose
	//moving the zinc will move the whole pose (our intention)
	core::kinematics::FoldTree tree(partner2.total_residue());
	tree.clear();
	int const p2_zinc( partner2.total_residue() ); // zinc is last residue of partner2, after it's been grafted
                                                 // note that 'partner2length' excludes the added zinc

	using core::kinematics::Edge;
	//Edges: PEPTIDE = -1, CHEMICAL = -2, jump = 1, 2...
	tree.add_edge(Edge(1, partner2length,  Edge::PEPTIDE));//entire partner2 protein
	tree.add_edge(Edge(p2_zinc, 1, 1));//zinc to arbitrary residue 1, 1st jump

	//define zinc as root of the fold tree
	tree.reorder(p2_zinc);
	partner2.fold_tree(tree);
	TR << "partner2 fold tree isolating " << p2_zinc << ": " << partner2.fold_tree() << std::flush;


	////////////////////////////////////////////////////////////////////////////////////
	//superimpose partner2's zinc onto the partner1's zinc
	//store current partner2 jumps
	using core::kinematics::Jump;
	Jump const j1(partner2.jump(1));

	//replace zinc on (moving) scaffold with zinc on (stationary) target
	partner2.replace_residue(p2_zinc, partner1.residue(metal_res_num/*zinc is last residue of partner 1*/), false);
	core::pose::remove_variant_type_from_pose_residue(partner2, core::chemical::LOWER_TERMINUS_VARIANT, p2_zinc);
 	core::pose::remove_variant_type_from_pose_residue(partner2, core::chemical::UPPER_TERMINUS_VARIANT, p2_zinc);

	//restore old jump transforms, which causes remainder of partner2 to align itself properly
 	partner2.set_jump(1, j1);

	////////////////////////////////////////////////////////////////////////////////////
	//build the combined pose
	Pose combined(partner1);
	combined.append_residue_by_jump(partner2.residue(1), 1);
	for (core::Size i=2; i<=partner2length; ++i){
		combined.append_residue_by_bond(partner2.residue(i));
	}
	//combined.dump_pdb("combined.pdb");

	//fix the pose into three chains
	//	combined.conformation().insert_chain_ending( partner1length ); // with this line, chains become A, B(zinc), C.  Without this line, chains become A, B
	combined.conformation().insert_chain_ending( metal_res_num );

	return combined;

}//build_combined_pose_with_zinc_overlay


void
MatchGrafter::ensure_proper_his_tautomers(
	Pose & combined,
	utility::vector1< core::Size > metalsite_seqpos
)
{
	/////////////////////////////////ensure proper tautomer for all histidines////////////////////////

	assert( combined.residue( metalsite_seqpos[1] ).is_ligand() );


	for( core::Size i(2) /*metal is 1*/; i <= metalsite_seqpos.size(); ++i){

		assert( combined.residue( metalsite_seqpos[i] ).is_protein() );

		core::conformation::Residue const old_rsd(combined.residue(metalsite_seqpos[i]));

		if( combined.residue_type(metalsite_seqpos[i]).aa() != core::chemical::aa_his) continue; //matters for HIS/HIS_D


		core::Vector const & metal_atom(combined.residue(metalsite_seqpos[1]).atom(1).xyz());


		core::Real dis_ND1 = (old_rsd.atom("ND1").xyz()).distance(metal_atom); //ND1
		core::Real dis_NE2 = (old_rsd.atom("NE2").xyz()).distance(metal_atom); //NE2
		TR << "line 4" << std::endl;

		bool const HIS_D(combined.residue_type(metalsite_seqpos[i]).name() == "HIS_D"); //true if HIS_D
		bool const ND1(dis_ND1 < dis_NE2); //true if should be HIS
		TR << "line 5" << std::endl;

		//if the two bools match, we have case 1 or 2 above
		if( HIS_D == ND1 ){

			std::string const new_name(HIS_D ? "HIS" : "HIS_D");

			//get type set from original residue; query it for a residue type of the other name
			core::chemical::ResidueType const & new_type(old_rsd.residue_type_set().name_map(new_name));
			TR << "mutating from " << old_rsd.name() << " to " << new_type.name() << " at combined position "
				 << metalsite_seqpos[i] << std::endl;
			core::conformation::ResidueOP new_rsd(core::conformation::ResidueFactory::create_residue(new_type, old_rsd, combined.conformation()));
			new_rsd->set_chi( 1, old_rsd.chi(1) );
			new_rsd->set_chi( 2, old_rsd.chi(2) );
			combined.replace_residue( metalsite_seqpos[i], *new_rsd, true );
		}
	}
}//ensure_proper_his_tautomers

} // namespace metal_interface
} // namespace protocols

