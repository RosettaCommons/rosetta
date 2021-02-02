// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
//
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/glycopeptide_docking/utils.cc
/// @brief util functions for peptide/glycopeptide sampling and refinement
/// @author Sai Pooja Mahajan (saipooja@gmail.com)
//
#include <protocols/glycopeptide_docking/utils.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <numeric/xyzVector.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <protocols/carbohydrates/SimpleGlycosylateMover.hh>
#include <protocols/carbohydrates/GlycanSampler.hh>
#include <protocols/glycopeptide_docking/GlycopeptideDockingFlags.hh>
#include <protocols/docking/util.hh>
#include <protocols/docking/metrics.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh> //needed for rmsd_no_super
#include <core/scoring/ScoreType.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/extra_pose_info_util.hh>


static basic::Tracer TR( "protocols.glycopeptide_docking.utils" );
namespace protocols {
namespace glycopeptide_docking {

/// @details Tracks distance between donor and acceptor.
/// If the distance exceeds a specified distance, the protocol
/// terminates sampling early.
core::Real
calculate_sampled_distance(core::pose::Pose const &pose, core::Size const glycosylation_residue, core::Size const donor_residue) {

	numeric::xyzVector<core::Real> vecdistance(0.0,0.0,0.0);
	//peptide
	if ( pose.residue( donor_residue ).name() == "pdb_HWU" || pose.residue( donor_residue ).name() == "pdb_UPG" ) {

		//distance between anomeric carbon on sugar and sidechain oxygen or nitrogen on the peptide substrate that gets modified
		if ( pose.residue(glycosylation_residue).name() == "ASN" ) { //N-linked
			vecdistance = (pose.residue( glycosylation_residue ).xyz("ND2") - pose.residue( donor_residue ).xyz("C1\'"));
		} else { //O-linked
			vecdistance = (pose.residue( glycosylation_residue ).xyz("OG1") - pose.residue( donor_residue ).xyz("C1\'"));
		}

	}
	//glycopeptide
	if ( pose.residue( donor_residue ).name() == "UDP:non-conjugated" ) {

		if ( pose.residue(glycosylation_residue).name() != "ASN" ) { //important distance for retaining glycosyltransferases
			vecdistance = (pose.residue(glycosylation_residue).xyz("N") - pose.residue( donor_residue ).xyz("3OPB"));
		}

	}
	core::Real distance(999.0);
	try {
		distance = static_cast<core::Real>(vecdistance.norm());
	} catch (...) {
		TR<<"distance too large to calculate"<<std::endl;
	}
	return distance;


}
/// @details For specific glycosylation cases when the sugar is donated from a UDP to a substrate
///  This function will be updated as more cases are added for sugar donors
///  For example, pdb_HWU is encountered in O-glycosylating enzymes GalNAcTs (20 enzymes)
///  TODO: Create a dictionary of possible donors and how to treat them.
std::map<std::string,core::Real>
calculate_additional_glycosylation_metrics(core::pose::Pose const &pose, core::Size const glycosylation_residue, core::Size const donor_residue ){
	/*  For specific glycosylation cases when the sugar is donated from a UDP to a substrate
	*  This function will be updated as more cases are added for sugar donors
	*  For example, pdb_HWU is encountered in O-glycosylating enzymes GalNAcTs (20 enzymes)
	*  */
	std::map<std::string,core::Real> values;
	if ( pose.residue( donor_residue ).name() == "pdb_HWU" || pose.residue( donor_residue ).name() == "pdb_UPG" ) {

		values[std::string("distance_catalysis")]=calculate_sampled_distance(pose,glycosylation_residue,donor_residue);

		numeric::xyzVector<core::Real> vec_catalysisdist_1(0.0,0.0,0.0);
		numeric::xyzVector<core::Real> vec_catalysisdist_2(0.0,0.0,0.0);

		if ( pose.residue(glycosylation_residue).name() == "SER" || pose.residue(glycosylation_residue).name() == "THR" ) {
			vec_catalysisdist_1 = (pose.residue(glycosylation_residue).xyz("N") - pose.residue(donor_residue).xyz("O1B"));
			vec_catalysisdist_2 = (pose.residue(glycosylation_residue).xyz("N") - pose.residue(donor_residue).xyz("O1\'"));
			values[std::string( "distance_catalysis_HWUO1B-THR7N")]=vec_catalysisdist_1.norm();
			values[std::string( "distance_catalysis_HWUO1-THR7N")]=vec_catalysisdist_2.norm();
		}
	}

	TR.Debug<<"Donor residue "<< pose.residue( donor_residue ).name() <<std::endl;
	/*  When UDP is not conjugated to sugar
	*  Can be used to study the glycosylated final product*/
	if ( pose.residue( donor_residue ).name() == "UDP:non-conjugated" ) {

		numeric::xyzVector<core::Real> vec_catalysisdist_1(0.0,0.0,0.0);
		numeric::xyzVector<core::Real> vec_catalysisdist_2(0.0,0.0,0.0);

		if ( pose.residue(glycosylation_residue).name() != "ASN" ) {
			vec_catalysisdist_1 = (pose.residue(glycosylation_residue).xyz("N") - pose.residue( donor_residue ).xyz("2OPB"));
			vec_catalysisdist_2 = (pose.residue(glycosylation_residue).xyz("N") - pose.residue( donor_residue ).xyz("3OPB"));
			values[std::string( "distance_catalysis_UDP2OPB-THR7N")]=vec_catalysisdist_1.norm();
			values[std::string( "distance_catalysis_UDP3OPB-THR7N")]=vec_catalysisdist_2.norm();
		}
	}


	return values;
}

/// @details Write additional pdb files during low resolution and
/// high resolution sampling. Useful for debugging and making movies.
void
write_debug_pdb(core::pose::Pose const &pose, core::Size const nstruct_max, core::Size const nstruct_index, std::string name)
{
	core::Size nstruct_width = 0;
	for ( core::Size i = 1; i <= nstruct_max || nstruct_width < 4; i *= 10 ) {
		nstruct_width += 1;
	}
	std::ostringstream ss;
	std::string prefix = basic::options::option[basic::options::OptionKeys::out::prefix]();
	std::string suffix = basic::options::option[basic::options::OptionKeys::out::suffix]();
	ss << prefix << "_" << name <<"_";
	ss << std::setfill('0') << std::setw(nstruct_width) << nstruct_index << "_" << suffix;
	TR << "Dumping metrics pose " << ss.str() << std::endl;
	pose.dump_pdb(ss.str() + ".pdb");


}

/// @details Calls simpleglycosylate mover to generate
/// pre-glycosylated peptides for sampling.
void
glycosylate_residues(core::pose::Pose &pose,utility::vector1<core::Size> const &sugar_residues,utility::vector1<std::string> &sugar_names ){
	using namespace protocols::carbohydrates;
	using namespace core::select::residue_selector;
	SimpleGlycosylateMover glycosylator;
	TR << "Positions, Sugars " << sugar_residues << " " << sugar_names<< std::endl;
	glycosylator.set_glycosylations(sugar_names);
	ResidueIndexSelectorCOP sugar_residue_selector =  ResidueIndexSelectorCOP (new ResidueIndexSelector(sugar_residues)) ;
	glycosylator.set_residue_selector(sugar_residue_selector);
	glycosylator.apply(pose);
}

/// @details setup specialized foldtree.
void
setup_glycosylation_foldtree( core::pose::Pose & pose, protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags, core::kinematics::FoldTreeOP ft_docking)
{
	using namespace core::kinematics;
	/*Glycosylation foldtree must always be setup for every decoy (overwrite previous decoy foldtree)
	//Overwrite any other foldtrees
	//Outside this function, the user can only specify the type of foldtree - docking or outward
	*/
	TR << "Setting up or Overwriting foldtree." << std::endl;
	// Setup jumps
	ft_docking = core::kinematics::FoldTreeOP(utility::pointer::make_shared<core::kinematics::FoldTree>());
	utility::vector1<Edge> jump_edges = pose.fold_tree().get_jump_edges();
	ft_docking->add_edge(1, flags->anchor_residue_substrate(), flags->jump_num_substrate()); // Jump to substrate anchor
	for ( core::Size jump = 1; jump <= jump_edges.size() - flags->jump_num_substrate(); jump++ ) {
		// Jumps to remaining ligands, cosubstrates, etc.
		ft_docking->add_edge(flags->first_residue_enzyme(), jump_edges[jump].stop(), jump + flags->jump_num_substrate());
	}

	// Setup backbone
	ft_docking->add_edge(flags->first_residue_enzyme(), flags->last_residue_enzyme(), -1); // Enzyme backbone
	if ( flags->substrate_type() == "peptide" and flags->tree_type() == "outward" ) {
		ft_docking->add_edge(flags->anchor_residue_substrate(), flags->first_residue_substrate(), -1); // Downstream substrate
		ft_docking->add_edge(flags->anchor_residue_substrate(), flags->last_residue_substrate(), -1);  // Upstream substrate
		// Setup branches
		utility::vector1<Edge> chemical_edges = pose.fold_tree().get_chemical_edges(); // Connections to glycan residues from anchor out
		for ( core::Size chemical = 1; chemical <= chemical_edges.size(); chemical++ ) {
			ft_docking->add_edge(chemical_edges[chemical].start(), chemical_edges[chemical].stop(), -2);
		}

		pose.fold_tree(*ft_docking);
	} else if ( flags->substrate_type() == "peptide" and flags->tree_type() == "docking" ) {
		using namespace docking;
		//do nothing because the chemical bonds will take care of the foldtree
		std::string const partners(flags->get_upstream_chain() + "_" + flags->get_downstream_chain());
		utility::vector1<int> movable_jumps(1, 1);
		protocols::docking::setup_foldtree(pose, partners, movable_jumps);
		*ft_docking = pose.fold_tree();
	} else if ( flags->substrate_type() == "lipid" ) {
		//not implemented
		TR << "Lipid substrate not implemented yet." << std::endl; //TODO: Change to error.
	} else {
		//not implemented
		TR << "Unknown substrate type. Not implemented yet." << std::endl; //TODO: Change to error.
	}
}


/// @details Record a collection of decoy metrics.
/// Includes distances, rmsds and interaction energies.
void
record_pose_metrics( core::pose::Pose & pose,
	protocols::glycopeptide_docking::GlycopeptideDockingFlagsOP flags,
	utility::vector1< int > const jumps,
	core::pose::PoseOP ref_pose)
{
	using namespace core::scoring;
	using namespace protocols::docking;

	ScoreFunctionOP sf = get_score_function();
	sf->set_weight( atom_pair_constraint, 1.0 ); /* TODO: Make this an option */

	TR << std::endl << " Metrics for this decoy:" << std::endl;
	TR << std::endl << " Final score: " << std::endl;
	sf->show( std::cout, pose );
	TR << std::endl;

	ScoreFunctionOP sf_minus_constraint = get_score_function();
	sf_minus_constraint->set_weight( atom_pair_constraint, 0.0 );

	/* TODO: Add a flag that indicates the glycan residues */
	std::list< core::Size > glycan_and_substrate_residues;
	std::list< core::Size > glycan_residues;
	core::Size glycan_residue = flags->last_residue_substrate() + 1;
	glycan_residues.push_back( glycan_residue );
	glycan_and_substrate_residues.push_back( glycan_residue );

	std::list< core::Size > substrate_residues;
	for ( core::Size res = flags->first_residue_substrate(); res <= flags->last_residue_substrate(); res++ ) {
		substrate_residues.push_back( res );
		glycan_and_substrate_residues.push_back( res );
	}

	std::list< core::Size > substrate_sequon_residues;

	core::Size res = flags->anchor_residue_substrate();
	core::Size sequon_res_upstream = 1; //Can be a flag
	core::Size sequon_res_downstream = 3; //Can be a flag
	substrate_sequon_residues.push_back( res - sequon_res_upstream );
	substrate_sequon_residues.push_back( res );
	substrate_sequon_residues.push_back( res + sequon_res_downstream );
	core::Real const interaction_energy( calc_interaction_energy( pose, sf, jumps ) );
	core::Real const atom_pair_constraint_sc( sf->score_by_scoretype( pose, atom_pair_constraint ) );
	core::Real const substrate_ca_rmsd( CA_rmsd( pose, *ref_pose, substrate_residues ) );
	core::Real const substrate_and_glycan_ca_rmsd( CA_rmsd(pose, *ref_pose, glycan_and_substrate_residues) );
	core::Real const fraction_native_contacts( calc_Fnat( pose, *ref_pose, sf, jumps ) );

	core::pose::setPoseExtraScore( pose, "interaction_energy", interaction_energy );
	core::pose::setPoseExtraScore( pose, "atom_pair_constraint", atom_pair_constraint_sc );
	core::pose::setPoseExtraScore( pose, "substrate_ca_rmsd", substrate_ca_rmsd);
	core::pose::setPoseExtraScore(pose, "substrate_and_glycan_ca_rmsd", substrate_and_glycan_ca_rmsd);
	core::pose::setPoseExtraScore(pose, "Fnat", fraction_native_contacts);

	/*no super - no rotation or translation - get rmsd as is - this is what we want.
	we dont want peptide to be moved to the native peptide position*/
	ObjexxFCL::FArray1D_bool substrate_residues_bool( pose.size(), false );
	ObjexxFCL::FArray1D_bool substrate_sequon_residues_bool( pose.size(), false );

	// Debugging block
	TR.Debug << "SubRes ";
	for ( auto &ires : substrate_residues ) {
		substrate_residues_bool(ires) = true;
		TR.Debug << ires << " ";
	}
	TR.Debug << std::endl << "SubSeqRes ";
	for ( auto &ires : substrate_sequon_residues ) {
		substrate_sequon_residues_bool(ires) = true;
		TR.Debug << ires << " ";
	}
	TR.Debug<<std::endl;

	//CA atoms
	core::Real const substrate_ca_no_super_rmsd(rmsd_no_super_subset( pose, *ref_pose, substrate_residues_bool, is_protein_CA ));
	core::Real const substrate_sequon_ca_no_super_rmsd(rmsd_no_super_subset( pose, *ref_pose, substrate_sequon_residues_bool, is_protein_CA ));
	core::pose::setPoseExtraScore(pose, "substrate_ca_no_super_rmsd", substrate_ca_no_super_rmsd);
	core::pose::setPoseExtraScore(pose, "substrate_sequon_ca_no_super_rmsd", substrate_sequon_ca_no_super_rmsd);

	TR << "  Interaction Energy:          " << interaction_energy << std::endl;
	TR << "  Atom Pair Constraint:        " << atom_pair_constraint_sc << std::endl;
	TR << "  Substrate RMSD:              " << substrate_ca_rmsd << std::endl;
	TR << "  Substrate and Glycan CA equivalent RMSD: " << substrate_and_glycan_ca_rmsd << std::endl;
	TR << "  Fraction of native contacts: " << fraction_native_contacts << std::endl;
	TR << "  Substrate CA No Super RMSD: " << substrate_ca_no_super_rmsd << std::endl;
	TR << "  Substrate Sequon CA No Super RMSD: " << substrate_sequon_ca_no_super_rmsd << std::endl;

	if ( flags->additional_metrics() ) {
		TR << " Additional metrics added for this run" << std::endl;
		TR << pose.residue( flags->get_sugar_donor() ).name() << std::endl;
		std::map<std::string,core::Real> additional_metrics = calculate_additional_glycosylation_metrics(pose,flags->glycosylation_residue_substrate(),flags->get_sugar_donor());
		for ( auto &entry : additional_metrics ) {
			// job.add_string_real_pair( entry.first ,entry.second );
			core::pose::setPoseExtraScore(pose, entry.first ,entry.second );
			TR << entry.first << " " << entry.second << std::endl;
		}
	}

}

} //glycopeptide_docking
} //protocols
