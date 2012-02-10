// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )


// Package headers

// project headers

static basic::Tracer TR("devel.matdes.util");

namespace devel {
namespace matdes {

PackerTaskOP
make_interface_design_packertask(core:pose:Pose & pose) {
}

void
add_native_bias_constraints(Pose & pose, Real cst_weight, const std::set<Size>& design_pos) {
 	utility::vector1<core::scoring::constraints::ConstraintOP> favor_native_constraints;
	for(std::set<Size>::iterator pos = design_pose.begin(); pos != design_pos.end() ++ pos) {
		core::scoring::constraints::ConstraintOP cst = new core::scoring::constraints::ResidueTypeConstraint(pose, *pos, weight);
		favor_native_constraints.push_back(cst);
		cst->show(TR);
		TR << std::endl;
	}
}

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius) {
  using core::id::AtomID; 
  utility::vector1<Real> rsd_sasa(pose.n_residue(),0.0);
  core::id::AtomID_Map<Real> atom_sasa;
  core::id::AtomID_Map<bool> atom_mask;
  core::pose::initialize_atomid_map(atom_sasa,pose,0.0);
  core::pose::initialize_atomid_map(atom_mask,pose,false);
  for(Size i = 1; i <= pose.n_residue(); i++) {
    for(Size j = 1; j <= pose.residue(i).nheavyatoms(); j++) {
      atom_mask[AtomID(j,i)] = true;
    }
  } 
  core::scoring::calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius, false, atom_mask );
  utility::vector1<Real> sc_sasa(pose.n_residue(),0.0);
  for(Size i = 1; i <= pose.n_residue(); i++) {
    // Use CA as the side chain for Glys
    if(pose.residue(i).name3()=="GLY") sc_sasa[i] += atom_sasa[AtomID(2,i)];
    for(Size j = 5; j <= pose.residue(i).nheavyatoms(); j++) {
      sc_sasa[i] += atom_sasa[AtomID(j,i)];
    }
  }
  return sc_sasa;
}


/// @brief Find out which positions are near the inter-subunit interfaces
/// These will be further screened below, then passed to design()
std::set<Size>
pick_design_position(core::pose::Pose & pose, Size nsub_bblock, Real contact_dist, Real bblock_dist, probe_radius) {
	
	Real bblock_dist_sq = bblock_dist * bblock_dist;
	vector1<Size> design_pos;
	// score the pose so the we check for clashes
	ScoreFunctionOP score12 = ScoreFunctionFactory::create_score_function("standard", "score12");
	score12->apply(pose);
	// get the accesible residues
	// THIS NEED TO BE FIXED, MONO SHOULD COME FROM SYMINFO
	vector1<Real> sc_sasa = sidechain_sasa(mono,probe_radius); 
	
	// get the residues in chain A that make contacts within the building block
	std::set<Size> bblock_interface;
	std::set<Size> bblock_interface_clashes;

	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		for (Size jr=ir+1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			// skip if residues are in same building block
   		if ( sym_info->subunit_index(jr) > nsub_bblock or sym_info->subunti_index(jr) == 1)  continue;
			// loop over the atoms in the residue 
			for (Size ia = 1; ia<=pose_for_design.residue(ir).nheavyatoms(); ia++) {
				bool residue_registered = false;
				for (Size ja = 1; ja<=pose_for_design.residue(jr).nheavyatoms(); ja++) {
					if (pose.residue(ir).xyz(ia).distance_squared(pose.residue(jr).xyz(ja)) <= bblock_dist_sq)	{
						// However, if the residue in question is clashing badly (usually with a
						// residue from another building block), it needs to be designed.
						core::scoring::EnergyMap em1 = pose_for_design.energies().residue_total_energies(ir);
						Real resi_fa_rep = em1[core::scoring::fa_rep];
						if (resi_fa_rep > 3.0) {
							bblock_interface_clashes.add(ir);
							residue_registered = true;
							break;
						} else {
							bblock_interface.add(ir);
							residue_registered = true;
							break;
						}
					}
				}
				if(residue_registered) break;
			}
		}
	}

	
	Real const contact_dist_sq = contact_dist * contact_dist;
	vector1<bool> indy_resis = sym_info->independent_residues();
	Sizes design_pos;

	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if (!indy_resis[ir]) continue;
		std::string atom_i = "";
		if (pose_for_design.residue(ir).name3() == "GLY") {
			atom_i = "CA";
		} else {
			atom_i = "CB";
		}
		for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			// skip residue jr if it is in the same building block
   		if ( sym_info->subunit_index(jr) <= nsub_bblock ) continue;
			std::string atom_j = "";
			if (pose_for_design.residue(jr).name3() == "GLY") {
				atom_j = "CA";
			} else {
				atom_j = "CB";
			}
			// Here we are filtering the residues that are in contact 
			if (pose_for_design.residue(ir).xyz(atom_i).distance_squared(pose_for_design.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
				if( ! bblock_interface.count(ir) && bblock_interface_clashes.count(ir) && sc_sasa[ir] > 0 )
					design_pos.insert(ir);
				break;
			}
		}
	}
	return design_pos;
}


} // devel
} // matdes
