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
#include <devel/matdes/util.hh>

// project headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/id/AtomID_Map.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>


static basic::Tracer TR("devel.matdes.util");

namespace devel {
namespace matdes {

core::pack::task::PackerTaskOP
make_interface_design_packertask(core::pose::Pose & /*pose*/) {
	return NULL;
}

void
add_native_bias_constraints(core::pose::Pose & pose, Real cst_weight, const std::set<Size>& design_pos) {
 	utility::vector1<core::scoring::constraints::ConstraintOP> favor_native_constraints;
	for(std::set<Size>::iterator pos = design_pos.begin(); pos != design_pos.end(); ++pos) {
		core::scoring::constraints::ConstraintOP cst = new core::scoring::constraints::ResidueTypeConstraint(pose, *pos, cst_weight);
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
pick_design_position(core::pose::Pose const & pose, Size nsub_bblock, Real contact_dist, Real bblock_dist,Real probe_radius) {
	TR.Debug << "Picking design positions" << std::endl;
	using namespace core::scoring;	
	using namespace core::pose;
	using utility::vector1;
	using core::conformation::symmetry::SymmetryInfoCOP;
	core::pose::Pose p(pose);
  SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(p);
	std::set<Size> design_pos;
	Real bblock_dist_sq = bblock_dist * bblock_dist;
	// score the pose so the we check for clashes
	ScoreFunctionOP sf = new  core::scoring::symmetry::SymmetricScoreFunction( ScoreFunctionFactory::create_score_function("standard", "score12") );
	core::scoring::methods::EnergyMethodOptions eo = sf->energy_method_options();
	eo.exclude_monomer_hack_elec(true);
	sf->set_energy_method_options(eo);
	sf->score(p);
	// get the accesible residues
	Pose mono(p, 1, sym_info->get_nres_subunit());
	vector1<Real> sc_sasa = sidechain_sasa(mono,probe_radius); 
	
	// get the residues in chain A that make contacts within the building block
	std::set<Size> intra_bblock_interface;
	std::set<Size> intra_bblock_interface_clashes;

	TR.Debug << "Looking for the residues in the interface or that are clashing" << std::endl;
	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		for (Size jr=ir+1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			TR.Debug << "res_i = " << ir << "\tres_j = " << jr << std::endl;
			// skip if residues are in same building block
   		if ( sym_info->subunit_index(jr) > nsub_bblock or sym_info->subunit_index(jr) == 1){
					TR.Debug << "skipping residue " << jr << ", same building block as residue " << ir << std::endl;
					continue;
			}
			// loop over the atoms in the residue 
			for (Size ia = 1; ia<=p.residue(ir).nheavyatoms(); ia++) {
				bool residue_registered = false;
				for (Size ja = 1; ja<=p.residue(jr).nheavyatoms(); ja++) {
					if (p.residue(ir).xyz(ia).distance_squared(p.residue(jr).xyz(ja)) <= bblock_dist_sq)	{
						// However, if the residue in question is clashing badly (usually with a
						// residue from another building block), it needs to be designed.
						core::scoring::EnergyMap em1 = p.energies().residue_total_energies(ir);
						Real resi_fa_rep = em1[core::scoring::fa_rep];
						if (resi_fa_rep > 3.0) {
							TR.Debug << "residue " << ir << " is clashing with residue " <<jr << std::endl;
							intra_bblock_interface_clashes.insert(ir);
							intra_bblock_interface_clashes.insert(jr);
							residue_registered = true;
							break;
						} else {
							TR.Debug << "residue " << ir << " is in the interface." << std::endl;
							intra_bblock_interface.insert(ir);
							intra_bblock_interface.insert(jr);
							residue_registered = true;
							break;
						}
					}
				}
				if(residue_registered) break;
			}
		}
	}

	TR.Debug  << "interface residues: " << stringify_iterable(intra_bblock_interface) << std::endl;
	TR.Debug  << "clashing  residues: " << stringify_iterable(intra_bblock_interface_clashes) << std::endl;

	Real const contact_dist_sq = contact_dist * contact_dist;
	vector1<bool> indy_resis = sym_info->independent_residues();

	TR.Debug << "filtering residues that are in contact with other residues on the interface (distance less than "<< contact_dist  <<")" << std::endl;
	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if (!indy_resis[ir]) continue;
		std::string atom_i = (p.residue(ir).name3() == "GLY") ? "CA" : "CB";

		for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
			// skip residue jr if it is in the same building block
   		if ( sym_info->subunit_index(jr) <= nsub_bblock ) continue;
			std::string atom_j = (p.residue(jr).name3() == "GLY") ? "CA" : "CB";
			// Here we are filtering the residues that are in contact 
			if (p.residue(ir).xyz(atom_i).distance_squared(p.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
				TR.Debug << "residue " << ir << " in contact with residue " << jr << std::endl;
				TR.Debug << "\tinterface = " <<  intra_bblock_interface.count(ir) << std::endl;
				TR.Debug << "\tclashes = " <<  intra_bblock_interface_clashes.count(ir) << std::endl;
				TR.Debug << "\tsasa = " << sc_sasa[ir] << std::endl;
				// Get only the exposed residues that are in the interface between building blocks or clashing in the building
				// between the components of the building block.
				if( (!intra_bblock_interface.count(ir) ||  intra_bblock_interface_clashes.count(ir)) && sc_sasa[ir] > 0.0 )
					design_pos.insert(ir);
				break;
			}
		}
	}
	return design_pos;
}


} // devel
} // matdes
