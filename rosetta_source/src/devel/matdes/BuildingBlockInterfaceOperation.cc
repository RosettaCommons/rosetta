// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/matdes/BuildingBlockInterfaceOperation.cc
/// @brief  Restrict design to only residues at inter-building block interfaces
/// @author Neil King (neilking@uw.edu) Rocco Moretti (rmoretti@u.washington.edu)

// Unit Headers
#include <devel/matdes/BuildingBlockInterfaceOperation.hh>
#include <devel/matdes/BuildingBlockInterfaceOperationCreator.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// C++ Headers

static basic::Tracer TR("devel.matdes.BuildingBlockInterfaceOperation" );

namespace devel {
namespace matdes {

core::pack::task::operation::TaskOperationOP
BuildingBlockInterfaceOperationCreator::create_task_operation() const
{
	return new BuildingBlockInterfaceOperation;
}


BuildingBlockInterfaceOperation::BuildingBlockInterfaceOperation( core::Size nsub_bblock /* = 1 */,
		core::Real contact_dist /* = 10*/, core::Real bblock_dist /*= 5 */, core::Real fa_rep_cut /* = 3.0 */ ):
	nsub_bblock_(nsub_bblock),
	contact_dist_(contact_dist),
	bblock_dist_(bblock_dist),
	fa_rep_cut_(fa_rep_cut)
{}

BuildingBlockInterfaceOperation::~BuildingBlockInterfaceOperation() {}

core::pack::task::operation::TaskOperationOP BuildingBlockInterfaceOperation::clone() const
{
	return new BuildingBlockInterfaceOperation( *this );
}

void
BuildingBlockInterfaceOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	core::pose::Pose scored_pose( pose );
	core::scoring::ScoreFunctionFactory::create_score_function( "standard", "score12" )->score(scored_pose);

	core::conformation::symmetry::SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
	utility::vector1<bool> indy_resis = sym_info->independent_residues();
	core::Real bblock_dist_sq = bblock_dist_ * bblock_dist_;
	core::Real const contact_dist_sq = contact_dist_ * contact_dist_;
	std::set<Size> design_pos;
  std::set<Size> intra_bblock_interface;

  // get the residues in chain A that make contacts within the building block
  TR.Debug << "Looking for the residues in the interface or that are clashing" << std::endl;
  for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
		if (!indy_resis[ir]) continue;
    for (Size jr=ir+1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
      TR.Debug << "res_i = " << ir << "\tres_j = " << jr << std::endl;
      // skip if residues are in same building block
      if ( sym_info->subunit_index(jr) > nsub_bblock_ || sym_info->subunit_index(jr) == 1){
          TR.Debug << "skipping residue " << jr << ", same building block as residue " << ir << std::endl;
          continue;
      }
      // loop over the atoms in the residue
      for (Size ia = 1; ia<=pose.residue(ir).nheavyatoms(); ia++) {
        bool residue_registered = false;
        for (Size ja = 1; ja<=pose.residue(jr).nheavyatoms(); ja++) {
          if (pose.residue(ir).xyz(ia).distance_squared(pose.residue(jr).xyz(ja)) <= bblock_dist_sq)  {
            // However, if the residue in question is clashing badly (usually with a
            // residue from another building block), it needs to be designed.
            core::scoring::EnergyMap em = scored_pose.energies().residue_total_energies(ir);
            core::Real resi_fa_rep = em[core::scoring::fa_rep];
            if (resi_fa_rep > 3.0) {
              TR.Debug << "residue " << ir << " is clashing with residue " <<jr << std::endl;
              residue_registered = true;
              break;
            } else {
              TR.Debug << "residue " << ir << " is in the interface." << std::endl;
              intra_bblock_interface.insert(ir);
              residue_registered = true;
              break;
            }
          }
        }
        if(residue_registered) break;
      }
    }
  }

//  TR.Debug  << "intra-building block interface residues: " << stringify_iterable(intra_bblock_interface) << std::endl;

	// Now get the resides that are at the inter-building block interface
  TR.Debug << "filtering residues that are in contact with other residues on the interface (distance less than "<< contact_dist_ <<")" << std::endl;
  for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
    if (!indy_resis[ir]) continue;
    std::string atom_i = (pose.residue(ir).name3() == "GLY") ? "CA" : "CB";
    for (Size jr=1; jr<=sym_info->num_total_residues_without_pseudo(); jr++) {
      // skip residue jr if it is in the same building block
      if ( sym_info->subunit_index(jr) <= nsub_bblock_ ) continue;
      std::string atom_j = (pose.residue(jr).name3() == "GLY") ? "CA" : "CB";
      // Here we are filtering the residues that are in contact
      if (pose.residue(ir).xyz(atom_i).distance_squared(pose.residue(jr).xyz(atom_j)) <= contact_dist_sq) {
        TR.Debug << "residue " << ir << " in contact with residue " << jr << std::endl;
        TR.Debug << "\tinterface = " <<  intra_bblock_interface.count(ir) << std::endl;
        design_pos.insert(ir);
        break;
      }
    }
  }

	// Now combine the above two filters, and prevent_repacking at all positions that are either:
	// a) not at the inter-building block interface, or
	// b) are, but also make intra-building block contacts and are not clashing
	std::string output = "design_pos ";
	for (Size ir=1; ir<=sym_info->num_total_residues_without_pseudo(); ir++) {
    if (!indy_resis[ir]) continue;
		if ((design_pos.find(ir) != design_pos.end()) && (intra_bblock_interface.find(ir) == intra_bblock_interface.end())) {
			output += ObjexxFCL::string_of(ir)+"+";
		} else {
			TR.Debug << "resi " << ir << " will not be designed" << std::endl;
			task.nonconst_residue_task(ir).prevent_repacking();
		}
	}
	TR.Debug << output << std::endl;

}

void
BuildingBlockInterfaceOperation::parse_tag( TagPtr tag )
{
  nsub_bblock_ = tag->getOption<core::Size>("nsub_bblock", 1);
	contact_dist_ = tag->getOption<core::Real>("contact_dist", 10.0);
	bblock_dist_ = tag->getOption<core::Real>("bblock_dist", 5.0);
	fa_rep_cut_ = tag->getOption<core::Real>("fa_rep_cut", 3.0);


}

} //namespace matdes
} //namespace devel
