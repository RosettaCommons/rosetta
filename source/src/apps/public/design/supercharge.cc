// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    apps/pilot/bder/supercharge.cc
/// @brief   This protocol supercharges the surface of an input pdb with either positive or negatively charged residues.
/// @details There are two modes for supercharging.  The first is called AvNAPSA developed by the David Liu lab at Harvard.  In this approach, surface residues are defined by the Average # Neighbor Atoms Per Sidechain Atom (AvNAPSA value), with a cutoff of 150.  I think 100 is a safer cutoff.  Arg, Lys, Asp, Glu, Asn, Gln are the only residues allowed to mutated.  Lys is always chosen for positive, Glu is always chosen for negative, unless the native is Asn, then Asp is chosen.  Thus, the sequence is deterministic.  If one desires a particular net charge, the residues would be sorted from low to high AvNAPSA value and mutated one at a time until the target charge is achieved - this ignores the ceiling of 150 or 100.  The second approach uses the Rosetta score function to guide the surface mutagenesis.  The user must specifiy if Arg, Lys, or Asp, Glu are desired, and what the reference weights are.  Alternatively, the  user can specify a target net charge, and the reference weights of the charged residues will be incremented/decremented until the net charge is reached.
/// @author Bryan Der


//AvNAPSA-mode, target charge
//1. Define surface.  sort NQ and RK/DE residues by AvNAPSA value (low to high)
//2. Next residue in sorted list: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
//3. If net charge = target net charge, output pdb

//AvNAPSA-mode, no target charge
//1. Define surface by AvNAPSA value (<100 default)
//2. For each NQ and DE/RK residue in the surface: Positive: mutate DENQ-->K, Negative: mutate RKQ-->E and N-->D
//3. Output pdb

//Rosetta-mode, target charge
//1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default
// or Define surface by AvNAPSA value (<100 default)
//2. Set design task
//   read user resfile, if provided
//   dont_mutate gly, pro, cys
//   dont_mutate h-bonded sidechains
//   dont_mutate correct charge residues
//3. Set reference energies for RK/DE, starting at user input values
//4. pack rotamers mover
//5. check net charge, increment/decrement reference energies (back to step 3.)
//6. Once a pack rotamers run results in the correct net charge, output pdb

//Rosetta-mode, no target charge
//1. Define surface.  Neighbor by distance calculator (CB dist.), <16 neighbors default
// or Define surface by AvNAPSA value (<100 default)
//2. Set design task
//   read user resfile, if provided
//   dont_mutate gly, pro, cys
//   dont_mutate h-bonded sidechains
//   dont_mutate correct charge residues
//3. Set reference energies for RK/DE, using the user input values
//4. pack rotamers mover
//5. Output pdb


#include <devel/init.hh>

#include <protocols/design_opt/Supercharge.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh> // used to create a resfile
#include <utility/excn/Exceptions.hh>
#include <sstream>
#include <string>

#include <utility/vector1.hh>
#include <utility/sort_predicates.hh>
#include <cmath>


//tracers
static basic::Tracer TR( "apps.public.design.supercharge" );

using namespace core;
using Pose = core::pose::Pose;
using SizeSet = std::set<Size>;

//local options
namespace local {

//AvNAPSA-mode
basic::options::BooleanOptionKey const AvNAPSA_positive("AvNAPSA_positive");
basic::options::BooleanOptionKey const AvNAPSA_negative("AvNAPSA_negative");

basic::options::BooleanOptionKey const target_net_charge_active("target_net_charge_active"); //ideally I'd use .user() to see if the option is active, but implementation on the ROSIE server requires this separate flag
basic::options::IntegerOptionKey const target_net_charge("target_net_charge");

//AvNAPSA-mode or Rosetta-mode
basic::options::IntegerOptionKey const surface_atom_cutoff("surface_atom_cutoff"); // if target_net_charge is specified, the AvNAPSA cutoff is ignored

//Rosetta-mode (these will be ignored if AvNAPSA mode is on via AvNAPSA_positive or AvNAPSA_negative)
basic::options::IntegerOptionKey const surface_residue_cutoff("surface_residue_cutoff"); //for choosing surface residues, cannot be done in AvNAPSA mode
basic::options::BooleanOptionKey const include_arg("include_arg");
basic::options::BooleanOptionKey const include_lys("include_lys");
basic::options::BooleanOptionKey const include_asp("include_asp");
basic::options::BooleanOptionKey const include_glu("include_glu");
basic::options::RealOptionKey const refweight_arg("refweight_arg");
basic::options::RealOptionKey const refweight_lys("refweight_lys");
basic::options::RealOptionKey const refweight_asp("refweight_asp");
basic::options::RealOptionKey const refweight_glu("refweight_glu");
basic::options::BooleanOptionKey const dont_mutate_glyprocys("dont_mutate_glyprocys"); // true by default
basic::options::BooleanOptionKey const dont_mutate_correct_charge("dont_mutate_correct_charge"); // true by default
basic::options::BooleanOptionKey const dont_mutate_hbonded_sidechains("dont_mutate_hbonded_sidechains"); // true by default
basic::options::BooleanOptionKey const pre_packminpack("pre_packminpack"); // true by default
basic::options::IntegerOptionKey const nstruct("nstruct"); // custom nstruct, not used in AvNAPSA mode bc that sequence is deterministic

//AvNAPSA-mode or Rosetta-mode
basic::options::BooleanOptionKey const compare_residue_energies_all("compare_residue_energies_all");
basic::options::BooleanOptionKey const compare_residue_energies_mut("compare_residue_energies_mut");

//Note: either mode will read a user input resfile, but be sure to use ALLAA as the default, because NATAA default will make the surface residues undesignable.  Either mode will make a second resfile with NATAA as default.

}//local

/// @brief Adds charged residues to a protein surface

int main( int argc, char* argv[] )
{
	try {

		using basic::options::option;
		option.add( local::AvNAPSA_positive, "AvNAPSA positive supercharge").def(false);
		option.add( local::AvNAPSA_negative, "AvNAPSA negative supercharge").def(false);

		option.add( local::target_net_charge_active, "target net charge active").def(false);
		option.add( local::target_net_charge, "target net charge").def(0);

		option.add( local::surface_atom_cutoff, "AvNAPSA neighbor atom cutoff").def(120); // this is how AvNAPSA defines surface, can be used in the Rosetta approach
		option.add( local::surface_residue_cutoff, "cutoff for surface residues ( <= # is surface)" ).def(16);

		option.add( local::include_arg, "include arg in supercharge design").def(false);
		option.add( local::include_lys, "include lys in supercharge design").def(false);
		option.add( local::include_asp, "include asp in supercharge design").def(false);
		option.add( local::include_glu, "include glu in supercharge design").def(false);
		option.add( local::refweight_arg, "reference energy for Arg").def(-0.14916);
		option.add( local::refweight_lys, "reference energy for Lys").def(-0.287374);
		option.add( local::refweight_asp, "reference energy for Asp").def(-1.28682);
		option.add( local::refweight_glu, "reference energy for Glu").def(-1.55374);
		option.add( local::dont_mutate_glyprocys, "don't mutate gly, pro, cys").def(true);
		option.add( local::dont_mutate_correct_charge, "don't mutate correct charge").def(true);
		option.add( local::dont_mutate_hbonded_sidechains, "don't mutate hbonded sidechains").def(true);
		option.add( local::pre_packminpack, "pack-min-pack before supercharging").def(false);

		option.add( local::nstruct, "local nstruct").def(1);

		option.add( local::compare_residue_energies_all, "compare energy terms for all residues").def(false);
		option.add( local::compare_residue_energies_mut, "compare energy terms for mutated residues only").def(true);

		devel::init(argc, argv);

		protocols::design_opt::SuperchargeOP mover = utility::pointer::make_shared< protocols::design_opt::Supercharge >();

		if ( option[ local::AvNAPSA_positive ].user() ) {
			mover->AvNAPSA_positive( option[ local::AvNAPSA_positive ] );
		}
		if ( option[ local::AvNAPSA_negative ].user() ) {
			mover->AvNAPSA_negative( option[ local::AvNAPSA_negative ] );
		}

		if ( option[ local::target_net_charge_active ].user() ) {
			mover->target_net_charge_active( option[ local::target_net_charge_active ] );
		}
		if ( option[ local::target_net_charge ].user() ) {
			mover->target_net_charge( option[ local::target_net_charge ] );
		}
		if ( option[ local::surface_atom_cutoff ].user() ) {
			mover->surface_atom_cutoff( option[ local::surface_atom_cutoff ] );
		}
		if ( option[ local::surface_residue_cutoff ].user() ) {
			mover->surface_residue_cutoff( option[ local::surface_residue_cutoff ] );
		}

		if ( option[ local::include_arg ].user() ) {
			mover->include_arg( option[ local::include_arg ] );
		}
		if ( option[ local::include_lys ].user() ) {
			mover->include_lys( option[ local::include_lys ] );
		}
		if ( option[ local::include_asp ].user() ) {
			mover->include_asp( option[ local::include_asp ] );
		}
		if ( option[ local::include_glu ].user() ) {
			mover->include_glu( option[ local::include_glu ] );
		}
		if ( option[ local::refweight_arg ].user() ) {
			mover->refweight_arg( option[ local::refweight_arg ] );
		}
		if ( option[ local::refweight_lys ].user() ) {
			mover->refweight_lys( option[ local::refweight_lys ] );
		}
		if ( option[ local::refweight_asp ].user() ) {
			mover->refweight_asp( option[ local::refweight_asp ] );
		}
		if ( option[ local::refweight_glu ].user() ) {
			mover->refweight_glu( option[ local::refweight_glu ] );
		}

		if ( option[ local::dont_mutate_glyprocys ].user() ) {
			mover->dont_mutate_glyprocys( option[ local::dont_mutate_glyprocys ] );
		}
		if ( option[ local::dont_mutate_correct_charge ].user() ) {
			mover->dont_mutate_correct_charge( option[ local::dont_mutate_correct_charge ] );
		}
		if ( option[ local::dont_mutate_hbonded_sidechains ].user() ) {
			mover->dont_mutate_hbonded_sidechains( option[ local::dont_mutate_hbonded_sidechains ] );
		}
		if ( option[ local::pre_packminpack ].user() ) {
			mover->pre_packminpack( option[ local::pre_packminpack ] );
		}

		if ( option[ local::nstruct ].user() ) {
			mover->local_nstruct( option[ local::nstruct ] );
		}

		if ( option[ local::compare_residue_energies_all ].user() ) {
			mover->compare_residue_energies_all( option[ local::compare_residue_energies_all ] );
		}
		if ( option[ local::compare_residue_energies_mut ].user() ) {
			mover->compare_residue_energies_mut( option[ local::compare_residue_energies_mut ] );
		}

		protocols::jd2::JobDistributor::get_instance()->go(mover);

		TR << "************************d**o**n**e**************************************" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
