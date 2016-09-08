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

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperation.hh>
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
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh> // used to create a resfile
#include <utility/excn/Exceptions.hh>
#include <sstream>
#include <string>
#include <fstream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/sort_predicates.hh>
#include <math.h>


//tracers
using basic::Error;
using basic::Warning;
using basic::T;
static THREAD_LOCAL basic::Tracer TR( "apps.public.design.supercharge" );

using namespace core;
typedef core::pose::Pose Pose;
typedef std::set< Size > SizeSet;

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
class supercharge : public protocols::moves::Mover {
public:
	supercharge(){}
	virtual ~supercharge(){};


	virtual
	void
	apply( Pose & pose ) {
		using namespace basic::options;
		out_path_ = basic::options::option[ OptionKeys::out::path::path ]();

		//check for chain ID
		char chain = pose.pdb_info()->chain(1);
		if ( chain == ' ' ) {
			TR << "chain is whitespace, setting chain ID to 'A' " << std::endl;
			for ( Size i=1; i<=pose.size(); ++i ) {
				pose.pdb_info()->chain(i, 'A');
			}
		}


		//If the target net charge is -10, current net charge is -4, need to perform positive-supercharging
		if ( option[local::target_net_charge_active] ) {

			int current_net_charge = get_net_charge( pose );
			int target_net_charge = option[local::target_net_charge];
			int delta_charge = target_net_charge - current_net_charge;

			if ( delta_charge < 0 ) {

				if ( !option[local::include_asp] && !option[local::include_glu] && !option[local::AvNAPSA_negative] ) {
					TR << "Current charge: " << current_net_charge << ".  Target charge: " << target_net_charge << ".  Incompatible user inputs.  Cannot add negative charge with current options.  Try using the flags include_asp include_glu (Rosetta-mode) or AvNAPSA_negative (AvNAPSA-mode)." << std::endl;
					set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
					return;
				}
			} else if ( delta_charge > 0 ) {
				if ( !option[local::include_arg] && !option[local::include_lys] && !option[local::AvNAPSA_positive] ) {
					TR << "Current charge: " << current_net_charge << ".  Target charge: " << target_net_charge << ".  Incompatible user inputs.  Cannot add positive charge with current options.  Try using the flags include_arg include_lys (Rosetta-mode) or AvNAPSA_positive (AvNAPSA-mode)." << std::endl;
					set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
					return;
				}
			} else if ( delta_charge == 0 ) {
				TR << "Current charge: " << current_net_charge << ".  Target charge: " << target_net_charge << ".  Current charge equals target charge, no supercharging necessary." << std::endl;
				set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
				return;
			}
		}


		option[OptionKeys::packing::use_input_sc](true); // always use_input_sc
		prepack_input_structure( pose );

		Pose const starting_pose( pose ); //save starting pose to list mutations
		Pose native( starting_pose );

		//AvNAPSA mode.  Does not use Rosetta energy calculation, chooses mutable positions solely by number of atom neighbors per sidechain
		if ( option[local::AvNAPSA_positive] || option[local::AvNAPSA_negative] ) {
			AvNAPSA_values( pose );
			set_resfile_AvNAPSA( pose );
			design_supercharge_AvNAPSA( starting_pose, pose );
		} else {
			//Rosetta mode.  Uses Rosetta energy function to choose which residues to mutate and to what residue type.
			//score pose for hbond detection
			using namespace core::scoring;
			ScoreFunctionOP scorefxn = get_score_function();
			scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
			energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
			scorefxn->set_energy_method_options( energymethodoptions );
			scorefxn->score( pose );

			set_surface( pose );
			set_resfile( pose ); // default=NATRO, same-charge=NATAA, allow native+desired charge
			design_supercharge( starting_pose, pose ); // sets reference energies, designs the surface, outputs with an informative name
		}

		if ( option[local::compare_residue_energies_all] || option[local::compare_residue_energies_mut] ) {
			energy_comparison( native, pose );
		}

		return;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////////         BEGIN AVNAPSA MODE (average number of neighboring atoms per sidechain atom)   //////////////////////////////////////
	//////////          doesn't consider surface energetics, operates only on surface accessibility   //////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	AvNAPSA_values( Pose const & pose ) {  //average number of neighboring atoms per sidechain atom

		for ( Size i=1; i <= pose.size(); ++i ) {
			Real avnapsa_value( 9999 ); //high value will mean don't mutate, don't mutate by default

			std::string name3 = pose.residue(i).name3();

			if ( name3 == "ASP" || name3 == "GLU" || name3 == "ARG" || name3 == "LYS" || name3 == "ASN" || name3 == "GLN" ) {

				//don't mutate correct charge
				if ( basic::options::option[local::AvNAPSA_positive] ) {
					if ( name3 == "ARG" || name3 == "LYS" ) {
						AvNAPSA_values_.push_back( avnapsa_value );
						TR << "residue " << i << " is already positive" << std::endl;
						continue;
					}
				} else if ( basic::options::option[local::AvNAPSA_negative] ) {
					if ( name3 == "ASP" || name3 == "GLU" ) {
						AvNAPSA_values_.push_back( avnapsa_value );
						TR << "residue " << i << " is already negative" << std::endl;
						continue;
					}
				}

				conformation::Residue const & i_rsd( pose.residue( i ) );
				Size total_atom_neighbors_of_sidechain( 1 );

				for ( Size ii = i_rsd.first_sidechain_atom(); ii <= i_rsd.nheavyatoms(); ++ii ) {
					conformation::Atom const & ii_atom( i_rsd.atom( ii ) );
					Vector const & ii_atom_xyz = ii_atom.xyz();

					for ( Size j=1; j <= pose.size(); ++j ) {

						conformation::Residue const & j_rsd( pose.residue( j ) );


						for ( Size jj = 1; jj <= j_rsd.nheavyatoms(); ++jj ) {

							//TR << i << " " << ii << " " << j << " " << jj << " total " << total_atom_neighbors_of_sidechain << std::endl;

							conformation::Atom const & jj_atom( j_rsd.atom( jj ) );
							Vector const & jj_atom_xyz = jj_atom.xyz();

							if ( ii_atom_xyz.distance( jj_atom_xyz ) < 10.0 ) {
								++total_atom_neighbors_of_sidechain;
							}
						}
					}
				}


				Size number_of_sidechain_atoms = i_rsd.nheavyatoms() - 4; // four backbone atoms, would equal 0 if residue is Glycine and this would cause a crash
				avnapsa_value = total_atom_neighbors_of_sidechain / number_of_sidechain_atoms;

				TR << "residue: " << i << "  heavy: " << i_rsd.nheavyatoms() << "  " << "sidechain: " << i_rsd.first_sidechain_atom() << "  AvNAPSA_value: " << avnapsa_value << std::endl;

			}

			AvNAPSA_values_.push_back( avnapsa_value ); //index must equal the residue number
		}

		TR << "there are " << pose.size() << " residues and " << AvNAPSA_values_.size() << " AvNAPSA values" << std::endl;

		return;
	}


	virtual
	void
	set_resfile_AvNAPSA( Pose const & pose ) {

		TR << "Creating a resfile, it will be saved as " << out_path_ << "resfile_output_Asc" << std::endl;
		utility::io::ozstream OutputResfile;
		OutputResfile.open( out_path_ + '/' + "resfile_output_Asc" );

		OutputResfile << "NATAA" << std::endl;
		OutputResfile << "start" << std::endl;

		std::stringstream pymol_avnapsa_residues;
		utility::vector1< Size > residues_to_mutate; //will be appended either to acheive correct charge or based on AvNAPSA value cutoff

		if ( ! basic::options::option[local::target_net_charge_active] ) {
			largest_mutated_AvNAPSA_ = (Size) basic::options::option[local::surface_atom_cutoff]; // no specified net charge, largest AvNAPSA allowed equals the cutoff.  This value is used to name output PDBs.
			for ( Size i(1); i <= AvNAPSA_values_.size(); ++i ) {
				if ( AvNAPSA_values_[i] < basic::options::option[local::surface_atom_cutoff] && AvNAPSA_values_[i] != 9999 ) {
					residues_to_mutate.push_back( i );
					TR << "Mutate " << i << std::endl;
				}
			}
		} else {
			// if a target net charge is included, the AvNAPSA cutoff is ignored
			//sort residues by their avnapsa value in ascending order
			utility::vector1<std::pair<Size,Real> > pair_residue_avnapsa;
			for ( Size i(1); i<= AvNAPSA_values_.size(); ++i ) {
				if ( AvNAPSA_values_[i] != 9999 ) {
					pair_residue_avnapsa.push_back( std::pair<Size,Real>(i,AvNAPSA_values_[i]) );
				}
			}
			std::sort(pair_residue_avnapsa.begin(), pair_residue_avnapsa.end(), utility::SortSecond< Size, Real >() );


			//append residues to mutate until the resulting net charge would be correct +- 1
			int net_charge = get_net_charge( pose );
			TR << "Starting net charge is " << net_charge << " and target_net_charge is " << basic::options::option[local::target_net_charge] << std::endl;

			for ( Size i(1); i <= pair_residue_avnapsa.size(); ++i ) {
				Size this_res = pair_residue_avnapsa[i].first;
				std::string name3 = pose.residue(this_res).name3();

				if ( basic::options::option[local::AvNAPSA_positive] && net_charge >= basic::options::option[local::target_net_charge] ) {
					break; // if positive enough
				} else if ( basic::options::option[local::AvNAPSA_negative] && net_charge <= basic::options::option[local::target_net_charge] ) {
					break; // if negative enough
				} else if ( pair_residue_avnapsa[i].second != 9999 ) { // if not charged enough, add another residue to mutate
					residues_to_mutate.push_back(this_res);
					TR << "Mutate " << this_res << std::endl;
				}

				if ( basic::options::option[local::AvNAPSA_positive] ) {
					if ( name3 == "ASP" || name3 == "GLU" ) {
						net_charge += 2;
					} else if ( name3 == "ARG" || name3 == "LYS" ) {
						net_charge += 0;
					} else {
						net_charge += 1;
					}
				}
				if ( basic::options::option[local::AvNAPSA_negative] ) {
					if ( name3 == "ARG" || name3 == "LYS" ) {
						net_charge -= 2;
					} else if ( name3 == "ASP" || name3 == "GLU" ) {
						net_charge -= 0;
					} else {
						net_charge -= 1;
					}
				}
				//net charge will only be correct plus/minus 1
				//pose net charge will actually change upon pack_rotamers->apply

			}

			Size number_to_mutate = residues_to_mutate.size();
			largest_mutated_AvNAPSA_ = (Size) pair_residue_avnapsa[number_to_mutate].second; //for naming output PDBs
		}


		////////// if no target net charge //////////////
		for ( Size i(1); i <= residues_to_mutate.size(); ++i ) {
			Size j = residues_to_mutate[i];
			pymol_avnapsa_residues << j << "+";
			if ( basic::options::option[local::AvNAPSA_positive] ) {
				OutputResfile << "   " << pose.pdb_info()->pose2pdb( j ) /*prints resnum and chain*/ << "  PIKAA  K" << std::endl; //lysine only, see SuppInfo of Lawrence et al. 2007
			} else if ( basic::options::option[local::AvNAPSA_negative] ) {
				std::string resname = pose.residue(j).name3();
				if ( resname == "ASN" ) {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb( j ) /*prints resnum and chain*/ << "  PIKAA  D" << std::endl; //Asp only if Asn is native
				} else {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb( j ) /*prints resnum and chain*/ << "  PIKAA  E" << std::endl; //Glu by default
				}
			}
		}

		TR << "PYMOL_SELECT AVNAPSA residues: " << pymol_avnapsa_residues.str() << std::endl;

		return;
	}


	virtual
	void
	design_supercharge_AvNAPSA( Pose const & starting_pose, Pose & pose ) { //there are no choices, either NATRO or the only PIKAA residue from the AvNAPSA resfile

		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace basic::options;
		TaskFactoryOP task_factory( new TaskFactory() );
		task_factory->push_back(TaskOperationCOP( new operation::InitializeFromCommandline() )); //ex1, ex1, minimize sidechains, use_input_sc

		// first, read in user resfile.  Intended to only contain NATAA or NATRO to specify additional residues to not mutation  MUST have ALLAA as default!!
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			task_factory->push_back( TaskOperationCOP( new operation::ReadResfile ) );
			TR << "Reading resfile from user input... make sure ALLAA is set as the default in your resfile!!" << std::endl;
		}

		task_factory->push_back( TaskOperationCOP( new operation::ReadResfile( out_path_ + '/' + "resfile_output_Asc" ) ) ); // reads the resfile previously created, adds to the user resfile (if provided)

		using namespace core::scoring;
		ScoreFunctionOP scorefxn = get_score_function();
		scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
		energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options( energymethodoptions );
		scorefxn_ = scorefxn;

		utility::file::FileName input_pdbname( pose.pdb_info()->name() );
		std::string input_namebase( input_pdbname.base() );

		protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover );
		packrot_mover->score_function( scorefxn );
		packrot_mover->task_factory( task_factory );
		packrot_mover->apply( pose ); /////////////////APPLY PACKROT MOVER
		packrot_mover->apply( pose ); /////////////////APPLY PACKROT MOVER
		packrot_mover->apply( pose ); /////////////////APPLY PACKROT MOVER


		std::string pos_or_neg;
		if ( option[local::AvNAPSA_positive] ) {
			pos_or_neg = "pos";
		} else {
			pos_or_neg = "neg";
		}

		int final_net_charge = get_net_charge( pose );
		std::stringstream name_info;
		name_info << "_AvNAPSA_" << pos_or_neg << "_" << largest_mutated_AvNAPSA_ << "_" << final_net_charge;

		outputname_ = input_namebase + name_info.str() + ".pdb";
		TR << "NEW NAME FOR SUPERCHARGED OUTPUT: " << outputname_ << std::endl;

		scorefxn->score( pose );

		pose.dump_scored_pdb( out_path_ + '/' + outputname_, *scorefxn );

		print_netcharge_and_mutations(starting_pose, pose );


		return;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////         BEGIN ROSETTA MODE               /////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	virtual
	void
	prepack_input_structure( Pose & pose ) {

		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace basic::options;
		TaskFactoryOP task_factory( new TaskFactory() );
		task_factory->push_back(TaskOperationCOP( new operation::InitializeFromCommandline() )); //use_input_sc

		using namespace core::scoring;
		ScoreFunctionOP scorefxn = get_score_function();
		scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
		energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options( energymethodoptions );
		scorefxn_ = scorefxn;

		protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover );
		packrot_mover->score_function( scorefxn_ );

		core::pack::task::operation::RestrictToRepackingOP restrict_to_repack( new core::pack::task::operation::RestrictToRepacking() );
		task_factory->push_back( restrict_to_repack );

		packrot_mover->task_factory( task_factory );

		kinematics::MoveMapOP movemap_sc( new kinematics::MoveMap );
		kinematics::MoveMapOP movemap_scbb( new kinematics::MoveMap );
		movemap_sc->set_chi( true );
		movemap_sc->set_bb( false );
		movemap_scbb->set_chi( true );
		movemap_scbb->set_bb( true );

		protocols::simple_moves::MinMoverOP min_sc( new protocols::simple_moves::MinMover( movemap_sc, scorefxn_, "lbfgs_armijo", 0.01, true ) );
		protocols::simple_moves::MinMoverOP min_scbb( new protocols::simple_moves::MinMover( movemap_scbb, scorefxn_, "lbfgs_armijo", 0.01, true ) );


		TR << "Packrotamers" << std::endl;
		packrot_mover->apply( pose );
		packrot_mover->apply( pose );
		packrot_mover->apply( pose );

		if ( basic::options::option[local::pre_packminpack] ) {
			TR << "Minimizing sidechains..." << std::endl;
			min_sc->apply( pose );
			TR << "Packrotamers" << std::endl;
			packrot_mover->apply( pose );
			packrot_mover->apply( pose );
			TR << "Minimizing sidechains and backbone..." << std::endl;
			min_sc->apply( pose );
			min_scbb->apply( pose );
			TR << "Packrotamers" << std::endl;
			packrot_mover->apply( pose );
			packrot_mover->apply( pose );
			TR << "Minimizing sidechains and backbone..." << std::endl;
			min_sc->apply( pose );
			min_scbb->apply( pose );
			TR << "Packrotamers" << std::endl;
			packrot_mover->apply( pose );
			packrot_mover->apply( pose );
		}

		return;
	}


	void set_surface( Pose const & pose ){

		//define surface by residue neighbors
		if ( ! basic::options::option[local::surface_atom_cutoff].user() ) {
			// registering the calculators (in this way will not allow nstruct > 1)
			Size biggest_calc(0);
			std::string const calc_stem("nbr_dist_calc_");
			std::ostringstream calcname;
			for ( Size res(1); res <= pose.size(); ++res ) {
				if ( biggest_calc < res ) { //create calculator
					calcname << calc_stem << res;
					if ( pose::metrics::CalculatorFactory::Instance().check_calculator_exists( calcname.str() ) ) {
						basic::Warning() << "Calculator " << calcname.str() << " already exists, this is hopefully correct for your purposes" << std::endl;
					} else {
						using pose::metrics::PoseMetricCalculatorOP;
						pose::metrics::CalculatorFactory::Instance().register_calculator( calcname.str(), PoseMetricCalculatorOP( new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(res) ) );
					}
					calcname.str("");
					++biggest_calc;
				}
			}
			//find surface residues (number of neighbors cutoff)


			basic::MetricValue< Size > num_n; // number of neighbors
			TR << pose.pdb_info()->name() << std::endl;
			for ( Size i(1); i<=pose.size(); ++i ) {
				calcname << calc_stem << i;
				pose.metric( calcname.str(), "num_neighbors", num_n);
				calcname.str("");

				//TR << "residue " << i << " num_neighbors " << num_n.value() << std::endl;

				if ( num_n.value() <= core::Size(basic::options::option[local::surface_residue_cutoff].value()) ) {
					TR << "adding " << i << " to surface set" << std::endl;
					surface_res_.insert(i);
				}
			}
		} else {
			//define surface by atom neighbors

			for ( Size i=1; i <= pose.size(); ++i ) {
				Real avnapsa_value( 9999 ); //high value will mean don't mutate, don't mutate by default

				std::string name3 = pose.residue(i).name3();

				if ( name3 == "GLY" ) {
					continue;
				}


				conformation::Residue const & i_rsd( pose.residue( i ) );
				Size total_atom_neighbors_of_sidechain( 1 );

				for ( Size ii = i_rsd.first_sidechain_atom(); ii <= i_rsd.nheavyatoms(); ++ii ) {
					conformation::Atom const & ii_atom( i_rsd.atom( ii ) );
					Vector const & ii_atom_xyz = ii_atom.xyz();

					for ( Size j=1; j <= pose.size(); ++j ) {

						conformation::Residue const & j_rsd( pose.residue( j ) );


						for ( Size jj = 1; jj <= j_rsd.nheavyatoms(); ++jj ) {

							conformation::Atom const & jj_atom( j_rsd.atom( jj ) );
							Vector const & jj_atom_xyz = jj_atom.xyz();

							if ( ii_atom_xyz.distance( jj_atom_xyz ) < 10.0 ) {
								++total_atom_neighbors_of_sidechain;
							}
						}
					}
				}


				Size number_of_sidechain_atoms = i_rsd.nheavyatoms() - 4; // four backbone atoms, would equal 0 if residue is Glycine and this would cause a crash
				avnapsa_value = total_atom_neighbors_of_sidechain / number_of_sidechain_atoms;

				//TR << "residue: " << i << "  heavy: " << i_rsd.nheavyatoms() << "  " << "sidechain: " << i_rsd.first_sidechain_atom() << "  AvNAPSA_value: " << avnapsa_value << std::endl;

				AvNAPSA_values_.push_back( avnapsa_value ); //index must equal the residue number
				TR << "AvNAPSA score for residue " << i << "  " << avnapsa_value << std::endl;

				for ( Size i(1); i<= AvNAPSA_values_.size(); ++i ) {
					if ( AvNAPSA_values_[i] <= basic::options::option[local::surface_atom_cutoff] ) { //every residue has an AvNAPSA value, so i equals residue number
						surface_res_.insert(i);
					}
				}
			}

		} //else


		return;
	}


	void set_resfile( Pose const & pose ){
		using namespace basic::options;

		TR << "Creating a resfile, it will be saved as " << out_path_ << "resfile_output_Rsc" << std::endl;
		utility::io::ozstream OutputResfile;
		OutputResfile.open( out_path_ + '/' + "resfile_output_Rsc" );

		OutputResfile << "NATAA" << std::endl;
		OutputResfile << "start" << std::endl;


		//compute set of hbonds
		core::scoring::hbonds::HBondSet hbond_set;
		hbond_set.setup_for_residue_pair_energies( pose, false, false );

		for ( SizeSet::const_iterator it(surface_res_.begin()), end(surface_res_.end()); it!=end; ++it ) {

			char NATAA_oneletter = pose.residue(*it).name1();

			//gly, pro, and cys are specialized residues that might be better of unchanged
			if ( option[local::dont_mutate_glyprocys] ) {
				if ( NATAA_oneletter == 'G' || NATAA_oneletter == 'P' || NATAA_oneletter == 'C' ) {
					continue;
				}
			}

			//if same-charge, leave as NATAA
			if ( option[local::dont_mutate_correct_charge] ) {
				if ( option[local::include_arg] && NATAA_oneletter == 'R' ) {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATAA  #same charge" << std::endl;
					continue;
				}
				if ( option[local::include_lys] && NATAA_oneletter == 'K' ) {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATAA  #same charge" << std::endl;
					continue;
				}
				if ( option[local::include_asp] && NATAA_oneletter == 'D' ) {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATAA  #same charge" << std::endl;
					continue;
				}
				if ( option[local::include_glu] && NATAA_oneletter == 'E' ) {
					OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATAA  #same charge" << std::endl;
					continue;
				}
			}

			//dont_mutate strong sidechain hbonds
			if ( option[local::dont_mutate_hbonded_sidechains] ) {
				//hbond detection
				bool found_sc_hbond( false );

				for ( Size i = 1; i<= hbond_set.nhbonds(); i++ ) {
					core::scoring::hbonds::HBondCOP hbond(hbond_set.hbond(i).get_self_ptr());
					if ( hbond->energy() > -0.5 ) { // a fun semi-arbitrary value for hbond strength cutoff
						continue;
					}
					if ( hbond->don_res() == *it && !hbond->don_hatm_is_backbone() ) {

						OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATRO  #has sc hbond energy=" << hbond->energy() << std::endl;

						found_sc_hbond = true;
						break;
					}
					if ( hbond->acc_res() == *it && !hbond->acc_atm_is_protein_backbone() ) {
						OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  NATRO  #has sc hbond energy=" << hbond->energy() << std::endl;
						found_sc_hbond = true;
						break;
					}
				}
				if ( found_sc_hbond ) {
					continue;
				}
			}


			//if there hasn't been a continue, it's a mutable position
			utility::vector1<char> PIKAA_residues;
			PIKAA_residues.push_back(NATAA_oneletter); // the native amino acid is allowed
			if ( option[local::include_arg].value() ) {
				PIKAA_residues.push_back('R');
			}
			if ( option[local::include_lys].value() ) {
				PIKAA_residues.push_back('K');
			}
			if ( option[local::include_asp].value() ) {
				PIKAA_residues.push_back('D');
			}
			if ( option[local::include_glu].value() ) {
				PIKAA_residues.push_back('E');
			}

			// print lines in the resfile uch as:  3  A  PIKAA NRK
			OutputResfile << "   " << pose.pdb_info()->pose2pdb(*it) /*prints resnum and chain*/ << "  PIKAA ";
			for ( core::Size j(1); j<=PIKAA_residues.size(); ++j ) {
				OutputResfile << PIKAA_residues[j];
			}
			OutputResfile << std::endl;

		}//iterate over surface residues

		return;
	}


	utility::vector1< Real >
	set_reference_energies(){
		using namespace basic::options;
		utility::vector1< Real > ref_weights;
		ref_weights.push_back(0.16);  //A 1
		ref_weights.push_back(1.7);   //C 2
		ref_weights.push_back( option[local::refweight_asp] ); //D 3 (default=-0.67)
		ref_weights.push_back( option[local::refweight_glu] ); //E 4 (default=-0.81)
		ref_weights.push_back(0.63);  //F 5
		ref_weights.push_back(-0.17); //G 6
		ref_weights.push_back(0.56);  //H 7
		ref_weights.push_back(0.24);  //I 8
		ref_weights.push_back( option[local::refweight_lys] ); //K 9 (default=-0.65)
		ref_weights.push_back(-0.1);  //L 10
		ref_weights.push_back(-0.34); //M 11
		ref_weights.push_back(-0.89); //N 12
		ref_weights.push_back(0.02);  //P 13
		ref_weights.push_back(-0.97); //Q 14
		ref_weights.push_back( option[local::refweight_arg] ); //R 15 (default=-0.98)
		ref_weights.push_back(-0.37); //S 16
		ref_weights.push_back(-0.27); //T 17
		ref_weights.push_back(0.29);  //V 18
		ref_weights.push_back(0.91);  //W 19
		ref_weights.push_back(0.51);  //Y 20
		//0.16 1.7 -0.67 -0.81 0.63 -0.17 0.56 0.24 -0.65 -0.1 -0.34 -0.89 0.02 -0.97 -0.98 -0.37 -0.27 0.29 0.91 0.51

		return ref_weights;
	}

	void
	design_supercharge( Pose const & starting_pose, Pose & pose ){

		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace basic::options;
		TaskFactoryOP task_factory( new TaskFactory() );
		task_factory->push_back(TaskOperationCOP( new operation::InitializeFromCommandline() )); //ex1, ex1, minimize sidechains, use_input_sc

		// first, read in user resfile.  Intended to only contain NATAA or NATRO to specify additional residues to not mutation  MUST have ALLAA as default!!
		if ( option[ OptionKeys::packing::resfile ].user() ) {
			task_factory->push_back( TaskOperationCOP( new operation::ReadResfile ) );
			TR << "Reading resfile from user input... make sure ALLAA is set as the default in your resfile!!" << std::endl;
		}

		task_factory->push_back( TaskOperationCOP( new operation::ReadResfile( out_path_ + '/' + "resfile_output_Rsc" ) ) ); // reads the resfile previously created, adds to the user resfile (if provided)

		using namespace core::scoring;
		ScoreFunctionOP scorefxn = get_score_function();
		scoring::methods::EnergyMethodOptions energymethodoptions( scorefxn->energy_method_options() );
		energymethodoptions.hbond_options().decompose_bb_hb_into_pair_energies(true);
		scorefxn->set_energy_method_options( energymethodoptions );

		ScoreFunctionOP customref_scorefxn = get_score_function();
		scoring::methods::EnergyMethodOptions energymethodoptions_customref( customref_scorefxn->energy_method_options() );
		energymethodoptions_customref.hbond_options().decompose_bb_hb_into_pair_energies(true);
		customref_scorefxn->set_energy_method_options( energymethodoptions_customref );
		utility::vector1< Real > custom_ref_weights = set_reference_energies(); // SET REFERENCE ENERGIES function
		customref_scorefxn->set_method_weights( ref, custom_ref_weights );
		TR << "Using SCORE12 with custom reference weights\n" << *customref_scorefxn << std::flush;
		scorefxn_ = scorefxn;

		protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover );
		packrot_mover->score_function( customref_scorefxn );
		packrot_mover->task_factory( task_factory );

		//Pose pre_design( pose );
		//native_ = pre_design;


		//if a target net charge is given as an option, iterate between packrot and incrementing refweights until target charge is acheived
		if ( option[local::target_net_charge_active] ) {

			int net_charge_target = option[local::target_net_charge];
			int charge_diff = std::abs( get_net_charge(pose) - net_charge_target );


			Real refweight_max_absvalue(3.2);
			bool refweight_under_max( true );
			Real refweight_increment(0.10);
			Size counter(0);
			while ( charge_diff > 1 && counter < 40 && refweight_under_max ) {

				int net_charge = get_net_charge( pose );

				//fine-tunes the changes to refweight depending on the charge difference
				if ( std::abs(net_charge - net_charge_target) > 10 )     { refweight_increment = 0.5; }
				else if ( std::abs(net_charge - net_charge_target) > 2 ) { refweight_increment = 0.1; }
				else if ( std::abs(net_charge - net_charge_target) < 2 ) { refweight_increment = 0.02;}


				//not positive enough
				if ( net_charge < net_charge_target && option[local::include_arg] ) {
					custom_ref_weights[15] = custom_ref_weights[15] - refweight_increment;
				}
				if ( net_charge < net_charge_target && option[local::include_lys] ) {
					custom_ref_weights[9] = custom_ref_weights[9] - refweight_increment;
				}
				//too positive
				if ( net_charge > net_charge_target && option[local::include_arg] ) {
					custom_ref_weights[15] = custom_ref_weights[15] + refweight_increment;
				}
				if ( net_charge > net_charge_target && option[local::include_lys] ) {
					custom_ref_weights[9] = custom_ref_weights[9] + refweight_increment;
				}

				//not negative enough
				if ( net_charge > net_charge_target && option[local::include_asp] ) {
					custom_ref_weights[3] = custom_ref_weights[3] - refweight_increment;
					//TR << "TEST " << refweight_increment << std::endl;
				}
				if ( net_charge > net_charge_target && option[local::include_glu] ) {
					custom_ref_weights[4] = custom_ref_weights[4] - refweight_increment;
				}
				//too negative
				if ( net_charge < net_charge_target && option[local::include_asp] ) {
					custom_ref_weights[3] = custom_ref_weights[3] + refweight_increment;
				}
				if ( net_charge < net_charge_target && option[local::include_glu] ) {
					custom_ref_weights[4] = custom_ref_weights[4] + refweight_increment;
				}

				customref_scorefxn->set_method_weights( ref, custom_ref_weights );
				packrot_mover->score_function( customref_scorefxn );

				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////  DESIGN STEP  //////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////

				TR << "Supercharging the protein surface... current charge: " << net_charge << "  target charge: " << net_charge_target << std::endl;
				TR << "Refweights D E K R " << custom_ref_weights[3] << " " << custom_ref_weights[4] << " " << custom_ref_weights[9] << " " << custom_ref_weights[15] << std::endl;
				packrot_mover->apply( pose );

				charge_diff = std::abs( get_net_charge(pose) - net_charge_target );

				counter++;

				if ( fabs(custom_ref_weights[15]) > refweight_max_absvalue || fabs(custom_ref_weights[9]) > refweight_max_absvalue || fabs(custom_ref_weights[3]) > refweight_max_absvalue || fabs(custom_ref_weights[4]) > refweight_max_absvalue ) {
					refweight_under_max = false;
				}

			}

		}


		/////////////// rename pdb for output ////////////////////////////////////////

		utility::file::FileName input_pdbname( pose.pdb_info()->name() );
		std::string input_namebase( input_pdbname.base() );
		std::string include_RKDE = "_";
		std::string weights_RKDE = "_";

		if ( option[local::include_arg] ) {
			include_RKDE = include_RKDE + "R";
			std::string s;
			std::stringstream out;
			out << std::setprecision(2) << std::fixed << custom_ref_weights[15];
			//out << option[local::refweight_arg];
			s = out.str();
			if ( custom_ref_weights[15] > 0.0 ) { s = "+" + s; }
			weights_RKDE = weights_RKDE + s + "_";
		}
		if ( option[local::include_lys] ) {
			include_RKDE = include_RKDE + "K";
			std::string s;
			std::stringstream out;
			out << std::setprecision(2) << std::fixed << custom_ref_weights[9];
			//out << option[local::refweight_lys];
			s = out.str();
			if ( custom_ref_weights[9] > 0.0 ) { s = "+" + s; }
			weights_RKDE = weights_RKDE + s + "_";
		}
		if ( option[local::include_asp] ) {
			include_RKDE = include_RKDE + "D";
			std::string s;
			std::stringstream out;
			out << std::setprecision(2) << std::fixed << custom_ref_weights[3];
			//out << option[local::refweight_asp];
			s = out.str();
			if ( custom_ref_weights[3] > 0.0 ) { s = "+" + s; }
			weights_RKDE = weights_RKDE + s + "_";
		}
		if ( option[local::include_glu] ) {
			include_RKDE = include_RKDE + "E";
			std::string s;
			std::stringstream out;
			out << std::setprecision(2) << std::fixed << custom_ref_weights[4];
			//out << option[local::refweight_glu];
			s = out.str();
			if ( custom_ref_weights[4] > 0.0 ) { s = "+" + s; }
			weights_RKDE = weights_RKDE + s + "_";
		}


		std::stringstream ss_i;
		std::string i_string;

		if ( ! option[local::target_net_charge_active] ) {

			Size nstruct = (Size) option[local::nstruct].value();
			for ( Size i=1; i <= nstruct; ++i ) {
				ss_i << i;
				if ( i < 10 ) { i_string = "000" + ss_i.str(); }
				else if ( i < 100 ) { i_string = "00" + ss_i.str(); }
				else if ( i < 1000 ) { i_string = "0" + ss_i.str(); }
				else { i_string = ss_i.str(); }
				ss_i.str(""); // clear stringstream

				//////////////////////////////////////////////////////////////////////////////
				/////////////////////////////  DESIGN STEP  //////////////////////////////////
				//////////////////////////////////////////////////////////////////////////////

				TR << "Supercharging the protein surface... nstruct = " << i << std::endl;
				packrot_mover->apply( pose );
				packrot_mover->apply( pose );
				packrot_mover->apply( pose );

				outputname_ = input_namebase + include_RKDE + weights_RKDE + i_string + ".pdb";
				TR << "NEW NAME FOR SUPERCHARGED OUTPUT: " << outputname_ << std::endl;

				scorefxn->score( pose );

				pose.dump_scored_pdb( out_path_ + '/' + outputname_, *scorefxn );  //score with regular-weighted scorefunction

				print_netcharge_and_mutations(starting_pose, pose );


			}//for nstruct

		} else { //don't do nstruct if specifying a target charge
			//apply packrot mover was already done until converging on target charge
			ss_i << get_net_charge( pose );
			i_string = ss_i.str();

			outputname_ = input_namebase + include_RKDE + weights_RKDE + i_string + ".pdb";
			TR << "NEW NAME FOR SUPERCHARGED OUTPUT: " << outputname_ << std::endl;

			scorefxn->score( pose );

			pose.dump_scored_pdb( out_path_ + '/' + outputname_, *scorefxn );  //score with regular-weighted scorefunction

			print_netcharge_and_mutations(starting_pose, pose );

		}


	}

	void
	print_netcharge_and_mutations( Pose const & starting_pose, Pose const & pose ){
		// net charge
		Size num_arg = 0;
		Size num_lys = 0;
		Size num_asp = 0;
		Size num_glu = 0;

		for ( Size i=1; i<=pose.size(); i++ ) {
			if ( pose.residue(i).name3() == "ARG" )      { num_arg++; }
			else if ( pose.residue(i).name3() == "LYS" ) { num_lys++; }
			else if ( pose.residue(i).name3() == "ASP" ) { num_asp++; }
			else if ( pose.residue(i).name3() == "GLU" ) { num_glu++; }
		}
		int net_charge = num_arg + num_lys - num_asp - num_glu;
		TR << outputname_ << "  Net Charge of input structure = " << get_net_charge( starting_pose ) << std::endl;
		TR << outputname_ << "  R=" << num_arg << ", K=" << num_lys << ", D=" << num_asp << ", E=" << num_glu << std::endl;
		TR << outputname_ << "  Net Charge = " << net_charge << std::endl;

		Size num_mutations = 0;
		std::string pymol_sel_mutations = "select ";
		TR << outputname_ << "  Mutations: ";
		for ( Size i=1; i<=pose.size(); i++ ) {
			//TR << "name1: " << starting_pose.residue(i).name1() << "  " << pose.residue(i).name1() << std::endl;
			if ( pose.residue(i).name1() != starting_pose.residue(i).name1() ) {
				num_mutations++;
				Size whitespace = pose.pdb_info()->pose2pdb(i).find(" ");
				std::string resnum_only = pose.pdb_info()->pose2pdb(i).substr(0, whitespace);
				TR << starting_pose.residue(i).name1() << resnum_only << pose.residue(i).name1() << ", ";
				pymol_sel_mutations = pymol_sel_mutations + resnum_only + "+";
			}
		}
		TR << std::endl;
		TR << outputname_ << "  # of mutations: " << num_mutations << std::endl;
		TR << outputname_ << "  " << pymol_sel_mutations << std::endl;

	}


	int
	get_net_charge( Pose const & pose ) {
		int net_charge(0);
		for ( Size i(1); i<=pose.size(); ++i ) {
			std::string name3 = pose.residue(i).name3();
			if ( name3 == "ARG" || name3 == "LYS" ) {
				net_charge++;
			} else if ( name3 == "ASP" || name3 == "GLU" ) {
				net_charge--;
			}
		}
		return net_charge;
	}


	void
	energy_comparison( Pose & native, Pose & pose ) {

		TR << "Comparing energies of native and designed" << std::endl;
		assert(native.size() == pose.size());

		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		using namespace basic::options;
		TaskFactoryOP task_factory( new TaskFactory() );
		task_factory->push_back(TaskOperationCOP( new operation::InitializeFromCommandline() )); //need for use_input_sc
		core::pack::task::operation::RestrictToRepackingOP restrict_to_repack( new core::pack::task::operation::RestrictToRepacking() );
		task_factory->push_back( restrict_to_repack );

		protocols::simple_moves::PackRotamersMoverOP packrot_mover( new protocols::simple_moves::PackRotamersMover );
		packrot_mover->score_function( scorefxn_ );
		packrot_mover->task_factory( task_factory );

		//kinematics::MoveMapOP movemap_sc = new kinematics::MoveMap;
		//kinematics::MoveMapOP movemap_scbb = new kinematics::MoveMap;
		//movemap_sc->set_chi( true );
		//movemap_scbb->set_chi( true );
		//movemap_scbb->set_bb( true );

		//protocols::simple_moves::MinMoverOP min_sc = new protocols::simple_moves::MinMover( movemap_sc, scorefxn_, "lbfgs_armijo", 0.01, true );
		//protocols::simple_moves::MinMoverOP min_scbb = new protocols::simple_moves::MinMover( movemap_scbb, scorefxn_, "lbfgs_armijo", 0.01, true );

		TR << "Packing designed" << std::endl;
		packrot_mover->apply( pose );
		packrot_mover->apply( pose );
		packrot_mover->apply( pose );
		TR << "Packing native" << std::endl;
		packrot_mover->apply( native );
		packrot_mover->apply( native );
		packrot_mover->apply( native );

		//min_sc->apply( pose );
		//min_scbb->apply( pose );
		//packrot_mover->apply( pose );


		TR << "Scoring native and design" << std::endl;
		scorefxn_->score(native);
		scorefxn_->score(pose);

		//native.dump_scored_pdb("native.pdb", *scorefxn_);
		//pose.dump_scored_pdb("pose.pdb", *scorefxn_);

		//talaris2013 score terms
		//total fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama omega fa_dun p_aa_pp ref

		utility::vector1< Real > total_native;
		utility::vector1< Real > fa_atr_native;
		utility::vector1< Real > fa_rep_native;
		utility::vector1< Real > fa_sol_native;
		utility::vector1< Real > fa_intra_rep_native;
		utility::vector1< Real > fa_elec_native;
		utility::vector1< Real > pro_close_native;
		utility::vector1< Real > hbond_sr_bb_native;
		utility::vector1< Real > hbond_lr_bb_native;
		utility::vector1< Real > hbond_bb_sc_native;
		utility::vector1< Real > hbond_sc_native;
		utility::vector1< Real > dslf_fa13_native;
		utility::vector1< Real > rama_native;
		utility::vector1< Real > omega_native;
		utility::vector1< Real > fa_dun_native;
		utility::vector1< Real > p_aa_pp_native;
		utility::vector1< Real > ref_native;

		utility::vector1< Real > total;
		utility::vector1< Real > fa_atr;
		utility::vector1< Real > fa_rep;
		utility::vector1< Real > fa_sol;
		utility::vector1< Real > fa_intra_rep;
		utility::vector1< Real > fa_elec;
		utility::vector1< Real > pro_close;
		utility::vector1< Real > hbond_sr_bb;
		utility::vector1< Real > hbond_lr_bb;
		utility::vector1< Real > hbond_bb_sc;
		utility::vector1< Real > hbond_sc;
		utility::vector1< Real > dslf_fa13;
		utility::vector1< Real > rama;
		utility::vector1< Real > omega;
		utility::vector1< Real > fa_dun;
		utility::vector1< Real > p_aa_pp;
		utility::vector1< Real > ref;

		utility::vector1< Real > diff_total;
		utility::vector1< Real > diff_fa_atr;
		utility::vector1< Real > diff_fa_rep;
		utility::vector1< Real > diff_fa_sol;
		utility::vector1< Real > diff_fa_intra_rep;
		utility::vector1< Real > diff_fa_elec;
		utility::vector1< Real > diff_pro_close;
		utility::vector1< Real > diff_hbond_sr_bb;
		utility::vector1< Real > diff_hbond_lr_bb;
		utility::vector1< Real > diff_hbond_bb_sc;
		utility::vector1< Real > diff_hbond_sc;
		utility::vector1< Real > diff_dslf_fa13;
		utility::vector1< Real > diff_rama;
		utility::vector1< Real > diff_omega;
		utility::vector1< Real > diff_fa_dun;
		utility::vector1< Real > diff_p_aa_pp;
		utility::vector1< Real > diff_ref;


		utility::vector1< Size > resnums;

		//talaris2013 score terms
		//total fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama omega fa_dun p_aa_pp ref

		TR << "Initialized vectors" << std::endl;
		TR << "Adding native energies" << std::endl;
		for ( Size i(1); i <= native.size(); ++i ) {

			total_native.push_back(        native.energies().residue_total_energy(i) );
			fa_atr_native.push_back(       native.energies().residue_total_energies(i)[scoring::fa_atr] );
			fa_rep_native.push_back(       native.energies().residue_total_energies(i)[scoring::fa_rep] );
			fa_sol_native.push_back(       native.energies().residue_total_energies(i)[scoring::fa_sol] );
			fa_intra_rep_native.push_back( native.energies().residue_total_energies(i)[scoring::fa_intra_rep] );
			fa_elec_native.push_back(    native.energies().residue_total_energies(i)[scoring::fa_elec] );
			pro_close_native.push_back(    native.energies().residue_total_energies(i)[scoring::pro_close] );
			hbond_sr_bb_native.push_back(  native.energies().residue_total_energies(i)[scoring::hbond_sr_bb] );
			hbond_lr_bb_native.push_back(  native.energies().residue_total_energies(i)[scoring::hbond_lr_bb] );
			hbond_bb_sc_native.push_back(  native.energies().residue_total_energies(i)[scoring::hbond_bb_sc] );
			hbond_sc_native.push_back(     native.energies().residue_total_energies(i)[scoring::hbond_sc] );
			dslf_fa13_native.push_back(    native.energies().residue_total_energies(i)[scoring::dslf_fa13] );
			rama_native.push_back(         native.energies().residue_total_energies(i)[scoring::rama] );
			omega_native.push_back(        native.energies().residue_total_energies(i)[scoring::omega] );
			fa_dun_native.push_back(       native.energies().residue_total_energies(i)[scoring::fa_dun] );
			p_aa_pp_native.push_back(      native.energies().residue_total_energies(i)[scoring::p_aa_pp] );
			ref_native.push_back(          native.energies().residue_total_energies(i)[scoring::ref] );
		}

		TR << "Adding designed energies" << std::endl;
		for ( Size i(1); i <= pose.size(); ++i ) {

			total.push_back(        pose.energies().residue_total_energy(i) );
			fa_atr.push_back(       pose.energies().residue_total_energies(i)[scoring::fa_atr] );
			fa_rep.push_back(       pose.energies().residue_total_energies(i)[scoring::fa_rep] );
			fa_sol.push_back(       pose.energies().residue_total_energies(i)[scoring::fa_sol] );
			fa_intra_rep.push_back( pose.energies().residue_total_energies(i)[scoring::fa_intra_rep] );
			fa_elec.push_back(    pose.energies().residue_total_energies(i)[scoring::fa_elec] );
			pro_close.push_back(    pose.energies().residue_total_energies(i)[scoring::pro_close] );
			hbond_sr_bb.push_back(  pose.energies().residue_total_energies(i)[scoring::hbond_sr_bb] );
			hbond_lr_bb.push_back(  pose.energies().residue_total_energies(i)[scoring::hbond_lr_bb] );
			hbond_bb_sc.push_back(  pose.energies().residue_total_energies(i)[scoring::hbond_bb_sc] );
			hbond_sc.push_back(     pose.energies().residue_total_energies(i)[scoring::hbond_sc] );
			dslf_fa13.push_back(    pose.energies().residue_total_energies(i)[scoring::dslf_fa13] );
			rama.push_back(         pose.energies().residue_total_energies(i)[scoring::rama] );
			omega.push_back(        pose.energies().residue_total_energies(i)[scoring::omega] );
			fa_dun.push_back(       pose.energies().residue_total_energies(i)[scoring::fa_dun] );
			p_aa_pp.push_back(      pose.energies().residue_total_energies(i)[scoring::p_aa_pp] );
			ref.push_back(          pose.energies().residue_total_energies(i)[scoring::ref] );

			resnums.push_back( i ); //only need for pose or native, one or the other

		}


		Real sum_total(0);
		Real sum_fa_atr(0);
		Real sum_fa_rep(0);
		Real sum_fa_sol(0);
		Real sum_fa_intra_rep(0);
		Real sum_fa_elec(0);
		Real sum_pro_close(0);
		Real sum_hbond_sr_bb(0);
		Real sum_hbond_lr_bb(0);
		Real sum_hbond_bb_sc(0);
		Real sum_hbond_sc(0);
		Real sum_dslf_fa13(0);
		Real sum_rama(0);
		Real sum_omega(0);
		Real sum_fa_dun(0);
		Real sum_p_aa_pp(0);
		Real sum_ref(0);


		TR << "Subtracting design from native energies" << std::endl;

		for ( Size jj(1); jj <= resnums.size(); ++jj ) {

			Size ii = resnums[jj];

			std::string native_name3 = native.residue(ii).name3();
			std::string pose_name3 = pose.residue(ii).name3();

			bool is_mutation( native_name3 != pose_name3);
			bool is_mutation_option( basic::options::option[local::compare_residue_energies_mut] && !basic::options::option[local::compare_residue_energies_all] );

			//TR << "comparing residue " << ii << " " << native_name3 << " " << pose_name3 << " mutation? " << is_mutation << " mutation_option " << is_mutation_option << std::endl;

			if ( is_mutation || !is_mutation_option ) {

				Real diff_total = total[ii] - total_native[ii];
				Real diff_fa_atr = fa_atr[ii] - fa_atr_native[ii];
				Real diff_fa_rep = fa_rep[ii] - fa_rep_native[ii];
				Real diff_fa_sol = fa_sol[ii] - fa_sol_native[ii];
				Real diff_fa_intra_rep = fa_intra_rep[ii] - fa_intra_rep_native[ii];
				Real diff_fa_elec = fa_elec[ii] - fa_elec_native[ii];
				Real diff_pro_close = pro_close[ii] - pro_close_native[ii];
				Real diff_hbond_sr_bb = hbond_sr_bb[ii] - hbond_sr_bb_native[ii];
				Real diff_hbond_lr_bb = hbond_lr_bb[ii] - hbond_lr_bb_native[ii];
				Real diff_hbond_bb_sc = hbond_bb_sc[ii] - hbond_bb_sc_native[ii];
				Real diff_hbond_sc = hbond_sc[ii] - hbond_sc_native[ii];
				Real diff_dslf_fa13 = dslf_fa13[ii] - dslf_fa13_native[ii];
				Real diff_rama = rama[ii] - rama_native[ii];
				Real diff_omega = omega[ii] - omega_native[ii];
				Real diff_fa_dun = fa_dun[ii] - fa_dun_native[ii];
				Real diff_p_aa_pp = p_aa_pp[ii] - p_aa_pp_native[ii];
				Real diff_ref = ref[ii] - ref_native[ii];

				sum_total        += diff_total        ;
				sum_fa_atr       += diff_fa_atr       ;
				sum_fa_rep       += diff_fa_rep       ;
				sum_fa_sol       += diff_fa_sol       ;
				sum_fa_intra_rep += diff_fa_intra_rep ;
				sum_fa_elec    += diff_fa_elec    ;
				sum_pro_close    += diff_pro_close    ;
				sum_hbond_sr_bb  += diff_hbond_sr_bb  ;
				sum_hbond_lr_bb  += diff_hbond_lr_bb  ;
				sum_hbond_bb_sc  += diff_hbond_bb_sc  ;
				sum_hbond_sc     += diff_hbond_sc     ;
				sum_dslf_fa13    += diff_dslf_fa13    ;
				sum_rama         += diff_rama         ;
				sum_omega        += diff_omega        ;
				sum_fa_dun       += diff_fa_dun       ;
				sum_p_aa_pp      += diff_p_aa_pp      ;
				sum_ref          += diff_ref          ;

				//will print in columns
				TR << pose.pdb_info()->name() << " " << ii << " ";
				TR << std::setprecision(4) << std::fixed << diff_total        << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_atr       << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_rep       << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_sol       << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_intra_rep << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_elec    << " ";
				TR << std::setprecision(4) << std::fixed << diff_pro_close    << " ";
				TR << std::setprecision(4) << std::fixed << diff_hbond_sr_bb  << " ";
				TR << std::setprecision(4) << std::fixed << diff_hbond_lr_bb  << " ";
				TR << std::setprecision(4) << std::fixed << diff_hbond_bb_sc  << " ";
				TR << std::setprecision(4) << std::fixed << diff_hbond_sc     << " ";
				TR << std::setprecision(4) << std::fixed << diff_dslf_fa13    << " ";
				TR << std::setprecision(4) << std::fixed << diff_rama         << " ";
				TR << std::setprecision(4) << std::fixed << diff_omega        << " ";
				TR << std::setprecision(4) << std::fixed << diff_fa_dun       << " ";
				TR << std::setprecision(4) << std::fixed << diff_p_aa_pp      << " ";
				TR << std::setprecision(4) << std::fixed << diff_ref          << " ";
				TR << std::endl;

			}

		}

		TR << pose.pdb_info()->name() << " SUM-DIFFS ";
		TR << sum_total       <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_atr      <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_rep      <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_sol      <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_intra_rep<<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_elec   <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_pro_close   <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_hbond_sr_bb <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_hbond_lr_bb <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_hbond_bb_sc <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_hbond_sc    <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_dslf_fa13   <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_rama        <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_omega       <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_fa_dun      <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_p_aa_pp     <<  std::setprecision(4) << std::fixed << " ";
		TR << sum_ref         <<  std::setprecision(4) << std::fixed << " ";
		TR << std::endl;

		TR << pose.pdb_info()->name() << " score_terms: total fa_atr fa_rep fa_sol fa_intra_rep fa_elec pro_close hbond_sr_bb hbond_lr_bb hbond_bb_sc hbond_sc dslf_fa13 rama omega fa_dun p_aa_pp ref" << std::endl;

		return;
	}


	virtual
	std::string
	get_name() const { return "supercharge"; }


private:
	SizeSet surface_res_;
	std::string outputname_;
	std::string out_path_;
	utility::vector1< Real > AvNAPSA_values_;
	Size largest_mutated_AvNAPSA_;
	Pose native_;

	core::scoring::ScoreFunctionOP scorefxn_;
};

typedef utility::pointer::shared_ptr< supercharge > superchargeOP;

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

		protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new supercharge ));

		TR << "************************d**o**n**e**************************************" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
