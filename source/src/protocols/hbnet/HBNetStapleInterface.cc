// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file <protocols/hbnet/HBNetStapleInterface.cc
/// @brief inherits from HBNet; protocol for designing h-bond networks to "staple" Prot/Prot interfaces;
/// @author Scott Boyken (sboyken@gmail.com)

// Header files
//#include <protocols/hbnet/HBNet.hh>
#include <protocols/hbnet/HBNet_util.hh>
#include <protocols/hbnet/HBNetStapleInterface.hh>
#include <protocols/hbnet/HBNetStapleInterfaceCreator.hh>
#include <core/select/util/SelectResiduesByLayer.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/database/open.hh>
#include <basic/database/sql_utils.hh>
#include <core/pack/rotamer_set/rotamer_building_functions.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/pose_creation/MakePolyXMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/simple_moves/FavorSequenceProfile.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/scoring/InterfaceInfo.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/LexicographicalIterator.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using namespace core;
using namespace pack;
using namespace task;
using namespace std;
using namespace utility;
using namespace core::scoring::hbonds;
using namespace interaction_graph;
using namespace core::pack::rotamer_set;
using namespace protocols::simple_moves;

namespace protocols {
namespace hbnet {

static basic::Tracer TR( "protocols.hbnet.HBNetStapleInterface" );

HBNetStapleInterface::HBNetStapleInterface( ) :
	HBNet( "HBNetStapleInterface" ),
	all_helices_(false),
	span_all_helices_(false),
	only_symm_interfaces_(false),
	allow_onebody_networks_(false),
	his_tyr_(false),
	pH_His_(false),
	boundary_his_must_to_hbond_pos_charge_(false),
	only_start_at_interface_pairs_(false),
	use_aa_dependent_weights_(false),
	min_networks_per_pose_(1),
	max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	min_asp_glu_hbonds_(3),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0),
	pair_lists_vec_(0)
{}

HBNetStapleInterface::HBNetStapleInterface( std::string const name ) :
	HBNet( name ),
	all_helices_(false),
	span_all_helices_(false),
	only_symm_interfaces_(false),
	allow_onebody_networks_(false),
	his_tyr_(false),
	pH_His_(false),
	boundary_his_must_to_hbond_pos_charge_(false),
	only_start_at_interface_pairs_(false),
	use_aa_dependent_weights_(false),
	min_networks_per_pose_(1),
	max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	min_asp_glu_hbonds_(3),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0),
	pair_lists_vec_(0)
{}

//Constructor from code
HBNetStapleInterface::HBNetStapleInterface(
	core::scoring::ScoreFunctionCOP scorefxn,
	Size max_unsat_Hpol,
	Size min_network_size, /* 3 */
	Real hb_threshold, /* -0.75 */
	Size max_network_size, /* 15 */
	std::string des_residues, /* "STRKHYWNQDE" */
	bool find_native, /*false*/
	bool only_native, /*false*/
	bool keep_existing, /*false*/
	bool extend_existing, /*false*/
	bool only_extend /*false*/
) :
	HBNet( std::move( scorefxn ), max_unsat_Hpol, min_network_size, hb_threshold, max_network_size, des_residues, find_native, only_native, keep_existing, extend_existing, only_extend ),
	all_helices_(false),
	span_all_helices_(false),
	only_symm_interfaces_(false),
	allow_onebody_networks_(false),
	his_tyr_(false),
	pH_His_(false),
	boundary_his_must_to_hbond_pos_charge_(false),
	only_start_at_interface_pairs_(false),
	use_aa_dependent_weights_(false),
	min_networks_per_pose_(1),
	max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	min_asp_glu_hbonds_(3),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0),
	pair_lists_vec_(0)
{}

//destructor
HBNetStapleInterface::~HBNetStapleInterface()= default;

moves::MoverOP
HBNetStapleInterface::clone() const {
	return pointer::make_shared< HBNetStapleInterface >( *this );
}

moves::MoverOP
HBNetStapleInterface::fresh_instance() const {
	return pointer::make_shared< HBNetStapleInterface>( );
}

void
HBNetStapleInterface::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &data,
	filters::Filters_map const &fmap,
	moves::Movers_map const &mmap,
	core::pose::Pose const & pose
){

	// for all options that inherit from base class HBNet
	HBNet::parse_my_tag( tag, data, fmap, mmap, pose);

	max_networks_per_pose_ = tag->getOption<Size>("max_networks_per_pose",1);
	min_networks_per_pose_ = tag->getOption<Size>("min_networks_per_pose",1);
	combos_ = tag->getOption<Size>("combos",1);
	min_intermolecular_hbonds_ = tag->getOption<Size>("min_intermolecular_hbonds",1);
	min_helices_contacted_by_network_ = tag->getOption<Size>("min_helices_contacted_by_network",0);
	min_asp_glu_hbonds_ = tag->getOption<Size>("min_asp_glu_hbonds",3);
	span_all_helices_ = tag->getOption<bool>( "span_all_helices", false );
	allow_onebody_networks_ = tag->getOption<bool>("allow_onebody_networks",false);
	his_tyr_ = tag->getOption< bool >( "his_tyr", false);
	pH_His_ = tag->getOption<bool>("pH_His", false);
	only_start_at_interface_pairs_ = tag->getOption<bool>("only_start_at_interface_pairs",false);
	use_aa_dependent_weights_ = tag->getOption<bool>("use_aa_dependent_weights",false);

	if ( pH_His_ ) {
		set_upweight_starting_twobody( 2.0 ); //need better way to implement this; need to update this for non-XML use
		boundary_his_must_to_hbond_pos_charge_ = tag->getOption<bool>("boundary_his_must_to_hbond_pos_charge", false);
	}

	all_helices_ = tag->getOption<bool>("all_helical_interfaces",false);

	//17/07/28 fixing bug for intramolecular, all helical interaces case
	if ( all_helices_ ) {
		min_intermolecular_hbonds_ = 0;
	}

	only_symm_interfaces_ = tag->getOption<bool>("only_symm_interfaces",false);
	if ( tag->hasOption("interface_distance") ) {
		interf_distance_ = tag->getOption<Real>("interface_distance",8.0);
	}

	if ( tag->hasOption("jump") ) {
		jump_nums_.clear();
		std::string str = tag->getOption< std::string >( "jump", "" );
		utility::vector1<std::string> const jumps( utility::string_split( str , ',' ) );
		for ( auto const & jump : jumps ) {
			jump_nums_.push_back(utility::string2int(jump));
		}
	}
}// end parse_my_tag

void
HBNetStapleInterface::setup_packer_task_and_starting_residues( core::pose::Pose const & pose )
{
	//17/07/28 fixing bug for intramolecular, all helical interaces case
	if ( all_helices_ ) {
		min_intermolecular_hbonds_ = 0;
	}

	Size total_ind_res( ( symmetric() ) ? get_symm_info()->num_independent_residues() : pose.total_residue() );

	pair_lists_vec_.resize(total_ind_res,std::list<Size>(0));

	bool no_init_taskfactory(false);
	if ( task_factory() == nullptr ) {
		no_init_taskfactory = true;
		task_factory( pointer::make_shared< TaskFactory>( ) );
	}

	if ( jump_nums_.empty() ) {
		if ( symmetric() && ( multi_component() || only_symm_interfaces_ ) ) {
			utility::vector1<std::string> sym_dof_name_list = core::pose::symmetry::sym_dof_names( pose );
			for ( Size i = 1; i <= sym_dof_name_list.size(); i++ ) {
				Size sym_aware_jump_id(core::pose::symmetry::sym_dof_jump_num(pose,sym_dof_name_list[i]));
				jump_nums_.push_back(sym_aware_jump_id);
			}
		} else {
			for ( Size i=1; i <= pose.num_jump(); ++i ) {
				jump_nums_.push_back(i);
			}
		}
	} else {
		for ( auto const & jump_num : jump_nums_ ) {
			if ( pose.num_jump() < jump_num ) {
				if ( TR.visible() ) {
					TR.flush();
					TR << " You have specied jump number '" << jump_num << "', ";
					TR << " however the pose only has '" << pose.num_jump() << "' ";
					TR << " jumps." << endl;
				}
				runtime_assert( pose.num_jump() < jump_num );
			}
		}
	}

	if ( all_helices_ || span_all_helices_ || min_helices_contacted_by_network_ ) {
		Pose dssp_pose = pose;
		core::scoring::dssp::Dssp new_ss(dssp_pose);
		new_ss.insert_ss_into_pose(dssp_pose);
		std::string ss = dssp_pose.secstruct();
		Size start(1), end(0), h_count(1);
		for ( Size pos = 1; pos <= dssp_pose.total_residue(); ++pos ) {
			char ss_pos = dssp_pose.secstruct(pos); //HSE
			if ( (ss_pos == 'H' && (pos == dssp_pose.total_residue() || dssp_pose.secstruct(pos+1) != 'H')) ||
					(pos < dssp_pose.total_residue() && ss_pos == 'H' && !(dssp_pose.residue(pos+1).is_bonded(dssp_pose.residue(pos)))) ) {
				end = pos;
				Size pos2(pos-1);
				while ( pos2 > 1 && dssp_pose.secstruct(pos2-1) == 'H' && dssp_pose.residue(pos2).is_bonded(dssp_pose.residue(pos2-1)) ) --pos2;
				start = pos2;
				std::pair<Size,Size> boundary(start,end);
				TR.Debug << "Helix #" << h_count << ": boundary = " << boundary.first << "," << boundary.second << std::endl;
				helix_boundaries_.push_back(boundary);
				h_count++;
			}
		}
		TR.Debug << "# helices detected = " << helix_boundaries_.size() << std::endl;
		if ( span_all_helices_ ) {
			min_helices_contacted_by_network_ = helix_boundaries_.size();
		}
	}

	if ( no_init_taskfactory ) {
		if ( get_start_res_vec().empty() ) {
			for ( auto const & jump_num : jump_nums_ ) {
				scoring::InterfaceOP interf = pointer::make_shared< scoring::Interface >( jump_num );
				interf->distance( interf_distance_ );
				interf->calculate( pose ); //selects residues: sq dist of nbr atoms < interf_dist_sq (8^2)
				if ( TR.visible() ) TR << " Storing interface residues for interface " << jump_num << ":" << std::endl;
				for ( Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
					if ( !(pose.residue(resnum).is_protein()) ) {
						continue;
					}
					for ( Size resnum2 = resnum+1; resnum2 <= pose.total_residue(); ++resnum2 ) {
						if ( !(pose.residue(resnum2).is_protein()) ) {
							continue;
						}
						if ( interf->is_pair(pose.residue(resnum), pose.residue(resnum2)) ) {
							// need pair_lists_vec to keep track of interface residues in symmetric poses
							//    at pose level, numbering is linear, but in other cases, numbering repeats for each symm subunit
							Size res1_ind(resnum), res2_ind(resnum2);
							if ( symmetric() ) {
								//get_ind_res( pose, resnum, resnum2, res1_ind, res2_ind );
								res1_ind = get_ind_res( pose, resnum);
								res2_ind = get_ind_res( pose, resnum2);
							}
							if ( verbose() && TR.visible() ) {
								TR << "res " << res1_ind << " and res " << res2_ind << " are pairs across the interface." << std::endl;
							}
							std::list< Size > templist1 = pair_lists_vec_[ res1_ind ];
							if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
								pair_lists_vec_[ res1_ind ].push_back( res2_ind );
							}
							std::list< Size > templist2 = pair_lists_vec_[ res2_ind ];
							if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
								pair_lists_vec_[ res2_ind ].push_back( res1_ind );
							}
							std::pair<Size,Size> key(res1_ind,res2_ind);
							//interface_jump_map[key] = *jit;

							add_start_res( res1_ind );
						}
					}
				}
			}
			if ( all_helices_ && get_start_res_vec().empty() ) {
				for ( Size resnum = 1; resnum <= pose.total_residue(); ++resnum ) {
					if ( !(pose.residue(resnum).is_protein()) ) {
						continue;
					}
					for ( Size resnum2 = resnum+1; resnum2 <= pose.total_residue(); ++resnum2 ) {
						if ( !(pose.residue(resnum2).is_protein()) ) {
							continue;
						}
						//if (interhelical_contact(helix_boundaries_,resnum,resnum2,*ala_pose_)){
						if ( interhelical_contact(helix_boundaries_,resnum,resnum2,pose) ) {
							//TR << "TRUE: res1 = " << resnum << " and res2 = " << resnum2 << std::endl;
							// need pair_lists_vec to keep track of interface residues in symmetric poses
							//    at pose level, numbering is linear, but in other cases, numbering repeats for each symm subunit
							Size res1_ind(resnum), res2_ind(resnum2);
							if ( symmetric() ) {
								//get_ind_res( pose, resnum, resnum2, res1_ind, res2_ind );
								res1_ind = get_ind_res( pose, resnum);
								res2_ind = get_ind_res( pose, resnum2);
							}
							if ( verbose() && TR.visible() ) {
								TR << "res " << res1_ind << " and res " << res2_ind << " are pairs across an inter-helical interface." << std::endl;
							}
							std::list< Size > templist1 = pair_lists_vec_[ res1_ind ];
							if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
								pair_lists_vec_[ res1_ind ].push_back( res2_ind );
							}
							std::list< Size > templist2 = pair_lists_vec_[ res2_ind ];
							if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
								pair_lists_vec_[ res2_ind ].push_back( res1_ind );
							}
							std::pair<Size,Size> key(res1_ind,res2_ind);
							//interface_jump_map[key] = jump_nums_.size() + 1;

							add_start_res( res1_ind );
						}
					}
				}
			}
		}
		//task_factory_ = new TaskFactory;
		toolbox::task_operations::DesignAroundOperationOP desaround = pointer::make_shared< toolbox::task_operations::DesignAroundOperation >();
		std::set< core::Size > temp_start_vec( get_start_res_vec() );
		for ( core::Size it : temp_start_vec ) {
			desaround->include_residue(it);
		}
		desaround->repack_shell(12.0);
		desaround->design_shell(8.0);
		task_factory()->push_back(desaround);

		set_task( create_ptask( pose ) );
	} else {
		set_task( create_ptask(pose) ); //temp task
		if ( get_start_res_vec().empty() ) {
			for ( auto const & jump_num : jump_nums_ ) {
				scoring::InterfaceOP interf = pointer::make_shared< scoring::Interface >( jump_num );
				interf->distance(interf_distance_);
				interf->calculate( pose ); //within 8 angstroms of Ala nbr residue
				//interf->calculate( *ala_pose_ ); //within 8 angstroms of Ala nbr residue

				//for ( std::set<Size>::iterator resnum = repack_and_des_residues_.begin(); resnum != repack_and_des_residues_.end(); ++resnum ){
				for ( Size resnum = 1; resnum < pose.total_residue(); ++resnum ) {
					for ( Size resnum2 = resnum+1; resnum2 <= pose.total_residue(); ++resnum2 ) {
						if ( !(pose.residue(resnum2).is_protein()) ) {
							continue;
						}
						if ( interf->is_pair(pose.residue(resnum),pose.residue(resnum2)) ) {
							// need pair_lists_vec to keep track of interface residues in symmetric poses
							//    at pose level, numbering is linear, but in other cases, numbering repeats for each symm subunit
							Size res1_ind(resnum), res2_ind(resnum2);
							if ( symmetric() ) {
								//get_ind_res( pose, *resnum, resnum2, res1_ind, res2_ind );
								res1_ind = get_ind_res( pose, resnum);
								res2_ind = get_ind_res( pose, resnum2);
							}
							if ( std::find(pair_lists_vec_[ res1_ind ].begin(),pair_lists_vec_[ res1_ind ].end(),res2_ind) == pair_lists_vec_[ res1_ind ].end() ) {
								pair_lists_vec_[ res1_ind ].push_back( res2_ind );
							}
							if ( std::find(pair_lists_vec_[ res2_ind ].begin(),pair_lists_vec_[ res2_ind ].end(),res1_ind) == pair_lists_vec_[ res2_ind ].end() ) {
								pair_lists_vec_[ res2_ind ].push_back( res1_ind );
							}
							//std::pair<Size,Size> key(res1_ind,res2_ind);
							//interface_jump_map[key] = *jit;
							if ( verbose() && TR.visible() ) {
								TR << "res " << res1_ind << " and res " << res2_ind << " are pairs across the interface." << std::endl;
							}
							add_start_res(res1_ind);
						}
					}
				}
			}
		}
		if ( all_helices_ && get_start_res_vec().empty() ) {
			for ( Size resnum = 1; resnum < pose.total_residue(); ++resnum ) {
				if ( !(pose.residue(resnum).is_protein()) ) {
					continue;
				}
				for ( Size resnum2 = resnum+1; resnum2 <= pose.total_residue(); ++resnum2 ) {
					if ( !(pose.residue(resnum2).is_protein()) ) {
						continue;
					}
					//if (interhelical_contact(helix_boundaries_,*resnum,resnum2,*ala_pose_))
					if ( interhelical_contact(helix_boundaries_,resnum,resnum2,pose) ) {
						// need pair_lists_vec to keep track of interface residues in symmetric poses
						//    at pose level, numbering is linear, but in other cases, numbering repeats for each symm subunit
						Size res1_ind(resnum), res2_ind(resnum2);
						if ( symmetric() ) {
							//get_ind_res( pose, *resnum, resnum2, res1_ind, res2_ind );
							res1_ind = get_ind_res( pose, resnum);
							res2_ind = get_ind_res( pose, resnum2);
						}
						if ( verbose() && TR.visible() ) {
							TR << "res " << res1_ind << " and res " << res2_ind << " are pairs across an inter-helical interface." << std::endl;
						}
						std::list< Size > templist1 = pair_lists_vec_[ res1_ind ];
						if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
							pair_lists_vec_[ res1_ind ].push_back( res2_ind );
						}
						std::list< Size > templist2 = pair_lists_vec_[ res2_ind ];
						if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
							pair_lists_vec_[ res2_ind ].push_back( res1_ind );
						}
						add_start_res(res1_ind);
					}
				}
			}
		}
	}
}//setup

//careful! only works if pose has not added or deleted residues; careful in symmetric bridging_waters case!
bool
HBNetStapleInterface::residues_are_interface_pairs( Size const res1, Size const res2 )
{
	Size const res1_ind( get_ind_res( get_orig_pose(), res1 ) );
	Size const res2_ind( get_ind_res( get_orig_pose(), res2 ) );
	if ( std::find(pair_lists_vec_[ res1_ind ].begin(),pair_lists_vec_[ res1_ind ].end(),res2_ind) != pair_lists_vec_[ res1_ind ].end()
			|| std::find(pair_lists_vec_[ res2_ind ].begin(),pair_lists_vec_[ res2_ind ].end(),res1_ind) != pair_lists_vec_[ res2_ind ].end() ) {
		return true;
	}
	return false;
}

bool
HBNetStapleInterface::has_pH_His( core::pose::Pose const & pose, hbond_net_struct & i )
{
	if ( !(network_contains_aa( 'H', i ) ) ) {
		return false;
	}
	runtime_assert( !(i.hbond_vec.empty() ) );

	bool found(false);
	std::vector<Size> his_acceptors(0);
	std::vector<Size> his_donors(0);
	for ( utility::vector1<HBondCOP>::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
		Size arsd((*h)->acc_res());
		Size drsd((*h)->don_res());
		char a_aa = pose.residue( arsd ).name1();
		char d_aa = pose.residue( drsd ).name1();
		if ( a_aa == 'H' && !((*h)->acc_atm_is_backbone()) ) {
			// ensure that his acceptor and its donor are different chains
			//    //for cages, need to check different building blocks rather than chains
			if ( multi_component() && ( get_symm_info()->get_component( arsd ) != get_symm_info()->get_component( drsd ) ) ) {
				his_acceptors.push_back( arsd );
			} else if ( pose.chain( arsd ) != pose.chain( drsd ) ) {
				his_acceptors.push_back( arsd );
			}
		} else if ( d_aa == 'H' && !((*h)->don_hatm_is_backbone()) ) {
			//   // ensure his donor and its acceptor are same chain (preorganized)
			//   if ( multi_component_ && ( symm_info_->get_component( arsd ) == symm_info_->get_component( drsd ) ) ) {
			//    his_donors.push_back( drsd );
			//   } else if ( pose.chain( drsd ) == pose.chain( arsd ) ) {
			//    his_donors.push_back( drsd );
			//   }
			his_donors.push_back( drsd );
		}
	}
	std::sort( his_acceptors.begin(), his_acceptors.end() );
	std::sort( his_donors.begin(), his_donors.end() );
	std::vector<Size> intersec(0);
	std::set_intersection( his_acceptors.begin(), his_acceptors.end(), his_donors.begin(), his_donors.end(), std::back_inserter(intersec) );
	if ( intersec.size() > 0 ) {
		found = true;
		if ( TR.visible() ) TR << "FOUND PH SENSITIVE HIS AT: ";
		for ( std::vector<Size>::const_iterator it = intersec.begin(); it != intersec.end(); ++it ) {
			if ( TR.visible() ) TR << *it << ", ";
		}
		if ( TR.visible() ) TR << std::endl;
	}
	//NEED TO ADD ADDITIONAL LOGIC HERE
	//    if ( boundary_his_must_to_hbond_pos_charge_ ){
	//
	//    }
	return found;
}

bool
HBNetStapleInterface::network_meets_initial_criteria( hbond_net_struct const & network )
{
	if ( !symmetric() && min_helices_contacted_by_network_ && ( num_helices_w_hbond( network ) < min_helices_contacted_by_network_ ) ) {
		//if (verbose_)
		//    TR << "min_helices_contacted_by_network_ = " << num_helices_w_hbond( i ) << std::endl;
		return false;
	}
	return true;
}

bool
HBNetStapleInterface::network_meets_final_criteria( Pose const & pose, hbond_net_struct & i )
{
	// update core residues
	//core::select::residue_selector::ResidueSubset temp_core_residues = get_core_residues();
	//core::select::residue_selector::ResidueSubset temp_boundary_residues = get_boundary_residues();
	//if ( bridging_waters_ ){
	//    update_core_and_boundary_residues( pose );
	//}
	if ( pH_His_ && !(has_pH_His( pose, i )) ) {
		return false;
	}
	if ( his_tyr_ && ( ( !(network_contains_aa('Y', i)) || !(network_contains_aa('H', i) )) || !(his_tyr_connectivity( pose, i )) ) ) { //need to generalize for user-specified connectivites
		return false;
	}
	if ( min_asp_glu_hbonds_ > 0 ) {
		// check the actual h-bonds rather than iterate through residues because res # can change in water cases (residues added/removed)
		for ( utility::vector1< HBondCOP >::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
			if ( res_is_core( (*h)->acc_res() ) && ( pose.residue((*h)->acc_res()).name1() == 'D' || pose.residue((*h)->acc_res()).name1() == 'E' ) ) {
				// network hbond_set does not include bb_bb h-bonds
				if ( i.hbond_set->nhbonds( (*h)->acc_res(), false /* include_only_allowed */ ) < min_asp_glu_hbonds_ ) {
					return false;
				}
			}
		}
	}

	if ( min_intermolecular_hbonds_ ) {
		//i.num_intermolecular_hbs = (symmetric_) ? symm_num_intermolecular_hbonds( i, pose ) : num_intermolecular_hbonds( i, pose );

		//i.num_intermolecular_hbs = num_intermolecular_hbonds( i, pose );
		if ( verbose() ) {
			TR << "i.num_intermolecular_hbs = " << num_intermolecular_hbonds( i, pose ) << std::endl;
		}
		if ( num_intermolecular_hbonds( i, pose ) < min_intermolecular_hbonds_ ) {
			return false;
		}
	}
	//for symmetric PPI interfaces, check if network is one connected network that spans entire symmetric interface (or a disjoint network repeated symmetrically)
	if ( symmetric() && !multi_component() && !( pose.total_residue() > 500 ) ) {
		core::conformation::symmetry::SymmetricConformation const & SymmConf(dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation()));
		bool network_spans_entire_symm_interface(false);
		for ( utility::vector1< HBondCOP >::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
			//i.show(in_pose, 1, TR);
			Size drsd((*h)->don_res());
			//Size don_hatm((*h)->don_hatm());
			Size arsd((*h)->acc_res());
			//Size aatm((*h)->acc_atm());
			//if ( arsd == get_ind_res( *orig_pose_, drsd) || drsd == get_ind_res( *orig_pose_, arsd) ) {
			utility::vector1<Size> acc_clones( SymmConf.Symmetry_Info()->bb_clones( arsd ) );
			utility::vector1<Size> don_clones( SymmConf.Symmetry_Info()->bb_clones( drsd ) );
			if ( std::find( acc_clones.begin(), acc_clones.end(), drsd) != acc_clones.end() || std::find( don_clones.begin(), don_clones.end(), arsd) != don_clones.end() ) {
				network_spans_entire_symm_interface = true;
				break;
			}
		}
		//utility::vector1< HBondResStructCOP > temp_residues( i.residues );
		if ( network_spans_entire_symm_interface ) {
			for ( HBondResStructCOP const & ir : i.residues ) {
				//for ( utility::vector1< HBondResStructCOP >::const_iterator ir = i.residues.begin(); ir != i.residues.end(); ++ir ) {
				if ( find_hbond_res_struct( i.asymm_residues, ir->resnum ) == i.asymm_residues.end() ) {
					i.asymm_residues.push_back( ir );
				}
				utility::vector1< Size > resi_clones( SymmConf.Symmetry_Info()->bb_clones( ir->resnum ) );
				for ( Size r : resi_clones ) {
					if ( find_hbond_res_struct( i.residues, r ) == i.residues.end() ) {
						i.asymm_residues.push_back(
							pointer::make_shared< hbond_res_struct >(
							r,
							0,
							ir->aa,
							pose.pdb_info()->chain(r),
							pose.residue(r).is_protein(),
							pose.residue(r).name1() == 'w',
							pose.residue(r).is_ligand()
							)
						);
					}
				}
			}
		}
	}
	//ALL CRITERIA AFTER THIS POINT DEPEND ON i.asymmetric residues being populated if symmetric_

	//TR << "num_helices_contacted = " << num_helices_w_hbond( i ) << std::endl;
	if ( ( min_helices_contacted_by_network_ ) && ( num_helices_w_hbond( i ) < min_helices_contacted_by_network_ ) ) {
		if ( verbose() && TR.visible() ) {
			TR << "min_helices_contacted_by_network_ = " << num_helices_w_hbond( i ) << std::endl;
		}
		return false;
	}

	//if ( bridging_waters_ ){
	//    set_core_residues( temp_core_residues );
	//    set_boundary_residues( temp_boundary_residues );
	//}

	return true;
} // network_meets_final_criteria

bool
HBNetStapleInterface::state_is_starting_aa_type( Size const res, Size const rot_id )
{

	//char const aa( rotamer_sets_->rotamer_set_for_residue( (platform::uint)(res) )->rotamer(rot_id)->name1() );
	char const aa( get_aa_for_state( res, rot_id ) );
	if ( pH_His_ && aa != 'H' ) {
		return false;
	} else if ( his_tyr_ && !( aa == 'H' || aa == 'Y' ) ) {
		return false;
	}
	return true;
}

bool
HBNetStapleInterface::pair_meets_starting_criteria( core::Size const res1, core::Size const rot1, core::Size const res2, core::Size const rot2 )
{
	if ( only_start_at_interface_pairs_ && !(residues_are_interface_pairs( res1, res2 )) ) return false;
	if ( his_tyr_ && !( get_aa_for_state( res1, rot1 ) == 'Y' && get_aa_for_state( res2, rot2 ) == 'H' )
			&& !( get_aa_for_state( res1, rot1 ) == 'H' && get_aa_for_state( res2, rot2 ) == 'H' ) ) {
		return false;
	}
	//    else if ( pH_His_ && boundary_his_must_to_hbond_pos_charge_ ){
	//        char c1( get_aa_for_state( res1, rot1 ) );
	//        char c2( get_aa_for_state( res2, rot2 ) );
	//        if ( !( c1 == 'H' && !(res_is_core( res1 )) && ( c2 == 'K' || c2 == 'R' || c2 == 'H' ) )
	//               && !( c2 == 'H' && !(res_is_core( res2 )) && ( c1 == 'K' || c1 == 'R' || c1 == 'H' ) ) ){
	//            return false;
	//        }
	//    }
	if ( pH_His_ ) {
		char c1( get_aa_for_state( res1, rot1 ) );
		char c2( get_aa_for_state( res2, rot2 ) );
		if ( c1 == 'H' || c2 == 'H' ) {
			return true;
		} else {
			return false;
		}
	}
	return true;
}

Real
HBNetStapleInterface::scale_twobody_energy( core::Real input_twobody_energy, char res1, char res2 )
{
	if ( use_aa_dependent_weights_ ) {
		Real scale_factor_1(1.0);
		Real scale_factor_2(1.0);

		if ( res1 == 'S' || res1 == 'T' || res1 == 'Y' ) scale_factor_1 = 1.3;
		else if ( res1 == 'N' || res1 == 'H' || res1 == 'W' ) scale_factor_1 = 1.15;
		else if ( res1 == 'Q' ) scale_factor_1 = 0.85;
		else if ( res1 == 'D' ) scale_factor_1 = 0.8;
		else if ( res1 == 'E' || res1 == 'K' || res1 == 'R' ) scale_factor_1 = 0.7;

		if ( res2 == 'S' || res2 == 'T' || res2 == 'Y' ) scale_factor_2 = 1.3;
		else if ( res2 == 'N' || res2 == 'H' || res2 == 'W' ) scale_factor_2 = 1.15;
		else if ( res2 == 'Q' ) scale_factor_2 = 0.85;
		else if ( res2 == 'D' ) scale_factor_2 = 0.8;
		else if ( res2 == 'E' || res2 == 'K' || res2 == 'R' ) scale_factor_2 = 0.7;

		return scale_factor_1*scale_factor_2*input_twobody_energy;
	}
	//TR << "input_twobody = " << input_twobody_energy << std::endl;
	//TR << "scale_factor_1*scale_factor_2*input_twobody_energy = " << scale_factor_1*scale_factor_2*input_twobody_energy << std::endl;

	return input_twobody_energy;
}

void
HBNetStapleInterface::prepare_output()
{
	output_networks( min_networks_per_pose_ == 1 ); //will add single networks to output vector

	if ( only_native() && ( get_native_vec().size() < min_networks_per_pose_ ) ) {
		get_native_vec().clear();
		return;
	}
	if ( get_extend_existing_networks() ) std::sort( get_net_vec().begin(), get_net_vec().end(), compare_net_vec() ); //sort all networks to put extended first
	if ( max_networks_per_pose_ > 1 ) { //add multiple network sets to output vector
		for ( std::vector< HBondNetStructOP >::const_iterator netit = get_net_vec().begin(); netit != get_net_vec().end(); ++netit ) {
			//std::string network( (pdb_numbering() ) ? ( print_list_to_string( get_orig_pose(), **netit) ) : (print_list_to_string( **netit) ) );
			//if ( TR.visible() ) TR << "combining networks " << (*netit)->id << ": " << network;
			//if ( TR.visible() ) TR << std::endl << "combining networks " << (*netit)->id << " ";

			HBondNetStructOP new_network = *netit;
			std::set< Size > net_ids;
			net_ids.insert( new_network->id );
			Size staple_count(1);

			if ( get_keep_existing_networks() ) {
				if ( (*netit)->is_extended ) {
					for ( auto j = netit; ++j != get_net_vec().end(); ) {
						if ( (*j)->is_extended ) {
							bool compatible( true );
							for ( auto const & net_id : net_ids ) {
								if ( net_clash( *(get_network_by_id(net_id)), **j ) ) {
									compatible = false;
									break;
								}
							}
							if ( compatible ) {
								net_ids.insert( (*j)->id );
								staple_count++;
							}
						}
					}
				} else {
					for ( auto const & k : get_net_vec() ) {
						if ( k->is_extended ) {
							bool compatible( true );
							for ( auto const & net_id : net_ids ) {
								if ( net_clash( *(get_network_by_id(net_id)), *k ) ) {
									compatible = false;
									break;
								}
							}
							if ( compatible ) {
								net_ids.insert( k->id );
								staple_count++;
							}
						}
					}

				}
			}
			rec_add_staple( netit, net_ids, 1 );
		}
	}
	TR.flush();
}

bool
HBNetStapleInterface::network_spans_all_helices( hbond_net_struct const & i )
{
	if ( num_helices_w_hbond( i ) == helix_boundaries_.size() ) {
		return true;
	}
	return false;
}

Size
HBNetStapleInterface::num_helices_w_hbond( hbond_net_struct const & i )
{
	return ( i.asymm_residues.empty() ) ? num_helices_w_hbond( i.residues ) : num_helices_w_hbond( i.asymm_residues );
}

Size
HBNetStapleInterface::num_helices_w_hbond( utility::vector1< HBondResStructCOP > const & residues )
{
	Size num_helices(0);
	utility::vector1< bool > helix_has_hbond_residue( helix_boundaries_.size(), false );
	for ( auto const & residue : residues ) {
		Size helix_id( get_helix_id( residue->resnum ) ); //returns 0 if residues is not part of a helix
		if ( helix_id > 0 ) {
			helix_has_hbond_residue[ helix_id ] = true;
		}
	}
	for ( auto && hel : helix_has_hbond_residue ) {
		if ( hel ) {
			num_helices++;
		}
	}
	return num_helices;
}

void
HBNetStapleInterface::rec_add_staple( std::vector< HBondNetStructOP >::const_iterator netit, std::set< Size > net_ids, Size staple_count )
{
	Size combo_count(1);
	auto next_netit(netit);
	while ( combo_count <= combos_ && ++next_netit != get_net_vec().end() )
			{ //number of combinations of multiple networks to try (default = 1)
		bool compatible( true );
		for ( auto const & net_id : net_ids ) {
			runtime_assert( get_network_by_id(net_id) != nullptr );
			bool branch(false);
			//            if ( !(is_sub_residues( (get_network_by_id(*net_id))->residues, (get_network_by_id(*net_id))->residues, branch ))
			//                && !branch && !(net_clash( *(get_network_by_id(*net_id)), **next_netit )) ) {
			if ( is_sub_residues( (get_network_by_id(net_id))->residues, (get_network_by_id((*next_netit)->id)->residues), branch )
					|| branch || net_clash( *(get_network_by_id(net_id)), **next_netit ) ) {
				compatible = false;
				break;
			}
		}
		if ( compatible ) {
			//Size net_index1 = next_netit - get_net_vec().begin();
			//std::string network( (pdb_numbering() ) ? ( print_list_to_string( get_orig_pose(), **next_netit ) ) : (print_list_to_string( **next_netit ) ) );
			//if ( TR.visible() ) TR << "; and " << (*next_netit)->id << ": " << network;
			//if ( TR.visible() ) TR << " and " << (*next_netit)->id;
			std::set< Size > new_net_ids( net_ids );
			new_net_ids.insert( (*next_netit)->id );
			staple_count++;

			if ( staple_count >= min_networks_per_pose_ ) {
				bool already_stored(false);
				for ( auto const & output_vec : get_output_vector() ) {
					if ( output_vec == new_net_ids ) {
						already_stored = true;
					}
				}
				if ( !already_stored ) {
					get_output_vector().push_back( new_net_ids );
					combo_count++;
				}
			}
			// < not <= here because will reach max with next rec_add_staple() call
			if ( staple_count < max_networks_per_pose_ ) { // < not <= here because will reach max with next rec_add_staple() call
				rec_add_staple(next_netit, new_net_ids, staple_count);
				staple_count--;
			}
		}
	}
}

bool
HBNetStapleInterface::same_helix( utility::vector1< std::pair<Size,Size> > const helix_boundaries, Size const r1, Size const r2)
{
	for ( auto const & helix_boundarie : helix_boundaries ) {
		if ( (r1 >= helix_boundarie.first && r1 <= helix_boundarie.second) && (r2 >= helix_boundarie.first && r2 <= helix_boundarie.second) ) {
			return true;
		}
	}
	return false;
}

///@details returns false if either of the residues are not Helix (H); returns true if residues contact eachother from different helices (interhelical contact)
bool
HBNetStapleInterface::interhelical_contact( utility::vector1< std::pair<Size,Size> > const helix_boundaries, Size const r1, Size const r2, Pose const & pose)
{
	if ( pose.residue(r1).nbr_atom_xyz().distance_squared(pose.residue(r2).nbr_atom_xyz()) < (interf_distance_*interf_distance_) ) {
		if ( !same_helix(helix_boundaries,r1,r2) ) {
			return true;
		}
	}
	return false;
}

///@details returns false if residues is not part of a helix
Size
HBNetStapleInterface::get_helix_id( Size r1 )
{
	runtime_assert( !(helix_boundaries_.empty()) );
	//for (utility::vector1< std::pair<Size,Size> >::iterator h = helix_boundaries.begin(); h != helix_boundaries.end(); ++h){
	for ( Size h = 1; h <= helix_boundaries_.size(); h++ ) {
		if ( r1 >= helix_boundaries_[h].first && r1 <= helix_boundaries_[h].second ) {
			return h;
		}
	}
	return 0;
}

std::string
HBNetStapleInterface::print_additional_headers()
{
	std::stringstream output;
	if ( min_intermolecular_hbonds_ > 0 ) {
		output << "interf_hbs \t";
	}
	if ( min_helices_contacted_by_network_ > 0 ) {
		output << "helics_contacted \t";
	}
	return output.str();
}

std::string
HBNetStapleInterface::print_additional_info_for_net( hbond_net_struct & i, Pose const & pose )
{
	std::stringstream output;
	if ( min_intermolecular_hbonds_ > 0 ) {
		output << num_intermolecular_hbonds( i, pose ) << "\t";
	}
	if ( min_helices_contacted_by_network_ > 0 ) {
		output << num_helices_w_hbond( i ) << "\t";
	}
	return output.str();
}

Size
HBNetStapleInterface::num_intermolecular_hbonds( hbond_net_struct & i, core::pose::Pose const & pose )
{
	// NEED TO SCORE FIRST IF GOING TO CALL get_hbond_atom_pairs()

	runtime_assert( !(i.hbond_vec.empty() ) );

	//return get_intermolecular_hbonds( i );

	Size num_intermol_hbs(0);
	for ( utility::vector1<HBondCOP>::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
		Size arsd((*h)->acc_res());
		Size drsd((*h)->don_res());
		//        if ( pose.chain(arsd) != pose.chain(drsd) && !( pose.residue(arsd).is_water() && pose.residue(drsd).is_water() ) ) {
		if ( pose.chain(arsd) != pose.chain(drsd) && !( pose.residue(arsd).name1() == 'w' && pose.residue(drsd).name1() == 'w' ) ) {
			//if ( pose.chain(arsd) != pose.chain(drsd) ) {
			num_intermol_hbs++;
		}
	}
	return num_intermol_hbs;
}

std::string HBNetStapleInterface::get_name() const {
	return mover_name();
}

std::string HBNetStapleInterface::mover_name() {
	return "HBNetStapleInterface";
}

void HBNetStapleInterface::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	HBNet::attributes_for_hbnet( attlist );
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "combos", xsct_non_negative_integer, "when combining multiple networks (max_networks_per_pose greater than 1) this is the max number of combinations tried for each network", "1")
		+ XMLSchemaAttribute::attribute_w_default( "max_networks_per_pose", xsct_non_negative_integer, "Maximum number of networks to make in a pose", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "min_networks_per_pose", xsct_non_negative_integer, "Minimum number of hydrogen bonds to generate per pose", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "min_intermolecular_hbonds", xsct_non_negative_integer, "Minimum number of intermolecular hydrogen bonds to form", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "min_helices_contacted_by_network", xsct_non_negative_integer, "Minimum number of helices the formed network should contact", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_asp_glu_hbonds", xsct_non_negative_integer, "Minimum number of h-bonds Asp or Glu must participate in", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "span_all_helices", xsct_rosetta_bool, "Must the network span all helices in the pose?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "allow_onebody_networks", xsct_rosetta_bool, "Allow networks within a pose", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "his_tyr", xsct_rosetta_bool, "Include histidine and tyrosine in the network", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "pH_His", xsct_rosetta_bool, "Identify and handle pH-sensitive histidines", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "only_start_at_interface_pairs", xsct_rosetta_bool, "Only start IG traversal with h-bonds that span across interface (different chains)", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_aa_dependent_weights", xsct_rosetta_bool, "weight twobody IG energies depending on aa type", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "all_helical_interfaces", xsct_rosetta_bool, "Interfaces must be composed of helices", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "only_symm_interfaces", xsct_rosetta_bool, "Only staple symmetric interfaces", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "interface_distance", xsct_real, "Desired distance between partners at interface", "8.0" )
		+ XMLSchemaAttribute( "jump", xs_string, "A comma-separated list of jumps across interface" );

	moves::xsd_type_definition_w_attributes( xsd, mover_name(), "\"Staples\" interface specificity with a hydrogen bond network.", attlist );
}

std::string HBNetStapleInterfaceCreator::keyname() const {
	return HBNetStapleInterface::mover_name();
}

moves::MoverOP
HBNetStapleInterfaceCreator::create_mover() const {
	return pointer::make_shared< HBNetStapleInterface >();
}

void HBNetStapleInterfaceCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HBNetStapleInterface::provide_xml_schema( xsd );
}


} //hbnet
} //protocols
