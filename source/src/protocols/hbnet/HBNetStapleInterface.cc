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
//#include <core/pack/rotamer_set/rotamer_building_functions.cc>

#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/toolbox/task_operations/PreventResiduesFromRepackingOperation.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <protocols/simple_moves/FavorSequenceProfile.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/scoring/InterfaceInfo.hh>

#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/symmetry/util.hh>
//#include <core/io/Remarks.hh>
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
//#include <core/scoring/BridgingWaterPotential.hh>
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
//#include <core/graph/Graph.hh>

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

//#include <devel/matdes/RestrictIdentitiesOperation.hh>


using namespace core;
using namespace pack;
using namespace task;
using namespace std;
//using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace interaction_graph;
using namespace core::pack::rotamer_set;
using namespace protocols::simple_moves;
//using namespace boost::numeric::ublas;

namespace protocols {
namespace hbnet {

static THREAD_LOCAL basic::Tracer TR( "protocols.hbnet.HBNetStapleInterface" );

std::string
HBNetStapleInterfaceCreator::keyname() const {
	return HBNetStapleInterfaceCreator::mover_name();
}

protocols::moves::MoverOP
HBNetStapleInterfaceCreator::create_mover() const {
	return protocols::moves::MoverOP( new HBNetStapleInterface );
}

std::string
HBNetStapleInterfaceCreator::mover_name() {
	return "HBNetStapleInterface";
}

HBNetStapleInterface::HBNetStapleInterface( ) :
	HBNet( "HBNetStapleInterface" ),
	all_helices_(0),
	span_all_helices_(0),
	only_symm_interfaces_(0),
	allow_onebody_networks_(0),
	his_tyr_(0),
	pH_His_(0),
	boundary_his_must_to_hbond_pos_charge_(0),
	runcount_(0),
	//min_staples_per_interface_(1),
    min_networks_per_pose_(1),
	//max_staples_per_interface_(1),
    max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0)
{}

HBNetStapleInterface::HBNetStapleInterface( std::string const name ) :
	HBNet( name ),
	all_helices_(0),
	span_all_helices_(0),
	only_symm_interfaces_(0),
	allow_onebody_networks_(0),
	his_tyr_(0),
	pH_His_(0),
	boundary_his_must_to_hbond_pos_charge_(0),
	runcount_(0),
    //min_staples_per_interface_(1),
    min_networks_per_pose_(1),
    //max_staples_per_interface_(1),
    max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0)
{}

//Constructor from code
HBNetStapleInterface::HBNetStapleInterface( core::scoring::ScoreFunctionCOP scorefxn,
	Size max_unsat,
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
	HBNet( scorefxn, max_unsat, min_network_size, hb_threshold, max_network_size, des_residues, find_native, only_native, keep_existing, extend_existing, only_extend ),
	all_helices_(0),
	span_all_helices_(0),
	only_symm_interfaces_(0),
	allow_onebody_networks_(0),
	his_tyr_(0),
	pH_His_(0),
	boundary_his_must_to_hbond_pos_charge_(0),
	runcount_(0),
	min_networks_per_pose_(1),
	max_networks_per_pose_(1),
	combos_(1),
	min_intermolecular_hbonds_(1),
	min_helices_contacted_by_network_(0),
	interf_distance_(8.0),
	jump_nums_(0),
	helix_boundaries_(0)
{}

//destructor
HBNetStapleInterface::~HBNetStapleInterface(){}

protocols::moves::MoverOP
HBNetStapleInterface::clone() const {
	return( protocols::moves::MoverOP( new HBNetStapleInterface( *this ) ) );
}

protocols::moves::MoverOP
HBNetStapleInterface::fresh_instance() const {
	return protocols::moves::MoverOP( new HBNetStapleInterface );
}

void
HBNetStapleInterface::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const &fmap, protocols::moves::Movers_map const &mmap, core::pose::Pose const & pose ){
	//TR << "verbose_ = " << verbose_ << std::endl;
	HBNet::parse_my_tag( tag, data, fmap, mmap, pose);

	max_networks_per_pose_ = tag->getOption<Size>("max_staples_per_interface",1); //legacy, should be deprecated
    max_networks_per_pose_ = tag->getOption<Size>("max_networks_per_pose",1);
	min_networks_per_pose_ = tag->getOption<Size>("min_staples_per_interface",1); //legacy, should be deprecated
    min_networks_per_pose_ = tag->getOption<Size>("min_networks_per_pose",1);
	combos_ = tag->getOption<Size>("combos",1);
	min_intermolecular_hbonds_ = tag->getOption<Size>("min_intermolecular_hbonds",1);
	min_helices_contacted_by_network_ = tag->getOption<Size>("min_helices_contacted_by_network",0);
	//span_interface_ = tag->getOption<bool>( "hbond_must_span_interface", 1 );
	span_all_helices_ = tag->getOption<bool>( "span_all_helices", 0 );
	//only_onebody_ = tag->getOption<bool>("only_onebody",0);
	allow_onebody_networks_ = tag->getOption<bool>("allow_onebody_networks",0);
	his_tyr_ = tag->getOption< bool >( "his_tyr", 0);
	pH_His_ = tag->getOption<bool>("pH_His", false);

	if ( pH_His_ ) {
		set_upweight_starting_twobody( 2.0 ); //need better way to implement this; need to update this for non-XML use
		boundary_his_must_to_hbond_pos_charge_ = tag->getOption<bool>("boundary_his_must_to_hbond_pos_charge", false);
	}

	all_helices_ = tag->getOption<bool>("all_helical_interfaces",0);
	only_symm_interfaces_ = tag->getOption<bool>("only_symm_interfaces",0);
	if ( tag->hasOption("interface_distance") ) {
		interf_distance_ = tag->getOption<Real>("interface_distance",8.0);
	}

	if ( tag->hasOption("jump") ) {
		jump_nums_.clear();
		std::string str = tag->getOption< std::string >( "jump", "" );
		utility::vector1<std::string> const jumps( utility::string_split( str , ',' ) );
		for ( utility::vector1<std::string>::const_iterator it = jumps.begin(); it != jumps.end(); ++it ) {
			jump_nums_.push_back(utility::string2int(*it));
		}
	}
}// end parse_my_tag

void
HBNetStapleInterface::setup( core::pose::Pose & pose )
{
	Size total_ind_res(pose.total_residue());
	if ( symmetric() ) {
		total_ind_res = get_symm_info()->num_independent_residues();
	}
	utility::vector1< std::list< Size > > pair_lists_vec(total_ind_res);

	bool no_init_taskfactory(false);
	if ( task_factory() == 0 ) {
		no_init_taskfactory = true;
		task_factory( TaskFactoryOP( new TaskFactory ) );
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
		for ( utility::vector1<Size>::iterator jit = jump_nums_.begin(); jit != jump_nums_.end(); ++jit ) {
			if ( pose.num_jump() < *jit ) {
				if ( TR.visible() ) {
					TR.flush();
					TR << " You have specied jump number '" << *jit << "', ";
					TR << " however the pose only has '" << pose.num_jump() << "' ";
					TR << " jumps." << endl;
				}
				runtime_assert( pose.num_jump() < *jit );
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
			for ( utility::vector1<Size>::iterator jit = jump_nums_.begin(); jit != jump_nums_.end(); ++jit ) {
				protocols::scoring::InterfaceOP interf = protocols::scoring::InterfaceOP( new protocols::scoring::Interface( *jit ) );
				interf->distance(interf_distance_);
				interf->calculate( pose ); //selects residues: sq dist of nbr atoms < interf_dist_sq (8^2)
				//interf->calculate( *ala_pose_ ); //selects residues: sq dist of nbr atoms < interf_dist_sq (8^2)
				if ( TR.visible() ) TR << " Storing interface residues for interface " << *jit << ":" << std::endl;
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
							std::list< Size > templist1 = pair_lists_vec[ res1_ind ];
							if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
								pair_lists_vec[ res1_ind ].push_back( res2_ind );
							}
							std::list< Size > templist2 = pair_lists_vec[ res2_ind ];
							if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
								pair_lists_vec[ res2_ind ].push_back( res1_ind );
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
							std::list< Size > templist1 = pair_lists_vec[ res1_ind ];
							if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
								pair_lists_vec[ res1_ind ].push_back( res2_ind );
							}
							std::list< Size > templist2 = pair_lists_vec[ res2_ind ];
							if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
								pair_lists_vec[ res2_ind ].push_back( res1_ind );
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
		protocols::toolbox::task_operations::DesignAroundOperationOP desaround = protocols::toolbox::task_operations::DesignAroundOperationOP( new protocols::toolbox::task_operations::DesignAroundOperation() );
		std::set< core::Size > temp_start_vec( get_start_res_vec() );
		for ( std::set< core::Size >::const_iterator it = temp_start_vec.begin(); it != temp_start_vec.end(); ++it ) {
			desaround->include_residue(*it);
		}
		desaround->repack_shell(12.0);
		desaround->design_shell(8.0);
		task_factory()->push_back(desaround);

		set_task( create_ptask( pose ) );
	} else {
		set_task( create_ptask(pose) ); //temp task
		if ( get_start_res_vec().empty() ) {
			for ( utility::vector1<Size>::iterator jit = jump_nums_.begin(); jit != jump_nums_.end(); ++jit ) {
				protocols::scoring::InterfaceOP interf = protocols::scoring::InterfaceOP( new protocols::scoring::Interface( *jit ) );
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
							if ( std::find(pair_lists_vec[ res1_ind ].begin(),pair_lists_vec[ res1_ind ].end(),res2_ind) == pair_lists_vec[ res1_ind ].end() ) {
								pair_lists_vec[ res1_ind ].push_back( res2_ind );
							}
							if ( std::find(pair_lists_vec[ res2_ind ].begin(),pair_lists_vec[ res2_ind ].end(),res1_ind) == pair_lists_vec[ res2_ind ].end() ) {
								pair_lists_vec[ res2_ind ].push_back( res1_ind );
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
						std::list< Size > templist1 = pair_lists_vec[ res1_ind ];
						if ( std::find(templist1.begin(),templist1.end(),res2_ind) == templist1.end() ) {
							pair_lists_vec[ res1_ind ].push_back( res2_ind );
						}
						std::list< Size > templist2 = pair_lists_vec[ res2_ind ];
						if ( std::find(templist2.begin(),templist2.end(),res1_ind) == templist2.end() ) {
							pair_lists_vec[ res2_ind ].push_back( res1_ind );
						}
						//std::pair<Size,Size> key(res1_ind,res2_ind);
						//interface_jump_map[key] = jump_nums_.size() + 1;

						add_start_res(res1_ind);
					}
					//else
					//    TR << "FALSE: res1 = " << *resnum << " and res2 = " << resnum2 << std::endl;
				}
			}
		}
		pose.update_residue_neighbors();
	}
	//jump_nums_.clear();
}//setup

bool
HBNetStapleInterface::has_pH_His( core::pose::Pose & pose, hbond_net_struct & i )
{
	if ( !(network_contains_aa( 'H', i ) ) ) {
		return false;
	}

	if ( i.hbond_vec.empty() ) {
		//i.hbond_vec = ( native_ ) ? get_hbond_atom_pairs( i.residues, pose ) : get_hbond_atom_pairs( i.residues, pose, get_packer_graph() );
		//i.hbond_vec = get_hbond_atom_pairs( i.residues, pose, get_packer_graph() );
		//i.hbond_vec = get_hbond_atom_pairs( i.residues, pose );
		get_hbond_atom_pairs( i, pose );
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
HBNetStapleInterface::network_meets_criteria( Pose & pose, hbond_net_struct & i )
{
	if ( pH_His_ && !(has_pH_His( pose, i )) ) {
		return false;
	}

	if ( his_tyr_ && ( ( !(network_contains_aa('Y', i)) || !(network_contains_aa('H', i) )) || !(his_tyr_connectivity( pose, i )) ) ) { //need to generalize for user-specified connectivites
		return false;
	}

	if ( min_intermolecular_hbonds_ ) {
		//i.num_intermolecular_hbs = (symmetric_) ? symm_num_intermolecular_hbonds( i, pose ) : num_intermolecular_hbonds( i, pose );
		i.num_intermolecular_hbs = num_intermolecular_hbonds( i, pose ); //NEED TO FIX THIS FOR FULL-NETWORK SYMMETRY CASES
		//if (verbose_)
		//    TR << "i.num_intermolecular_hbs = " << i.num_intermolecular_hbs << std::endl;
		if ( i.num_intermolecular_hbs < min_intermolecular_hbonds_ ) {
			return false;
		}
	}

	//TR << "num_helices_contacted = " << num_helices_w_hbond( i ) << std::endl;
	if ( ( min_helices_contacted_by_network_ ) && ( num_helices_w_hbond( i ) < min_helices_contacted_by_network_ ) ) {
		//if (verbose_)
		//    TR << "min_helices_contacted_by_network_ = " << num_helices_w_hbond( i ) << std::endl;
		return false;
	}
	return true;
}

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

void
HBNetStapleInterface::prepare_output()
{
	output_networks((min_networks_per_pose_ == 1)); //will add single networks to output vector
    if ( get_extend_existing_networks() ) std::sort( get_net_vec().begin(), get_net_vec().end(), compare_net_vec() ); //sort all networks to put extended first
	if ( max_networks_per_pose_ > 1 ) { //add multiple network sets to output vector
		//use_jd2_out_num_=0;
		for ( std::vector< HBondNetStructOP >::const_iterator netit = get_net_vec().begin(); netit != get_net_vec().end(); ++netit ) {
			std::string network( (pdb_numbering() ) ? ( print_list_to_string( get_orig_pose(), (*netit)->residues) ) : (print_list_to_string( (*netit)->residues) ) );
			if ( TR.visible() ) TR << "combining networks " << (*netit)->id << ": " << network;

			HBondNetStructOP new_network = *netit;
			//rec_add_staple(netit, new_network, 1);
            std::vector< Size > net_ids(0);
            net_ids.push_back( new_network->id );
            rec_add_staple( netit, net_ids, 1 );
		}
	}
}

bool
HBNetStapleInterface::network_spans_all_helices( hbond_net_struct & i )
{
	if ( num_helices_w_hbond( i ) == helix_boundaries_.size() ) {
		return true;
	}
	return false;
}

Size
HBNetStapleInterface::num_helices_w_hbond( hbond_net_struct & i )
{
	return ( i.asymm_residues.empty() ) ? num_helices_w_hbond( i.residues ) : num_helices_w_hbond( i.asymm_residues );
}

Size
HBNetStapleInterface::num_helices_w_hbond( utility::vector1< HBondResStructCOP > const & residues )
{
	Size num_helices(0);
	utility::vector1< bool > helix_has_hbond_residue( helix_boundaries_.size(), false );
	for ( utility::vector1< HBondResStructCOP >::const_iterator r = residues.begin(); r != residues.end(); ++r ) {
		Size helix_id( get_helix_id( (*r)->resnum ) ); //returns 0 if residues is not part of a helix
		if ( helix_id > 0 ) {
			helix_has_hbond_residue[ helix_id ] = true;
		}
	}
	for ( utility::vector1< bool >::iterator hel = helix_has_hbond_residue.begin(); hel != helix_has_hbond_residue.end(); ++hel ) {
		if ( *hel ) {
			num_helices++;
		}
	}
	return num_helices;
}

void
//HBNetStapleInterface::rec_add_staple( std::vector< HBondNetStructOP >::const_iterator netit, HBondNetStructOP new_network, Size staple_count )
HBNetStapleInterface::rec_add_staple( std::vector< HBondNetStructOP >::const_iterator netit, std::vector< Size > net_ids, Size staple_count )
{
    //runtime_assert( get_network_by_id(new_network_id) != nullptr );
	//numeric::xyzVector<core::Real> const & first_begin_coordinates = (*netit)->rotlist.front()->atom("CA").xyz();
	//numeric::xyzVector<core::Real> const & first_end_coordinates = (*netit)->rotlist.back()->atom("CA").xyz();
	Size combo_count(1);
	std::vector< HBondNetStructOP >::const_iterator next_netit(netit);
	while ( combo_count <= combos_ && ++next_netit != get_net_vec().end() )
    { //number of combinations of multiple networks to try (default = 1)
        bool compatible( true );
        for ( auto net_id = net_ids.begin(); net_id != net_ids.end(); ++ net_id ){
            runtime_assert( get_network_by_id(*net_id) != nullptr );
            bool branch(false);
//            if ( !(is_sub_residues( (get_network_by_id(*net_id))->residues, (get_network_by_id(*net_id))->residues, branch ))
//                && !branch && !(net_clash( *(get_network_by_id(*net_id)), **next_netit )) ) {
            if ( is_sub_residues( (get_network_by_id(*net_id))->residues, (get_network_by_id(*net_id))->residues, branch )
                || branch || net_clash( *(get_network_by_id(*net_id)), **next_netit ) ) {
                compatible = false;
                break;
            }
        }
        if (compatible ){
            //Size net_index1 = next_netit - get_net_vec().begin();
            std::string network( (pdb_numbering() ) ? ( print_list_to_string( get_orig_pose(), (*next_netit)->residues) ) : (print_list_to_string( (*next_netit)->residues) ) );
            if ( TR.visible() ) TR << "; and " << (*next_netit)->id << ": " << network;
            net_ids.push_back( (*next_netit)->id );
            std::vector< Size > new_net_ids( net_ids );
            
            //			HBondNetStructOP merged_nets( new hbond_net_struct() );
            //            //TODO NEED BETTER SOLUTION FOR COMBINING NETWORKS THAN MERGING INTO SINGLE NETWORK
            //			merge_2_networks( *new_network, **next_netit, merged_nets );
            //			merged_nets->score = (new_network->score + (*next_netit)->score)/2;
            //			//Size net_index2 = next_netit - network_vector_.begin();
            //			HBondNetStructOP new_out_struct( new hbond_net_struct(*merged_nets) );
            
            staple_count++;
            
            if ( staple_count >= min_networks_per_pose_ ) {
                //get_output_net_vec().push_back( new_out_struct );
                get_output_vector().push_back( new_net_ids );
            }
            
            if ( staple_count < max_networks_per_pose_ ) {
                rec_add_staple(next_netit, net_ids, staple_count);
                staple_count--;
            }
            combo_count++;
        }
	}
	TR.flush();
}

bool
HBNetStapleInterface::same_helix(utility::vector1< std::pair<Size,Size> > helix_boundaries, Size r1, Size r2)
{
	for ( utility::vector1< std::pair<Size,Size> >::iterator h = helix_boundaries.begin(); h != helix_boundaries.end(); ++h ) {
		if ( (r1 >= h->first && r1 <= h->second) && (r2 >= h->first && r2 <= h->second) ) {
			return true;
		}
	}
	return false;
}

//returns false if either of the residues are not Helix (H)
//returns true if residues contact eachother from different helices (interhelical contact)
bool
HBNetStapleInterface::interhelical_contact(utility::vector1< std::pair<Size,Size> > helix_boundaries, Size r1, Size r2, Pose & pose)
{
	if ( pose.residue(r1).nbr_atom_xyz().distance_squared(pose.residue(r2).nbr_atom_xyz()) < (interf_distance_*interf_distance_) ) {
		if ( !same_helix(helix_boundaries,r1,r2) ) {
			return true;
		}
	}
	return false;
}

//returns 0 if residues is not part of a helix
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
HBNetStapleInterface::print_additional_info_for_net( hbond_net_struct & i )
{
	std::stringstream output;
	if ( min_intermolecular_hbonds_ > 0 ) {
		output << i.num_intermolecular_hbs << "\t";
	}
	if ( min_helices_contacted_by_network_ > 0 ) {
		output << num_helices_w_hbond( i ) << "\t";
	}
	return output.str();
}

Size
HBNetStapleInterface::num_intermolecular_hbonds( hbond_net_struct & i, core::pose::Pose & pose )
{
	if ( i.hbond_vec.empty() ) {
		get_hbond_atom_pairs( i, pose );
		i.total_hbonds = i.hbond_vec.size();
	}
	//return get_intermolecular_hbonds( i );

	Size num_intermol_hbs(0);
	for ( utility::vector1<HBondCOP>::const_iterator h = i.hbond_vec.begin(); h != i.hbond_vec.end(); ++h ) {
		Size arsd((*h)->acc_res());
		Size drsd((*h)->don_res());
		if ( pose.chain(arsd) != pose.chain(drsd) ) {
			num_intermol_hbs++;
		}
	}
	return num_intermol_hbs;
}

} //hbnet
} //protocols
