// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/AtomicContactCountFilter.cc
/// @brief
/// @author Alex Ford (fordas@uw.edu)

#include <algorithm>
#include <iterator>

#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>
#include <utility/string_util.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/make_symmetric_task.hh>

#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/protein_interface_design/filters/AtomicContactCountFilter.hh>
#include <protocols/protein_interface_design/filters/AtomicContactCountFilterCreator.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.AtomicContactCountFilter" );


AtomicContactCountFilter::AtomicContactCountFilter() :
	protocols::filters::Filter( "AtomicContactCount" ),
	distance_cutoff_(4.5)
{
	initialize_all_atoms(NULL);
}

AtomicContactCountFilter::AtomicContactCountFilter(core::Real distance_cutoff) :
	Filter( "AtomicContactCount" ),
	distance_cutoff_(distance_cutoff)
{
	initialize_all_atoms(NULL);
}

AtomicContactCountFilter::AtomicContactCountFilter( AtomicContactCountFilter const & src ) :
	Filter(src),
	task_factoryA_(src.task_factoryA_),
	task_factoryB_(src.task_factoryB_),
	distance_cutoff_(src.distance_cutoff_),
	filter_mode_(src.filter_mode_),
	normalize_by_sasa_(src.normalize_by_sasa_),
	ss_only_(src.ss_only_),
	normalize_by_carbon_count_(src.normalize_by_carbon_count_),
	jump_(src.jump_),
	sym_dof_name_(src.sym_dof_name_)
{
}

AtomicContactCountFilter::~AtomicContactCountFilter() {}

protocols::filters::FilterOP AtomicContactCountFilter::clone() const { return protocols::filters::FilterOP( new AtomicContactCountFilter( *this ) ); }
protocols::filters::FilterOP AtomicContactCountFilter::fresh_instance() const { return protocols::filters::FilterOP( new AtomicContactCountFilter() ); }

void AtomicContactCountFilter::initialize_all_atoms( core::pack::task::TaskFactoryOP task_factoryA, bool individual_tasks, core::pack::task::TaskFactoryOP task_factoryB, bool normalize_by_carbon_count)
{
	task_factoryA_ = task_factoryA;
	individual_tasks_ = individual_tasks;
	task_factoryB_ = task_factoryB;
	normalize_by_carbon_count_ = normalize_by_carbon_count;
	jump_ = 0;
	sym_dof_name_ = "";
	normalize_by_sasa_ = false;
	ss_only_ = false;

	filter_mode_ = ALL;
}

void AtomicContactCountFilter::initialize_cross_jump(core::Size jump, std::string sym_dof_name, core::pack::task::TaskFactoryOP task_factoryA, bool normalize_by_sasa, bool individual_tasks, core::pack::task::TaskFactoryOP task_factoryB, bool normalize_by_carbon_count)
{
	jump_ = jump;
	sym_dof_name_ = sym_dof_name;
	task_factoryA_ = task_factoryA;
	normalize_by_sasa_ = normalize_by_sasa;
	individual_tasks_ = individual_tasks;
	task_factoryB_ = task_factoryB;
	normalize_by_carbon_count_ = normalize_by_carbon_count;

	filter_mode_ = CROSS_JUMP;
}
void AtomicContactCountFilter::initialize_cross_chain( core::pack::task::TaskFactoryOP task_factoryA, bool normalize_by_sasa, bool detect_chains_for_interface, bool individual_tasks, core::pack::task::TaskFactoryOP task_factoryB, bool normalize_by_carbon_count)
{
	task_factoryA_ = task_factoryA;
	normalize_by_sasa_ = normalize_by_sasa;
	individual_tasks_ = individual_tasks;
	task_factoryB_ = task_factoryB;
	jump_ = 0;
	sym_dof_name_ = "";
	normalize_by_carbon_count_ = normalize_by_carbon_count;

	filter_mode_ = detect_chains_for_interface ? CROSS_CHAIN_DETECTED : CROSS_CHAIN_ALL;
}

void AtomicContactCountFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	distance_cutoff_ = tag->getOption< core::Real >( "distance", 4.5 );

	std::string specified_mode = tag->getOption< std::string >( "partition", "none" );
	std::string specified_normalized_by_sasa = tag->getOption< std::string >( "normalize_by_sasa", "0" );
	bool normalize_by_sasa = tag->getOption< bool >( "normalize_by_sasa", false );
	bool normalize_by_carbon_count = tag->getOption< bool >( "normalize_by_carbon_count", false );
	if ( normalize_by_sasa && normalize_by_carbon_count ) {
		utility_exit_with_message("Can specify normalize_by_sasa or normalize_by_carbon_count, but not both.");
	}

	ss_only_ = tag->getOption< bool >( "ss_only", false );

	bool individual_tasks = false;
	core::pack::task::TaskFactoryOP task_factoryA( new core::pack::task::TaskFactory );
	core::pack::task::TaskFactoryOP task_factoryB( new core::pack::task::TaskFactory );
	if ( tag->hasOption( "taskA" ) ) {
		utility::vector1< std::string > taskA_names = utility::string_split( tag->getOption< std::string >( "taskA" ), ',' );
		for ( Size i = 1; i <= taskA_names.size(); i++ ) {
			task_factoryA->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", taskA_names[i] ) );
		}
		if ( tag->hasOption( "taskB" ) ) {
			utility::vector1< std::string > taskB_names = utility::string_split( tag->getOption< std::string >( "taskB" ), ',' );
			for ( Size i = 1; i <= taskB_names.size(); i++ ) {
				task_factoryB->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", taskB_names[i] ) );
			}
		} else {
			utility_exit_with_message("Must specify both TaskA and TaskB if using indivual tasks.");
		}
		if ( tag->hasOption( "task_operations" ) ) {
			utility_exit_with_message("Cannot specify task_operations together with individual tasks TaskA and TaskB.");
		}
		individual_tasks = true;
	} else if ( tag->hasOption( "taskB" ) ) {
		utility_exit_with_message("Must specify both TaskA and TaskB if using individual tasks.");
	} else if ( tag->hasOption( "task_operations" ) ) {
		task_factoryA = protocols::rosetta_scripts::parse_task_operations( tag, data );
		task_factoryB = protocols::rosetta_scripts::parse_task_operations( tag, data );
	}

	if ( specified_mode == "none" ) {
		initialize_all_atoms(
			task_factoryA,
			individual_tasks,
			task_factoryB,
			normalize_by_carbon_count);
	} else if ( specified_mode == "jump" ) {
		initialize_cross_jump(
			tag->getOption< core::Size >( "jump", 1 ),
			tag->getOption< std::string >( "sym_dof_name", "" ),
			task_factoryA,
			specified_normalized_by_sasa != "0",
			individual_tasks,
			task_factoryB,
			normalize_by_carbon_count);
	} else if ( specified_mode == "chain" ) {

		initialize_cross_chain(
			task_factoryA,
			specified_normalized_by_sasa != "0",
			specified_normalized_by_sasa == "detect_by_task",
			individual_tasks,
			task_factoryB,
			normalize_by_carbon_count);
	}

	if ( filter_mode_ == ALL && specified_normalized_by_sasa != "0" ) {
		TR.Error << "Must specify jump or chain partition mode in AtomicContactFilter with normalize_by_sasa: " << tag << std::endl;
		utility_exit_with_message("Must specify jump or chain partition mode in AtomicContactFilter with normalize_by_sasa.");
	}

	TR.Debug << "Parsed AtomicContactCount filter: <AtomicContactCount" <<
		" distance=" << distance_cutoff_ <<
		" jump=" << jump_ <<
		" normalize_by_sasa=" << (normalize_by_sasa_ ? (filter_mode_ != CROSS_CHAIN_DETECTED ? "1" : "detect_by_task") : "0") <<
		" />" << std::endl;
}

core::Real AtomicContactCountFilter::compute(core::pose::Pose const & pose) const
{
	// Create map of taskA and taskB using taskoperations
	core::pack::task::PackerTaskOP taskA = core::pack::task::TaskFactory::create_packer_task( pose );
	core::pack::task::PackerTaskOP taskB = core::pack::task::TaskFactory::create_packer_task( pose );

	if ( task_factoryA_ != 0 ) {
		taskA = task_factoryA_->create_task_and_apply_taskoperations( pose );
		TR << "Initializing taskA from packer task." << std::endl;
		//TR.Debug << "Initializing taskA from packer task." << std::endl;
	} else {
		TR << "No packer taskA specified, using default task." << std::endl;
		//TR.Debug << "No packer taskA specified, using default task." << std::endl;
	}
	if ( task_factoryB_ != 0 ) {
		taskB = task_factoryB_->create_task_and_apply_taskoperations( pose );
		TR << "Initializing taskB from packer task." << std::endl;
		//TR.Debug << "Initializing taskB from packer task." << std::endl;
	} else {
		TR << "No packer taskB specified, using default task." << std::endl;
		//TR.Debug << "No packer taskB specified, using default task." << std::endl;
	}

	bool symmetric = core::pose::symmetry::is_symmetric( pose );

	if ( symmetric ) {
		taskA = core::pack::make_new_symmetric_PackerTask_by_requested_method(pose, taskA);
		taskB = core::pack::make_new_symmetric_PackerTask_by_requested_method(pose, taskB);
	}

	// Create lookup of setA and setB residues
	utility::vector1<core::Size> setA, setB;
	for ( core::Size resi = 1; resi <= pose.n_residue(); resi++ ) {
		if ( taskA->pack_residue(resi) ) {
			setA.push_back(resi);
		}
		if ( taskB->pack_residue(resi) ) {
			setB.push_back(resi);
		}
	}
	if ( !individual_tasks_ ) {
		setB = setA;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "SetA residues from task: ";
		std::copy(setA.begin(), setA.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
		TR.Debug << std::endl;
		TR.Debug << "SetB residues from task: ";
		std::copy(setB.begin(), setB.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
		TR.Debug << std::endl;
	}

	// Divide residues into partitions based on filter mode
	utility::vector1<core::Size> residue_partition;
	utility::vector1<core::Size> target_jumps;

	if ( filter_mode_ == ALL ) {
		TR.Debug << "Partitioning by residue number." << std::endl;
		// Each residue is a separate partition
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			residue_partition.push_back(i);
		}
	} else if ( filter_mode_ == CROSS_CHAIN_DETECTED || filter_mode_ == CROSS_CHAIN_ALL ) {
		TR.Debug << "Partitioning by residue chain." << std::endl;

		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			residue_partition.push_back(pose.chain(i));
		}
	} else if ( filter_mode_ == CROSS_JUMP ) {
		TR << "Partitioning by jump." << std::endl;
		//TR.Debug << "Partitioning by jump." << std::endl;

		// Lookup symmetry-aware jump identifier
		if ( sym_dof_name_ != "" ) {
			target_jumps.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_ ) );
		} else if ( !symmetric ) {
			target_jumps.push_back( jump_ );
		} else {
			// all slidable jumps
			Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
			for ( Size j = 1; j <= nslidedofs; j++ ) {
				target_jumps.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
			}
		}

		// Partition pose by jump
		ObjexxFCL::FArray1D<bool> jump_partition ( pose.total_residue(), false );
		if ( !symmetric ) {
			pose.fold_tree().partition_by_jump( target_jumps[1], jump_partition );
		} else {
			core::pose::symmetry::partition_by_symm_jumps( target_jumps, pose.fold_tree(),
				core::pose::symmetry::symmetry_info(pose), jump_partition );
		}

		for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
			residue_partition.push_back(jump_partition(i));
		}
	}

	if ( TR.Trace.visible() ) {
		TR.Trace << "Residue partitions from task: ";
		for ( core::Size i = 1; i <= residue_partition.size(); ++i ) {
			TR.Trace << i << ":" << residue_partition[i] << ",";
		}
		TR.Trace << std::endl;
	}

	// Count all cross-partition contacts
	core::Size contact_count = 0;
	core::Size carbon_count = 0;
	utility::vector1<bool>  indy_resis;
	if ( symmetric ) {
		core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		indy_resis = symm_info->independent_residues();
	}

	std::string pose_ss;
	if ( ss_only_ ) {
		// get secstruct
		core::scoring::dssp::Dssp dssp( pose );
		pose_ss = dssp.get_dssp_reduced_IG_as_L_secstruct();
	}
	TR << "Entering outer loop" << std::endl; // Remove
	for ( core::Size i = 1; i <= setA.size(); i++ ) {
		if ( symmetric && !indy_resis[setA[i]] ) {
			continue;
		}
		core::conformation::Residue const & residue_i = pose.residue(setA[i]);
		if ( normalize_by_carbon_count_ ) {
			for ( core::Size atom_i = residue_i.first_sidechain_atom(); atom_i <= residue_i.nheavyatoms(); atom_i++ ) {
				if ( residue_i.atom_type(atom_i).element() == "C" ) {
					carbon_count += 1;
				}
			}
		}
		core::Size start_index = 1;
		if ( !individual_tasks_ ) {
			start_index = i+1;
		}
		TR << "Entering inner loop" << std::endl;
		for ( core::Size j = start_index; j <= setB.size(); j++ ) {
			//fpd ss filter
			//jbb shouldn't this be target[i]-1 and target[j]-1 not i-1 and j-1?
			//if (ss_only_ && (pose_ss[i-1] == 'L' || pose_ss[j-1] == 'L') ) continue;
			if ( ss_only_ && (pose_ss[setA[i]-1] == 'L' || pose_ss[setB[j]-1] == 'L') ) continue;

			//if (target[i] != j)
			//if (residue_partition[target[i]] != residue_partition[target[j]])
			if ( residue_partition[setA[i]] != residue_partition[setB[j]] ) {
				core::conformation::Residue const & residue_j = pose.residue(setB[j]);

				for ( core::Size atom_i = residue_i.first_sidechain_atom(); atom_i <= residue_i.nheavyatoms(); atom_i++ ) {
					if ( residue_i.atom_type(atom_i).element() != "C" ) {
						continue;
					}

					for ( core::Size atom_j = residue_j.first_sidechain_atom(); atom_j <= residue_j.nheavyatoms(); atom_j++ ) {
						if ( residue_j.atom_type(atom_j).element() != "C" ) {
							continue;
						}

						if ( residue_i.xyz(atom_i).distance(residue_j.xyz(atom_j)) <= distance_cutoff_ ) {
							TR << "select (resi " << setA[i] << " and name " << residue_i.atom_name(atom_i) << ") + (resi " << setB[j] << " and name " << residue_j.atom_name(atom_j) << ")" << std::endl;
							//TR.Debug << "select (resi " << target[i] << " and name " << residue_i.atom_name(atom_i) << ") + (resi " << target[j] << "and name " << residue_j.atom_name(atom_j) << ")" << std::endl;
							contact_count += 1;
						}
					}
				}
			}
		}
	}

	TR.Debug << "Found cross partition contacts:" << contact_count << std::endl;

	if ( normalize_by_sasa_ ) {
		TR.Debug << "Normalizing cross partition contacts by sasa." << std::endl;

		core::Real interface_sasa = 0.0;

		if ( filter_mode_ == ALL ) {
			utility_exit_with_message("AtomicContactCount filter unable to normalize by sasa in all-contact mode.");
		} else if ( filter_mode_ == CROSS_CHAIN_DETECTED ||
				filter_mode_ == CROSS_CHAIN_ALL ) {
			std::set<core::Size> interface_chains;

			if ( filter_mode_ == CROSS_CHAIN_DETECTED ) {
				// Detect chains containing target residues.
				for ( core::Size i = 1; i <= setA.size(); ++i ) {
					interface_chains.insert(pose.chain(setA[i]));
				}
			} else if ( filter_mode_ == CROSS_CHAIN_ALL ) {
				// Add all pose chains
				for ( core::Size i = 1; i <= pose.conformation().num_chains(); ++i ) {
					interface_chains.insert(i);
				}
			}

			TR.Debug << "Identified interface chains:";
			std::copy(interface_chains.begin(), interface_chains.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
			TR.Debug << std::endl;

			if ( interface_chains.size() != pose.conformation().num_chains() ) {
				TR.Debug << "Pruning target pose to interface chains.";
				core::pose::Pose interface_pose(pose);

				// Delete chains not in the chain list in reverse order in order to preserve chain number during deletion.
				for ( core::Size i = interface_pose.conformation().num_chains(); i >= 1; --i ) {
					if ( interface_chains.count(i) == 0 ) {
						TR.Debug << "Pruning target pose chain: " << i <<
							"[" << interface_pose.conformation().chain_begin( i ) <<
							"]" << interface_pose.conformation().chain_end( i ) << std::endl;

						interface_pose.conformation().delete_residue_range_slow(
							interface_pose.conformation().chain_begin( i ),
							interface_pose.conformation().chain_end( i ));
					}
				}

				if ( TR.Trace.visible() ) {
					TR.Trace << "Pruned pose:" << std::endl << interface_pose << std::endl;
				}

				// Calculate sasa across all remaining jumps
				protocols::simple_filters::InterfaceSasaFilter sasa_filter;
				utility::vector1<core::Size> sasa_jumps;
				for ( core::Size i = 1; i <= interface_pose.num_jump(); ++i ) {
					TR.Debug << "Adding jump to sasa filter: " << i << std::endl;
					sasa_jumps.push_back(i);
				}
				sasa_filter.jumps(sasa_jumps);

				interface_sasa = sasa_filter.compute(interface_pose);
			} else {
				TR.Debug << "Using original pose, interface chains include entire pose." << std::endl;

				// Calculate sasa across all jumps
				protocols::simple_filters::InterfaceSasaFilter sasa_filter;
				utility::vector1<core::Size> sasa_jumps;
				for ( core::Size i = 1; i <= pose.num_jump(); ++i ) {
					TR.Debug << "Adding jump to sasa filter: " << i << std::endl;
					sasa_jumps.push_back(i);
				}
				sasa_filter.jumps(sasa_jumps);
				interface_sasa = sasa_filter.compute(pose);
			}
		} else if ( filter_mode_ == CROSS_JUMP ) {
			TR.Debug << "Normalizing on jump." << std::endl;

			// Calculate sasa across the specified jump
			protocols::simple_filters::InterfaceSasaFilter sasa_filter;

			TR.Debug << "Adding jump to sasa filter: " << target_jumps[1] << std::endl;
			if ( symmetric && (sym_dof_name_ != "") ) {
				sasa_filter.sym_dof_names(sym_dof_name_);
			} else {
				sasa_filter.jumps(target_jumps);
			}
			interface_sasa = sasa_filter.compute(pose);
		}

		TR.Debug << "Calculated interface sasa: " << interface_sasa << std::endl;

		if ( interface_sasa != 0 ) {
			return contact_count / interface_sasa;
		} else {
			return 0;
		}
	} else if ( normalize_by_carbon_count_ ) {
		return (core::Real)(contact_count) / (core::Real)(carbon_count);
	} else {
		return contact_count;
	}
}

protocols::filters::FilterOP
AtomicContactCountFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AtomicContactCountFilter ); }

std::string
AtomicContactCountFilterCreator::keyname() const { return "AtomicContactCount"; }

}
}
}
