// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
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

static basic::Tracer TR("protocols.protein_interface_design.filters.AtomicContactCountFilter");


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

AtomicContactCountFilter::AtomicContactCountFilter( AtomicContactCountFilter const & src ) : Filter(src),
		task_factory_(src.task_factory_),
		distance_cutoff_(src.distance_cutoff_),
		filter_mode_(src.filter_mode_),
		normalize_by_sasa_(src.normalize_by_sasa_),
		ss_only_(src.ss_only_),
		jump_(src.jump_),
		sym_dof_name_(src.sym_dof_name_)
{}

AtomicContactCountFilter::~AtomicContactCountFilter() {}

protocols::filters::FilterOP AtomicContactCountFilter::clone() const { return new AtomicContactCountFilter( *this ); }
protocols::filters::FilterOP AtomicContactCountFilter::fresh_instance() const { return new AtomicContactCountFilter(); }

void AtomicContactCountFilter::initialize_all_atoms(core::pack::task::TaskFactoryOP task_factory)
{
	task_factory_ = task_factory;
	jump_ = 0;
	sym_dof_name_ = "";
	normalize_by_sasa_ = false;
	ss_only_ = false;

	filter_mode_ = ALL;
}

void AtomicContactCountFilter::initialize_cross_jump(core::Size jump, std::string sym_dof_name, core::pack::task::TaskFactoryOP task_factory, bool normalize_by_sasa)
{
	task_factory_ = task_factory;
	jump_ = jump;
	sym_dof_name_ = sym_dof_name;
	normalize_by_sasa_ = normalize_by_sasa;

	filter_mode_ = CROSS_JUMP;
}

void AtomicContactCountFilter::initialize_cross_chain(core::pack::task::TaskFactoryOP task_factory, bool normalize_by_sasa, bool detect_chains_for_interface)
{
	task_factory_ = task_factory;
	jump_ = 0;
	sym_dof_name_ = "";
	normalize_by_sasa_ = normalize_by_sasa;

	filter_mode_ = detect_chains_for_interface ? CROSS_CHAIN_DETECTED : CROSS_CHAIN_ALL;
}

void AtomicContactCountFilter::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	distance_cutoff_ = tag->getOption< core::Real >( "distance", 4.5 );

	std::string specified_mode = tag->getOption< std::string >( "partition", "none" );
	std::string specified_normalized_by_sasa = tag->getOption< std::string >( "normalize_by_sasa", "0" );

	ss_only_ = tag->getOption< bool >( "ss_only", false );

	if (specified_mode == "none")
	{
		initialize_all_atoms(protocols::rosetta_scripts::parse_task_operations( tag, data ));
	}
	else if (specified_mode == "jump")
	{
		initialize_cross_jump(
				tag->getOption< core::Size >( "jump", 1 ),
				tag->getOption< std::string >( "sym_dof_name", "" ),
				protocols::rosetta_scripts::parse_task_operations( tag, data ),
				specified_normalized_by_sasa != "0");
	}
	else if (specified_mode == "chain")
	{

		initialize_cross_chain(
				protocols::rosetta_scripts::parse_task_operations( tag, data ),
				specified_normalized_by_sasa != "0",
				specified_normalized_by_sasa == "detect_by_task");
	}

	if (filter_mode_ == ALL && specified_normalized_by_sasa != "0")
	{
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
	// Create map of target residues using taskoperation
  core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

  if ( task_factory_ != 0 )
	{
    task = task_factory_->create_task_and_apply_taskoperations( pose );
		TR.Debug << "Initializing from packer task." << std::endl;
  }
	else
	{
		TR.Debug << "No packer task specified, using default task." << std::endl;
  }

	bool symmetric = core::pose::symmetry::is_symmetric( pose );

	if ( symmetric )
	{
		task = core::pack::make_new_symmetric_PackerTask_by_requested_method(pose, task);
	}

	// Create lookup of target residues
	utility::vector1<core::Size> target;
	for (core::Size resi = 1; resi <= pose.n_residue(); resi++)
	{
		if( task->pack_residue(resi) )
		{
			target.push_back(resi);
		}
	}

	if (TR.Debug.visible())
	{
		TR.Debug << "Targets from task: ";
		std::copy(target.begin(), target.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
		TR.Debug << std::endl;
	}

	// Divide residues into partitions based on filter mode
	utility::vector1<core::Size> residue_partition;
	utility::vector1<core::Size> target_jumps;

	if (filter_mode_ == ALL)
	{
		TR.Debug << "Partitioning by residue number." << std::endl;
		// Each residue is a separate partition
		for (core::Size i = 1; i <= pose.total_residue(); ++i)
		{
			residue_partition.push_back(i);
		}
	}
	else if ( filter_mode_ == CROSS_CHAIN_DETECTED || filter_mode_ == CROSS_CHAIN_ALL )
	{
		TR.Debug << "Partitioning by residue chain." << std::endl;

		for (core::Size i = 1; i <= pose.total_residue(); ++i)
		{
			residue_partition.push_back(pose.chain(i));
		}
	}
	else if (filter_mode_ == CROSS_JUMP)
	{
		TR.Debug << "Partitioning by jump." << std::endl;

		// Lookup symmetry-aware jump identifier
		if ( sym_dof_name_ != "" ) {
			target_jumps.push_back( core::pose::symmetry::sym_dof_jump_num( pose, sym_dof_name_ ) );
		} else if (!symmetric) {
			target_jumps.push_back( jump_ );
		} else {
			// all slidable jumps
			Size nslidedofs = core::pose::symmetry::symmetry_info(pose)->num_slidablejumps();
			for (Size j = 1; j <= nslidedofs; j++)
				target_jumps.push_back( core::pose::symmetry::get_sym_aware_jump_num(pose, j ) );
		}

		// Partition pose by jump
		ObjexxFCL::FArray1D<bool> jump_partition ( pose.total_residue(), false );
		if (!symmetric) {
			pose.fold_tree().partition_by_jump( target_jumps[1], jump_partition );
		} else {
			core::pose::symmetry::partition_by_symm_jumps( target_jumps, pose.fold_tree(),
				core::pose::symmetry::symmetry_info(pose), jump_partition );
		}

		for (core::Size i = 1; i <= pose.total_residue(); ++i)
		{
			residue_partition.push_back(jump_partition(i));
		}
	}

	if (TR.Trace.visible())
	{
		TR.Trace << "Residue partitions from task: ";
		for (core::Size i = 1; i <= residue_partition.size(); ++i)
		{
			TR.Trace << i << ":" << residue_partition[i] << ",";
		}
		TR.Trace << std::endl;
	}

	// Count all cross-partition contacts
	core::Size contact_count = 0;
  utility::vector1<bool>  indy_resis;
	if ( symmetric )
	{
  	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
		indy_resis = symm_info->independent_residues();
	}

	std::string pose_ss;
	if (ss_only_) {
		// get secstruct
		core::scoring::dssp::Dssp dssp( pose );
		pose_ss = dssp.get_dssp_reduced_IG_as_L_secstruct();
	}

	for (core::Size i = 1; i <= target.size(); i++)
	{
		if ( symmetric && (filter_mode_ == CROSS_JUMP) )
		{
			// TODO: ALEX FORD: Fix so that can take multiple tasks.  One that specifies the residues for which counts are analyzed and the other task specifies all of the other residues with which to look for interactions.
			// The residue_partition logic works in this case, but may not make sense for symmetric assemblies in cross_chain mode.
			// For multicomponent systems we only want to count contacts from residues in the primary subunit corresponding to the user-specified symdof
			if (!indy_resis[target[i]] || residue_partition[target[i]]) continue;
		}
		for (core::Size j = i+1; j <= target.size(); j++)
		{
			//fpd ss filter
			if (ss_only_ && (pose_ss[i-1] == 'L' || pose_ss[j-1] == 'L') ) continue;

			if (residue_partition[target[i]] != residue_partition[target[j]])
			{
				core::conformation::Residue const & residue_i = pose.residue(target[i]);
				core::conformation::Residue const & residue_j = pose.residue(target[j]);

				for (core::Size atom_i = residue_i.first_sidechain_atom(); atom_i <= residue_i.nheavyatoms(); atom_i++)
				{
					if (residue_i.atom_type(atom_i).element() != "C")
					{
						continue;
					}

					for (core::Size atom_j = residue_j.first_sidechain_atom(); atom_j <= residue_j.nheavyatoms(); atom_j++)
					{
						if (residue_j.atom_type(atom_j).element() != "C")
						{
							continue;
						}

						if (residue_i.xyz(atom_i).distance(residue_j.xyz(atom_j)) <= distance_cutoff_)
						{
							TR.Debug << "select (resi " << target[i] << " and name " << residue_i.atom_name(atom_i) << ") + (resi " << target[j] << "and name " << residue_j.atom_name(atom_j) << ")" << std::endl;
							contact_count += 1;
						}
					}
				}
			}
		}
	}

	TR.Debug << "Found cross partition contacts:" << contact_count << std::endl;

	if (normalize_by_sasa_)
	{
		TR.Debug << "Normalizing cross partition contacts by sasa." << std::endl;

		core::Real interface_sasa = 0.0;

		if (filter_mode_ == ALL)
		{
			utility_exit_with_message("AtomicContactCount filter unable to normalize by sasa in all-contact mode.");
		}
		else if (filter_mode_ == CROSS_CHAIN_DETECTED ||
				filter_mode_ == CROSS_CHAIN_ALL)
		{
			std::set<core::Size> interface_chains;

			if (filter_mode_ == CROSS_CHAIN_DETECTED)
			{
				// Detect chains containing target residues.
				for (core::Size i = 1; i <= target.size(); ++i)
				{
					interface_chains.insert(pose.chain(target[i]));
				}
			}
			else if (filter_mode_ == CROSS_CHAIN_ALL)
			{
				// Add all pose chains
				for (core::Size i = 1; i <= pose.conformation().num_chains(); ++i)
				{
					interface_chains.insert(i);
				}
			}

			TR.Debug << "Identified interface chains:";
			std::copy(interface_chains.begin(), interface_chains.end(), std::ostream_iterator<core::Size>(TR.Debug, ","));
			TR.Debug << std::endl;

			if (interface_chains.size() != pose.conformation().num_chains())
			{
				TR.Debug << "Pruning target pose to interface chains.";
				core::pose::Pose interface_pose(pose);

				// Delete chains not in the chain list in reverse order in order to preserve chain number during deletion.
				for (core::Size i = interface_pose.conformation().num_chains(); i >= 1; --i)
				{
					if (interface_chains.count(i) == 0)
					{
						TR.Debug << "Pruning target pose chain: " << i <<
							"[" << interface_pose.conformation().chain_begin( i ) <<
							"]" << interface_pose.conformation().chain_end( i ) << std::endl;

						interface_pose.conformation().delete_residue_range_slow(
							interface_pose.conformation().chain_begin( i ),
							interface_pose.conformation().chain_end( i ));
					}
				}

				if (TR.Trace.visible())
				{
					TR.Trace << "Pruned pose:" << std::endl << interface_pose << std::endl;
				}

				// Calculate sasa across all remaining jumps
				protocols::simple_filters::InterfaceSasaFilter sasa_filter;
				utility::vector1<core::Size> sasa_jumps;
				for (core::Size i = 1; i <= interface_pose.num_jump(); ++i)
				{
					TR.Debug << "Adding jump to sasa filter: " << i << std::endl;
					sasa_jumps.push_back(i);
				}
				sasa_filter.jumps(sasa_jumps);

				interface_sasa = sasa_filter.compute(interface_pose);
			}
			else
			{
				TR.Debug << "Using original pose, interface chains include entire pose." << std::endl;

				// Calculate sasa across all jumps
				protocols::simple_filters::InterfaceSasaFilter sasa_filter;
				utility::vector1<core::Size> sasa_jumps;
				for (core::Size i = 1; i <= pose.num_jump(); ++i)
				{
					TR.Debug << "Adding jump to sasa filter: " << i << std::endl;
					sasa_jumps.push_back(i);
				}
				sasa_filter.jumps(sasa_jumps);
				interface_sasa = sasa_filter.compute(pose);
			}
		}
		else if (filter_mode_ == CROSS_JUMP)
		{
			TR.Debug << "Normalizing on jump." << std::endl;

			// Calculate sasa across the specified jump
			protocols::simple_filters::InterfaceSasaFilter sasa_filter;

			TR.Debug << "Adding jump to sasa filter: " << target_jumps[1] << std::endl;
			if ( symmetric && (sym_dof_name_ != ""))
			{
				sasa_filter.sym_dof_names(sym_dof_name_);
			}
			else
			{
				sasa_filter.jumps(target_jumps);
			}
			interface_sasa = sasa_filter.compute(pose);
		}

		TR.Debug << "Calculated interface sasa: " << interface_sasa << std::endl;

		if (interface_sasa != 0)
		{
			return contact_count / interface_sasa;
		}
		else
		{
			return 0;
		}
	}
	else
	{
		return contact_count;
	}
}

protocols::filters::FilterOP
AtomicContactCountFilterCreator::create_filter() const { return new AtomicContactCountFilter; }

std::string
AtomicContactCountFilterCreator::keyname() const { return "AtomicContactCount"; }

}
}
}
