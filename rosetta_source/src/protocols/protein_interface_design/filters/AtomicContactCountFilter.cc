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

#include <ObjexxFCL/FArray1D.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/scoring/ScoreType.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/make_symmetric_task.hh>

#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/protein_interface_design/filters/AtomicContactCountFilter.hh>
#include <protocols/protein_interface_design/filters/AtomicContactCountFilterCreator.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {


AtomicContactCountFilter::AtomicContactCountFilter() :
	protocols::filters::Filter( "AtomicContactCount" ),
	task_factory_(NULL),
	jump_(1),
	distance_cutoff_(4.5),
	normalize_by_sasa_(false)
{
}

AtomicContactCountFilter::AtomicContactCountFilter(
		core::pack::task::TaskFactoryOP task_factory,
		core::Size const jump,
		core::Real const distance_cutoff,
		bool normalize_by_sasa) :
	Filter( "AtomicContactCount" ),
	task_factory_(task_factory),
	jump_(jump),
	distance_cutoff_(distance_cutoff),
	normalize_by_sasa_(normalize_by_sasa)
{
}

AtomicContactCountFilter::AtomicContactCountFilter( AtomicContactCountFilter const & src ) :
	Filter(src),
	task_factory_(src.task_factory_),
	jump_(src.jump_),
	distance_cutoff_(src.distance_cutoff_),
	normalize_by_sasa_(src.normalize_by_sasa_)
{
}

AtomicContactCountFilter::~AtomicContactCountFilter() {}

protocols::filters::FilterOP AtomicContactCountFilter::clone() const { return new AtomicContactCountFilter( *this ); }
protocols::filters::FilterOP AtomicContactCountFilter::fresh_instance() const { return new AtomicContactCountFilter(); }

void AtomicContactCountFilter::parse_my_tag(
	utility::tag::TagPtr const tag,
	protocols::moves::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
	jump_ = tag->getOption< core::Size >( "jump", 1 );
	distance_cutoff_ = tag->getOption< core::Real >( "distance", 4.5 );
	normalize_by_sasa_ = tag->getOption< bool >( "normalize_by_sasa", false );
}

core::Real AtomicContactCountFilter::compute(core::pose::Pose const & pose) const
{
	// Lookup symmetry-aware jump identifier
	core::Size target_jump = core::pose::symmetry::get_sym_aware_jump_num( pose, jump_ );
	
	// Create map of target residues using taskoperation
  core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );

  if ( task_factory_ != 0 )
	{
    task = task_factory_->create_task_and_apply_taskoperations( pose );
  }
	else
	{
		//TODO alexford Error on init, die.
  }

	if (core::pose::symmetry::is_symmetric( pose ))
	{
		task = core::pack::make_new_symmetric_PackerTask_by_requested_method(pose, task);
	}
	
	// Create lookup of target residues
	utility::vector1<core::Size> target;
	for (core::Size resi = 1; resi <= pose.n_residue(); resi++)
	{
		if( pose.residue(resi).is_protein() && task->pack_residue(resi) )
		{
			target.push_back(resi);
		}
	}

	// Partition pose by jump
	ObjexxFCL::FArray1D<bool> jump_partition ( pose.total_residue(), false );
	pose.fold_tree().partition_by_jump( target_jump, jump_partition );

	core::Size contact_count = 0;

	for (core::Size i = 1; i <= target.size(); i++)
	{
		for (core::Size j = i+1; j <= target.size(); j++)
		{
			if (jump_partition(target[i]) != jump_partition(target[j]))
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
							contact_count += 1;
						}
					}
				}
			}
		}
	}

	if (normalize_by_sasa_)
	{
		protocols::simple_filters::InterfaceSasaFilter sasa_filter;
		sasa_filter.jump(target_jump);

		core::Real interface_sasa = sasa_filter.compute(pose);

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
