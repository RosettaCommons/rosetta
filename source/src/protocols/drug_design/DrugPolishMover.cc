// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/DrugPolishMover.cc
/// @brief MonteCarlo protocol to design drugs in a protein context
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/drug_design/DrugPolishMover.hh>
#include <protocols/drug_design/DrugPolishMoverCreator.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <protocols/chemistries/Chemistry.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/AtomRefMapping.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/min_pack.hh>
#include <core/pack/rtmin.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>
#include <utility/tag/Tag.hh>

// Debugging output
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <utility/io/ozstream.hh>

// External headers
#include <boost/foreach.hpp>

// C/C++ headers
#include <string>

#define foreach BOOST_FOREACH

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.DrugPolishMover");

using namespace core;
using namespace protocols::moves;
using namespace protocols::filters;

std::string
DrugPolishMoverCreator::keyname() const
{
	return DrugPolishMoverCreator::mover_name();
}

protocols::moves::MoverOP
DrugPolishMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DrugPolishMover );
}

std::string
DrugPolishMoverCreator::mover_name()
{
	return "DrugPolishMover";
}

/// @brief default constructor
DrugPolishMover::DrugPolishMover() {}

/// @brief destructor
DrugPolishMover::~DrugPolishMover() {}

/// @brief clone this object
MoverOP
DrugPolishMover::clone() const
{
	return MoverOP( new DrugPolishMover( *this ) );
}

/// @brief create this type of object
MoverOP
DrugPolishMover::fresh_instance() const
{
	return MoverOP( new DrugPolishMover() );
}

void
DrugPolishMover::apply( Pose & pose )
{
	using namespace core::chemical;
	using namespace core::pose;

	if ( chemistries().size() != 1 ) {
		utility_exit_with_message("Something's wrong - chemistry not found for DrugPolishMover.");
	}
	if ( !redocker() ) {
		TR.Warning << "Redocking mover is empty! " << std::endl;
		return;
	}
	if ( !scorer() ) {
		TR.Warning << "No scoring Filter specified! " << std::endl;
		return;
	}

	Real initial_score = scorer()->report_sm( pose );
	core::Size res_pos( find_design_position( pose ) );
	TR << "Starting score for ligand " << pose.residue( res_pos ).name() << " is " << initial_score << std::endl;

	// Repeated application until we can't find an applicable reaction to use.
	for ( Size ii=1; ii<=maxtrials(); ii++ ) {
		TR << "Cycle number: " << ii <<std::endl;

		utility::vector1< core::pose::PoseOP > poses;
		utility::vector1< core::Real > scores;

		// Make a working copy of the restype of interest.
		res_pos = find_design_position( pose );
		MutableResidueTypeOP restype( new MutableResidueType( pose.residue_type( res_pos ) ) );
		IndexVDMapping starting_index_vd_mapping( combine( IndexNameMapping( pose.residue_type( res_pos ) ), NameVDMapping(*restype) ) ); // Starting index to current vds.

		std::string original_name( restype->name() ); // Need this as chemistries aren't necessarily all that good with naming.

		if ( pre_process_restype(restype, starting_index_vd_mapping, pose ) ) {
			TR.Warning << "Issue pre-processing ResidueType " << restype->name() << " in cycle " << ii << std::endl;
			// ResidueType is unusable for futher processing.
			break;
		}

		// Apply the chemistry
		debug_assert( chemistries()[1] );
		protocols::chemistries::Chemistry & chemistry( *chemistries()[1] );

		dump_molecule( *restype, "before_mod" );
		chemistry.apply( *restype, pose ); // Allow context sensitive
		dump_molecule( *restype, "after_mod" );

		core::Size subiteration( 0 );

		while ( restype && chemistry.get_last_status() == core::chemical::modifications::SUCCESS ) {
			++subiteration;

			IndexVDMapping index_vd_mapping( combine( starting_index_vd_mapping, chemistry.get_mapping() ) );

			std::string new_name( find_new_res_name( original_name, ii, subiteration ) );

			if ( post_process_restype(restype, index_vd_mapping, new_name, pose ) ) {
				// ResidueType is unusable
				continue;
			}

			core::pose::PoseOP working_pose( pose.clone() );

			if ( ! emplace_residue_type(*working_pose, restype, index_vd_mapping ) ) {
				// No errors on emplacement

				// Because the position may change during docking
				res_pos = find_design_position( *working_pose );
				dump_molecule( working_pose->residue( res_pos ), "tested" );

				Real score = scorer()->report_sm( *working_pose );
				debug_assert( poses.size() == scores.size() );
				poses.push_back( working_pose );
				scores.push_back( score );
				TR << "For cycle " << ii << " output " << poses.size() << " ligand " << working_pose->residue( res_pos ).name() << " score " << score << std::endl;
			}

			// Try next ResidueType
			restype = chemistry.get_additional_output();
		}

		// Figure out the best pose.
		if ( poses.size() == 0 ) {
			TR << "After " << ii << " cycles, no more applicable reactions." << std::endl;
			break;
		}

		debug_assert( poses.size() == scores.size() );
		core::Size selected( 0 );
		core::Real best_score( 999999999 );
		for ( core::Size jj(1); jj <= scores.size(); ++jj ) {
			if ( scores[jj] - bonus_ >= initial_score ) {
				continue;
			}
			if ( scores[jj] < best_score ) {
				best_score = scores[jj];
				selected = jj;
			}
		}

		if ( selected == 0 ) {
			TR << "After " << ii << " cycles, no products passed the score cutoffs. " << std::endl;
			break;
		}

		// Substitute in the best pose, and repeat the reactions.
		pose = *poses[ selected ];

	} // i<=maxtrials_

	TR << "DrugPolishMover finished." << std::endl;

	// Output the final docked ligand as sdf, for convenience
	res_pos = find_design_position( pose );
	dump_molecule( pose.residue( res_pos ), "final" );

}// apply

std::string
DrugPolishMover::get_name() const {
	return DrugPolishMoverCreator::mover_name();
}

/// @brief parse xml file
void
DrugPolishMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{
	maxtrials( tag->getOption< core::Size >( "maxcycles", 10 ) );
	chain( tag->getOption< char >( "chain", 'X' ) );

	if ( ! tag->hasOption( "scorer" ) ) {
		utility_exit_with_message("You must provide a scorer option to the DrugPolishMover for scoring.");
	}
	std::string const filter_name( tag->getOption< std::string >( "scorer" ) );
	scorer( protocols::rosetta_scripts::parse_filter(filter_name, data) );
	std::string const mover_name( tag->getOption< std::string >( "redocker", "null_mover" ) );
	redocker( protocols::rosetta_scripts::parse_mover(mover_name, data) );

	bonus( tag->getOption< core::Real >( "bonus" ) );

	std::string const prefilter_name( tag->getOption< std::string >( "prefilter", "" ) );
	if ( prefilter_name.size() ) {
		prefilter( protocols::rosetta_scripts::parse_filter(prefilter_name, data) );
	} else {
		prefilter( nullptr );
	}

	std::string const postfilter_name( tag->getOption< std::string >( "postfilter", "" ) );
	if ( postfilter_name.size() ) {
		postfilter( protocols::rosetta_scripts::parse_filter(postfilter_name, data) );
	} else {
		postfilter( nullptr );
	}

	utility::vector1< TagCOP > const sub_tags( tag->getTags() );
	BOOST_FOREACH ( TagCOP const subtag, sub_tags ) {
		if ( subtag->getName() == "Before" ) {
			add_before_chemistry( chemistry_from_subtag( subtag, data ) );
		} else if ( subtag->getName() == "After" ) {
			add_after_chemistry( chemistry_from_subtag( subtag, data ) );
		} else {
			utility_exit_with_message( "tag name " + subtag->getName() + " unrecognized." );
		}
	}//foreach subtag
}

} // ns drug_design
} // ns protocols
