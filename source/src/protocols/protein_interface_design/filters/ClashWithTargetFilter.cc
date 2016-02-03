// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/ClashWithTargetFilter.cc
/// @brief ClashWithTarget filtering
/// @author Lei Shi (shileiustc@gmail.com)
#include <protocols/protein_interface_design/filters/ClashWithTargetFilter.hh>
#include <protocols/protein_interface_design/filters/ClashWithTargetFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/protein_interface_design/movers/BuildAlaPose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/simple_filters/DdgFilter.hh>
//#include <protocols/motif_grafting/movers/MotifGraftMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <algorithm>
#include <list>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

using namespace ObjexxFCL;
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;


ClashWithTargetFilter::ClashWithTargetFilter() :
	protocols::filters::Filter( "ClashWithTarget" ) { }

ClashWithTargetFilter::~ClashWithTargetFilter() {}

protocols::filters::FilterOP
ClashWithTargetFilter::clone() const {
	return protocols::filters::FilterOP( new ClashWithTargetFilter( *this ) );
}

static THREAD_LOCAL basic::Tracer TR( "protocols.protein_interface_design.filters.ClashWithTargetFilter" );

core::Size ClashWithTargetFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core;

	//superposition
	core::pose::Pose pose_copy=pose;

	//pose_copy.dump_pdb("input.pdb");
	if ( !align_to_pdbname_.empty() ) {
		core::pose::PoseOP ref_pose= core::import_pose::pose_from_file( align_to_pdbname_, false , core::import_pose::PDB_file);
		protocols::simple_moves::SuperimposeMoverOP SuperimposeMoverOP( new protocols::simple_moves::SuperimposeMover(*ref_pose, ref_start_, ref_end_, pose_start_, pose_end_, true) );
		SuperimposeMoverOP->apply(pose_copy);
	}
	//pose_copy.dump_pdb("inputAlign.pdb");

	//Append target pose
	core::pose::PoseOP context_pose= core::import_pose::pose_from_file( context_pdbname_, false , core::import_pose::PDB_file);
	Size chain1end(pose_copy.conformation().chain_end( 1 ) );
	pose_copy.append_pose_by_jump(*context_pose,chain1end);
	//pose_copy.dump_pdb("target_inputAlign.pdb");

	//mutate interface Ala
	protocols::protein_interface_design::movers::BuildAlaPose toAla(1,2,20,clash_residues_);
	toAla.apply(pose_copy);
	//pose_copy.dump_pdb("target_Ala_inputAlign.pdb");

	//Clash check
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function("empty") );
	scorefxn->set_weight( core::scoring::fa_atr, 0.8);
	scorefxn->set_weight( core::scoring::fa_rep, 0.44);
	protocols::simple_filters::DdgFilter ddg = protocols::simple_filters::DdgFilter( clash_score_cutoff_, scorefxn, 1, 1);

	core::Real score = ddg.compute( pose_copy );
	return score;

	//typedef utility::pointer::owning_ptr< MotifGraftMover >  MotifGraftMoverOP;
	//bool b_full_motif_bb_alignment=true;
	//MotifGraftMoverOP motif_grafting_moverOP = new protocols::motif_grafting::movers::MotifGraftMover();
	//motif_grafting_moverOP->parse_my_string_arguments_and_cast_to_globalPrivateSpaceVariables(context_pdbname_,align_to_pdbname_,align_RMSD_tolerance_,align_RMSD_tolerance_,clash_score_cutoff_,utility::to_string("0:0"),utility::to_string("0:0"),clash_residues_,hotspots_,b_full_motif_bb_alignment,false,false,false,false,false,false);

	//core::pose::Pose pose_copy=pose;
	//core::Size grafting=1;
	//try{
	// motif_grafting_moverOP->apply(pose_copy);
	//  pose_copy.dump_pdb("temp.pdb");
	//} catch ( utility::excn::EXCN_Base& excn ) {
	// grafting=0;
	//}
	//return grafting;
}

bool
ClashWithTargetFilter::apply( core::pose::Pose const & pose ) const {

	core::Size const ddG( compute( pose ));
	if ( ddG < clash_score_cutoff_ ) {
		TR<<"Passing: ddG " << ddG << " < clash_score_cutoff_: " << clash_score_cutoff_  << std::endl;
		return( true );
	} else {
		TR<<"Failure: ddG " << ddG << " >= clash_score_cutoff_: " << clash_score_cutoff_  << std::endl;
		return( false );
	}
}

void
ClashWithTargetFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const ddG( compute( pose ));
	if ( ddG < clash_score_cutoff_ ) {
		out<<"Passing: ddG " << ddG << " < clash_score_cutoff_: " << clash_score_cutoff_  << std::endl;
	} else {
		out<<"Failure: ddG " << ddG << " >= clash_score_cutoff_: " << clash_score_cutoff_  << std::endl;
	}
}

core::Real
ClashWithTargetFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const ddG( compute( pose ));
	return( ddG);
}

void
ClashWithTargetFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & , protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	if ( tag->hasOption("align_to_pdbname") ) {
		align_to_pdbname_=tag->getOption<std::string>("align_to_pdbname");
	}

	if ( tag->hasOption("context_pdbname") ) {
		context_pdbname_=tag->getOption<std::string>("context_pdbname");
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption("Must specify the context_pdbname, required");
	}

	clash_score_cutoff_= tag->getOption<core::Size>( "clash_score_cutoff", 100 );
	clash_residues_= tag->getOption<std::string>( "clash_residues", "ALA" );

	if ( tag->hasOption("ref_start") && tag->hasOption("ref_end") && tag->hasOption("pose_start") && tag->hasOption("pose_end")  ) {
		ref_start_ = tag->getOption< Size >("ref_start");
		ref_end_ = tag->getOption< Size >("ref_end");
		pose_start_ = tag->getOption< Size >("pose_start");
		pose_end_ = tag->getOption< Size >("pose_end");
		runtime_assert(ref_start_ > 0 && ref_start_ < ref_end_);
		runtime_assert(pose_start_ > 0 && pose_start_ < pose_end_);
	} else if ( !tag->hasOption("ref_start") && !tag->hasOption("ref_end") && !tag->hasOption("pose_start") && !tag->hasOption("pose_end") ) {
		TR<<"No residues for alignment is specified. Use TMalign"  << std::endl;
	} else  {
		throw utility::excn::EXCN_RosettaScriptsOption("Must specify all ref_start/ref_end (start/end residue number in the align_to_pdbname for alignment) and pose_start/pose_end, start/end residue number in the current pose for alignment");
	}

}

protocols::filters::FilterOP
ClashWithTargetFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ClashWithTargetFilter() ); }

std::string
ClashWithTargetFilterCreator::keyname() const { return "ClashWithTarget"; }

} // filters
} // protein_interface_design
} // devel


