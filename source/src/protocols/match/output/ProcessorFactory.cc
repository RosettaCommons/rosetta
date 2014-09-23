// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/ProcessorFactory.hh
/// @brief  Implementation for a factory class responsible for instantiating
///         MatchProcessor classes.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini
/// @author Roland A Pache

// Unit headers
#include <protocols/match/output/ProcessorFactory.hh>

// Package headers
#include <protocols/match/Matcher.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/output/DownstreamRMSEvaluator.hh>
#include <protocols/match/output/MatchConsolidator.hh>
#include <protocols/match/output/MatchOutputter.hh>
#include <protocols/match/output/MatchProcessor.hh>
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchFilter.hh>
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/output/OutputWriter.hh>
#include <protocols/match/output/SameChiBinComboGrouper.hh>
#include <protocols/match/output/SameRotamerComboGrouper.hh>
#include <protocols/match/output/SameSequenceGrouper.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>
#include <protocols/match/output/PDBWriter.hh>
#include <protocols/match/output/UpstreamCollisionFilter.hh>
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.hh>
#include <protocols/match/output/WriteKinemageOutputter.hh>
#include <protocols/match/output/WriteUpstreamCoordinateKineamge.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>
//#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// C++ headers
#include <list>

#include <utility/vector0.hh>


namespace protocols {
namespace match {
namespace output {

static thread_local basic::Tracer TR( "protocols.match.output.ProcessorFactory" );

ProcessorFactory::ProcessorFactory() {}

ProcessorFactory::~ProcessorFactory() {}

MatchProcessorOP
ProcessorFactory::create_processor(
	MatcherCOP matcher,
	MatcherTaskCOP mtask
)
{
	MatchProcessorOP processor;

	UpstreamHitCacherOP cacher( new UpstreamHitCacher( matcher ) );

	if ( mtask->consolidate_matches() ) {
		TR << "Matches will be consolidated before output." << std::endl;
		if( mtask->output_writer_name()  == "CloudPDB" ) TR << "NOTICE: match consolidation and CloudPDB writing are both active. This is fine but somewhat redundant. In this case, the -output_matches_per_group parameter should be set to a higher number than without cloud pdb writing, say 100 or so." << std::endl;
		MatchConsolidatorOP consolidator( new MatchConsolidator );
		consolidator->set_grouper( create_grouper( matcher, mtask, cacher ) );
		consolidator->set_output_writer( create_output_writer( matcher, mtask, cacher ) );
		consolidator->set_n_to_output_per_group( mtask->n_to_output_per_group() );

		processor = consolidator;
	} else {
		MatchOutputterOP outputter( new MatchOutputter );

		outputter->set_output_writer( create_output_writer( matcher, mtask, cacher ) );

		processor = outputter;
	}

	// kylebarlow - Evaluator is now stored in the MatchProcessor parent class to allow match score output
	// Right now this isn't that helpful, but could be more helpful if DownstreamRMSEvaluator is changed to be able to score match_dpos1
	processor->set_evaluator( create_evaluator( matcher, mtask, cacher ) );
	MatchScoreWriterOP match_score_writer = create_match_score_writer();
	if( mtask->output_scores() ) {
		match_score_writer->set_output_filename( mtask->score_output_file_name() );
	}
	processor->set_match_score_writer( match_score_writer );

	std::list< MatchFilterOP > filters( create_filters( matcher, mtask, cacher ) );
	for ( std::list< MatchFilterOP >::const_iterator
					iter = filters.begin(), iter_end = filters.end();
				iter != iter_end; ++iter ) {
		processor->add_filter( *iter );
	}
	return processor;

}

MatchGrouperOP
ProcessorFactory::create_grouper(
	MatcherCOP matcher,
	MatcherTaskCOP mtask,
	UpstreamHitCacherOP cacher
)
{
	if (mtask->grouper_name() == "SameChiBinComboGrouper" ) {
		SameChiBinComboGrouperOP chibin_grouper( new SameChiBinComboGrouper );
		chibin_grouper->set_n_geometric_constraints( matcher->n_geometric_constraints() );
		chibin_grouper->set_hit_cacher( cacher );
		return chibin_grouper;
	} else if ( mtask->grouper_name() == "SameSequenceGrouper" ) {
		SameSequenceGrouperOP seq_grouper( new SameSequenceGrouper );
		seq_grouper->set_n_geometric_constraints( matcher->n_geometric_constraints() );
		seq_grouper->set_hit_cacher( cacher );
		return seq_grouper;
	}else if ( mtask->grouper_name() == "SameSequenceAndDSPositionGrouper" ) {
		SameSequenceAndDSPositionGrouperOP seq_ds_grouper( new SameSequenceAndDSPositionGrouper );
		seq_ds_grouper->set_n_geometric_constraints( matcher->n_geometric_constraints() );
		seq_ds_grouper->set_hit_cacher( cacher );
		for ( Size ii = 1; ii <= matcher->n_geometric_constraints() ; ++ii ) {
			seq_ds_grouper->set_downstream_builder( ii, matcher->downstream_builder( ii ) );
		}
		seq_ds_grouper->set_relevant_atom_ids( mtask->relevant_downstream_atoms()  );
		seq_ds_grouper->set_rms_group_cutoff( mtask->grouper_ds_rmsd() );
		return seq_ds_grouper;
	} else if ( mtask->grouper_name() == "SameRotamerComboGrouper" ) {
		SameRotamerComboGrouperOP rot_grouper( new SameRotamerComboGrouper );
		rot_grouper->set_n_geometric_constraints( matcher->n_geometric_constraints() );
		return rot_grouper;
	} else {
		utility_exit_with_message( "Could not recognize requested MatchGrouper named: " + mtask->grouper_name() );
		return 0;
	}
}



MatchEvaluatorOP
ProcessorFactory::create_evaluator(
	MatcherCOP matcher,
	MatcherTaskCOP mtask,
	UpstreamHitCacherOP /*cacher*/
)
{
	if ( mtask->evaluator_name() == "DownstreamRMSEvaluator" ) {
		DownstreamRMSEvaluatorOP rms_eval( new DownstreamRMSEvaluator );
		rms_eval->set_n_geometric_constraints( matcher->n_geometric_constraints() );
        for ( Size ii = 1; ii <= matcher->n_geometric_constraints(); ++ii ) {
			/// HACK -- all RigidLigandBuilders are equivalent -- FIX THIS!
			rms_eval->set_downstream_builder( ii, matcher->downstream_builder( ii ) );
		}
        rms_eval->set_downstream_pose(matcher->downstream_pose());
		return rms_eval;
	} else {
		utility_exit_with_message( "Could not recognize requested MatchEvaluator named: " + mtask->evaluator_name() );
		return 0;
	}
}

std::list< MatchFilterOP >
ProcessorFactory::create_filters(
	MatcherCOP matcher,
	MatcherTaskCOP mtask,
	UpstreamHitCacherOP cacher
)
{
	if ( ! mtask->filter_names().empty() ) {
		std::cerr << "ERROR: match::output::ProcessorFactory currently lacks logic to instantiate any of the desired filters" << std::endl;
		for ( std::list< std::string >::const_iterator
				filtiter = mtask->filter_names().begin(),
				filtiter_end = mtask->filter_names().end();
				filtiter != filtiter_end; ++filtiter ) {
			std::cerr << "Requested filter '" << *filtiter << "' cannot be instantiated" << std::endl;
		}
		utility_exit_with_message( "Processor factory cannot create requested filter(s)" );

	}

	std::list< MatchFilterOP > filter_list;
	if ( mtask->filter_upstream_residue_collisions() ) {
		output::UpstreamCollisionFilterOP collfilt( new output::UpstreamCollisionFilter( "UpstreamCollisionFilter", cacher ) );
		if ( mtask->filter_upstream_collisions_by_score() ) {
			collfilt->set_filter_by_lj( true );
			collfilt->set_lj_cutoff( mtask->upstream_residue_collision_score_cutoff() );
			collfilt->set_lj_atr_weight( mtask->upstream_residue_collision_Wfa_atr() );
			collfilt->set_lj_rep_weight( mtask->upstream_residue_collision_Wfa_rep() );
			collfilt->set_lj_sol_weight( mtask->upstream_residue_collision_Wfa_sol() );
		} else {
			collfilt->set_tolerated_overlap( mtask->upstream_residue_collision_tolerance() );
		}
		filter_list.push_back( collfilt );
	}

	if ( mtask->filter_upstream_downstream_collisions() || matcher->has_upstream_only_geomcsts() ) {
		output::UpstreamDownstreamCollisionFilterOP collfilt( new output::UpstreamDownstreamCollisionFilter( "UpstreamDownstreamCollisionFilter", cacher ) );
		if ( mtask->filter_upstream_downstream_collisions_by_score() ) {
			collfilt->set_filter_by_lj( true );
			collfilt->set_lj_cutoff( mtask->upstream_downstream_residue_collision_score_cutoff() );
			collfilt->set_lj_atr_weight( mtask->upstream_downstream_residue_collision_Wfa_atr() );
			collfilt->set_lj_rep_weight( mtask->upstream_downstream_residue_collision_Wfa_rep() );
			collfilt->set_lj_sol_weight( mtask->upstream_downstream_residue_collision_Wfa_sol() );
		} else {
			collfilt->set_tolerated_overlap( mtask->upstream_downstream_atom_collision_tolerance() );
		}
		collfilt->set_downstream_pose( * matcher->downstream_pose() );
		collfilt->set_num_geometric_constraints( matcher->n_geometric_constraints() );
		for ( Size ii = 1; ii <= matcher->n_geometric_constraints(); ++ii ) {
			collfilt->set_downstream_builder( ii, matcher->downstream_builder( ii ) );
			if ( mtask->enz_input_data()->mcfi_list( ii )->mcfi(1)->is_covalent() ) {
				collfilt->set_chemical_bond_from_upstream_to_downstream( ii );
			}
		}

		filter_list.push_back( collfilt );
	}

	return filter_list;
}


OutputWriterOP
ProcessorFactory::create_output_writer(
	MatcherCOP matcher,
	MatcherTaskCOP mtask,
	UpstreamHitCacherOP cacher
)
{
	using namespace basic::options;


	if ( mtask->output_writer_name() == "KinWriter" ) {


		output::WriteKinemageOutputterOP kin_match_writer( new output::WriteKinemageOutputter );
		WriteUpstreamHitKinemageOP kin_hit_writer( new WriteUpstreamHitKinemage( mtask->output_file_name() ) );
		kin_match_writer->set_n_geomcst( matcher->n_geometric_constraints() );
		kin_match_writer->set_kin_writer( kin_hit_writer );

		for ( Size ii = 1; ii <= matcher->n_geometric_constraints() ; ++ii ) {
			/// DIRTY ASSUMPTION -- SINGLE RESIDUE IN DOWNSTREAM PARTNER
			if ( matcher->downstream_builder( ii ) ) {
				SingleDownstreamResidueWriterOP downstream_writer( new SingleDownstreamResidueWriter );
				downstream_writer->set_restype( mtask->downstream_pose()->residue(1).type().get_self_ptr() );
				downstream_writer->set_downstream_builder( matcher->downstream_builder( ii ) );
				downstream_writer->set_downstream_master( mtask->downstream_pose()->residue(1).name3() );

				kin_match_writer->set_downstream_writer( ii, downstream_writer );
			}
		}
		kin_match_writer->set_coordinate_cacher( cacher );
		return kin_match_writer;

	} else if ( (mtask->output_writer_name()  == "PDB" ) || (mtask->output_writer_name()  == "pdb" )  || (mtask->output_writer_name()  == "CloudPDB" ) || ( mtask->output_writer_name() == "PoseMatchOutputWriter" ) ) {
		output::PDBWriterOP pdb_writer;
		if ( (mtask->output_writer_name()  == "PDB" ) || (mtask->output_writer_name()  == "pdb" ) ){
			pdb_writer = output::PDBWriterOP( new output::PDBWriter );
		}
		else if( mtask->output_writer_name() == "PoseMatchOutputWriter" ){
			pdb_writer = output::PDBWriterOP( new output::PoseMatchOutputWriter( create_grouper( matcher, mtask, cacher) ) );
		}
		else{
			pdb_writer = output::PDBWriterOP( new output::CloudPDBWriter( create_grouper( matcher, mtask, cacher) ) );
		}
		pdb_writer->set_coordinate_cacher( cacher );
		pdb_writer->initialize_from_matcher_task( mtask );

		if ( option[ OptionKeys::match::ligand_rotamer_index ].user() ) {
			pdb_writer->set_prefix( "UM_LIGROT_" +
				utility::to_string( option[ OptionKeys::match::ligand_rotamer_index ]() ) +
				"_" );
		}

		runtime_assert( matcher->downstream_builder( 1 /* wrong */ ) != 0 );
		for ( Size ii = 1; ii <= matcher->n_geometric_constraints() ; ++ii ) {
			pdb_writer->set_downstream_builder( ii, matcher->downstream_builder( ii ) );
		}
		return pdb_writer;
	}

	else {
		utility_exit_with_message( "Could not recognize requested OutputWriter named: '" + mtask->output_writer_name() + "'" );
		return 0;

	}

}

MatchScoreWriterOP
ProcessorFactory::create_match_score_writer() {
	return MatchScoreWriterOP( new MatchScoreWriter() );
}


}
}
}
