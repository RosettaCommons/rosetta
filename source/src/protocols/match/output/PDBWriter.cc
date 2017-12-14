// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/output/PDBWriter.cc
/// @brief  Forward declaration of class to write output matches that have (presumably) passed the output filters.
/// @author Florian Richter (floric@u.washington.edu), june 09

// Unit headers
#include <protocols/match/output/PDBWriter.hh>

// Package Headers
#include <protocols/match/Hit.hh>
#include <protocols/match/MatcherTask.hh>
#include <protocols/match/output/MatchGrouper.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/output/MatchScoreWriter.hh>

// Project Headers
#include <core/conformation/Residue.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/io/Remarks.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/func/Func.hh>
#include <basic/Tracer.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>

// Utility headers
#include <utility>
#include <utility/string_util.hh>
#include <utility/pointer/ReferenceCount.hh>

//numeric headers
#include <numeric/random/random.hh>

#include <fstream>

#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

static basic::Tracer TR( "protocols.match.output.PDBWriter" );

PDBWriter::PDBWriter() :
	coordinate_cacher_(/* NULL */),
	write_matchres_only_(false),
	scaf_name_(""),
	cstfile_name_(""),
	prefix_("UM"),
	num_written_(0),
	num_geom_cst_(0)
{
	signature_map_.clear();
	upstream_only_geom_cst_.clear();
}

PDBWriter::~PDBWriter() = default;


void
PDBWriter::prepare_for_output_writing()
{
	signature_map_.clear();
	num_written_ = 0;
}


void
PDBWriter::record_match( match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer )
{
	utility::vector1< core::conformation::ResidueCOP > upstream_matchres;
	utility::vector1< core::pose::PoseCOP > downstream_poses;
	utility::vector1< core::Size > ex_geom_ids_for_upstream_res;

	//TR << "NPW recording match w/ hits ";
	//for( core::Size ii =1; ii <= m.size(); ++ii) TR << &(m[ii]) << ", ";
	//TR << std::endl;


	std::map< core::Size, core::Size > redundant_upstream_res;
	determine_redundant_upstream_matchres( match_dspos1( m, 1 ), redundant_upstream_res );

	for ( Size ii = 1; ii <= m.size(); ++ii ) {
		core::conformation::ResidueCOP conf =
			coordinate_cacher_->upstream_conformation_for_hit( ii, m[ ii ] );
		upstream_matchres.push_back( conf );
		ex_geom_ids_for_upstream_res.push_back( m[ ii ].external_geom_id() );
		if ( ! output_dsgeom_for_geomcst_[ ii ] ) continue;
		if ( dsbuilders_[ ii ] ) {
			downstream_poses.push_back( dsbuilders_[ ii ]->downstream_pose_from_hit( m[ii] ) );
		} else {
			utility_exit_with_message( "Cannot output downstream pose for geomcst id " +
				utility::to_string( ii ) + " which does not have a downstream builder!" );
		}
	}

	//utility::vector1< core::Size > matchres_seqpos;
	core::pose::PoseCOP up_outpose = create_output_upstream_pose(
		upstream_matchres, redundant_upstream_res, ex_geom_ids_for_upstream_res );
	std::string outtag = assemble_outtag( upstream_matchres );
	core::Size outcounter(1);

	for ( Size ii = 1; ii <= m.size(); ++ii ) {
		if ( ! output_dsgeom_for_geomcst_[ ii ] ) continue;
		std::string this_tag = outtag + "_" + utility::to_string( ii );
		match_score_writer->add_match( this_tag , evaluator->score(m) );
		this_tag += ".pdb";

		//every match in its own file for now
		std::ofstream file_out( this_tag.c_str() );
		//Size atom_counter(0);

		up_outpose->dump_pdb( file_out );

		downstream_poses[ outcounter ]->dump_pdb( file_out );
		file_out.close();
		outcounter++;
	}

	num_written_++;

	//if( num_written_ % 100 == 0 ) std::cout << "recored match function called " << num_written_ << "times" <<std::endl;
}

void
PDBWriter::record_match( match_dspos1 const & m )
{

	//TR << "NPW recording dspos1 match w/ upstream hits ";
	//for( core::Size ii =1; ii <= m.upstream_hits.size(); ++ii) TR << &(m.upstream_hits[ii]) << ", ";
	//TR << std::endl;


	/// What are we doing here?  The Matcher should never have sent us this match.
	if ( ! output_dsgeom_for_geomcst_[ m.originating_geom_cst_for_dspos ] ) return;

	if ( ! dsbuilders_[ m.originating_geom_cst_for_dspos ] ) {
		utility_exit_with_message( "Cannot output downstream pose for geomcst id " +
			utility::to_string( m.originating_geom_cst_for_dspos ) +
			" which does not have a downstream builder!" );
	}

	utility::vector1< core::conformation::ResidueCOP > upstream_matchres;
	utility::vector1< core::Size > ex_geom_ids_for_upstream_res;

	Hit fullhit = full_hit( m );

	core::pose::PoseCOP downstream_pose = dsbuilders_[ m.originating_geom_cst_for_dspos ]->
		downstream_pose_from_hit( fullhit );

	/// Does any single residue do double-duty? e.g. is there a residue whose
	/// sidechain is used in satisfying one geometric constraint, while it's
	/// backbone is used in another?
	std::map< core::Size, core::Size > redundant_upstream_res;
	determine_redundant_upstream_matchres( m, redundant_upstream_res );

	for ( Size ii = 1; ii <= m.upstream_hits.size(); ++ii ) {
		core::conformation::ResidueCOP conf =
			coordinate_cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] ) );
		upstream_matchres.push_back( conf );
		ex_geom_ids_for_upstream_res.push_back( fullhit.external_geom_id() );
	}

	//utility::vector1< core::Size > matchres_seqpos;

	core::pose::PoseCOP up_outpose = create_output_upstream_pose(
		upstream_matchres, redundant_upstream_res, ex_geom_ids_for_upstream_res );

	std::string outtag = assemble_outtag( upstream_matchres );

	std::string this_tag = outtag + "_" + utility::to_string( m.originating_geom_cst_for_dspos );
	this_tag += ".pdb";
	std::ofstream file_out( this_tag.c_str() );
	up_outpose->dump_pdb( file_out );
	downstream_pose->dump_pdb( file_out );

	num_written_++;

}

//void
//PDBWriter::create_output_pose_from_dspos1_match(
// match_dspos1 const & m )
//{
//}

void
PDBWriter::set_coordinate_cacher( UpstreamHitCacherOP cacher )
{
	coordinate_cacher_ = cacher;
}

void
PDBWriter::set_prefix( std::string const & prefix )
{
	prefix_ = prefix;
}

void
PDBWriter::initialize_from_matcher_task(
	MatcherTaskCOP mtask
){
	Parent::initialize_from_matcher_task( mtask );

	orig_upstream_pose_ = mtask->upstream_pose();
	downstream_pose_from_task_ = mtask->downstream_pose();
	write_matchres_only_ = mtask->output_matchres_only();
	num_geom_cst_ = mtask->enz_input_data()->mcfi_lists_size();
	output_dsgeom_for_geomcst_.resize( num_geom_cst_ );
	std::fill( output_dsgeom_for_geomcst_.begin(), output_dsgeom_for_geomcst_.end(), false );
	dsbuilders_.resize( num_geom_cst_ );
	for ( core::Size i = 1; i <= mtask->geom_csts_downstream_output().size(); ++i ) {
		output_dsgeom_for_geomcst_[ mtask->geom_csts_downstream_output()[i] ] = true;
	}
	scaf_name_ = mtask->upstream_pose_name();
	utility::vector1< std::string > cst_path = utility::string_split( mtask->cstfile_name(), '/' );
	cstfile_name_ = (utility::string_split( cst_path[ cst_path.size() ] , '.' ) )[1];
	upstream_only_geom_cst_ = mtask->upstream_only_geom_cst();
}

void
PDBWriter::set_downstream_builder(
	Size geomcst_id,
	downstream::DownstreamBuilderCOP dsbuilder
){
	dsbuilders_[ geomcst_id ] = dsbuilder;
}


void
PDBWriter::assemble_remark_lines(
	core::pose::Pose & outpose,
	utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres,
	std::map< core::Size, core::Size > const & redundant_upstream_res,
	utility::vector1< core::Size > const & ex_geom_ids_for_upstream_res
) const
{

	using namespace core::io::pdb;

	core::io::Remarks rems;

	std::string ds_resname = downstream_pose_from_task_->residue(1).name3();

	for ( core::Size i = 1; i <= upstream_matchres.size(); ++i ) {

		core::io::RemarkInfo ri;
		ri.num = 666; /// really now?
		std::string upname3( upstream_matchres[ i ]->name3() );
		auto red_it = redundant_upstream_res.find( i );
		if ( red_it != redundant_upstream_res.end() ) {
			upname3 = upstream_matchres[ red_it->second ]->name3();
		}
		std::string upres_chain( utility::to_string( orig_upstream_pose_->pdb_info()->chain( upstream_matchres[ i ]->seqpos() ) ) );

		std::string targ_resname( ds_resname );
		std::string targ_chain("X");
		int targ_seqpos(0);
		auto upstream_only_it(upstream_only_geom_cst_.find( i ) );
		if ( upstream_only_it != upstream_only_geom_cst_.end() ) {
			targ_resname = upstream_matchres[ upstream_only_it->second ]->name3();
			targ_chain =  utility::to_string( orig_upstream_pose_->pdb_info()->chain( upstream_matchres[ upstream_only_it->second ]->seqpos() ) );
			targ_seqpos = orig_upstream_pose_->pdb_info()->number( upstream_matchres[ upstream_only_it->second ]->seqpos() );
		}
		ri.value =  toolbox::match_enzdes_util::assemble_remark_line(
			targ_chain, targ_resname, targ_seqpos, upres_chain, upname3, orig_upstream_pose_->pdb_info()->number( upstream_matchres[ i ]->seqpos() ), i , ex_geom_ids_for_upstream_res[ i ] );
		rems.push_back( ri );

	}
	for ( auto const & r : orig_upstream_pose_->pdb_info()->remarks() ) {
		if ( r.num != 666 ) {
			rems.push_back( r );
		}
	}
	outpose.pdb_info()->remarks( rems );

}


core::pose::PoseCOP
PDBWriter::create_output_upstream_pose(
	utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres,
	std::map< core::Size, core::Size > const & redundant_upstream_res,
	utility::vector1< core::Size > const & ex_geom_ids_for_upstream_res
)
{
	runtime_assert( ex_geom_ids_for_upstream_res.size() == upstream_matchres.size() );
	//matchres_seqpos.clear();

	core::pose::PoseOP outpose;

	if ( write_matchres_only_ ) {
		outpose = core::pose::PoseOP( new core::pose::Pose() );

		core::pose::PDBInfoOP pdbinf( new core::pose::PDBInfo( *outpose ) );
		outpose->pdb_info( pdbinf );
	} else outpose = core::pose::PoseOP( new core::pose::Pose( *orig_upstream_pose_ ) );

	for ( core::Size i = 1, count_non_redundant = 1; i <= upstream_matchres.size(); ++i ) {
		if ( redundant_upstream_res.find( i ) != redundant_upstream_res.end() ) continue;
		//matchres_seqpos.push_back( upstream_matchres[i]->seqpos() );
		if ( write_matchres_only_ ) {
			if ( i == 1 ) {
				outpose->append_residue_by_jump( *(upstream_matchres[i]), 1 );
			} else {
				outpose->append_residue_by_bond( *(upstream_matchres[i]) );
			}
			outpose->pdb_info()->chain( count_non_redundant, orig_upstream_pose_->pdb_info()->chain( upstream_matchres[i]->seqpos() ) );
			outpose->pdb_info()->number( count_non_redundant, orig_upstream_pose_->pdb_info()->number( upstream_matchres[i]->seqpos() ) );
		} else {
			outpose->replace_residue( upstream_matchres[i]->seqpos(), *(upstream_matchres[i]), false );
		}
		++count_non_redundant; // only increment if we did not encounter a redundant residue on this iteration
	}

	assemble_remark_lines( *outpose, upstream_matchres, redundant_upstream_res, ex_geom_ids_for_upstream_res );

	//looks like we have to set the pdbinfo obsolete to false
	outpose->pdb_info()->obsolete( false );

	return outpose;
}

std::string
PDBWriter::signature_string(
	utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres
) const
{
	std::string signature("");
	for ( core::Size i =1; i <= upstream_matchres.size(); ++i ) {
		signature = signature + upstream_matchres[i]->name1() +
			utility::to_string( orig_upstream_pose_->pdb_info()->number( upstream_matchres[i]->seqpos() ) );
	}
	return signature;
}

core::Size
PDBWriter::num_geom_cst() const
{
	return num_geom_cst_;
}

std::string
PDBWriter::scaf_name() const
{
	return scaf_name_;
}

std::string
PDBWriter::cstfile_name() const
{
	return cstfile_name_;
}

std::string
PDBWriter::prefix() const
{
	return prefix_;
}


UpstreamHitCacherOP
PDBWriter::coordinate_cacher()
{
	return coordinate_cacher_;
}

utility::vector1< downstream::DownstreamBuilderCOP > const &
PDBWriter::dsbuilders()
{
	return dsbuilders_;
}


std::string
PDBWriter::assemble_outtag(
	utility::vector1< core::conformation::ResidueCOP > const & upstream_matchres
)
{

	std::string signature( signature_string( upstream_matchres ));
	//for ( core::Size i =1; i <= upstream_matchres.size(); ++i ){
	// signature = signature + upstream_matchres[i]->name1() +
	//  utility::to_string( upstream_matchres[i]->seqpos() );
	//}

	auto map_it = signature_map_.find( signature );
	if ( map_it == signature_map_.end() ) {
		core::Size num_this_um = signature_map_.size() + 1;
		signature_map_[ signature ] = SizePair( num_this_um, 1 );
		map_it = signature_map_.find( signature );
	} else {
		map_it->second.second++;
	}

	std::string unique_string( utility::to_string( map_it->second.first ) + "_" +
		signature + "_" + utility::to_string( map_it->second.second ) );

	return prefix_ + "_" + unique_string + "_" + scaf_name_ + "_" + cstfile_name_;
}


CloudPDBWriter::CloudPDBWriter( MatchGrouperOP grouper ) :
	PDBWriter(),
	grouper_(std::move(grouper))
{
	clear_match_data();
}

CloudPDBWriter::~CloudPDBWriter() = default;

void
CloudPDBWriter::prepare_for_output_writing()
{
	parent::prepare_for_output_writing();
	clear_match_data();
	grouper_->reset();
}


/// @details this class writes output in a form that has no link
/// between a downstream conformation and the upstream conf that
/// it came from.
void
CloudPDBWriter::end_output_writing()
{
	TR << "CloudPDBWriter beginning to output matches for " << match_groups_ushits_.size() << " match groups." << std::endl;
	core::Size total_ushit_count(0), total_dshit_count(0);
	for ( core::Size ii = 1; ii <= match_groups_ushits_.size(); ++ii ) {
		for ( core::Size jj =1; jj <= num_geom_cst(); ++jj ) {
			total_ushit_count += match_groups_ushits_[ii][jj].size();
			total_dshit_count += match_groups_dshits_[ii][jj].size();
		}
	}
	TR << "A total of " << total_ushit_count << " upstream hits and " << total_dshit_count << " downstream hits were recorded." << std::endl;
	write_match_groups();

	//let's save some memory
	clear_match_data();
}

void
CloudPDBWriter::record_match(  match const & m , MatchEvaluatorOP evaluator, MatchScoreWriterOP match_score_writer )
{
	runtime_assert( m.size() == num_geom_cst() );
	core::Size mgroup( grouper_->assign_group_for_match( m ) );
	//TR << "CPW recording match for assigned group " << mgroup << " ... " << std::endl;

	// check if this is a new unqiue match
	if ( mgroup > match_groups_ushits_.size() ) {
		runtime_assert( mgroup == match_groups_ushits_.size() + 1 );
		match_groups_ushits_.push_back( UpstreamHitSets( num_geom_cst() ) );
		match_groups_dshits_.push_back( DownstreamHitSets( num_geom_cst() ) );
		representative_group_matches_.push_back( match_dspos1( m, 1 ) );
		//TR << "CPW pushed back first match for group " << mgroup << " match_group_ushits_.size is " << match_groups_ushits_.size() << ", ex geom id for hit2 is " << m[2].external_geom_id() << std::endl;

		utility::vector1< core::conformation::ResidueCOP > upstream_matchres;
		for ( Size ii = 1; ii <= m.size(); ++ii ) {
			core::conformation::ResidueCOP conf = coordinate_cacher()->upstream_conformation_for_hit( ii, m[ ii ] );
			upstream_matchres.push_back( conf );

			//note the special case here: we don't push back the first because
			//that will be written out separately
			if ( (ii != 1) && (dsbuilders()[ii]) ) {
				match_groups_dshits_[mgroup][ii].insert( downstream_hit( m[ii] ) );
			}
		}
		std::string unique_match_name = prefix() + "_" + utility::to_string( mgroup ) + "_" + signature_string( upstream_matchres ) + "_" + scaf_name() + "_" + cstfile_name();
		unique_match_names_.push_back( unique_match_name );
		match_score_writer->add_match( unique_match_name , evaluator->score(m) );
	} else {
		//TR << "CPW group " << mgroup << " has been previously encountered.. ";
		for ( core::Size ii = 1; ii <= num_geom_cst(); ++ii ) {
			upstream_hit uhit( m[ii] );
			if ( match_groups_ushits_[mgroup][ii].find( uhit ) == match_groups_ushits_[mgroup][ii].end() ) {
				match_groups_ushits_[mgroup][ii].insert( uhit );
				//TR << " got new uphit for geomcst " << ii << ";   ";
			}
			if ( dsbuilders()[ii] ) {
				downstream_hit dhit( m[ii] );
				if ( match_groups_dshits_[mgroup][ii].find( dhit ) == match_groups_dshits_[mgroup][ii].end() ) {
					match_groups_dshits_[mgroup][ii].insert( dhit );
					//TR << " got new downhit for geomcst " << ii << ";   ";
				}
			}
		}
		//TR << std::endl;
	}
}


void
CloudPDBWriter::record_match( match_dspos1 const & /*m*/ )
{
	utility_exit_with_message( "not implemented yet");
}

void
CloudPDBWriter::write_match_groups()
{
	runtime_assert( match_groups_ushits_.size() == unique_match_names_.size() );
	for ( core::Size ii = 1; ii <= match_groups_ushits_.size(); ++ii ) {

		TR << "beginning writing cloud for group " << ii << std::endl;
		for ( core::Size jj =1; jj <= num_geom_cst(); ++jj ) {
			TR << match_groups_ushits_[ii][jj].size() << " upstream hits and "  << match_groups_dshits_[ii][jj].size() << " downstream hits " << " for geom cst " << jj << std::endl;
		}
		TR.flush();

		//1. for each group, we have to write out the whole pose once
		std::map< core::Size, core::Size > redundant_upstream_res;
		utility::vector1< core::conformation::ResidueCOP > upstream_matchres;
		utility::vector1< core::Size > ex_geom_ids_for_upstream_res;
		match_dspos1 const & rep_match( representative_group_matches_[ii] );

		determine_redundant_upstream_matchres( rep_match, redundant_upstream_res );
		setup_hitset_iterators_for_group( ii );

		for ( Size jj = 1; jj <= rep_match.upstream_hits.size(); ++jj ) {
			core::conformation::ResidueCOP conf =
				coordinate_cacher()->upstream_conformation_for_hit( jj, fake_hit( rep_match.upstream_hits[ jj ] ) );
			upstream_matchres.push_back( conf );
			ex_geom_ids_for_upstream_res.push_back( rep_match.upstream_hits[jj].external_geom_id() );

			if ( redundant_upstream_res.find( jj ) != redundant_upstream_res.end() ) {
				us_hitset_its_[jj] = us_hitset_end_its_[jj];
			}
		}
		core::pose::PoseCOP up_outpose = create_output_upstream_pose(
			upstream_matchres, redundant_upstream_res, ex_geom_ids_for_upstream_res );
		core::pose::PoseCOP downstream_pose = dsbuilders()[ rep_match.originating_geom_cst_for_dspos ]->downstream_pose_from_hit(  full_hit( rep_match ) );

		std::ofstream file_out( (unique_match_names_[ii] + "_1.pdb").c_str() );
		if ( !file_out.is_open() ) utility_exit_with_message("Could not open file with name " + unique_match_names_[ii] + "_1.pdb for outputting matches.");

		file_out << "MODEL    1\n";
		up_outpose->dump_pdb( file_out );
		downstream_pose->dump_pdb( file_out );
		file_out << "ENDMDL \n";

		//now spit out the remaining hits in this group
		core::Size allmodelcount(0), num_files_this_group(1), modelcount(2), atomcounter(0);
		bool all_hitset_iterators_at_end( false );
		//UpstreamHitSets const & us_hitset_ref( match_groups_ushits_[ ii ] );

		core::io::StructFileRepOptionsOP pdb_writer_options( new core::io::StructFileRepOptions );
		pdb_writer_options->set_skip_connect_info( true );

		while ( ! all_hitset_iterators_at_end ) {
			all_hitset_iterators_at_end = true;
			bool modeltag_written( false );

			for ( core::Size jj = 1; jj <= num_geom_cst(); ++jj ) {
				if ( us_hitset_its_[jj] == us_hitset_end_its_[jj] ) continue;
				all_hitset_iterators_at_end = false;

				if ( !modeltag_written ) {
					file_out << "MODEL    "+utility::to_string( modelcount++ )+"\n";
					modeltag_written = true;
				}
				core::conformation::ResidueCOP conf = coordinate_cacher()->upstream_conformation_for_hit( jj, fake_hit(*us_hitset_its_[jj] ) );
				core::io::pdb::dump_pdb_residue( *conf, atomcounter, file_out, pdb_writer_options );
				//TR <<"nm upstream res for cstid " << jj << "; ";
				us_hitset_its_[jj]++;
			} //jj loop over all geom csts

			//write out one ds pose
			for ( core::Size jj = 1; jj <= num_geom_cst(); ++jj ) {
				if ( ds_hitset_its_[jj] == ds_hitset_end_its_[jj] ) continue;
				all_hitset_iterators_at_end = false;
				if ( !modeltag_written ) {
					file_out << "MODEL    "+utility::to_string( modelcount++ )+"\n";
					modeltag_written = true;
				}
				dsbuilders()[jj]->downstream_pose_from_hit( fake_hit(*(ds_hitset_its_[jj])) )->dump_pdb( file_out );
				//TR <<"nm downstream res for cstid " << jj << "; ";
				ds_hitset_its_[jj]++;
				break;
			}
			if ( modeltag_written ) file_out << "ENDMDL \n";

			//note: we don't want to go overboard w/ the number of different models,
			//bc the packer can only handle so much. so if there are more than 100,
			//let's open a new file
			if ( modelcount == 101 && (! all_hitset_iterators_at_end ) ) {
				file_out.close();
				atomcounter = 0;
				allmodelcount += (modelcount - 2);
				modelcount = 2;
				num_files_this_group++;
				file_out.open( (unique_match_names_[ii] + "_"+utility::to_string( num_files_this_group ) +".pdb").c_str() );
				file_out << "MODEL    1\n";
				up_outpose->dump_pdb( file_out );
				downstream_pose->dump_pdb( file_out );
				file_out << "ENDMDL \n";
			}

		} //while there are still hits
		allmodelcount += modelcount;
		TR << "A total of " << allmodelcount - 1 << " models were written for " << " match group " << ii << " in " << num_files_this_group << " files." << std::endl;
		file_out.close();

		//at the end, let's clear the stuff of this group to relieve some memory
		match_groups_ushits_[ ii ].clear();
		match_groups_dshits_[ ii ].clear();
	} //loop over all match groups
} //write_match_groups


utility::vector1< CloudPDBWriter::UpstreamHitSets > const &
CloudPDBWriter::match_groups_ushits() const
{
	return match_groups_ushits_;
}

utility::vector1< CloudPDBWriter::DownstreamHitSets > const &
CloudPDBWriter::match_groups_dshits() const
{
	return match_groups_dshits_;
}

utility::vector1< match_dspos1 > const &
CloudPDBWriter::representative_group_matches() const
{
	return representative_group_matches_;
}

utility::vector1< std::set< downstream_hit >::const_iterator > const &
CloudPDBWriter::ds_hitset_its() const
{
	return ds_hitset_its_;
}

utility::vector1< std::set< downstream_hit >::const_iterator > const &
CloudPDBWriter::ds_hitset_end_its() const
{
	return ds_hitset_end_its_;
}

void
CloudPDBWriter::clear_match_data()
{
	match_groups_ushits_.clear();
	match_groups_dshits_.clear();
	unique_match_names_.clear();
	representative_group_matches_.clear();
	us_hitset_its_.clear();
	us_hitset_end_its_.clear();
	ds_hitset_its_.clear();
	ds_hitset_end_its_.clear();
}

void
CloudPDBWriter::setup_hitset_iterators_for_group(
	core::Size const group )
{
	us_hitset_its_.clear();
	us_hitset_end_its_.clear();
	ds_hitset_its_.clear();
	ds_hitset_end_its_.clear();

	UpstreamHitSets const & us_hitset_ref( match_groups_ushits_[ group ] );
	DownstreamHitSets const & ds_hitset_ref( match_groups_dshits_[ group ] );
	runtime_assert( us_hitset_ref.size() == ds_hitset_ref.size() );
	for ( core::Size ii = 1; ii <= us_hitset_ref.size(); ++ii ) {
		//runtime_assert( us_hitset_ref[ii].size() > 0 );
		us_hitset_its_.push_back( us_hitset_ref[ ii ].begin() );
		us_hitset_end_its_.push_back( us_hitset_ref[ii].end() );
		ds_hitset_its_.push_back( ds_hitset_ref[ ii ].begin() );
		ds_hitset_end_its_.push_back( ds_hitset_ref[ii].end() );
	}
}


PoseMatchOutputWriter::PoseMatchOutputWriter( MatchGrouperOP grouper )
: parent( grouper )
{}

PoseMatchOutputWriter::~PoseMatchOutputWriter()= default;

/// @brief this function doesn't do anything, but needs
/// to exist to suppress the unwanted call to the base
/// class function
void
PoseMatchOutputWriter::end_output_writing(){}

void
PoseMatchOutputWriter::insert_match_into_pose(
	core::pose::Pose & pose,
	core::Size match_group
)
{

	std::map< core::Size, core::Size > redundant_upstream_res;
	utility::vector1< core::conformation::ResidueCOP > upstream_matchres;
	utility::vector1< core::Size > ex_geom_ids_for_upstream_res;
	match_dspos1 const & rep_match( representative_group_matches()[match_group] );

	determine_redundant_upstream_matchres( rep_match, redundant_upstream_res );
	//setup_hitset_iterators_for_group( match_group );

	//first put in the upstream residues
	for ( Size jj = 1; jj <= rep_match.upstream_hits.size(); ++jj ) {
		core::conformation::ResidueCOP conf =
			coordinate_cacher()->upstream_conformation_for_hit( jj, fake_hit( rep_match.upstream_hits[ jj ] ) );
		upstream_matchres.push_back( conf );
		ex_geom_ids_for_upstream_res.push_back( full_hit( rep_match).external_geom_id() );
		if ( redundant_upstream_res.find( jj ) == redundant_upstream_res.end() ) {
			pose.replace_residue( conf->seqpos(), *conf, false );
		}
	}

	//assemble remarks
	assemble_remark_lines( pose, upstream_matchres, redundant_upstream_res, ex_geom_ids_for_upstream_res );

	//then the downstream pose
	core::pose::PoseCOP downstream_pose = dsbuilders()[ rep_match.originating_geom_cst_for_dspos ]->downstream_pose_from_hit(  full_hit( rep_match ) );

	if ( downstream_pose->size() == 1 ) {
		pose.append_residue_by_jump( downstream_pose->residue(1), 1, "", "", true );
		pose.pdb_info()->chain( pose.size(), 'X' );
		pose.pdb_info()->number( pose.size(), 1 );
		pose.pdb_info()->obsolete( false );
		//pose.append_residue_by_jump( downstream_pose->residue(1), pose.size(), "CA", downstream_pose->residue(1).atom_name(1), true );

		utility::vector1< core::conformation::ResidueCOP > additional_lig_confs;
		for ( core::Size i = 1; i <= num_geom_cst(); ++i ) {
			for ( auto ds_hits_it( match_groups_dshits()[match_group][i].begin()), ds_hits_end( match_groups_dshits()[match_group][i].end() ); ds_hits_it != ds_hits_end; ++ds_hits_it ) {
				additional_lig_confs.push_back( dsbuilders()[i]->downstream_pose_from_hit( fake_hit(*(ds_hits_it)) )->residue(1).get_self_ptr() );
				if ( additional_lig_confs.size() > 99 ) break;
			}
			if ( additional_lig_confs.size() > 99 ) break;
		}
		if ( additional_lig_confs.size() > 1 ) {
			protocols::toolbox::match_enzdes_util::get_enzdes_observer( pose )->set_rigid_body_confs_for_lig( pose.size(), additional_lig_confs );
		}
	} else utility_exit_with_message("PoseMatchOutputWriter not set up to put a downstream pose containing more than one ligand into the upstream pose");
}

void
PoseMatchOutputWriter::insert_match_into_pose(
	core::pose::Pose & pose
)
{
	core::Size num_match_groups( match_groups_ushits().size() );
	if ( num_match_groups > 1 ) {
		core::Size mgroup( numeric::random::random_range( 1, num_match_groups ) );
		TR.Warning << "Matcher produced " << num_match_groups << " unique match groups, randomly picked " << mgroup << " to be inserted into the pose." << std::endl;
		insert_match_into_pose( pose, mgroup );
	} else if ( num_match_groups == 1 ) insert_match_into_pose( pose, 1 );
	else {
		TR << "Apparently no matches were found, so it's impossible to insert a match into the pose." << std::endl;
	}
}


}
}
}
