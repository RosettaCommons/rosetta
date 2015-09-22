// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/RestrictRegion.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <devel/denovo_design/RestrictRegion.hh>
#include <devel/denovo_design/RestrictRegionCreator.hh>

// Project Headers
#include <devel/denovo_design/task_operations/HighestEnergyRegion.hh>
#include <devel/denovo_design/task_operations/DesignBySecondaryStructure.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/picking_old/vall/util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

//C++ Headers


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif

static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.RestrictRegion" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {

std::string
RestrictRegionCreator::keyname() const
{
	return RestrictRegionCreator::mover_name();
}

protocols::moves::MoverOP
RestrictRegionCreator::create_mover() const {
	return protocols::moves::MoverOP( new RestrictRegion() );
}

std::string
RestrictRegionCreator::mover_name()
{
	return "RestrictRegion";
}

///  ---------------------------------------------------------------------------------
///  RestrictRegion main code:
///  ---------------------------------------------------------------------------------

/// @brief initialize static member variable
utility::vector1< core::Size > RestrictRegion::last_residues_restricted_;
utility::vector1< std::string > RestrictRegion::permanently_restricted_residues_;
core::pose::PoseOP RestrictRegion::previous_pose_ = NULL;

/// @brief default constructor
RestrictRegion::RestrictRegion() :
	Mover( "RestrictRegion" ),
	type_( "score" ),
	permanent_restriction_( false ),
	resfile_( "" ),
	blueprint_file_( "" ),
	psipred_cmd_( "" ),
	enable_max_trp_( false ),
	max_trp_( 1 ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	regions_to_mutate_( 1 ),
	last_type_( "" )
{
	metric_stats_.clear();
	restricted_residues_.clear();
}

/// @brief copy constructor
RestrictRegion::RestrictRegion( RestrictRegion const & rval ) :
	Mover( rval ),
	type_( rval.type_ ),
	permanent_restriction_( rval.permanent_restriction_ ),
	resfile_( rval.resfile_ ),
	blueprint_file_( rval.blueprint_file_ ),
	psipred_cmd_( rval.psipred_cmd_ ),
	enable_max_trp_( rval.enable_max_trp_ ),
	max_trp_( rval.max_trp_ ),
	task_factory_( rval.task_factory_->clone() ),
	restricted_residues_( rval.restricted_residues_ ),
	scorefxn_( rval.scorefxn_->clone() ),
	regions_to_mutate_( rval.regions_to_mutate_ ),
	last_type_( rval.last_type_ ),
	metric_stats_( rval.metric_stats_ ),
	highestEnergyRegionOperation_ops_ (rval.highestEnergyRegionOperation_ops_)
{}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
RestrictRegion::~RestrictRegion()
{}


/// Return a copy of ourselves
protocols::moves::MoverOP
RestrictRegion::clone() const {
	return protocols::moves::MoverOP( new RestrictRegion(*this) );
}

void
RestrictRegion::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	type_ = tag->getOption< std::string >( "type", type_ );
	initialize_resfile( tag->getOption< std::string >( "resfile", "" ) );
	psipred_cmd_ = tag->getOption< std::string >( "psipred_cmd", psipred_cmd_ );
	blueprint_file_ = tag->getOption< std::string >( "blueprint", blueprint_file_ );
	task_factory_ = protocols::rosetta_scripts::parse_task_operations( tag, data );
	regions_to_mutate_ = tag->getOption< core::Size >( "num_regions", regions_to_mutate_ );
	permanent_restriction_ = tag->getOption< core::Size >( "permanent_restriction", permanent_restriction_ );
	if ( tag->hasOption( "max_trp" ) ) {
		enable_max_trp_ = true;
		max_trp_ = tag->getOption< core::Size >( "max_trp", max_trp_ );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		scorefxn_ = data.get_ptr< core::scoring::ScoreFunction >( "scorefxns",
			tag->getOption< std::string >( "scorefxn" ) )->clone();
		assert( scorefxn_ );
	}
	bool tmpOpSet = false;
	std::string type( type_ );
	if ( (type == "psipred")||( type == "random" ) ) {
		task_operations::HighestEnergyRegionOperationOP op( new task_operations::DesignBySecondaryStructureOperation(
			blueprint_file_,psipred_cmd_,false,false ) );
		highestEnergyRegionOperation_ops_.push_back(op);
		tmpOpSet = true;
	}
	if ( type == "score" ) {
		task_operations::HighestEnergyRegionOperationOP op( new task_operations::HighestEnergyRegionOperation() );
		op->set_scorefxn( scorefxn_ );
		highestEnergyRegionOperation_ops_.push_back(op);
		tmpOpSet = true;
	}
	if ( (type == "packstat") ||(type == "random" ) ) {
		task_operations::HighestEnergyRegionOperationOP op( new task_operations::DesignByPackStatOperation() );
		highestEnergyRegionOperation_ops_.push_back(op);
		tmpOpSet = true;
	}
	if ( (type == "random_mutation")||(type == "random" ) ) {
		task_operations::HighestEnergyRegionOperationOP op( new task_operations::DesignRandomRegionOperation() );
		highestEnergyRegionOperation_ops_.push_back(op);
		tmpOpSet = true;
	}
	if ( tmpOpSet == false ) {
		utility_exit_with_message( "Bad type specified to RestrictRegion op: " + type  );
	}
	//initialize metric calculator. maybe not the best location for this?
	for ( core::Size ii=1; ii<=highestEnergyRegionOperation_ops_.size(); ++ii ) {
		std::string opName = highestEnergyRegionOperation_ops_[ii]->get_name();
		metric_stats_[opName] = std::pair< core::Size, core::Size >( 0, 0 );
	}
}

std::string
RestrictRegion::get_name() const
{
	return "RestrictRegion";
}

/// @brief Does the RestrictRegion moves
void
RestrictRegion::apply( core::pose::Pose & pose )
{
	if ( resfile_ == "" ) {
		utility_exit_with_message( "A resfile was not specified to the RestrictRegion mover." );
	}

	if ( permanently_restricted_residues_.size() == 0 ) {
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			permanently_restricted_residues_.push_back( "" );
		}
	}

	if ( permanently_restricted_residues_.size() != pose.total_residue() ) {
		utility_exit_with_message( "WARNING: pose size changed! Pose size=" + boost::lexical_cast< std::string >( pose.total_residue() ) + ", cache size=" + boost::lexical_cast< std::string >( permanently_restricted_residues_.size() ) );
	}

	// initialize the restricted residue list if necessary
	restricted_residues_.clear();
	TR << "Clearing list of restricted residues" << std::endl;
	for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
		restricted_residues_.push_back( permanently_restricted_residues_[i] );
	}

	if ( restricted_residues_.size() != pose.total_residue() ) {
		utility_exit_with_message( "WARNING: pose size changed! Pose size=" + boost::lexical_cast< std::string >( pose.total_residue() ) + ", cache size=" + boost::lexical_cast< std::string >( restricted_residues_.size() ) );
	}

	// set up the packer task
	core::pack::task::PackerTaskOP task;
	if ( task_factory_ ) {
		task = task_factory_->create_task_and_apply_taskoperations( pose );
	} else {
		task = core::pack::task::TaskFactory::create_packer_task( pose );
	}

	// check the current pose against the previous pose -- if it has changed, make the restrictions permanent
	if ( previous_pose() ) {
		bool changed( false );
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( previous_pose()->residue( i ).aa() != pose.residue( i ).aa() ) {
				TR << "Poses not identical for " << previous_pose()->residue( i ).name1() << i << pose.residue( i ).name1();
				// this means there has been a change
				changed = true;
				// if a change is found in the last set of residues restricted, permanently restrict it
				if ( permanent_restriction_ ) {
					bool found( false );
					for ( core::Size j=1; j<=last_residues_restricted().size(); ++j ) {
						TR << "-- saving changes!" << std::endl;
						permanently_restrict_aa( (*previous_pose_), task, last_residues_restricted()[j] );
						found = true;
					}
					if ( ! found ) {
						TR << "skipping" << std::endl;
					}
				}
			}
		}
		// increment counter if sequence changes occurred in the pose
		if ( changed ) {
			++( metric_stats_[ last_type_ ].first );
		}
	}
	previous_pose_ = core::pose::PoseOP( new core::pose::Pose( pose ) );
	// set operation to find worst region
	task_operations::HighestEnergyRegionOperationOP op = highestEnergyRegionOperation_ops_[numeric::random::random_range(1,highestEnergyRegionOperation_ops_.size())];
	type_ = op->get_name();
	//output stats on types
	TR << "type=" << type_ << std::endl;
	std::map< std::string, std::pair< core::Size, core::Size > >::iterator itr;
	for ( itr = metric_stats_.begin(); itr!= metric_stats_.end(); ++itr ) {
		TR << itr->first << ": " << itr->second.first << " / " << itr->second.second << " ; ";
	}
	TR << std::endl;
	last_type_ = type_;
	++(metric_stats_[ type_ ].second);
	op->set_regions_to_design( pose.total_residue() );
	utility::vector1< core::Size > residues( op->get_residues_to_design( pose ) );
	core::Size restrict_count( 0 );
	last_residues_restricted_.clear();
	for ( core::Size i=1; i<=residues.size(); ++i ) {
		core::Real random( numeric::random::uniform() );
		core::Size residue_idx( 0 );
		for ( core::Size j=1; j<=residues.size(); ++j ) {
			random = random * 2;
			if ( random > 1 ) {
				residue_idx = j;
				break;
			}
		}
		if ( residue_idx == 0 ) {
			residue_idx = residues.size();
		}
		TR << "Trying to restrict residue " << pose.residue( residues[residue_idx] ).name() << residues[residue_idx] << std::endl;
		if ( restrict_aa( pose, task, residues[residue_idx] ) ) {
			last_residues_restricted_.push_back( residues[residue_idx] );
			// mutate residue to something that's allowed
			std::string const allowed_types( residues_allowed( task, residues[residue_idx] ) );
			char const target_aa( allowed_types[ numeric::random::random_range( 1, allowed_types.size() ) - 1 ] );
			TR << "Mutating " << pose.residue( residues[residue_idx] ).name1() << residues[residue_idx] << target_aa << std::endl;
			protocols::simple_moves::MutateResidue mut_res( residues[residue_idx], target_aa );
			mut_res.apply( pose );
			++restrict_count;
		}

		// if we are done
		if ( restrict_count >= regions_to_mutate_ ||
				( i >= 1 && i == residues.size() ) ) {
			// re-create task after mutation to make sure task ops haven't changed and write resfile
			core::pack::task::PackerTaskOP task;
			if ( task_factory_ ) {
				task = task_factory_->create_task_and_apply_taskoperations( pose );
			} else {
				task = core::pack::task::TaskFactory::create_packer_task( pose );
			}

			// restrict new trps if the max has been reached
			if ( enable_max_trp_ ) {
				utility::vector1< bool > allowed_aa( 20, true );
				core::Size trp_count( 0 );
				for ( core::Size j=1; j<=pose.total_residue(); ++j ) {
					if ( pose.residue( j ).name1() == 'W' ) {
						++trp_count;
					}
				}
				if ( trp_count >= max_trp_ ) {
					for ( core::Size j=1; j<=pose.total_residue(); ++j ) {
						if ( pose.residue( j ).name1() != 'W' ) {
							allowed_aa[ core::chemical::aa_from_oneletter_code( 'W' ) ] = false;
							task->nonconst_residue_task( j ).restrict_absent_canonical_aas( allowed_aa );
							allowed_aa[ core::chemical::aa_from_oneletter_code( 'W' ) ] = true;
						} else {
							TR << "Not restricting TRP at position " << pose.residue( j ).name() << j << std::endl;
						}
					}
				}
			}

			write_resfile( pose, task );
			//TR << *task;
			TR << "Success; type=" << type_ << std::endl;
			return;
		}

		TR << "Removing" << std::endl;
		utility::vector1< core::Size > new_residues;
		for ( core::Size j=1; j<=residues.size(); ++j ) {
			if ( j != residue_idx ) {
				new_residues.push_back( residues[j] );
			}
		}
		residues = new_residues;
	}
	TR.Warning << "No mutations possible; residues=" << residues << std::endl;

	//preserve_fragment( pose, worst );

}

/// @brief initialize the resfile -- takes an input resfile, but converts it to a resfile in the output directory
void
RestrictRegion::initialize_resfile( std::string const & orig_resfile )
{
	if ( orig_resfile == "" ) {
		resfile_ = orig_resfile;
		return;
	}
	core::Size const pos( orig_resfile.find_last_of( '/' ) );
	if ( pos == std::string::npos ) {
		resfile_ = orig_resfile;
	} else {
		resfile_ = orig_resfile.substr( pos+1, orig_resfile.size() );
	}
	/*if ( resfile_ == orig_resfile ) {
	resfile_ += "_out";
	}*/
	TR << "Resfile=" << resfile_ << std::endl;
	std::string resfile_string;
	utility::io::izstream file( orig_resfile );
	if ( !file ) {
		resfile_string = "start\n";
	} else {
		utility::slurp( file, resfile_string );
	}
	file.close();
	utility::io::ozstream outfile( resfile_ );
	if ( !outfile ) {
		utility_exit_with_message( "Could not write to " + resfile_ );
	}
	outfile << resfile_string;
	outfile.close();

}

/// @brief writes a resfile from the task
void
RestrictRegion::write_resfile( core::pose::Pose const & pose,  core::pack::task::PackerTaskCOP task ) const
{
	core::Real const distance_cutoff( 8.0 );
	utility::io::ozstream outfile( resfile_ );
	if ( !outfile ) {
		utility_exit_with_message( "Could not write to " + resfile_ );
	}
	outfile << "USE_INPUT_SC" << '\n' << "start" << '\n';
	for ( core::Size i=1; i<=task->total_residue(); ++i ) {
		if ( ! pose.residue( i ).is_protein() ) {
			continue;
		}

		// check to see if the existing residue is allowed; if not, add it to the PIKAA string or remove from NOTAA
		// to be allowed, the native AA must be 1) allowed by the task operations, and 2) not in the restricted list
		// note that allowing an amino acid in this way will remove it from the restricted list, which means that mutation
		// may be restricted again in a later call to this mover.
		bool const native_allowed( residue_is_allowed( task, i, pose.residue( i ).aa() ) &&
			! is_restricted( pose.residue( i ).name1(), i ) );

		// find the minimum distance to a residue that is being restricted. If it is far, it may not fall into the design shell and should be set to repack, which means the native AA must be allowed
		core::Real min_dist( 100000 ); // ridiculously huge distance
		for ( core::Size ii=1; ii<=last_residues_restricted().size(); ++ii ) {
			core::Real const distance( pose.residue( i ).xyz( pose.residue( i ).nbr_atom() ).distance( pose.residue( last_residues_restricted()[ii] ).xyz( pose.residue( last_residues_restricted()[ii] ).nbr_atom() ) ) );
			if ( distance < min_dist ) {
				min_dist = distance;
			}
		}
		std::string const allowed_resis( residues_allowed( task, i ) );
		bool const repacking_nonnative( !native_allowed &&
			(  min_dist >= distance_cutoff || allowed_resis.size() == 0 ) );
		std::string cmd( " PIKAA " );
		cmd += allowed_resis;
		std::string notaa( restricted_residues_[i] );
		if ( repacking_nonnative ) {
			TR << pose.residue(i).name1() << " NOT allowed at pos " << i << "; " << cmd << "; dist=" << min_dist << "; fixing!" << std::endl;
			cmd = allow_in_resfile_line( cmd, pose.residue( i ).name1() );
			notaa = "";
			for ( core::Size j=0; j<restricted_residues_[i].size(); ++j ) {
				if ( restricted_residues_[i][j] != pose.residue( i ).name1() ) {
					notaa += restricted_residues_[i][j];
				}
			}
			//restricted_residues_[i] = tmp_str;
		}

		// if we have previously restricted any of these residues, add them here
		if ( notaa.size() ) {
			cmd += " NOTAA " + notaa;
		}

		// if we are greater than 8 A away, set to repack only
		if ( min_dist >= distance_cutoff ) {
			cmd += " NATRO";
		}

		outfile << i << ' ' << pose.pdb_info()->chain( i ) << cmd << '\n';
	}
	outfile.close();
}


/// @brief restricts design at the specified position such that only the existing amino acid is allowed, and writes this to a resfile
void
RestrictRegion::preserve_residue( core::pose::Pose const & /*pose*/, core::pack::task::PackerTaskOP task, core::Size const seqpos ) const
{
	// just restrict to repacking and we're done
	task->nonconst_residue_task( seqpos ).restrict_to_repacking();
}

/// @brief restricts design at the specified position such that the existing amino acid is not allowed, and writes this action to a resfile so that other movers can use it. Also checks for allowed amino acids based on task factory to ensure there is at least one allowed residue type. Returns true if successful, false if not successful.
bool
RestrictRegion::restrict_aa( core::pose::Pose const & pose, core::pack::task::PackerTaskOP task, core::Size const seqpos )
{
	// don't restrict this residue if it is the only possibility at this position
	if ( residues_allowed( task, seqpos ).size() <= 1 ) {
		return false;
	}

	// don't restrict this residue if it's already restricted
	if ( is_restricted( pose.residue( seqpos ).name1(), seqpos ) ) {
		TR << "Restricted already." << std::endl;
		return true;
	}

	restricted_residues_[ seqpos ] += pose.residue( seqpos ).name1();
	return true;
}

/// @brief Permanently restricts design at the specified position such that the existing amino acid is not allowed.
void
RestrictRegion::permanently_restrict_aa( core::pose::Pose const & pose, core::pack::task::PackerTaskOP /*task*/, core::Size const seqpos )
{
	// if the amino acid is not found in the string of permanently restricted residues, add it in there
	if ( permanently_restricted_residues_[ seqpos ].find( pose.residue( seqpos ).name1() ) == std::string::npos ) {
		permanently_restricted_residues_[ seqpos ] += pose.residue( seqpos ).name1();
	}
}


/// @brief static method that tells the last residue affected by any RestrictRegion mover. This method also makes sure that the pose length hasn't changed. If this is called before any instance of the mover::apply() is called, returns 0.
utility::vector1< core::Size > const &
RestrictRegion::last_residues_restricted()
{
	return last_residues_restricted_;
}

/// @brief static method that tells the last pose returned by any RestrictRegion mover. This method also makes sure that the pose length hasn't changed. If this is called before any instance of the mover::apply() is called, returns 0.
core::pose::PoseCOP
RestrictRegion::previous_pose()
{
	return previous_pose_;
}

/// @brief tells whether a given fragment is compatible with the task
/// basically, this function counts the number of amino acids in each fragment that are compatible with the task
core::Size
RestrictRegion::compatible_with_task( core::pack::task::PackerTaskOP task,
	core::Size const frag_id,
	core::fragment::FrameOP frame ) const
{
	core::Size const task_start_seqpos( frame->start() );
	std::string const frag_seq( frame->fragment_ptr( frag_id )->sequence() );
	core::Size count( 0 );
	for ( core::Size i=0; i<frag_seq.size(); ++i ) {
		if ( residue_is_allowed( task, task_start_seqpos + i, core::chemical::aa_from_oneletter_code( frag_seq[i] ) ) ) {
			++count;
		}
	}
	return count;
}

/// @brief function that tells whether a given amino acid is allowed at a certain position
bool
residue_is_allowed( core::pack::task::PackerTaskCOP task, core::Size const seqpos, core::chemical::AA const aa )
{
	core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter res_iter( task->residue_task( seqpos ).allowed_residue_types_begin() );
	while ( res_iter != task->residue_task( seqpos ).allowed_residue_types_end() ) {
		if ( (*res_iter)->aa() == aa ) {
			return true;
		}
		++res_iter;
	}
	return false;
}

/// @brief function that checks to see whether an amino acid is restricted at a certain position
bool
RestrictRegion::is_restricted( char const aa, core::Size const seqpos ) const
{
	runtime_assert( seqpos <= restricted_residues_.size() );
	// if this amino acid is found in the restricted list for seqpos, will return true
	return ( restricted_residues_[seqpos].find( aa ) != std::string::npos );
}

/// @brief function that tells how many residues are allowed at a certain position
std::string
RestrictRegion::residues_allowed( core::pack::task::PackerTaskCOP task, core::Size const seqpos ) const
{
	std::string aas;
	core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter res_iter( task->residue_task( seqpos ).allowed_residue_types_begin() );
	while ( res_iter != task->residue_task( seqpos ).allowed_residue_types_end() ) {
		// don't do anything if this residue type is restricted
		if ( ! is_restricted( (*res_iter)->name1(), seqpos ) ) {
			// add the amino acid code if it's not in there already
			if ( aas.find( (*res_iter)->name1() ) == std::string::npos ) {
				aas += (*res_iter)->name1();
			}
		}
		++res_iter;
	}
	return aas;
}

/// @brief alter the resfile line such that the amino acid given is listed as allowed and return the modified version
std::string
allow_in_resfile_line( std::string const & cmd_orig, char const aa )
{
	utility::vector1< std::string > cmds( utility::string_split( cmd_orig, ' ' ) );
	TR << "looking at " << cmd_orig << std::endl;
	std::string ret_cmd;
	// advance through the list of tokens
	for ( core::Size i=1; i<=cmds.size(); ++i ) {
		core::Size pos( cmds[i].find( "NOTAA" ) );
		if ( pos != std::string::npos ) {
			TR << "NOTAA line = " << cmds[i] << std::endl;
			++i;
			std::string newstr;
			for ( core::Size j=0; j<cmds[i].size(); ++j ) {
				if ( cmds[i][j] != aa ) {
					newstr += cmds[i][j];
				}
			}
			if ( newstr.size() ) {
				ret_cmd += cmds[i-1] + ' ' + newstr + ' ';
			}
		}
		pos = cmds[i].find( "PIKAA" );
		if ( pos != std::string::npos ) {
			ret_cmd += cmds[i] + ' ';
			++i;
			if ( cmds[i].find( aa ) == std::string::npos ) {
				cmds[i] += aa;
			}
		}
		ret_cmd += cmds[i] + ' ';
	}

	return ret_cmd;
}


} // namespace denovo_design
} // namespace devel


