// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/MatchPositionModifiers.cc
/// @brief implementations for MatchPositionModifiers
/// @author Florian Richter (floric@u.washington.edu ), may 2010

// unit headers
#include <protocols/match/MatchPositionModifiers.hh>

//package headers
#include <protocols/match/MatcherTask.hh>

//project headers
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/MatchConstraintFileInfo.hh>

//Numeric headers
#include <numeric/constants.hh>
#include <numeric/xyzVector.hh>

//utility headers
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {

static THREAD_LOCAL basic::Tracer tr( "protocols.match.MatchPositionModifiers" );

/// @brief "factory" function to create the match position modifiers
MatchPositionModifierCOP
create_match_position_modifier(
	std::string const & mpm_name,
	core::Size geom_cst,
	utility::vector1< std::string > const & input_tokens )
{
	if ( mpm_name == "ss" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new SecondaryStructureMPM( input_tokens ) ) );
	else if ( mpm_name == "num_neighbors" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new NumNeighborsMPM( input_tokens ) ) );
	else if ( mpm_name == "bfactor" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new BfactorMPM( input_tokens ) ) );
	else if ( mpm_name == "all" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new AddAllPositionsMPM() ) );
	else if ( mpm_name == "no_c_n_term" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new RemoveNorCTermMPM( input_tokens ) ) );
	else if ( mpm_name == "task_operation" ) return MatchPositionModifierCOP( MatchPositionModifierOP( new TaskOperationMPM( geom_cst, input_tokens ) ) );
	return NULL;
}

MatchPositionModifier::MatchPositionModifier(){}

MatchPositionModifier::~MatchPositionModifier(){}

SecondaryStructureMPM::SecondaryStructureMPM( utility::vector1< std::string > const & input_tokens )
: MatchPositionModifier()
{
	if ( input_tokens.size() < 3 ) utility_exit_with_message("Not enough information given to initialize SecondaryStructureMPM");
	for ( core::Size i = 2; i < input_tokens.size(); ++i ) {
		if ( input_tokens[ i ] == "ss_char" ) {
			desired_ss_chars_.insert( input_tokens[i+1][0] );
			tr << "SecondaryStructureMPM requires positions to have ss char " << input_tokens[i+1][0] << "." << std::endl;
			i += 2;
		} else if ( input_tokens[ i ] == "ss_motif" ) {
			ss_motifs_.push_back( input_tokens[ i+1 ] );
			tr << "SecondaryStructureMPM requires positions to have motif " << input_tokens[i+1] << "." << std::endl;
			i += 2;
		} else tr << "Token " << input_tokens[ i ] << " could not be understood by SecondaryStructureMPM and will be ignored." << std::endl;
	}
}

SecondaryStructureMPM::~SecondaryStructureMPM(){}

utility::vector1< core::Size >
SecondaryStructureMPM::modified_match_positions(
	utility::vector1< core::Size > const & original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP //mtask
) const
{
	utility::vector1< core::Size > to_return;
	std::string remove_string("");

	core::scoring::dssp::Dssp pose_ss( match_pose );

	for ( core::Size i =1; i <= original_positions.size(); ++i ) {

		bool position_passes(false);
		core::Size seqpos( original_positions[i] );
		//first we'll check if any excat ss chars have been specified
		if ( desired_ss_chars_.size() != 0 ) {
			if ( desired_ss_chars_.find( pose_ss.get_dssp_secstruct( seqpos ) ) != desired_ss_chars_.end() ) {
				position_passes = true;
			}
		} //if( desired_ss_chars_.size() != 0 )

		//then we'll check in any more complex motifs have been specified
		for ( core::Size motif = 1; motif <= ss_motifs_.size(); ++motif ) {

			if ( position_passes ) break;

			/// @details helix_nterm logic: if seqpos has ss char H, it needs to be in
			/// the beginning of the helix, i.e. there needs to be another ss char within
			/// three positions upstream of seqpos
			/// is seqpos doesn't have ss char H, a helix needs to commence within 3 positions
			/// downstream of seqpos
			if ( ss_motifs_[ motif ] == "helix_nterm" ) {
				if ( pose_ss.get_dssp_secstruct( seqpos ) == 'H' ) {
					for ( core::Size j = seqpos - 1; (j >= seqpos - 3) && (j > 0 ); --j ) {
						if ( pose_ss.get_dssp_secstruct( j ) != 'H' ) {
							position_passes = true;
							break;
						}
					}
				} else { //if seqpos is in helix
					for ( core::Size j = seqpos; (j <= seqpos + 3) && ( j <= match_pose.total_residue() ); ++j ) {
						if ( pose_ss.get_dssp_secstruct( j ) == 'H' ) {
							position_passes = true;
							break;
						}
					}
				}
			} else { //helix_nterm motif
				tr << "WARNING: SecondaryStructureMPM doesn't know how to interpret motif '" << ss_motifs_[ motif ] << "'." << std::endl;
			}
		}//for( core::Size motif = 1; motif <= ss_motifs_.size(); ++motif)
		if ( position_passes ) {
			to_return.push_back( seqpos );
		} else remove_string += utility::to_string( seqpos ) + "+";
	} //loop over all original positions
	tr << "SecondaryStructureMPM removed the following match positions " << remove_string << "." << std::endl;
	return to_return;
}


NumNeighborsMPM::NumNeighborsMPM( utility::vector1< std::string > const & input_tokens )
: MatchPositionModifier(), min_neighbors_(0), max_neighbors_(0),
	com_vector_criterion_(false), both_criteria_needed_to_pass_(false),
	min_com_vector_ang_cos_( 1.0 ), max_com_vector_ang_cos_( -1.0 )
{
	if ( input_tokens.size() < 3 ) utility_exit_with_message("Not enough information given to initialize NumNeighborsMPM");
	for ( core::Size i = 2; i <= input_tokens.size(); ++i ) {
		if ( input_tokens[ i ] == "min_neighbors" ) {
			min_neighbors_ =  (core::Size) atoi( input_tokens[i+1].c_str() );
			tr << "NumNeighborsMPM will only allow positions that have at least " << min_neighbors_ << " 10A neighbors." << std::endl;
			i++;
		} else if ( input_tokens[ i ] == "max_neighbors" ) {
			max_neighbors_ =  (core::Size) atoi( input_tokens[i+1].c_str() );
			tr << "NumNeighborsMPM will only allow positions that have no more than " << max_neighbors_ << " 10A neighbors." << std::endl;
			i++;
		} else if ( input_tokens[ i ] == "min_com_vector_ang" ) {
			com_vector_criterion_ = true;
			min_com_vector_ang_cos_ = cos( (( core::Real) atof(input_tokens[i+1].c_str() ) )* numeric::constants::f::degrees_to_radians);
			i++;
		} else if ( input_tokens[ i ] == "max_com_vector_ang" ) {
			com_vector_criterion_ = true;
			max_com_vector_ang_cos_ = cos( (( core::Real) atof(input_tokens[i+1].c_str() ) )* numeric::constants::f::degrees_to_radians );;
			i++;
		} else if ( input_tokens[ i ] == "both_criteria_needed_to_pass" ) {
			both_criteria_needed_to_pass_ = true;
		} else tr << "Token " << input_tokens[ i ] << " could not be understood by NumNeighborsMPM and will be ignored." << std::endl;
	}
}

NumNeighborsMPM::~NumNeighborsMPM(){}

utility::vector1< core::Size >
NumNeighborsMPM::modified_match_positions(
	utility::vector1< core::Size > const & original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP //mtask
) const
{
	core::scoring::TenANeighborGraph const & cur_graph = match_pose.energies().tenA_neighbor_graph();
	utility::vector1< core::Size > to_return;
	std::string remove_string("");

	//if we need to calculate com
	core::Vector center_of_mass(0,0,0);
	if ( com_vector_criterion_ ) {
		for ( Size i = 1; i <= match_pose.total_residue(); ++i ) center_of_mass += match_pose.residue(i).nbr_atom_xyz();
		center_of_mass /= match_pose.total_residue();
		//tr << "Center of mass is " << center_of_mass.x() << " " << center_of_mass.y() << " " << center_of_mass.z() << std::endl;
	}

	for ( core::Size i =1; i <= original_positions.size(); ++i ) {

		core::Size neighbors( cur_graph.get_node( original_positions[i] )->num_neighbors_counting_self() - 1 );

		bool neighborpass( true ), com_vect_pass(true);
		if ( com_vector_criterion_ ) com_vect_pass = passes_com_vector_criterion( original_positions[i], match_pose, center_of_mass );

		if ( (min_neighbors_ != 0 ) && (neighbors < min_neighbors_) ) neighborpass = false;
		if ( (max_neighbors_ != 0 ) && (neighbors > max_neighbors_) ) neighborpass = false;
		bool pass(neighborpass);

		if ( com_vector_criterion_ ) {
			if ( both_criteria_needed_to_pass_ ) {
				pass = (neighborpass && com_vect_pass);
				//tr << "mpf resi " << original_positions[i] << " neighborpass is " << neighborpass << " compass is " << com_vect_pass << " pass is " << pass << std::endl;
			} else pass = (neighborpass || com_vect_pass);
		}
		if ( pass ) to_return.push_back( original_positions[i] );
		else remove_string += utility::to_string( original_positions[ i ] ) + "+";
	}
	tr << "NumNeighborsMPM removed the following match positions " << remove_string << "." << std::endl;
	return to_return;
}

bool
NumNeighborsMPM::passes_com_vector_criterion(
	core::Size seqpos,
	core::pose::Pose const & pose,
	core::Vector const & com
) const
{
	core::Vector seqpos_to_com( com - pose.residue( seqpos ).atom("CA").xyz() );
	core::Vector ca_cb( (pose.residue( seqpos ).aa() == core::chemical::aa_gly ? pose.residue( seqpos ).atom("2HA") : pose.residue( seqpos ).atom("CB")).xyz() - pose.residue( seqpos ).atom("CA").xyz() );

	//core::Real com_cos( numeric::cos_of( seqpos_to_com, pose.residue( seqpos ).nbr_atom_xyz() ) );
	core::Real com_cos( seqpos_to_com.dot( ca_cb ) / (seqpos_to_com.length() * ca_cb.length() ) );
	//std::cerr << seqpos << " com_cos is " << com_cos ;
	if ( ( com_cos < min_com_vector_ang_cos_ ) && (com_cos > max_com_vector_ang_cos_ ) ) {
		//std::cerr << " passing " << std::endl;
		return true;
	}
	//std::cerr << " not passing " << std::endl;
	return false;
}


BfactorMPM::BfactorMPM( utility::vector1< std::string > const & input_tokens )
: MatchPositionModifier(), use_relative_bfactors_(false), all_bfactors_zero_(false), max_bfactor_(0.0)
{
	if ( input_tokens.size() < 3 ) utility_exit_with_message("Not enough information given to initialize BfactorMPM");
	for ( core::Size i =2; i < input_tokens.size(); ++i ) {
		if ( input_tokens[ i ] == "relative" ) {
			use_relative_bfactors_ = true;
			max_bfactor_ = (core::Real) atof( input_tokens[i+1].c_str() );
			tr << "BfactorMPM will only allow positions that have a relative B-factor of not more than" << max_bfactor_ <<"." << std::endl;
		}
		if ( input_tokens[ i ] == "absolute" ) {
			use_relative_bfactors_ = false;
			max_bfactor_ = (core::Real) atof( input_tokens[i+1].c_str() );
			tr << "BfactorMPM will only allow positions that have an absolute B-factor of not more than" << max_bfactor_ <<"." << std::endl;
		} else tr << "Token " << input_tokens[ i ] << " could not be understood by BfactorMPM and will be ignored." << std::endl;
	}
}

BfactorMPM::~BfactorMPM(){}

utility::vector1< core::Size >
BfactorMPM::modified_match_positions(
	utility::vector1< core::Size > const & original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP //mtask
) const
{
	utility::vector1< core::Size > to_return;
	std::string remove_string("");
	utility::vector1< core::Real > bfactors( this->get_ca_bfactors( match_pose ) );

	if ( all_bfactors_zero_ ) {
		tr << "Warning: all bfactors in the pose were 0, meaning they were probably wiped. BfactorMPM will not modify match positions." << std::endl;
		to_return = original_positions;
		return to_return;
	}

	for ( core::Size i =1; i <= original_positions.size(); ++i ) {
		if ( bfactors[ original_positions[ i ] ] <= max_bfactor_ ) to_return.push_back( original_positions[i] );
		else remove_string += utility::to_string( original_positions[ i ] ) + "+";
	}
	tr << "BfactorMPM removed the following match positions " << remove_string << "." << std::endl;
	return to_return;
}

utility::vector1< core::Real >
BfactorMPM::get_ca_bfactors( core::pose::Pose const & pose ) const
{
	utility::vector1< core::Real > bfactors;
	core::pose::PDBInfo const & pdb_info( *(pose.pdb_info()) );
	all_bfactors_zero_ = true;

	if ( use_relative_bfactors_ ) {
		core::Real max_bfactor(0.0);
		for ( core::Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) {
			if ( ! pose.residue_type( seqpos ).is_protein() ) bfactors.push_back( pdb_info.temperature( seqpos, 1 ) );
			else bfactors.push_back( pdb_info.temperature( seqpos, pose.residue( seqpos ).atom_index( "CA" ) ) );
			if ( bfactors[ seqpos ] > max_bfactor ) max_bfactor = bfactors[ seqpos ];
			if ( bfactors[ bfactors.size() ] > 0.0 ) all_bfactors_zero_ = false;
		}
		if ( max_bfactor == 0.0 ) { //in this case we return
			all_bfactors_zero_ = true;
			return bfactors;
		}

		for ( core::Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) bfactors[ seqpos ] /= max_bfactor;
	} else {
		for ( core::Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ) {
			if ( ! pose.residue( seqpos ).is_protein() ) bfactors.push_back( pdb_info.temperature( seqpos, 1 ) );
			else bfactors.push_back( pdb_info.temperature( seqpos, pose.residue( seqpos ).atom_index( "CA" ) ) );

			if ( bfactors[ bfactors.size() ] > 0.0 ) all_bfactors_zero_ = false;
		}
	}
	return bfactors;
}

AddAllPositionsMPM::AddAllPositionsMPM()
: MatchPositionModifier()
{}

AddAllPositionsMPM::~AddAllPositionsMPM(){}

utility::vector1< core::Size >
AddAllPositionsMPM::modified_match_positions(
	utility::vector1< core::Size > const &, //original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP //mtask
) const
{
	utility::vector1< core::Size > to_return;
	for ( core::Size i = 1; i <= match_pose.total_residue(); ++i ) {
		if ( match_pose.residue_type( i ).is_protein() ) to_return.push_back( i );
	}
	return to_return;
}

RemoveNorCTermMPM::RemoveNorCTermMPM( utility::vector1< std::string > const & input_tokens )
: MatchPositionModifier(), cterm_length_(0), nterm_length_(0)
{
	for ( core::Size i =2; i < input_tokens.size(); ++i ) {
		if ( input_tokens[ i ] == "cterm" ) {
			runtime_assert( input_tokens.size() > i );
			cterm_length_ = (core::Size) atoi( input_tokens[i+1].c_str() );
		} else if ( input_tokens[ i ] == "nterm" ) {
			runtime_assert( input_tokens.size() > i );
			nterm_length_ = (core::Size) atoi( input_tokens[i+1].c_str() );
		}
	}
}

RemoveNorCTermMPM::~RemoveNorCTermMPM(){}

utility::vector1< core::Size >
RemoveNorCTermMPM::modified_match_positions(
	utility::vector1< core::Size > const & original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP //mtask
) const
{
	utility::vector1< core::Size > to_return;
	core::Size cterm = match_pose.total_residue();

	if ( cterm_length_ != 0 ) { //we have to determine the cterminus of the protein
		for ( core::Size i = match_pose.total_residue(); i > 0; --i ) {
			if ( match_pose.residue_type( i ).is_protein() ) {
				cterm = i;
				break;
			}
		}
		cterm = cterm - cterm_length_;
	}

	for ( core::Size i = 1; i <= original_positions.size(); ++i ) {
		if ( (original_positions[i] >= nterm_length_) && (original_positions[i] <= cterm) ) to_return.push_back( original_positions[i] );

	} // loop over all original positions

	return to_return;
}

/// @details
/// this is a little tricky, we have to reassemble a tag out of the tokens
/// and then call the TaskOperationFactory::init and some other stuff like that...
TaskOperationMPM::TaskOperationMPM(
	core::Size which_geom_cst,
	utility::vector1< std::string > const & input_tokens
)
: MatchPositionModifier(), which_geom_cst_(which_geom_cst), task_op_(/* NULL */)
{
	//1. reassemble tag components into string
	std::string tagstring(input_tokens[2]);
	for ( core::Size i = 3; i <= input_tokens.size(); ++i ) tagstring = tagstring + " " + input_tokens[i];

	//let's make sure the dumb user has actually supplied a proper tag
	runtime_assert(input_tokens[2].substr(0,1) == "<");
	runtime_assert(input_tokens[input_tokens.size()].substr( input_tokens[input_tokens.size()].length() - 1,1) == ">" );
	std::string task_op_name = utility::trim( input_tokens[2], "<");
	tr << "TaskOperationMPM getting task_op of type " << task_op_name << " with tag '" << tagstring << "'." << std::endl;

	//2. instantiate a tag object from the string
	utility::tag::TagOP tag( new utility::tag::Tag() );
	std::istringstream tagstream(tagstring);
	tag->read(tagstream);

	//3. make task op
	basic::datacache::DataMap datamap;
	core::pack::task::operation::TaskOperationFactory const * topfac( core::pack::task::operation::TaskOperationFactory::get_instance() );
	task_op_ = topfac->newTaskOperation( task_op_name, datamap, tag );

	//4. yay :)
}

TaskOperationMPM::~TaskOperationMPM(){}

/// @details
/// generates a task based on the parsed task operation,
/// and then allows matching at all positions where the
/// upstream residues specified in the cstfile are allowed
/// in the position's residue_task.
utility::vector1< core::Size >
TaskOperationMPM::modified_match_positions(
	utility::vector1< core::Size > const & original_positions,
	core::pose::Pose const & match_pose,
	protocols::match::MatcherTaskCOP mtask
) const
{
	using namespace core::pack::task;
	TaskFactory tfactory;
	tfactory.push_back( task_op_);
	PackerTaskOP ptask(tfactory.create_task_and_apply_taskoperations( match_pose ) );
	utility::vector1< core::Size > to_return;

	//the following code is not the most efficient, searching through vectors/lists, etc
	//but it'll only be done once per matcher run, and the vectors are usually not very large
	//so it shouldn't matter that much
	utility::vector1< core::chemical::ResidueTypeCOP> const & upstream_restypes( mtask->enz_input_data()->mcfi_list( which_geom_cst_ )->upstream_restypes() );

	for ( core::Size i = 1; i <= original_positions.size(); ++i ) {
		ResidueLevelTask const & restask( ptask->residue_task( original_positions[i] ));
		if ( restask.being_packed() ) { //let's only match repackable positions

			for ( ResidueLevelTask::ResidueTypeCOPListConstIter restype_it( restask.allowed_residue_types_begin()), restype_it_end( restask.allowed_residue_types_end() ); restype_it != restype_it_end; ++restype_it ) {

				//let's be somewhat generous: if any of the upstream residues specified in the cstfile
				//is allowed at this residue, we consider it good for matching
				//note: we're doing name3 comparison instead of pointer comparison here because
				//of variant type uncertainties
				bool name3_found( false );
				for ( utility::vector1< core::chemical::ResidueTypeCOP>::const_iterator upres_it( upstream_restypes.begin() ), upres_end(upstream_restypes.end()); upres_it != upres_end; ++ upres_it ) {
					if ( (*restype_it)->name3() == (*upres_it)->name3() ) {
						//tr << "ARRG residue " << (*restype_it)->name3() << " allowed at pos " << i << ", set to matching " << std::endl;
						to_return.push_back( original_positions[i] );
						name3_found = true;
						break;
					}
				}
				if ( name3_found ) break;
			} //loop over packer residue types at this position
		} //if being repacked
	} // loop over original positions
	return to_return;
}


}
}
