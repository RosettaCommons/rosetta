// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file devel/flxbb/filters/DesignBySecondaryStructure.cc
/// @brief Design residues that don't match the predicted secondary structure.
/// @author Tom Linsky (tlinsky@uw.edu)

// unit headers
#include <devel/denovo_design/task_operations/DesignBySecondaryStructure.hh>
#include <devel/denovo_design/task_operations/DesignBySecondaryStructureCreator.hh>

// package headers
#include <devel/denovo_design/filters/PsiPredInterface.hh>
#include <protocols/flxbb/utility.hh>
#include <protocols/moves/DsspMover.hh>

// project headers
#include <core/conformation/Residue.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/ss_prediction/SS_predictor.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>

// utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static thread_local basic::Tracer TR( "devel.denovo_design.task_operations.DesignBySecondaryStructure" );

namespace devel {
namespace denovo_design {
namespace task_operations {
	// helper functions

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Looks for unknown amino acids in the pose and returns their indices
utility::vector1<core::Size>
find_ligands( core::pose::Pose const & pose )
{
	utility::vector1<core::Size> retval;
	// look at each amino acid to see if it is of unknown type
	for( core::Size i = 1; i <= pose.total_residue(); i++ ) {
		if ( ! pose.residue( i ).is_protein() ) {
			TR << "Residue " << i << " (type=" << pose.residue(i).name3() << ") is probably a ligand" << std::endl;
			retval.push_back( i );
		}
	}
	// WARNING: This doesn't rearrange residue numbers to put the ligands at the end
	// However, this behavior IS important for many functions
	return retval;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Returns a secondary structure string including "L" characters for ligand residues
std::string
secstruct_with_ligand( std::string const & secstruct, utility::vector1<core::Size> const & liglist )
{
	utility::vector1<core::Size> ligands = liglist;
	std::string newstr = secstruct;
	// sorting the ligand list before this operation is mandatory
	sort( ligands.begin(), ligands.end() );
	for( core::Size i = 1; i <= ligands.size(); i++ ) {
		//tr << "Trying to add ligand at position " << ligands[i] << " into " << newstr << std::endl;
		runtime_assert( ligands[i] <= newstr.size() + 1 );
		runtime_assert( ligands[i] >= 1 );
		newstr = newstr.substr( 0, ligands[i]-1 ) + 'L' + newstr.substr( ligands[i]-1, (newstr.size()+1)-ligands[i] );
		//tr << "new string is " << newstr << std::endl;
	}
	return newstr;
}

  // default constructor
DesignBySecondaryStructureOperation::DesignBySecondaryStructureOperation()
	: HighestEnergyRegionOperation(),
		blueprint_ss_( "" ),
		pred_ss_( "" ),
		prevent_native_aa_( false ),
		prevent_bad_point_mutants_( false ),
		psipred_interface_( NULL ),
		ss_predictor_( NULL )
{}

// copy constructor
DesignBySecondaryStructureOperation::DesignBySecondaryStructureOperation( DesignBySecondaryStructureOperation const & rval )
  : HighestEnergyRegionOperation( rval ),
		blueprint_ss_( rval.blueprint_ss_ ),
		pred_ss_( "" ),
		prevent_native_aa_( rval.prevent_native_aa_ ),
		prevent_bad_point_mutants_( rval.prevent_bad_point_mutants_ ),
		psipred_interface_( rval.psipred_interface_ )
{}

// value constructor
DesignBySecondaryStructureOperation::DesignBySecondaryStructureOperation( std::string const bp_file,
																																					std::string const cmd,
																																					bool const prevent_native,
																																					bool const prevent_bad_point_mutants )
	: HighestEnergyRegionOperation(),
		pred_ss_( "" ),
		prevent_native_aa_( prevent_native ),
  	prevent_bad_point_mutants_( prevent_bad_point_mutants ),
	  psipred_interface_( NULL )
{
	if ( cmd != "" ) {
		psipred_interface_ = new denovo_design::filters::PsiPredInterface( cmd );
	} else {
		ss_predictor_ = new protocols::ss_prediction::SS_predictor( "HLE" );
	}
	initialize_blueprint_ss( bp_file );
}

// destructor
DesignBySecondaryStructureOperation::~DesignBySecondaryStructureOperation()
{}

/// @brief make clone
core::pack::task::operation::TaskOperationOP
DesignBySecondaryStructureOperation::clone() const {
  return new DesignBySecondaryStructureOperation( *this );
}

/// @brief utility function that compares two resid-probability pairs and returns true of the probability of the first is greater than probability of the second
bool compare_prob( std::pair< core::Size, core::Real > const & p1,
									 std::pair< core::Size, core::Real > const & p2 ) {
	return ( p1.second < p2.second );
}

/// @brief Gets a list of residues for design
utility::vector1< core::Size >
DesignBySecondaryStructureOperation::get_residues_to_design( core::pose::Pose const & pose ) const
{
	std::string wanted_ss( blueprint_ss_ );
	// if a blueprint is not specified, use DSSP to determine pose secondary structure.
	if ( wanted_ss == "" ) {
		protocols::moves::DsspMover dssp;
		core::pose::Pose posecopy( pose );
		dssp.apply( posecopy );
		wanted_ss = posecopy.secstruct();
		std::string pruned_ss( "" );
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue(i).is_protein() ) {
				pruned_ss += wanted_ss[i-1];
			}
		}
		wanted_ss = pruned_ss;
	}
	std::string bp_ss( wanted_ss );
	//utility::vector1<core::Size> const & ligands( find_ligands( pose ) );
	//std::string bp_ss( secstruct_with_ligand( wanted_ss, ligands ) );

	utility::vector1< core::Size > residues_to_design;
	TR << "Calculating predicted SS" << std::endl;
	std::string pred_ss( "" );
	utility::vector1< core::Real > psipred_prob;
	if ( ! psipred_interface_ ) {
		runtime_assert( ss_predictor_ );
		std::string sequence( "" );
		for ( core::Size i=1; i<=pose.total_residue(); ++i ) {
			if ( pose.residue( i ).is_protein() )
				sequence += pose.residue( i ).name1();
		}
		utility::vector1< utility::vector1< core::Real > > svm_pred( ss_predictor_->predict_ss( sequence ) );
		for ( core::Size i=1; i<=wanted_ss.size(); ++i ) {
			psipred_prob.push_back( protocols::ss_prediction::get_prob( wanted_ss[i-1], svm_pred[i] ) );
			pred_ss += protocols::ss_prediction::get_label( svm_pred[i] );
		}
	} else {
		denovo_design::filters::PsiPredResult const & psipred_result( psipred_interface_->run_psipred( pose, wanted_ss ) );
		psipred_prob = psipred_result.psipred_prob;
		pred_ss = psipred_result.pred_ss;
	}

	if ( regions_to_design() == 0 ) {
		TR << "Going to design all residues with psipred that doesn't match blueprint secondary structure." << std::endl;
		//pred_ss = secstruct_with_ligand( pred_ss, ligands );
		TR << "Blueprint SS = " << bp_ss << std::endl;
		TR << "Predicted SS = " << pred_ss << std::endl;
		residues_to_design = denovo_design::filters::nonmatching_residues( bp_ss, pred_ss );
	} else {
		TR << "Going to design the worst " << regions_to_design() << " residues." << std::endl;
		// create a map and insert pairs for residue and probabilities
		utility::vector1< std::pair< core::Size, core::Real > > res_to_prob;
		for ( core::Size j=1; j<=psipred_prob.size(); ++j ) {
			res_to_prob.push_back( std::pair< core::Size, core::Real >( j, psipred_prob[j] ) );
		}
		// sort the map based on the psipred probability
		std::sort( res_to_prob.begin(), res_to_prob.end(), compare_prob );
		TR << "Top member of sorted probability list is " << res_to_prob[1].second << std::endl;
		for ( core::Size j=1; ( j<=res_to_prob.size() && j<=regions_to_design() ); ++j ) {
			residues_to_design.push_back( res_to_prob[j].first );
		}
	}
	return residues_to_design;
}

/// @brief apply
void
DesignBySecondaryStructureOperation::apply( Pose const & pose, core::pack::task::PackerTask & task ) const
{
	utility::vector1< core::Size > resids_to_design( residues_to_design( pose ) );

	// for each bad residue, we want to prevent the original amino acid
	utility::vector1<bool> restrict_to_aa( core::chemical::num_canonical_aas, true );
	if ( prevent_native_aa_ ) {
		for ( core::Size i=1; i<=resids_to_design.size(); ++i ) {
			char aa;
			if ( cached_pose() ) {
				TR << "Using cached pose amino acid for pos " << resids_to_design[i] << " - " << cached_pose()->residue( resids_to_design[i] ).name1() << std::endl;
				aa = cached_pose()->residue( resids_to_design[i] ).name1();
			} else {
			  aa = pose.residue( resids_to_design[i] ).name1();
			}
			// retrieve restriction map
			restrict_to_aa[core::chemical::aa_from_oneletter_code( aa )] = false;
			task.nonconst_residue_task( resids_to_design[i] ).restrict_absent_canonical_aas( restrict_to_aa );
			restrict_to_aa[core::chemical::aa_from_oneletter_code( aa )] = true;
		}
	}

	std::string wanted_ss( blueprint_ss_ );
	if ( wanted_ss  == "" ) {
		core::pose::Pose posecopy( pose );
		protocols::moves::DsspMover dssp;
		dssp.apply( posecopy );
		std::string dssp_ss( posecopy.secstruct() );
		// look for ligands and don't include them in the SS
		for ( core::Size i=0; i<dssp_ss.size(); ++i ) {
				if ( pose.residue(i+1).is_protein() ) {
						wanted_ss += dssp_ss[i];
				}
		}
	}

	// if we are told to, scan for point mutants which would have a bad effect on psipred
	if ( prevent_bad_point_mutants_ ) {
		core::Size const baseline( psipred_interface_->run_psipred( pose, wanted_ss ).nres );
		for ( core::Size i=1; i<=resids_to_design.size(); ++i ) {
			// mutate the residue to all allowable positions and check psipred
			core::pack::task::ResidueLevelTask const & res_task( task.nonconst_residue_task( resids_to_design[i] ) );
			utility::vector1< char > types;
			for ( core::pack::task::ResidueLevelTask::ResidueTypeCOPListConstIter it=res_task.allowed_residue_types_begin(); it != res_task.allowed_residue_types_end(); ++it ) {
				types.push_back( (*it)->name1() );
			}
			for ( core::Size j=1; j<=types.size(); ++j ) {
				char const aa( types[j] );
				core::pose::Pose posecopy( pose );
				//TR << "after copy; residues to design=" << resids_to_design[i] << ", name=" << aa << std::endl;
				protocols::simple_moves::MutateResidue mutate( resids_to_design[i], aa );
				//TR << "after mut" << std::endl;
				mutate.apply( posecopy );
				//TR << "after apply" << std::endl;
				core::Size const newvalue( psipred_interface_->run_psipred( posecopy, wanted_ss ).nres );
				// the 0.1/total_residue is just to offset the value by a small amount to account for floating point precision
				// I didn't want to hard code a number because the precision of the psipred result is related to the size of the
				// pose
				//TR << "Mutating residue " << resids_to_design[i] << " to " << aa << " :  value=" << newvalue << "; baseline=" << baseline << std::endl;
				TR << pose.residue( i ).name1() << i << posecopy.residue( i ).name1() << " : " << newvalue << "(baseline=" << baseline << " ) -- " << "AA=" << aa;
				if ( newvalue < baseline ) {
					// failed mutation; prevent it
					TR << " Restricted." << std::endl;
					//TR << "AFTER char assign : " << aa << std::endl;
					restrict_to_aa[core::chemical::aa_from_oneletter_code( aa )] = false;
					core::pack::task::ResidueLevelTask & res_task_nonconst( task.nonconst_residue_task( resids_to_design[i] ) );
					//TR << "Successfully made nonconst res task." << std::endl;
					res_task_nonconst.restrict_absent_canonical_aas( restrict_to_aa );
					//TR << "successfully restricted aas" << std::endl;
					restrict_to_aa[core::chemical::aa_from_oneletter_code( aa )] = true;
				} else {
						TR << " Allowed." << std::endl;
				}
			}
		}
	}

  // now we can just apply a DesignAround operation using the residues that don't match
  protocols::toolbox::task_operations::DesignAroundOperation design_around;
  design_around.design_shell( region_shell() );
	if ( repack_non_selected() ) {
		design_around.repack_shell( 1000.0 );
	} else {
		design_around.repack_shell( region_shell() );
	}
  TR << "Residues to design are: ";
  for ( core::Size i=1; i<=resids_to_design.size(); i++) {
    TR << resids_to_design[i] << " ";
    design_around.include_residue( resids_to_design[i] );
  }
  TR << std::endl;
  design_around.apply( pose, task );
}

void
DesignBySecondaryStructureOperation::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) {
	HighestEnergyRegionOperation::parse_tag( tag, datamap );
	initialize_blueprint_ss( tag->getOption< std::string >( "blueprint", "" ) );
	std::string const cmd( tag->getOption< std::string >( "cmd", "" ) );
	if ( cmd == "" ) {
		utility_exit_with_message( "cmd must be set in DesignBySecondaryStructure." );
	}
	// now that we have a command, we can create the psipred interface object
	psipred_interface_ = new denovo_design::filters::PsiPredInterface( cmd );
	prevent_bad_point_mutants_ = tag->getOption< bool >( "prevent_bad_point_mutations", prevent_bad_point_mutants_ );
}

void
DesignBySecondaryStructureOperation::initialize_blueprint_ss( std::string const blueprint_file ) {
	if ( blueprint_file == "" ) {
		//utility_exit_with_message( "A blueprint filename must be provided to DesignBySecondaryStructure." );
		return;
	}
	protocols::jd2::parser::BluePrintOP bp = new protocols::jd2::parser::BluePrint( blueprint_file );
	if ( ! bp ) {
		utility_exit_with_message( "Error initializing the blueprint." );
	}
	blueprint_ss_ = bp->secstruct();
}

core::pack::task::operation::TaskOperationOP
DesignBySecondaryStructureOperationCreator::create_task_operation() const {
	return new DesignBySecondaryStructureOperation;
}

std::string
DesignBySecondaryStructureOperationCreator::keyname() const {
	return "DesignBySecondaryStructure";
}

//namespaces
}
}
}
