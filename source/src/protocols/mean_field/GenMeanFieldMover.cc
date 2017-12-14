// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Aliza Rubenstein (aliza.rubenstein@gmail.com)

// Unit headers
#include <protocols/mean_field/GenMeanFieldMover.hh>
#include <protocols/mean_field/GenMeanFieldMoverCreator.hh>

// Package headers
#include <protocols/mean_field/MeanField.hh>
#include <protocols/mean_field/DesignMeanField.hh>
#include <protocols/mean_field/FlexBBDesignMeanField.hh>
#include <protocols/mean_field/RotMatrix.hh>
#include <protocols/mean_field/MeanFieldFactory.hh>
#include <protocols/mean_field/AAProb.hh>
#include <protocols/mean_field/RotProb.hh>
#include <protocols/mean_field/jagged_array.functions.hh>
#include <protocols/mean_field/AAMatrix.hh>
#include <protocols/mean_field/ResHashMap.hh>

// Project Headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// util
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <cmath>
#include <utility/file/file_sys_util.hh> // file_exists
#include <utility/file/FileName.hh>
#include <utility/vector0.hh>
#include <utility/io/izstream.hh>

// option keys
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mean_field.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/prof.hh>


namespace protocols {
namespace mean_field {

static basic::Tracer TR( "protocols.mean_field.GenMeanFieldMover" );

using namespace core;
using namespace utility;

///@brief default ctor
GenMeanFieldMover::GenMeanFieldMover() :
	parent( "GenMeanFieldMover" ),
	task_factory_( /* NULL */ ),
	scorefxn_( /* NULL */ ),
	lambda_memory_( default_value_for_lambda_memory() ),
	tolerance_( default_value_for_tolerance() ),
	temperature_( default_value_for_temperature() ),
	init_option_( default_value_for_init_option() ),
	threshold_( default_value_for_threshold() ),
	unbound_( default_value_for_unbound() ),
	aa_matrix_()
{}

/// @details calls other methods to:
/// @details 1. read input AAMatrix and pdbs
/// @details 2. prepare task and poses
/// @details 3. call the mean-field algorithm
/// @details 4. report probabilities of rotamers, specificity profile (AAMatrix), BB Boltzmann probabilities, and distances
void
GenMeanFieldMover::apply( pose::Pose & pose )
{
	basic::prof_reset();
	PROF_START(basic::MEAN_FIELD);
	read_aa_matrix();

	read_input_pdbs();

	prepare_task_poses( pose );
	calc_mean_field();

	report_rot_prob( );
	if ( task_->design_any() ) {

		report_aa_matrix( );
	}

	PROF_STOP(basic::MEAN_FIELD);

	basic::prof_show();
}

/// @details reports AAMatrix to output stream and to file (if dump_transfac is set)
/// @details if gold standard AAMatrix has been set, outputs distances between AAMatrix and gold standard
/// @details if FlexBB, reports bb_boltz_probs, AAMatrixs and list of pdbs to output/files
void
GenMeanFieldMover::report_aa_matrix() const
{
	//only report if task is designable
	if ( task_->design_any() ) {

		//if MeanField has a FlexBB, report per-bb information (bb_boltz_probs and AAMatrixs)
		if ( dynamic_cast< FlexBBDesignMeanField * > ( mean_field_.get() ) ) {
			FlexBBDesignMeanFieldOP dmf =
				pointer::static_pointer_cast< protocols::mean_field::FlexBBDesignMeanField > ( mean_field_ );

			//   vector1< jagged_array< Real > > ems = dmf->energy_matrices();
			//   vector1< RotMatrix > rms = dmf->rot_matrices();

			vector1 < AAMatrix > sps = dmf->aa_matrices();
			jagged_array < Real > bb_boltz_probs = dmf->bb_boltz_probs();
			vector1 < jagged_array < Real > > bb_boltz_probs_per_aa = dmf->bb_boltz_probs_per_aa();
			//TODO: AR: change back to TR.Debug
			// output to TR.Debug if TR.Debug is visible
			if ( TR.Debug.visible() ) {
				TR.Debug << "Reporting BB Boltz Probs" << std::endl;

				for ( Size bb = 1; bb <= sps.size(); ++bb ) {
					TR.Debug << "Backbone " << bb << " ------------" << std::endl;
					TR.Debug << "Per position: " << " -------------" << std::endl;
					for ( Size pos = 1; pos <= bb_boltz_probs.size(); ++pos ) {
						TR.Debug << bb_boltz_probs[pos][bb] << "\t";
					}

					TR.Debug << std::endl;

					TR.Debug << "Per amino acid at each position: " << " --------" << std::endl;
					bb_boltz_probs_per_aa[bb].show( TR.Info );

					TR.Debug << "Reporting Specificity Profile (AAMatrix) " << bb << " out of " << sps.size() << std::endl;
					sps[bb].show( TR.Debug );
					if ( !aa_matrix_.empty() ) {
						report_dist( sps[bb] );
					}

				}

			}

			//if dump_transfac is active output bb_boltz_probs, AAMatrixs, and a list of the pdbs that converged to files
			if ( basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ].active() ) {
				//output AAMatrices
				for ( Size bb = 1; bb <= sps.size(); ++bb ) {
					std::string filename = basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ];
					filename += "_";
					filename += boost::lexical_cast<std::string>(bb);
					filename += ".transfac";
					sps[bb].dump_transfac( filename );
				}

				std::string boltz_prob_filename = basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ];
				boltz_prob_filename += "_boltz.txt";
				std::ofstream output( boltz_prob_filename.c_str() );

				std::string list_filename = basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ];
				list_filename += "_list.txt";
				std::ofstream list_output( list_filename.c_str() );

				for ( Size ii = 1 ; ii <= bb_boltz_probs.max_size_col() ; ++ii ) {
					output << poses_[ii]->pdb_info()->name() << "\t";
					list_output << poses_[ii]->pdb_info()->name() << std::endl;
					for ( Size jj = 1 ; jj <= bb_boltz_probs.size() ; ++jj ) {
						if ( ii <= bb_boltz_probs[jj].size() ) {
							output << bb_boltz_probs[jj][ii] << "\t";
						} else {
							output << "\t\t";
						}
					}
					output << std::endl;
				}
			}

		}

		//show predicted AAMatrix and the distances between it and the gold standard AAMatrix
		pred_aa_matrix_.show( TR.Info);

		if ( ! aa_matrix_.empty() ) {
			report_dist( pred_aa_matrix_ );
		}

		//if dump_transfac is on, dump predicted AAMatrix
		if ( basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ].active() ) {
			std::string filename = basic::options::option[ basic::options::OptionKeys::mean_field::dump_transfac ];
			filename += ".transfac";
			pred_aa_matrix_.dump_transfac( filename );

		}

	}
}

void
GenMeanFieldMover::report_rot_prob( ) const
{
	mean_field_->rot_matrix()->show( TR.Info );
}

/// @details only runs if Design is set for the MeanField and if there is a "gold standard" AAMatrix to compare to
void
GenMeanFieldMover::report_dist( AAMatrix const & am ) const
{
	if ( task_->design_any() && ! aa_matrix_.empty() ) {

		vector1 < Real > cosine_dist = aa_matrix_.cosine_distance( am );
		vector1 < Real > frob_dist = aa_matrix_.frob_distance( am );
		vector1 < Real > ave_abs_diff = aa_matrix_.ave_abs_diff( am );

		vector1 < vector1 < Real > > vecs;
		vecs.push_back( cosine_dist );
		vecs.push_back( frob_dist );
		vecs.push_back( ave_abs_diff );

		for ( vector1 < vector1 < Real > >::const_iterator vec = vecs.begin(); vec != vecs.end();
				++vec ) {
			vector1 < Real > v = *vec;
			for ( vector1 < Real >::const_iterator dist = v.begin(); dist != v.end();
					++dist ) {
				TR.Info << *dist << "\t";
			}
			TR.Info << std::endl;
		}
	}
}

/// @details reads the AAMatrix from aa_matrix
/// @details Error if file cannot be opened
void
GenMeanFieldMover::read_aa_matrix()
{

	std::istringstream am;
	if ( basic::options::option[ basic::options::OptionKeys::mean_field::spec_profile ].user() ) {
		std::string filename = basic::options::option[ basic::options::OptionKeys::mean_field::spec_profile];
		std::string am_str;
		utility::io::izstream file( filename );

		if ( ! file ) {
			TR.Error << "File: " << filename << " not found!" << std::endl;
		} else {
			slurp( file, am_str );
			am.str( am_str );
			aa_matrix_.build_aa_matrix( am );
			aa_matrix_.show( TR.Info );
		}

	}
}

/// @details read input pdbs into pdb_list_ from bb_list or s if bb_list is empty
/// @details Error if bb_list specified but there are no elements in the list
void
GenMeanFieldMover::read_input_pdbs()
{

	using Filenames = vector1<file::FileName>;


	if ( basic::options::option[ basic::options::OptionKeys::mean_field::bb_list ].active() ) {
		Filenames listnames( basic::options::option[ basic::options::OptionKeys::mean_field::bb_list ]().vector() );
		for ( Filenames::const_iterator filename( listnames.begin() );
				filename != listnames.end(); ++filename ) {

			utility::io::izstream list( (*filename).name().c_str() );

			while ( list ) {
				std::string pdbname;
				list >> pdbname;
				if ( pdbname != "" ) pdb_list_.push_back( pdbname );
			}
		}
		if ( pdb_list_.empty() ) {
			std::stringstream error_message;
			error_message
				<< "bb_list specified but no files found." << std::endl;
			throw CREATE_EXCEPTION(excn::Exception,  error_message.str() );
		}

	} else if ( basic::options::option[ basic::options::OptionKeys::in::file::s ].active() ) {
		pdb_list_ = basic::options::option[ basic::options::OptionKeys::in::file::s ]().vector();

	}

}

/// @details initalize poses_ from list of pdbs
/// @details creates a task from the task_factory
/// @details unbinds the poses by the last jump if unbound is set to true
/// @remarks Error if no residues are set as packable by the task_
void
GenMeanFieldMover::prepare_task_poses( pose::Pose const & pose )
{
	poses_.reserve( pdb_list_.size() );

	poses_ = import_pose::poseOPs_from_files( pdb_list_ );

	task_factory()->push_back( pack::task::operation::TaskOperationCOP( new pack::task::operation::InitializeFromCommandline ) );

	//uses input pose (in -s option) to create task
	task_ = task_factory()->create_task_and_apply_taskoperations( pose );
	task_->set_bump_check( true );

	if ( task_->num_to_be_packed() == 0 ) {
		std::stringstream error_message;
		error_message
			<< "Task does not pack anything." << std::endl;
		throw CREATE_EXCEPTION(excn::Exception,  error_message.str() );
	}


	//TODO AR: added for release which builds task on -s and has all new innovations (1task) in the future this should be implemented as an option or tag
	//tasks_.push_back( task_ );

	//unbind poses by last jump if unbound is on
	for ( pose::PoseOPs::const_iterator pose = poses_.begin();
			pose != poses_.end(); ++pose ) {
		//commented out for 1task
		task_ = task_factory()->create_task_and_apply_taskoperations( **pose );
		task_->set_bump_check( true );
		tasks_.push_back( task_ );

		if ( unbound() ) {
			protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( **pose, (*pose)->num_jump() ) ); // JB 120420
			translate->step_size( 1000.0 );
			translate->apply( **pose );
		}

	}
}

/// @details uses MeanFieldFactory to create a MeanField class, calls MeanField::process
/// @details since AAMatrix is used in DesignMeanField and FlexBBDesignMeanField but not in MeanField or FlexBBMeanField
/// @details must perform a cast of the mean_field_ to retrieve the predicted AAMatrix
void
GenMeanFieldMover::calc_mean_field()
{
	mean_field_ = protocols::mean_field::MeanFieldFactory::create_mean_field(
		init_option(), poses_, tasks_, scorefxn(), lambda_memory(), tolerance(), temperature(), threshold() );
	mean_field_->process();
	if ( task_->design_any() ) {
		if ( dynamic_cast< DesignMeanField * > ( mean_field_.get() ) ) {

			DesignMeanFieldOP dmf =
				pointer::static_pointer_cast< protocols::mean_field::DesignMeanField > ( mean_field_ );

			pred_aa_matrix_ = * dmf->aa_matrix();
		} else {
			assert( dynamic_cast< FlexBBDesignMeanField * > ( mean_field_.get() ) );

			FlexBBDesignMeanFieldOP fbbdmf =
				pointer::static_pointer_cast< protocols::mean_field::FlexBBDesignMeanField > ( mean_field_ );

			pred_aa_matrix_ = * fbbdmf->aa_matrix();
		}

	}

	//added 6/23/15
	//this is a method that subtracts bound-unbound energies and may be used again the future.  leaving it in for now.
	// if ( dynamic_cast< FlexBBMeanField * > ( mean_field_.get() ) )
	// {
	//
	//  FlexBBMeanFieldOP mf =
	//   pointer::static_pointer_cast< protocols::mean_field::FlexBBMeanField > ( mean_field_ );
	//
	//  protocols::mean_field::RotMatrix bound_rot_matrix = * mf->rot_matrix();
	//  protocols::mean_field::jagged_array< Real > bound_energy_matrix = mf->energy_matrix();
	//  for ( pose::PoseOPs::const_iterator pose = poses_.begin();
	//     pose != poses_.end(); ++pose )
	//  {
	//   protocols::rigid::RigidBodyTransMoverOP translate( new protocols::rigid::RigidBodyTransMover( **pose, (*pose)->num_jump() ) ); // JB 120420
	//   translate->step_size( 1000.0 );
	//   translate->apply( **pose );
	//  }
	//  mean_field_ = protocols::mean_field::MeanFieldFactory::create_mean_field(
	//    init_option(), poses_, tasks_, scorefxn(), lambda_memory(), tolerance(), temperature(), threshold() );
	//  mean_field_->process();
	//
	//  protocols::mean_field::RotMatrix unbound_rot_matrix = * mf->rot_matrix();
	//  protocols::mean_field::jagged_array< Real > unbound_energy_matrix = mf->energy_matrix();
	//
	//  protocols::mean_field::RotMatrix subt_rot_matrix;
	//  protocols::mean_field::jagged_array< Real > subt_energy_matrix( bound_rot_matrix.size(), utility::vector1< Real >() );
	//
	//
	//  for ( Size pos = 1; pos <= unbound_rot_matrix.size(); ++pos )
	//  {
	//
	//   subt_rot_matrix.push_back( utility::vector1< RotProb >() );
	//
	//      //insert unbound rotamers into ResHashMap, update unbound_rot_matrix rot_ind to correct numbers
	//   ResHashMap unbound_rothash;
	//
	//   for ( Size unbound_rot = 1; unbound_rot <= unbound_rot_matrix[ pos ].size(); unbound_rot++ )
	//   {
	//    SSize rot_ind = unbound_rothash.attempt_insert( unbound_rot_matrix[ pos ][ unbound_rot ].res() );
	//
	//    unbound_rot_matrix[ pos ][ unbound_rot ].rot_ind( rot_ind );
	//
	//   }
	//
	//   //loop through bound rotamers
	//   for ( Size bound_rot = 1; bound_rot <= bound_rot_matrix[ pos ].size(); bound_rot++ )
	//   {
	//    //get rot_ind of the bound_rotamer (to compare to unbound_rotamers)
	//    SSize rot_ind = unbound_rothash.get_rot_ind( unbound_rot_matrix[ pos ][ bound_rot ].res() );
	//
	//    //if bound_rotamer was found in hashmap (if it's not we don't care about it)
	//    if ( rot_ind != SSize( -1 ) )
	//    {
	//     //update rot_ind of bound_rot_matrix
	//     bound_rot_matrix[ pos ][ bound_rot ].rot_ind( rot_ind );
	//     //push back bound_rot into subt_rot_matrix and energy of bound_rot into subt_energy_matrix
	//     subt_rot_matrix[ pos ].push_back( bound_rot_matrix[ pos ][ bound_rot ] );
	//     subt_energy_matrix[ pos ].push_back( bound_energy_matrix[ pos ][ bound_rot ] );
	//
	//     //loop through some of unbound rotamers to look for equivalent unbound rotamer (== rot_ind)
	//     //begin at rot_ind because index will always be >= rot_ind (I think...)
	//     for ( Size unbound_rot = rot_ind; unbound_rot <= unbound_rot_matrix[ pos ].size(); unbound_rot++ )
	//     {
	//      //breaks once it finds the correct rotamer
	//      if ( unbound_rot_matrix[ pos ][ unbound_rot ].rot_ind() == Size( rot_ind ) )
	//      {
	//       subt_energy_matrix[ pos ][ bound_rot ] -= unbound_energy_matrix[ pos ][ unbound_rot ];
	//       break;
	//      }
	//     }
	//    }
	//
	//   }
	//
	//  }
	//
	//  //doesn't use jagged_array get_totals_columns because this is faster
	//  utility::vector1 < Real > totals( subt_energy_matrix.size(), core::Real( 0.0 ) );
	//
	//  for ( Size pos = 1 ; pos <= subt_energy_matrix.size() ; ++pos )
	//  {
	//   for ( Size rot = 1 ; rot <= subt_energy_matrix[ pos ].size() ; ++rot )
	//   {
	//    Real energy = subt_energy_matrix[ pos ][ rot ] < temperature() ? subt_energy_matrix[ pos ][ rot ] : threshold();
	//
	//    subt_rot_matrix[ pos ][ rot ].probability( exp( ( energy * Real( -1.0 ) ) / temperature() ) );
	//
	//    totals[ pos ] += subt_rot_matrix[ pos ][ rot ].probability();
	//   }
	//  }
	//
	//  subt_rot_matrix /= totals;
	//
	//  subt_rot_matrix.is_designed( bound_rot_matrix.is_designed() );
	//
	//  pred_aa_matrix_.build_aa_matrix( subt_rot_matrix, subt_energy_matrix, temperature() );
	//
	//  //TODO: AR: if I end up using this method - turn the EM into a RM
	//  //use energy_matrix + rotamer_matrix to sum over energies for AA's
	//
	//  //save resulting aamatrix as pred_aa_matrix_
	// }
}

void
GenMeanFieldMover::parse_my_tag( tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	pose::Pose const & )
{
	using namespace tag;

	TR << "GenMeanFieldMover"<<std::endl;
	task_factory( protocols::rosetta_scripts::parse_task_operations( tag, data ) );
	scorefxn( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	lambda_memory( tag->getOption< Real >( "lambda_memory", 0.5 ) );
	tolerance( tag->getOption< Real >( "tolerance", 0.0001 ) );
	temperature( tag->getOption< Real >( "temperature", 0.6 ) );
	init_option( tag->getOption< Size >( "init_option", 1 ) );
	threshold( tag->getOption< Real >( "threshold", 10.0 ) );
	unbound( tag->getOption< bool >( "unbound", false ) );

	TR<<"with options lambda (memory): "<<lambda_memory()<<" tolerance: "<<tolerance()<<" temperature: "<<temperature()<<" threshold: "<<threshold()<<" init option: "<<init_option()<<" unbound: "<<unbound()<<std::endl;
}

void GenMeanFieldMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "lambda_memory", xs_decimal, "Used in updating the original mean-field matrix to decrease oscillation.", utility::to_string( default_value_for_lambda_memory() ) )
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xs_decimal, "Tolerance used to determine if convergence is achieved.", utility::to_string( default_value_for_tolerance() ) )
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xs_decimal, "kT used in Boltzmann weighting.", utility::to_string( default_value_for_temperature() ) )
		+ XMLSchemaAttribute::attribute_w_default( "init_option", xs_decimal, "Determines initial values of mean-field matrix. A value of 1 forces initial values of 1/nrot and a value of 0 forces initial values of 0.  This has not been tested with a value of 0.", utility::to_string( default_value_for_init_option() ) )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xs_decimal, "Threshold for truncation of energy matrix values.", utility::to_string( default_value_for_threshold() ) )
		+ XMLSchemaAttribute::attribute_w_default( "unbound", xsct_rosetta_bool, "Unbind protein from substrate before running?", utility::to_string( default_value_for_unbound() ) );

	rosetta_scripts::attributes_for_parse_score_function( attlist );
	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Run mean-field on pose", attlist );

}

void GenMeanFieldMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	GenMeanFieldMover::provide_xml_schema( xsd );
}

protocols::moves::MoverOP
GenMeanFieldMover::fresh_instance() const{
	return protocols::moves::MoverOP( new GenMeanFieldMover );
}

GenMeanFieldMover::~GenMeanFieldMover()= default;

protocols::moves::MoverOP
GenMeanFieldMover::clone() const{
	return protocols::moves::MoverOP( new GenMeanFieldMover( *this ) );
}

protocols::moves::MoverOP
GenMeanFieldMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new GenMeanFieldMover );
}

std::string
GenMeanFieldMoverCreator::keyname() const
{
	return GenMeanFieldMover::mover_name();
}

std::string GenMeanFieldMover::get_name() const {
	return mover_name();
}

std::string
GenMeanFieldMover::mover_name()
{
	return "GenMeanFieldMover";
}

} // mean_field
} // protocols
