// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/SecondaryMatchProtocol.cc
/// @author Florian Richter

#include <protocols/enzdes/SecondaryMatchProtocol.hh>
#include <protocols/enzdes/EnzdesBaseProtocol.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzCstTemplateRes.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintParameters.hh>
#include <protocols/toolbox/match_enzdes_util/util_functions.hh>
#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/enzdes/EnzdesTaskOperations.hh>
#include <protocols/enzdes/EnzdesTaskOperations.fwd.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/io/pdb/file_data.hh> //reading remarks

#include <basic/options/option.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh> //reading remarks
//#include <core/pose/PDBPoseMap.hh> //for PDB-info-to-resid functionality
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/Func.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>

// AUTO-REMOVED #include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>

#include <basic/Tracer.hh>
#include <utility/string_util.hh>


// option key includes

#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.hh>


namespace protocols{
namespace enzdes{

static basic::Tracer tr("protocols.enzdes.SecondaryMatchProtocol");

SecondaryMatchProtocol::SecondaryMatchProtocol() :
	EnzdesBaseProtocol(),
	found_res_compatibility_determined_( false ),
	cut1_(0),
	cut2_(0),
	cut3_( basic::options::option[basic::options::OptionKeys::enzdes::cut3] ),
	cut4_( basic::options::option[basic::options::OptionKeys::enzdes::cut4] )
{

	reduced_scofx_ =  core::scoring::ScoreFunctionFactory::create_score_function( "enzdes_polyA_min" );
	match_params_.clear();
	found_res_compatibility_.clear();
}

SecondaryMatchProtocol::~SecondaryMatchProtocol(){}

void
SecondaryMatchProtocol::apply(
	core::pose::Pose & start_pose
){


	using namespace basic::datacache;
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;

	//in case we want to swich ligands, let's do it now
	if( basic::options::option[basic::options::OptionKeys::enzdes::change_lig].user() ){
		bool lig_switch = exchange_ligands_in_pose( start_pose, true, reduced_scofx_ );
		if( !lig_switch ){
			tr << "Warning: could not perform the requested ligand switch, aborting protocol... " << std::endl;
			return;
		}
	}


	//set up constraints (read cstfile, do mapping, etc, then add to pose)
	if( basic::options::option[basic::options::OptionKeys::enzdes::cstfile].user() ){
		enable_constraint_scoreterms();
		setup_enzdes_constraints( start_pose, true );
	}
	//if there is no cstfile but a change of lig requested, we need to dump the pose now before returning
	else if( basic::options::option[basic::options::OptionKeys::enzdes::change_lig].user() ){

		std::string outtag = start_pose.data().get_ptr< CacheableString >( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG )->str() + "_" + "lx.pdb";
		start_pose.dump_pdb( outtag );

		return;
	}

	else utility_exit_with_message("This protocol either needs a cstfile or change_lig file.");

	//TO DO: check whether the newly placed ligand clashes with the catalytic residues

	find_all_allowed_positions( start_pose );
	do_matching( start_pose );

}

std::string
SecondaryMatchProtocol::get_name() const {
	return "SecondaryMatchProtocol";
}

bool
SecondaryMatchProtocol::do_matching(
	core::pose::Pose & start_pose
){

	using namespace basic::datacache;
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;


	utility::vector1< core::pose::PoseOP > poses_to_process;
	poses_to_process.clear();

	//this might be a little iffy, come back to this if it causes problems
	poses_to_process.push_back( & start_pose );

	core::pose::PoseCOP ref_pose( new core::pose::Pose( start_pose ) );
	//core::pose::Pose ref_pose = start_pose;

	//clear_catalytic_res(); //make sure this gets done first thing

	found_resis_.clear();
	match_params_.clear();

	toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( enzutil::get_enzcst_io( start_pose ) );
	if( !cstio ){
		tr << "Weird. No constraints seem to have been setup in the pose." << std::endl;
		return false;
	}

	utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP > tmp_params = cstio->enz_cst_params_missing_in_pose( start_pose );
	//utility::vector1< EnzConstraintParametersCOP > match_params;

	//first we'll clone the missing params
	for( utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP >::const_iterator tmp_it = tmp_params.begin();
			 tmp_it != tmp_params.end(); ++tmp_it) {
		match_params_.push_back( new toolbox::match_enzdes_util::EnzConstraintParameters( **tmp_it) );
	}

	core::Size num_res_to_match = match_params_.size();

	if( num_res_to_match == 0 ){
		tr.Info << "The pose already has all the required constraints present, aborting protocol.\n";
		return true;
	}
	else tr.Info << "There are " << num_res_to_match << " interactions that need to be found.\n" << std::endl;

	for( utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP >::const_iterator param_it = match_params_.begin();
			 param_it != match_params_.end(); ++param_it) {

		toolbox::match_enzdes_util::EnzCstTemplateResCOP missing_template = (*param_it)->get_missing_template_res( start_pose );
		toolbox::match_enzdes_util::EnzCstTemplateResCOP present_template = (*param_it)->get_missing_template_other_res( start_pose );

		for( utility::vector1< core::pose::PoseOP >::iterator pose_it = poses_to_process.begin();
				 pose_it != poses_to_process.end(); ++pose_it){

			add_enz_cst_interaction_to_pose( **pose_it, *param_it, missing_template, present_template, cstio);

		} // loop over all poses to process

		if( found_resis_.size() == 0 ) {
			tr.Info << "Bummer :( Could not find the desired interaction... " << std::endl;

			//in case a switch of ligands was requested, we dump the pose(s) now
			if( basic::options::option[basic::options::OptionKeys::enzdes::change_lig].user() ){

				std::string outtag = poses_to_process[1]->data().get_ptr< CacheableString >( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG )->str() + "_" + "lx.pdb";
				poses_to_process[1]->dump_pdb( outtag );
			}

			return false;
		}

	} //loop over all params to process


	return generate_and_dump_pose_found_residues_combinations( ref_pose );


} //do matching function


void
SecondaryMatchProtocol::add_enz_cst_interaction_to_pose(
	core::pose::Pose & pose,
	toolbox::match_enzdes_util::EnzConstraintParametersCOP params,
	toolbox::match_enzdes_util::EnzCstTemplateResCOP missing_template,
	toolbox::match_enzdes_util::EnzCstTemplateResCOP present_template,
	toolbox::match_enzdes_util::EnzConstraintIOCOP cstio)
{
	using namespace toolbox::match_enzdes_util;

	utility::vector1< core::conformation::ResidueOP > found_res_this_param;

	//1. we need to find out which residues in the pose we are looking at
	utility::vector1< core::Size > target_residues;

	tr.Info << "Trying to add interaction of pose residue(s) ";
	for( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator
			pos_it = get_enzdes_observer( pose )->cst_cache()->param_cache( params->cst_block() )->template_res_cache( present_template->param_index() )->seqpos_map_begin(),
			pos_end = get_enzdes_observer( pose )->cst_cache()->param_cache( params->cst_block() )->template_res_cache( present_template->param_index() )->seqpos_map_end();
			pos_it != pos_end; ++pos_it ){
		target_residues.push_back( pos_it->first );
		tr.Info << pose.residue_type( pos_it->first ).name() << " " << pos_it->first << ", ";
	}

	if( target_residues.size() != 1 ){
		utility_exit_with_message("Error: there are more (or less) than 1 target residue. "
				"This protocol isn't setup to handle dealing with multiple target residues yet( and it shouldn't necessarily be.");
	}
	core::Size target_residue = target_residues[1];

	tr.Info << "with a newly placed residue of type(s) ";

	utility::vector1< std::string > trial_restypes = missing_template->allowed_res_types();

	if( missing_template->is_backbone() ){
		trial_restypes.clear();
		trial_restypes.push_back( "GLY" );
	}

	for( utility::vector1< std::string >::const_iterator resi_it = trial_restypes.begin(); resi_it != trial_restypes.end(); ++resi_it ){
		tr.Info << *resi_it << ", ";
	}
	tr.Info << std::endl;


	//create poly A pose of neighbors
	protocols::toolbox::pose_manipulation::construct_poly_ala_pose( pose, trial_positions_, true, true, true );

	core::conformation::Residue ala_res( restype_set()->name_map("ALA"), true );

	for( utility::vector1< std::string >::const_iterator resi_it = trial_restypes.begin(); resi_it != trial_restypes.end(); ++resi_it ){

		core::conformation::Residue trial_res( restype_set()->name_map( *resi_it ), true );

		tr << "starting search for restype " << *resi_it << "... ";

		for( utility::vector1< core::Size >::const_iterator pos_try_it = trial_positions_.begin(); pos_try_it != trial_positions_.end(); ++pos_try_it){

			tr << "searching position " << *pos_try_it << "... ";

			if( ! restype_possible_at_position( pose, & trial_res.type(), & pose.residue( target_residue ), *pos_try_it ) ) continue;

			utility::vector1< std::string > current_variants;
			bool match = variants_match( pose.residue_type( *pos_try_it ), trial_res.type() );

			if( ! match ){
				current_variants = pose.residue_type( *pos_try_it ).properties().get_list_of_variants();
			}

			pose.replace_residue( *pos_try_it, trial_res, true);

			if( ! match ){
				for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
					core::pose::add_variant_type_to_pose_residue(
							pose,
							core::chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ),
							*pos_try_it );
				}
			}

			//gotta mark the new residue in the constraint object
			cstio->set_position_for_missing_res_in_parameter_block( pose, params->cst_block(), *pos_try_it );

			cstio()->add_constraints_to_pose_for_block_without_clearing_and_header_processing( pose, reduced_scofx_, params->cst_block() );


			core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
			//task->initialize_from_command_line();
			for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i) {
				if( ( i == *pos_try_it ) || ( i == target_residue ) ){
					task->nonconst_residue_task( i ).restrict_to_repacking();
					task->nonconst_residue_task( i ).initialize_from_command_line();
				}
				else task->nonconst_residue_task( i ).prevent_repacking();
			}

			protocols::simple_moves::PackRotamersMoverOP trial_packer = new protocols::simple_moves::PackRotamersMover(reduced_scofx_, task);
			trial_packer->apply( pose );

			(*reduced_scofx_)( pose );

			//make a bump check and constraint scoreterm check
			//bump check necessary in case each of the rotamers clashed
			core::Real try_targ_clash(0.0);

			core::scoring::EnergyEdge const * eedge = pose.energies().energy_graph().find_energy_edge( *pos_try_it, target_residue );

			if( eedge != 0 ){
				// apl -- removing TwoBodyEnergyMap + making this call a little more efficient.
				try_targ_clash = eedge->dot( reduced_scofx_->weights() );
			}

			core::Real cst_sum = enzutil::sum_constraint_scoreterms( pose, *pos_try_it );


			//little hacky, not the real bbclash energy, but in this case probably a better approximation because it also includes the clash with ALAs, etc
			core::Real try_bb_clash = pose.energies().residue_total_energies( *pos_try_it )[ core::scoring::fa_rep ];

		tr.Info << "For restype " << pose.residue_type( *pos_try_it ).name() << " at pos " << *pos_try_it << " the cst_sum is " << cst_sum << " and the clashE is " << try_targ_clash << "." << std::endl;

			if( ( cst_sum < basic::options::option[basic::options::OptionKeys::enzdes::secmatch_Ecutoff] )
				&& ( try_targ_clash < basic::options::option[basic::options::OptionKeys::enzdes::sc_sc_bump_cutoff] )
				&& try_bb_clash < basic::options::option[basic::options::OptionKeys::enzdes::bb_bump_cutoff]){

				core::conformation::ResidueOP fres = new core::conformation::Residue( pose.residue( *pos_try_it ) );

				found_res_this_param.push_back( fres );

				tr.Info << "Success: " << *resi_it << " at position " << *pos_try_it << " satisfies constraints." << std::endl;
			}

			cstio->remove_constraints_from_pose_for_block( pose, params->cst_block(), true );

			cstio->remove_position_from_template_res_for_block( pose, *pos_try_it, params->cst_block() );

			cstio->clear_active_pose_constraints_for_block( pose, params->cst_block() );

			pose.replace_residue( *pos_try_it, ala_res, true);

			if( ! match ){
				for ( core::Size var = 1; var <= current_variants.size(); ++var ) {
					core::pose::add_variant_type_to_pose_residue(
							pose,
							core::chemical::ResidueProperties::get_variant_from_string( current_variants[ var ] ),
							*pos_try_it );
				}
			}

		} //iterator over all trial positions

	}// iterator over all trial residue types

	found_resis_.push_back( found_res_this_param );


} //add_enz_cst_interaction_to_pose


void
SecondaryMatchProtocol::find_all_allowed_positions(
	core::pose::Pose const & pose
)
{
	trial_positions_.clear();

	//still to be developed properly. for now we assume that we only allow attaching
	//residues at the design or repack positions as specified by the one ligand in the pose
	utility::vector1< bool > dummy( pose.total_residue() );
	utility::vector1< bool > allowed_res( pose.total_residue() );
	std::set< core::Size > interface_res;
	interface_res.insert( pose.fold_tree().downstream_jump_residue( pose.num_jump() ) );
	DetectProteinLigandInterfaceOP lig_prot_interface = new DetectProteinLigandInterface();
	lig_prot_interface->find_design_interface(pose, interface_res, cut1_, cut2_, cut3_, cut4_, allowed_res, dummy );

	//utility::vector1< core::Size > cat_res = catalytic_res();

	for(core::Size i = 1, i_end = pose.total_residue(); i <= i_end; ++i){

		if( allowed_res[ i ] == true ){


			if( pose.residue_type( i ).is_ligand() ) continue;

			if( ( pose.residue_type( i ).name3() == "GLY" )
				|| ( pose.residue_type( i ).name3() == "PRO" ) ) continue;

			//utility::vector1< Size >::iterator find_res = find( cat_res.begin(), cat_res.end(), i );

			//if( find_res == cat_res.end() ) trial_positions_.push_back( i );
			if( !is_catalytic_position( pose, i ) ) trial_positions_.push_back( i );

		}
	}

	tr << "trying at a maximum of " << trial_positions_.size() << " pose positions." << std::endl;

} //find_neighbors function



bool
SecondaryMatchProtocol::generate_and_dump_pose_found_residues_combinations( core::pose::PoseCOP ref_poseCOP )
{


	//utility::vector1< core::pose::PoseCOP > process_poses ;
	utility::vector1< PoseFoundResiduesCombinationOP > process_combos;

	//process_poses.push_back( ref_poseCOP );

	process_combos.push_back( new PoseFoundResiduesCombination( ref_poseCOP, this ) );


	for( utility::vector1< utility::vector1< core::conformation::ResidueOP > >::iterator param_it = found_resis_.begin();
			 param_it != found_resis_.end(); ++param_it )
		{

			//utility::vector1< core::pose::PoseCOP > temp_poses;
			utility::vector1< PoseFoundResiduesCombinationOP > temp_combos;

			for( utility::vector1< core::conformation::ResidueOP >::iterator res_it = param_it->begin();
					 res_it != param_it->end(); ++res_it )
				{

					for( utility::vector1< PoseFoundResiduesCombinationOP >::iterator pp_it = process_combos.begin();
					 pp_it != process_combos.end(); ++pp_it )
						{
							//core::pose::PoseOP success_pose = new core::pose::Pose(**pp_it);
							//success_pose->replace_residue( (*res_it)->seqpos(), **res_it, true);
							PoseFoundResiduesCombinationOP success_combo = new PoseFoundResiduesCombination( **pp_it );
							success_combo->add_residue( *res_it );
							temp_combos.push_back( success_combo );
						}

				}

			process_combos = temp_combos;

		}

	//now hopefully all the output poses are in the process combos array

	//before we can begin outputting, we first have to figure out which residues are compatible with which
	determine_found_residues_compatibility( ref_poseCOP );

	bool successful = false;

	for( utility::vector1< PoseFoundResiduesCombinationOP >::iterator pp_it = process_combos.begin();
					 pp_it != process_combos.end(); ++pp_it )
		{
			if( (*pp_it)->construct_and_dump_outpose( match_params_ ) ) successful = true;
		}

	return successful;

} //generate_output_names_for_found_poses



/// @brief rough check whether the two residues in question are close to each other
bool
SecondaryMatchProtocol::restype_possible_at_position(
	core::pose::Pose const & pose,
	core::chemical::ResidueTypeCOP restype,
	core::conformation::ResidueCOP target_residue,
	core::Size const trial_pos
)
{

	core::Size target_nbr_atom = target_residue->type().nbr_atom();
	//core::Size trial_pos_nbr_atom = pose.residue_type( trial_pos ).nbr_atom();

	core::PointPosition trial_xyz = pose.residue( trial_pos).xyz( "CB" );
	core::PointPosition targ_neighbor_xyz = target_residue->xyz( target_nbr_atom );

	core::Real dist = trial_xyz.distance( targ_neighbor_xyz );

	if( dist < ( 4.0 + restype->nbr_radius() + target_residue->type().nbr_radius() ) ) {
		return true;
	}
	else{
		tr.Info << "Can't place " << restype->name3() << " at position " << trial_pos << "because it's too far from target "<< target_residue->type().name3() << target_residue->seqpos() << "." << std::endl;
	}

	return false;

} //restype possible at position


/// @brief lookup function to determine whether to residues are compatible
core::Size
SecondaryMatchProtocol::residues_compatible(
	core::conformation::ResidueCOP res1,
	core::conformation::ResidueCOP res2
) const
{

	using namespace core::conformation;

	if( !found_res_compatibility_determined_ ){
		utility_exit_with_message( "Error: trying to lookup residue compatibility without having determined it previously." );
	}


	std::map< ResidueCOP, std::map< ResidueCOP, core::Size > >::const_iterator map_this_res = found_res_compatibility_.find( res1 );

	if( map_this_res == found_res_compatibility_.end() ){
		utility_exit_with_message( "Error: no residue compatibility map found for residue "+utility::to_string( res1->seqpos() )+"." );
	}

	std::map< ResidueCOP, core::Size>::const_iterator res2_it = map_this_res->second.find( res2 );

	if( res2_it == map_this_res->second.end() ){
			utility_exit_with_message( "Error: no residue compatibility info found for res1 "+utility::to_string( res1->seqpos() )+" and res2 "+utility::to_string( res2->seqpos() ) +"." );
	}

	return res2_it->second;


} //residues compatible



/// @brief this function does clash checks between all residues that were found,
/// and saves the information in a map
void
SecondaryMatchProtocol::determine_found_residues_compatibility(
	core::pose::PoseCOP ref_poseCOP )
{

	core::Real Ecutoff = basic::options::option[ basic::options::OptionKeys::enzdes::sc_sc_bump_cutoff ];
	found_res_compatibility_.clear();

	for( core::Size i = 1; i <= found_resis_.size(); ++i ){


		for( utility::vector1< core::conformation::ResidueOP >::const_iterator res1_it = found_resis_[i].begin();
				 res1_it != found_resis_[i].end(); ++res1_it){

			std::map< core::conformation::ResidueCOP, core::Size > map_this_res;

			for( core::Size j  = i+1; j<= found_resis_.size(); ++j ){

				for( utility::vector1< core::conformation::ResidueOP>::const_iterator res2_it = found_resis_[j].begin();
						 res2_it != found_resis_[j].end(); ++res2_it){

					core::Size compatible(0);
					if( (*res1_it)->seqpos() != (*res2_it)->seqpos() ){

						//perform the clash check between res1 and res2
						core::scoring::EnergyMap emap;
						reduced_scofx_->bump_check_full( **res1_it, **res2_it, *ref_poseCOP, emap );

						if( (reduced_scofx_->weights().dot( emap )) < Ecutoff ) compatible = 1;
						//std::cerr << "just checked compatibility for i=" << i << " and j=" << j;
						//std::cerr << ", res i is a " << (*res1_it)->name3() << " at " << (*res1_it)->seqpos() << ", and res j is a " << (*res2_it)->name3() << " at " << (*res2_it)->seqpos() <<  std::endl;
					}

					map_this_res.insert( std::pair<core::conformation::ResidueCOP, core::Size> (*res2_it, compatible ) );
				}
			}

				//save the map for this residue
			found_res_compatibility_.insert( std::pair< core::conformation::ResidueCOP, std::map< core::conformation::ResidueCOP, core::Size > > ( *res1_it, map_this_res) );

		}

	}

	found_res_compatibility_determined_ = true;

} // determine_found_residues_compatibility



PoseFoundResiduesCombination::~PoseFoundResiduesCombination() {}

PoseFoundResiduesCombination::PoseFoundResiduesCombination(
	core::pose::PoseCOP ref_pose_in,
	SecondaryMatchProtocolCAP secmatch_in
	) : ref_pose_( ref_pose_in ),
			secmatch_prot_( secmatch_in )
{
	combine_resis_.clear();
}

void
PoseFoundResiduesCombination::add_residue( core::conformation::ResidueOP res_in ){
	combine_resis_.push_back( res_in );
}


bool
PoseFoundResiduesCombination::construct_and_dump_outpose(
	utility::vector1< toolbox::match_enzdes_util::EnzConstraintParametersCOP > match_params )
{
	using namespace basic::datacache;
	//using core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG;

	//1.first we have to copy the refpose
	core::pose::Pose outpose( *ref_pose_ );

	CacheableStringOP outcache = outpose.data().get_ptr< CacheableString >( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG );
	core::pose::PDBInfoOP pose_pdbinfo = outpose.pdb_info();
	std::string outtag = outcache->str() + "_sm";

	//2. then we need to go through the found residues, figure
	// out if they are all compatible with each other and
	// generate the new header line

  bool all_residues_compatible = true;

	core::Size num_resis_to_combine = combine_resis_.size();

	if( num_resis_to_combine != match_params.size() ){
		utility_exit_with_message("ERROR: number of matching parameters doesn't fit number of residues.\n");
	}

	for( core::Size i = 1; i <= num_resis_to_combine; ++i){

		//generate the header line for this interaction
		std::string mis_chain(""), mis_name3("");
		std::string other_chain(""), other_name3("");
		core::Size mispos(0), otherpos(0);

		mispos = outpose.pdb_info()->number( combine_resis_[i]->seqpos() );

		//for traditional reasons, if residue is a ligand, it will be set to 0
		if( outpose.residue_type( combine_resis_[i]->seqpos() ).is_ligand() ) mispos = 0;

		mis_name3 = combine_resis_[i]->name3();
		mis_chain =  outpose.pdb_info()->chain( combine_resis_[i]->seqpos() );

		utility::vector1< core::Size > other_positions;
		toolbox::match_enzdes_util::EnzCstTemplateResCOP other_template( match_params[i]->get_missing_template_other_res( outpose ) );
		toolbox::match_enzdes_util::EnzCstTemplateResCacheCOP other_template_cache( toolbox::match_enzdes_util::get_enzdes_observer( outpose )->cst_cache()->param_cache( match_params[i]->cst_block() )->template_res_cache( other_template->param_index() ) );

		for( std::map< Size, toolbox::match_enzdes_util::EnzCstTemplateResAtomsOP >::const_iterator pos_it = other_template_cache->seqpos_map_begin(), pos_end = other_template_cache->seqpos_map_end();
				 pos_it != pos_end; ++pos_it ){
		//for( std::map< Size, EnzCstTemplateResAtomsOP >::const_iterator pos_it = match_params[i]->get_missing_template_other_res( outpose )->respos_map_begin();
				 //pos_it != match_params[i]->get_missing_template_other_res( outpose )->respos_map_end(); ++pos_it ){
			other_positions.push_back( pos_it->first );
		}
		if( other_positions.size() != 1 ){
			utility_exit_with_message("Impossible error just happened...");
		}

		otherpos = outpose.pdb_info()->number( other_positions[1] );

		//for traditional reasons, if residue is a ligand, it will be set to 0
		if( outpose.residue_type( other_positions[1] ).is_ligand() ) otherpos = 0;

		other_name3 = outpose.residue_type( other_positions[1] ).name3();
		other_chain = outpose.pdb_info()->chain( other_positions[1] );

		core::pose::RemarkInfo ri;

		if( match_params[i]->resA() == match_params[i]->get_missing_template_res( outpose ) ){
			//ri.value = "BONE TEMPLATE "+mis_remark+" MATCH MOTIF "+other_remark;
			ri.value = toolbox::match_enzdes_util::assemble_remark_line( mis_chain, mis_name3, mispos, other_chain, other_name3, otherpos, match_params[i]->cst_block() );
		}
		else if( match_params[i]->resB() == match_params[i]->get_missing_template_res( outpose ) ){
			//			ri.value = "BONE TEMPLATE "+other_remark+" MATCH MOTIF "+mis_remark;
			ri.value = toolbox::match_enzdes_util::assemble_remark_line( other_chain, other_name3, otherpos, mis_chain, mis_name3, mispos, match_params[i]->cst_block() );
		}
		else{
			utility_exit_with_message("Weird. Suddenly no residue is missing anymore.\n");
		}

		pose_pdbinfo->remarks().push_back( ri );

		//then determine if it's compatible with all other residues

		for( core::Size j = i + 1; j <= num_resis_to_combine; ++j){

			//if( ! residues_compatible_without_checking( combine_resis_[i], combine_resis_[j] ) ){
			if( secmatch_prot_->residues_compatible( combine_resis_[i], combine_resis_[j] ) == 0 ) {
				all_residues_compatible = false;
				return false;
			}
		}

		//finally put the residue into the pose
		outpose.replace_residue( combine_resis_[i]->seqpos(), *(combine_resis_[i]), true);

		outtag = outtag + combine_resis_[i]->name1() + utility::to_string( combine_resis_[i]->seqpos() ) ;


	} // loop over all residues

	outtag = outtag + ".pdb";

	//3. now we can output
	//if all the residues are compatible, we just dump
	//else we have to repack/minimize the pose, and only dump if
	//the cst_score is still below the cutoff

	//note: the all_residues_compatible check is not really necessary, but just to make sure...
	if( all_residues_compatible ){
		outpose.dump_pdb( outtag );
		return true;
	}
	return false;

} //construct and dump outpose



} //namespace enzdes
} //namespace protocols


