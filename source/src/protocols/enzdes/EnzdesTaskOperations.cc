// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/enzdes/DesignProteinLigandInterface.hh
/// @brief various task operations used in enzyme design
/// @author Florian Richter (floric@u.washington.edu), Sinisa Bjelic (sbjelic@u.washington.edu), Rocco Moretti (rmoretti@u.washington.edu)

#include <protocols/motifs/BuildPosition.hh>  // REQUIRED FOR WINDOWS
#include <protocols/motifs/Motif.hh>  // REQUIRED FOR WINDOWS

#include <protocols/enzdes/EnzdesTaskOperations.hh>
#include <protocols/enzdes/DetectProteinLigandInterfaceOperationCreator.hh>
#include <protocols/enzdes/ProteinLigandInterfaceUpweighterOperationCreator.hh>
#include <protocols/enzdes/AddRigidBodyLigandConfsCreator.hh>
#include <protocols/enzdes/SetCatalyticResPackBehaviorCreator.hh>
#include <protocols/enzdes/AddLigandMotifRotamersOperationCreator.hh>

#include <protocols/enzdes/enzdes_util.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCacheableObserver.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>
#include <protocols/toolbox/match_enzdes_util/EnzdesCstCache.hh>

#include <protocols/dna/util.hh> // Needed for arginine sweep interface detection

#include <protocols/motifs/LigandMotifSearch.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/toolbox/IGEdgeReweighters.hh>
#include <protocols/toolbox/rotamer_set_operations/RigidBodyMoveRotSetOps.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/IGEdgeReweightContainer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <basic/options/option.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <core/pack/task/operation/task_op_schemas.hh>
#include <utility/pointer/access_ptr.hh>

// option key includes

// utility includes
#include <utility/string_util.hh>

#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/chemical/VariantType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace enzdes {

using namespace core::pack::task::operation;
using namespace utility::tag;

static THREAD_LOCAL basic::Tracer tr( "protocols.enzdes.EnzdesTaskOperations" );

SetCatalyticResPackBehavior::SetCatalyticResPackBehavior():
	fix_catalytic_aa_(basic::options::option[basic::options::OptionKeys::enzdes::fix_catalytic_aa]),
	behavior_non_catalytic_("")
{}

SetCatalyticResPackBehavior::~SetCatalyticResPackBehavior()= default;

core::pack::task::operation::TaskOperationOP
SetCatalyticResPackBehavior::clone() const
{
	SetCatalyticResPackBehaviorOP to_return( new SetCatalyticResPackBehavior() );
	to_return->set_fix_catalytic_aa( this->fix_catalytic_aa_);
	to_return->behavior_non_catalytic_ = this->behavior_non_catalytic_;
	return to_return;
}

void
SetCatalyticResPackBehavior::apply(
	Pose const & pose,
	PackerTask & task ) const
{

	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return;
	toolbox::match_enzdes_util::EnzdesCstCacheCOP cst_cache( enz_obs->cst_cache() );
	if ( !cst_cache ) return;
	toolbox::match_enzdes_util::EnzConstraintIOCOP cstio( cst_cache->enzcst_io() );
	std::string catinfo("SetCatalyticResPackBehavior task operation touched the following catalytic residues: ");

	bool also_touch_noncatalytic( behavior_non_catalytic_ != "" );
	core::pack::task::ResfileCommandOP resf_command;

	//in case the user requested it, we'll use the resfile reading machinery
	//to also set the behavior for the noncatalytic residues
	if ( also_touch_noncatalytic ) {
		std::map< std::string, core::pack::task::ResfileCommandOP > command_map = core::pack::task::create_command_map();
		resf_command = command_map[ behavior_non_catalytic_ ]->clone();
	}

	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {

		if ( cst_cache->contains_position( i ) ) {
			catinfo = catinfo + utility::to_string( i ) + ", ";

			if ( fix_catalytic_aa_ ) {
				task.nonconst_residue_task( i ).prevent_repacking();
			} else {
				if ( cstio->is_backbone_only_cst( pose, i ) ) {
					utility::vector1< bool > allowed_aas( core::chemical::num_canonical_aas, false );
					utility::vector1< std::string > allowed_name3 = cstio->allowed_res_name3_at_position( pose, i );
					for ( core::Size j =1; j <= allowed_name3.size(); ++j ) {
						allowed_aas[ core::chemical::aa_from_name( allowed_name3[ j ] ) ] = true;
					}
					task.nonconst_residue_task(i).restrict_absent_canonical_aas( allowed_aas );
				} else {
					task.nonconst_residue_task(i).restrict_to_repacking();
					//is catalytic oversampling requested?
					if ( basic::options::option[basic::options::OptionKeys::enzdes::ex_catalytic_rot].user() &&
							pose.residue_type( i ).is_protein() ) {

						core::pack::task::ExtraRotSample cat_rot_sample( static_cast<core::pack::task::ExtraRotSample>(basic::options::option[basic::options::OptionKeys::enzdes::ex_catalytic_rot].value() ) );
						if ( pose.residue_type( i ).nchi() >= 1 ) {
							task.nonconst_residue_task(i).or_ex1( true );
							task.nonconst_residue_task(i).or_ex1_sample_level( cat_rot_sample );
						}
						if ( pose.residue_type( i ).nchi() >= 2 ) {
							task.nonconst_residue_task(i).or_ex2( true );
							task.nonconst_residue_task(i).or_ex2_sample_level( cat_rot_sample );
						}
						if ( pose.residue_type( i ).nchi() >= 3 ) {
							task.nonconst_residue_task(i).or_ex3( true );
							task.nonconst_residue_task(i).or_ex3_sample_level( cat_rot_sample );
						}
						if ( pose.residue_type( i ).nchi() >= 4 ) {
							task.nonconst_residue_task(i).or_ex4( true );
							task.nonconst_residue_task(i).or_ex4_sample_level( cat_rot_sample );
						}
					}//if  catalytic oversampling
				} //if( backbone only cst ) else
			} //if fix_catalytic_aa else
		} else if ( also_touch_noncatalytic ) resf_command->residue_action(task, i ); //if position is catalytic
	} //loop over all pose residues
	tr << catinfo << std::endl;
}

void
SetCatalyticResPackBehavior::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("fix_catalytic_aa") ) fix_catalytic_aa_ = tag->getOption<bool>( "fix_catalytic_aa", 1);
	if ( tag->hasOption("behavior_non_catalytic") ) behavior_non_catalytic_ = tag->getOption<std::string>("behavior_non_catalytic" );

}

void SetCatalyticResPackBehavior::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes
		+ XMLSchemaAttribute::attribute_w_default(  "fix_catalytic_aa", xsct_rosetta_bool, "XRW TO DO",  "true"  )
		+ XMLSchemaAttribute( "behavior_non_catalytic", xs_string , "XRW TO DO" );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

core::pack::task::operation::TaskOperationOP
SetCatalyticResPackBehaviorCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new SetCatalyticResPackBehavior );
}

void SetCatalyticResPackBehaviorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetCatalyticResPackBehavior::provide_xml_schema( xsd );
}

std::string SetCatalyticResPackBehaviorCreator::keyname() const
{
	return SetCatalyticResPackBehavior::keyname();
}


core::pack::task::operation::TaskOperationOP
DetectProteinLigandInterfaceOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new DetectProteinLigandInterface );
}

void DetectProteinLigandInterfaceOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DetectProteinLigandInterface::provide_xml_schema( xsd );
}

std::string DetectProteinLigandInterfaceOperationCreator::keyname() const
{
	return DetectProteinLigandInterface::keyname();
}

core::pack::task::operation::TaskOperationOP
ProteinLigandInterfaceUpweighterOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinLigandInterfaceUpweighter );
}

void ProteinLigandInterfaceUpweighterOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ProteinLigandInterfaceUpweighter::provide_xml_schema( xsd );
}

std::string ProteinLigandInterfaceUpweighterOperationCreator::keyname() const
{
	return ProteinLigandInterfaceUpweighter::keyname();
}

core::pack::task::operation::TaskOperationOP
AddLigandMotifRotamersOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new AddLigandMotifRotamers );
}

void AddLigandMotifRotamersOperationCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddLigandMotifRotamers::provide_xml_schema( xsd );
}

std::string AddLigandMotifRotamersOperationCreator::keyname() const
{
	return AddLigandMotifRotamers::keyname();
}

DetectProteinLigandInterface::DetectProteinLigandInterface():
	detect_design_interface_(true), catalytic_res_part_of_interface_(false), design_(true),
	repack_only_(false), score_only_(false),
	resfilename_(), //Empty string
	add_observer_cache_segs_to_interface_(false),
	no_design_cys_(true),
	catres_only_(false),
	use_cstid_list_(false),
	cstid_list_()
{
	init_from_options();
	design_target_res_.clear();
}

DetectProteinLigandInterface::~DetectProteinLigandInterface() = default;

core::pack::task::operation::TaskOperationOP
DetectProteinLigandInterface::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new DetectProteinLigandInterface(*this) );
}

/// @brief Initialize the class based on the command line options.
void
DetectProteinLigandInterface::init_from_options()
{
	detect_design_interface_ = basic::options::option[basic::options::OptionKeys::enzdes::detect_design_interface];
	arg_sweep_interface_ = basic::options::option[basic::options::OptionKeys::enzdes::arg_sweep_interface];
	arg_sweep_cutoff_ = basic::options::option[basic::options::OptionKeys::enzdes::arg_sweep_cutoff];
	score_only_ = basic::options::option[basic::options::OptionKeys::enzdes::enz_score];
	repack_only_ = basic::options::option[basic::options::OptionKeys::enzdes::enz_repack];
	cut1_ = basic::options::option[basic::options::OptionKeys::enzdes::cut1];
	cut2_ = basic::options::option[basic::options::OptionKeys::enzdes::cut2];
	cut3_ = basic::options::option[basic::options::OptionKeys::enzdes::cut3];
	cut4_ = basic::options::option[basic::options::OptionKeys::enzdes::cut4];
	if ( basic::options::option[basic::options::OptionKeys::packing::resfile].user() ) {
		resfilename_ = basic::options::option[basic::options::OptionKeys::packing::resfile].value().at(1);
	}
}

void
DetectProteinLigandInterface::parse_tag( TagCOP tag , DataMap & )
{
	detect_design_interface_ = true; //Using DetectProteinLigandInterface implies you want to.
	score_only_ = false; //Using DetectProteinLigandInterface implies you're doing more than scoring.
	if ( tag->hasOption("repack_only") ) repack_only_ = tag->getOption< bool >( "repack_only", false );
	cut1_ = tag->getOption< core::Real >( "cut1", 6.0 );
	cut2_ = tag->getOption< core::Real >( "cut2", 8.0 );
	cut3_ = tag->getOption< core::Real >( "cut3", 10.0 );
	cut4_ = tag->getOption< core::Real >( "cut4", 12.0 );
	arg_sweep_cutoff_ = tag->getOption< core::Real >( "arg_sweep_cutoff", 3.7 );
	design_ = tag->getOption< bool >( "design", 1 );
	resfilename_ =  tag->getOption< std::string >( "resfile", "");
	no_design_cys_ = ! ( tag->getOption< bool >( "design_to_cys", 0 ) );
	if ( tag->hasOption("segment_interface") ) add_observer_cache_segs_to_interface_ = tag->getOption< bool >( "segment_interface", true );

	if ( tag->hasOption("catres_interface") )  catalytic_res_part_of_interface_ = tag->getOption< bool >( "catres_interface", true );
	if ( tag->hasOption("arg_sweep_interface") )  arg_sweep_interface_ = tag->getOption< bool >( "arg_sweep_interface", true );
	if ( tag->hasOption("catres_only_interface") ) catres_only_ = tag->getOption< bool >( "catres_only_interface", true );
	if ( tag->hasOption("target_cstids") ) {
		use_cstid_list_ = true;
		cstid_list_ = tag->getOption< std::string >( "target_cstids", "" );
	}

}

void DetectProteinLigandInterface::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "repack_only", xsct_rosetta_bool, "XRW TO DO", "false" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "cut1", xsct_real, "XRW TO DO", "6.0" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "cut2", xsct_real, "XRW TO DO", "8.0" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "cut3", xsct_real, "XRW TO DO", "10.0" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "cut4", xsct_real, "XRW TO DO", "12.0" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "arg_sweep_cutoff", xsct_real, "XRW TO DO", "3.7" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "XRW TO DO", "true" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "resfile", xs_string, "XRW TO DO", "" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "design_to_cys", xsct_rosetta_bool, "XRW TO DO", "false" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "segment_interface", xsct_rosetta_bool, "XRW TO DO", "true" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "catres_interface", xsct_rosetta_bool, "XRW TO DO", "true" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "arg_sweep_interface", xsct_rosetta_bool, "XRW TO DO", "true" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "catres_only_interface", xsct_rosetta_bool, "XRW TO DO", "true" ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default( "target_cstids", xs_string, "XRW TO DO", "" ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

/// @brief Change a packer task in some way.  The input pose is the one to which the input
/// task will be later applied.
void DetectProteinLigandInterface::apply(
	Pose const & pose,
	PackerTask & task) const
{

	// First read the resfile, setting the appropriate operations in the task
	// Note that "AUTO" tags will be processed later by detect_design_interface
	if ( ! resfilename_.empty() ) {
		tr.Info << "Reading resfile input from: " << resfilename_ <<  std::endl;
		core::pack::task::operation::ReadResfileAndObeyLengthEvents resfile_read( resfilename_ );
		resfile_read.apply( pose, task );
		if ( !design_ ) {
			for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
				if ( task.design_residue(i) ) task.nonconst_residue_task(i).restrict_to_repacking();
			}
		}
	} // end (if resfile)

	// detect design interface, only at positions marked "AUTO" in resfile if there is a resfile
	if ( detect_design_interface_ ) {
		utility::vector1< bool > repack_res( pose.size(), false );
		utility::vector1< bool > design_res( pose.size(), false );
		utility::vector1< bool > detect_res( pose.size() );
		core::Real cut1(cut1_), cut2(cut2_), cut3(cut3_), cut4(cut4_);

		if ( cut1 < 0.0 || cut2 < cut1 || cut3 < cut2 || cut4 < cut3 ) {
			if ( cut1 < 0.0 ) { cut1 = 0.0; }
			if ( cut2 < cut1 ) { cut2 = cut1; }
			if ( cut3 < cut2 ) { cut3 = cut2; }
			if ( cut4 < cut3 ) { cut4 = cut3; }
			tr.Warning << "WARNING: detect design interface cutpoints should be in ascending order. Was " << cut1_ << " " << cut2_ << " " << cut3_ << " " << cut4_ << "; Reset to " << cut1 << " " << cut2 << " " << cut3 << " " << cut4 << std::endl;
		}
		if ( !design_ ) {
			cut1 = 0.0;
			cut2 = 0.0;
		}

		std::set< core::Size > interface_target_res = design_target_res_;

		// if picking neighbors of particular cstids only
		if ( ( interface_target_res.empty() /*size() == 0 */) && use_cstid_list_ ) {
			utility::vector1< core::Size > trg_res;
			protocols::enzdes::enzutil::get_resnum_from_cstid_list( cstid_list_, pose, trg_res );
			auto last = std::unique( trg_res.begin(), trg_res.end() );
			std::set< core::Size > trg_set( trg_res.begin(), last ); //trg_res.end() );
			interface_target_res = trg_set;
		}

		if ( ( interface_target_res.empty() /*size() == 0 */ ) && (catalytic_res_part_of_interface_ || catres_only_) ) {
			for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
				if ( enzutil::is_catalytic_seqpos( pose, i ) ) {
					if ( pose.residue_type( i ).is_ligand() && !catres_only_ ) interface_target_res.insert( i );
					if ( !pose.residue_type( i ).is_ligand() ) interface_target_res.insert( i );
				}
			}
		}

		if ( add_observer_cache_segs_to_interface_ ) add_observer_cache_segments_to_set( pose, interface_target_res );

		if ( interface_target_res.empty() /*size() == 0 */&& !catres_only_ ) {
			if ( core::pose::symmetry::is_symmetric(pose) ) {
				core::kinematics::FoldTree asymm_ft( core::conformation::symmetry::get_asymm_unit_fold_tree( pose.conformation() ) );
				interface_target_res.insert( asymm_ft.downstream_jump_residue( asymm_ft.num_jump() ) );
			} else interface_target_res.insert( pose.fold_tree().downstream_jump_residue( pose.num_jump() ) );
		}

		tr.Info << "Choosing the following residues as targets for detecting interface: ";
		for ( core::Size interface_target_re : interface_target_res ) {
			tr.Info << interface_target_re << "+";
		}
		tr.Info <<  std::endl;

		// initialize detect_res vector, specifies whether the designability of a residue should be decided by find_design_interface
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( ! resfilename_.empty() ) {
				if ( task.residue_task( i ).has_behavior("AUTO") ) {
					detect_res[i] = true;
					// tr.Info << "detectable:" << i << std::endl;
				} else {
					detect_res[i] = false;
					// tr.Info << "not detectable:" << i << std::endl;
				}
			} else { // end if (resfile)
				// if no resfile input, set all to automatic detection
				detect_res[i] = true;
			}
		}
		//Now we have two ways to find interface: typical sphere method, and new arg sweep method (if arginine rotamers can interact with ligand, position is designable)
		if ( arg_sweep_interface_ == true ) {
			core::Real arg_sweep_cutoff(arg_sweep_cutoff_);
			find_design_interface_arg_sweep( pose, interface_target_res, cut1, cut2, cut3, cut4, arg_sweep_cutoff, repack_res, design_res );
		} else {
			find_design_interface( pose, interface_target_res, cut1, cut2, cut3, cut4, repack_res, design_res );
		}

		//setup the task accordingly
		//default behavior for design will be everything except cys, and of course we want to leave disulfide bonds untouched
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			// only do auto-detection if the residue was set to AUTO above in detect_res initialization
			if ( detect_res[i] == true ) {
				if ( design_res[i] == true ) {
					if ( pose.residue( i ).aa() == core::chemical::aa_cys && pose.residue( i ).has_variant_type( core::chemical::DISULFIDE ) ) {
						task.nonconst_residue_task( i ).restrict_to_repacking();
					} else if ( no_design_cys_ ) {
						utility::vector1< bool > keep_aas( core::chemical::num_canonical_aas, true );
						if ( pose.residue( i ).aa() != core::chemical::aa_cys ) keep_aas[ core::chemical::aa_cys ] = false;
						task.nonconst_residue_task(i).restrict_absent_canonical_aas( keep_aas );
					}
				} else if ( repack_res[i] == true ) {
					task.nonconst_residue_task(i).restrict_to_repacking();
				} else if ( pose.residue(i).is_ligand() && ( ! core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( pose.residue_type(i) )) ) {
					task.nonconst_residue_task(i).restrict_to_repacking();
				} else {
					task.nonconst_residue_task( i ).prevent_repacking();
				}
			}
		}  // end pack/design assignment loop
	} // detect design interface


	//in case we are only interested in scoring
	if ( score_only_ ) {
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			task.nonconst_residue_task(i).prevent_repacking();
		}
	} else if ( repack_only_ ) {
		//in case we are only interested in repacking
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( task.design_residue(i) ) {
				task.nonconst_residue_task(i).restrict_to_repacking();
			}
		}
	}

	// As a final check print out everything designed
	tr.Info << "Final Design Shell Residues: ";
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( task.design_residue(i) ) {
			tr.Info << i << ", ";
		}
	}
	tr.Info << std::endl  << "Final Repack Shell Residues: ";
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( task.pack_residue(i) ) {
			tr.Info << i << ", ";
		}
	}
	tr.Info << std::endl;

	//lastly, in case of symmetry the task should be symmetrized by union
	if ( core::pose::symmetry::is_symmetric(pose) ) task.request_symmetrize_by_union();
} //apply

void
DetectProteinLigandInterface::find_design_interface(
	core::pose::Pose const & pose,
	std::set< core::Size > const & interface_target_res,
	core::Real cut1,
	core::Real cut2,
	core::Real cut3,
	core::Real cut4,
	utility::vector1< bool > & repack_res,
	utility::vector1< bool > & design_res
) const
{

	core::Real cut1_sq = cut1 * cut1;
	core::Real cut2_sq = cut2 * cut2;
	core::Real cut3_sq = cut3 * cut3;
	core::Real cut4_sq = cut4 * cut4;

	for ( auto targ_it( interface_target_res.begin()),targ_end(interface_target_res.end());
			targ_it != targ_end; ++targ_it ) {

		repack_res[ *targ_it ] = true;
		// on protein side, have to do distance check
		core::conformation::Residue const & targ_rsd = pose.residue( *targ_it );
		core::Size targ_res_atom_start = 1;
		if ( targ_rsd.is_protein() ) {
			design_res[ *targ_it ] = true; //might be designable
			repack_res[ *targ_it ] = false;
			targ_res_atom_start = targ_rsd.first_sidechain_atom();
		}

		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( design_res[i] ) continue; //in case this is already set to design, we don't have to loop over it again
			if ( interface_target_res.find( i ) != interface_target_res.end() ) continue;
			core::conformation::Residue const & prot_rsd = pose.residue(i);
			//tr << "Rsd has " << targ_rsd.nheavyatoms() << " heavy atoms " << std::endl;
			for ( core::Size k = targ_res_atom_start, k_end = targ_rsd.nheavyatoms(); k <= k_end; ++k ) {
				//tr << "Current heavy is " << targ_rsd.atom_name( k ) << std::endl;
				core::Vector prot_cb, prot_ca;
				if ( prot_rsd.has("CB") ) prot_cb = prot_rsd.xyz("CB");
				if ( prot_rsd.has("CA") ) prot_ca = prot_rsd.xyz("CA"); // GLY
				core::Real ca_dist2 = targ_rsd.xyz(k).distance_squared( prot_ca );
				if ( ca_dist2 <= cut4_sq ) {
					if ( ca_dist2 <= cut3_sq ) {
						if ( ca_dist2 <= cut2_sq ) {
							if ( ca_dist2 <= cut1_sq ) {
								design_res[i] = true;
								repack_res[i] = false;
								break;
							} else if ( prot_rsd.has("CB") ) { // cut1
								core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
								//                tr.Info << "cb_dist2 is " << cb_dist2 << "; ";
								if ( cb_dist2 < ca_dist2 ) {
									design_res[i] = true;
									repack_res[i] = false;
									break;
								} else {
									repack_res[i] = true;
								}
							} else if ( prot_rsd.has("2HA") ) {  // end of non-gly residues   //glycine doesn't have a CB, so use 2HA to get position where CB would be
								// use the name "cb" to describe the 2HA atom; design if 2HA < CA
								prot_cb = prot_rsd.xyz("2HA");
								core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
								if ( cb_dist2 < ca_dist2 ) {   // 2HA is closer than CA
									design_res[i] = true;
									repack_res[i] = false;
									break;
								} else {   // 2HA is further than CA
									repack_res[i] = true;
								}
							} else {  // end of gly residues  // Exception handling case for residue without CB or 2HA
								tr.Info << "Weird residue without CB or 2HA. Watch out! Residue:" << i << std::endl;
								design_res[i] = false;
								repack_res[i] = true;
								break;
							} // end of exception catching for neither CB nor 2HA
						} else { //cut2
							repack_res[i] = true;
						}
					} else if ( prot_rsd.has("CB") ) { //cut3
						core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
						if ( cb_dist2 < ca_dist2 ) {
							repack_res[i] = true;
						}
					}
				} //cut4
			} //loop over target res atoms
		} //loop over protein residues
	} //loop over target residues

	std::string repackres_string(""), designres_string;
	core::Size num_repack(0), num_design(0);
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( repack_res[i] ) {
			num_repack++;
			repackres_string = repackres_string + utility::to_string( i ) + "+";
		}
		if ( design_res[i] ) {
			num_design++;
			designres_string = designres_string + utility::to_string( i ) + "+";
		}
		if ( ( repack_res[i] == true)  && ( design_res[i] == true ) ) { tr.Info << "Huh? this should not happen. " << std::endl; }
	}

	tr.Info << "Design Interface: detected " << num_design << " design-shell residues and " << num_repack << " repack-shell residues, shell sizes cut1-4 used were " << cut1 << " " << cut2 << " " << cut3 << " " << cut4 << std::endl << "Design-shell Residues(pose-numbering): " << designres_string;
	tr.Info << std::endl << "Repack-shell Residues(pose-numbering): " << repackres_string << std::endl;

} //find_design_interface

void
DetectProteinLigandInterface::find_design_interface_arg_sweep(
	core::pose::Pose const & pose,
	std::set< core::Size > const & interface_target_res,
	core::Real cut1,
	core::Real cut2,
	core::Real cut3,
	core::Real cut4,
	core::Real arg_sweep_cutoff,
	utility::vector1< bool > & repack_res,
	utility::vector1< bool > & design_res
) const
{
	core::Real cut1_sq = cut1 * cut1;
	core::Real cut2_sq = cut2 * cut2;
	core::Real cut3_sq = cut3 * cut3;
	core::Real cut4_sq = cut4 * cut4;

	for ( auto targ_it( interface_target_res.begin()),targ_end(interface_target_res.end());
			targ_it != targ_end; ++targ_it ) {

		repack_res[ *targ_it ] = true;
		// on protein side, have to do distance check
		core::conformation::Residue const & targ_rsd = pose.residue( *targ_it );
		core::Size targ_res_atom_start = 1;
		if ( targ_rsd.is_protein() ) {
			design_res[ *targ_it ] = true; //might be designable
			repack_res[ *targ_it ] = false;
			targ_res_atom_start = targ_rsd.first_sidechain_atom();
		}

		//  tr.Info << "Design residues by arg_sweep are: ";
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( design_res[i] ) continue;  //in case this is already set to design, we don't have to loop over it again
			if ( interface_target_res.find( i ) != interface_target_res.end() ) continue;
			core::conformation::Residue const & prot_rsd = pose.residue(i);
			// Defining private members as magic numbers, shouldn't do this...
			core::Real contact_threshold_ = arg_sweep_cutoff * arg_sweep_cutoff ; // Originally 3.7*3.7 (13.7)
			core::Real close_threshold_ =  10. * 10. ;
			bool base_only_ = false; //Doesn't matter here, chooses between all atoms in residue and sidechain atoms for iterator

			//Super simple arg sweep interface detection routinue
			bool close_to_lig(  protocols::dna::close_to_dna( prot_rsd, targ_rsd, close_threshold_, base_only_ ) );
			if ( close_to_lig ) {
				core::Real dis2 =  protocols::dna::argrot_dna_dis2( pose, i, prot_rsd, targ_rsd, contact_threshold_, base_only_ ); // i is protein residue id
				if ( dis2 <= contact_threshold_ ) {
					design_res[i] = true;
					repack_res[i] = false;
					// tr.Info << i << ", ";
					continue; //If arg_sweep finds design_res, go onto next resi
				}
			}
			// tr.Info << "on protein resi number " << i << " , close_to_lig is " << close_to_lig << " and dis2 from arginine rotamer sweep is " << dis2 << std::endl ;
			for ( core::Size k = targ_res_atom_start, k_end = targ_rsd.nheavyatoms(); k <= k_end; ++k ) { //Now looping over atoms in ligand, I don't need to do that
				core::Vector prot_cb, prot_ca;
				if ( prot_rsd.has("CB") ) prot_cb = prot_rsd.xyz("CB");
				if ( prot_rsd.has("CA") ) prot_ca = prot_rsd.xyz("CA"); // GLY
				core::Real ca_dist2 = targ_rsd.xyz(k).distance_squared( prot_ca );
				if ( ca_dist2 <= cut4_sq ) {
					if ( ca_dist2 <= cut3_sq ) {
						if ( ca_dist2 <= cut2_sq ) {
							if ( ca_dist2 <= cut1_sq ) {
								design_res[i] = false;
								repack_res[i] = true;
								break;
							} else if ( prot_rsd.has("CB") ) { // cut1
								core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
								//                tr.Info << "cb_dist2 is " << cb_dist2 << "; ";
								if ( cb_dist2 < ca_dist2 ) {
									design_res[i] = false;
									repack_res[i] = true;
									break;
								} else {
									repack_res[i] = true;
								}
							} else if ( prot_rsd.has("2HA") ) {  // end of non-gly residues   //glycine doesn't have a CB, so use 2HA to get position where CB would be
								// use the name "cb" to describe the 2HA atom; design if 2HA < CA
								prot_cb = prot_rsd.xyz("2HA");
								core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
								if ( cb_dist2 < ca_dist2 ) {   // 2HA is closer than CA
									design_res[i] = false;
									repack_res[i] = true;
									break;
								} else {   // 2HA is further than CA
									repack_res[i] = true;
								}
							} else {  // end of gly residues  // Exception handling case for residue without CB or 2HA
								tr.Info << "Weird residue without CB or 2HA. Watch out! Residue:" << i << std::endl;
								design_res[i] = false;
								repack_res[i] = true;
								break;
							} // end of exception catching for neither CB nor 2HA
						} else { //cut2
							repack_res[i] = true;
						}
					} else if ( prot_rsd.has("CB") ) { //cut3
						core::Real cb_dist2 = targ_rsd.xyz(k).distance_squared( prot_cb );
						if ( cb_dist2 < ca_dist2 ) {
							repack_res[i] = true;
						}
					}
				} //cut4
			} //loop over target res atoms
		} //loop over protein residues
	} //loop over target residues

	std::string repackres_string(""), designres_string;
	core::Size num_repack(0), num_design(0);
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( repack_res[i] ) {
			num_repack++;
			repackres_string = repackres_string + utility::to_string( i ) + "+";
		}
		if ( design_res[i] ) {
			num_design++;
			designres_string = designres_string + utility::to_string( i ) + "+";
		}
		if ( ( repack_res[i] == true)  && ( design_res[i] == true ) ) { tr.Info << "Huh? this should not happen. " << std::endl; }
	}

	tr.Info << "Design Interface: detected " << num_design << " design-shell residues and " << num_repack << " repack-shell residues, shell sizes cut1-4 used were " << cut1 << " " << cut2 << " " << cut3 << " " << cut4 << std::endl << "Design-shell Residues(pose-numbering): " << designres_string;
	tr.Info << std::endl << "Repack-shell Residues(pose-numbering): " << repackres_string << std::endl;

} //find_design_interface_arg_sweep

void
DetectProteinLigandInterface::register_options()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( OptionKeys::packing::resfile  );

	protocols::simple_moves::MinMover::register_options();
	protocols::simple_moves::PackRotamersMover::register_options();

	option.add_relevant( OptionKeys::enzdes::detect_design_interface );
	option.add_relevant( OptionKeys::enzdes::cut1 );
	option.add_relevant( OptionKeys::enzdes::cut2 );
	option.add_relevant( OptionKeys::enzdes::cut3 );
	option.add_relevant( OptionKeys::enzdes::cut4 );
	option.add_relevant( OptionKeys::enzdes::enz_score );
	option.add_relevant( OptionKeys::enzdes::enz_repack );

}

void
DetectProteinLigandInterface::add_observer_cache_segments_to_set(
	core::pose::Pose const & pose,
	std::set< core::Size > & set
)
{
	using namespace basic::datacache;
	using namespace core::pose::datacache;

	if ( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER) ) {

		utility::vector1< std::pair< core::Size, core::Size > > const & segments = utility::pointer::static_pointer_cast< core::pose::datacache::SpecialSegmentsObserver const >(pose.observer_cache().get_const_ptr( core::pose::datacache::CacheableObserverType::SPECIAL_SEGMENTS_OBSERVER ) )->segments();
		for ( core::Size i = 1; i <= segments.size(); ++i ) {
			for ( core::Size j = segments[i].first; j < segments[i].second; ++j ) set.insert( j );
		}
	}
}


ProteinLigandInterfaceUpweighter::ProteinLigandInterfaceUpweighter()
{
	init_from_options();
	catres_packer_weight_=1.0;
}

ProteinLigandInterfaceUpweighter::~ProteinLigandInterfaceUpweighter() = default;

core::pack::task::operation::TaskOperationOP
ProteinLigandInterfaceUpweighter::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new ProteinLigandInterfaceUpweighter(*this) );
}

/// @brief Initialize the class based on the command line options.
void
ProteinLigandInterfaceUpweighter::init_from_options()
{
	lig_packer_weight_ = basic::options::option[basic::options::OptionKeys::enzdes::lig_packer_weight];
}

void
ProteinLigandInterfaceUpweighter::parse_tag( TagCOP tag , DataMap & )
{
	if ( tag->hasOption("interface_weight") ) lig_packer_weight_ = tag->getOption< core::Real >( "interface_weight", 1.0 );
	if ( tag->hasOption("catres_interface_weight") ) catres_packer_weight_ = tag->getOption< core::Real >( "catres_interface_weight", 1.0 );
}

void ProteinLigandInterfaceUpweighter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	AttributeList attributes;

	attributes.push_back( XMLSchemaAttribute::attribute_w_default(  "interface_weight", xsct_real, "XRW TO DO",  "1.0"  ) );
	attributes.push_back( XMLSchemaAttribute::attribute_w_default(  "catres_interface_weight", xsct_real, "XRW TO DO",  "1.0"  ) );

	task_op_schema_w_attributes( xsd, keyname(), attributes );
}

/// @brief Change a packer task in some way.  The input pose is the one to which the input
/// task will be later applied.
void ProteinLigandInterfaceUpweighter::apply(
	Pose const & pose,
	PackerTask & task) const
{
	//If applicable, set the ligand weigths to the specified value
	if ( lig_packer_weight_ != 1.0 ) {
		core::pack::task::IGEdgeReweighterOP lig_up( new protocols::toolbox::IGLigandDesignEdgeUpweighter( lig_packer_weight_ ) );
		core::pack::task::IGEdgeReweightContainerOP IGreweight = task.set_IGEdgeReweights();
		IGreweight->add_reweighter( lig_up );

		tr.Info << "Packer Energies between ligand and design residues are upweighted by factor " << lig_packer_weight_ << "." << std::endl;

	} //if different ligand weights are asked for
	if ( catres_packer_weight_ !=1.0 ) {
		utility::vector1< core::Size > catres;
		utility::vector1< core::Size > design_residues;
		for ( core::Size ii= 1; ii<=pose.size(); ++ii ) {
			if ( enzutil::is_catalytic_seqpos( pose, ii) && pose.residue( ii ).is_protein() ) catres.push_back( ii );
			else if ( task.design_residue( ii ) ) design_residues.push_back( ii );
		}

		core::pack::task::IGEdgeReweighterOP catres_up( new protocols::toolbox::ResidueGroupIGEdgeUpweighter( catres_packer_weight_, catres , design_residues ) );
		task.set_IGEdgeReweights()->add_reweighter( catres_up );

		tr.Info << "Packer Energies between catalytic residues and design residues are upweighted by factor " << catres_packer_weight_ << "." << std::endl;
	} //if different catalytic-residue weights are asked for

} //apply

void
ProteinLigandInterfaceUpweighter::register_options()
{
	basic::options::option.add_relevant( basic::options::OptionKeys::enzdes::lig_packer_weight );
}

core::pack::task::operation::TaskOperationOP
AddRigidBodyLigandConfsCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new AddRigidBodyLigandConfs );
}

void AddRigidBodyLigandConfsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddRigidBodyLigandConfs::provide_xml_schema( xsd );
}

std::string AddRigidBodyLigandConfsCreator::keyname() const
{
	return AddRigidBodyLigandConfs::keyname();
}

AddRigidBodyLigandConfs::AddRigidBodyLigandConfs(){}

AddRigidBodyLigandConfs::~AddRigidBodyLigandConfs()= default;

core::pack::task::operation::TaskOperationOP
AddRigidBodyLigandConfs::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new AddRigidBodyLigandConfs() );
}

/// @details doesn't do anything atm, because the operation right now
/// puts lig conformers it finds in the pose at apply time into the task
/// this could be expanded to actually having this operation create additional
/// conformers
void
AddRigidBodyLigandConfs::parse_tag( TagCOP /*tag*/ , DataMap & )
{
}

void AddRigidBodyLigandConfs::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

void
AddRigidBodyLigandConfs::apply(
	Pose const & pose,
	PackerTask & task
) const
{

	//std::cout << "starting apply func of AddRigidBodyLigandConfs" << std::endl;
	toolbox::match_enzdes_util::EnzdesCacheableObserverCOP enz_obs( toolbox::match_enzdes_util::get_enzdes_observer( pose ) );
	if ( !enz_obs ) return;
	//std::cout << "non zero observer given ";

	std::map< core::Size, utility::vector1< core::conformation::ResidueCOP > > const & rb_map = enz_obs->lig_rigid_body_confs();

	for ( auto lig_it = rb_map.begin(); lig_it != rb_map.end(); ++lig_it ) {

		if ( !task.being_packed( lig_it->first ) ) continue;
		if ( lig_it->second.size() == 0 ) continue;
		runtime_assert( pose.residue_type( lig_it->first ).name3() == lig_it->second[1]->type().name3() );

		protocols::toolbox::rotamer_set_operations::RigidBodyMoveRSOOP rb_rotsetop( new protocols::toolbox::rotamer_set_operations::RigidBodyMoveRSO( lig_it->first ) );
		rb_rotsetop->set_rigid_body_confs( lig_it->second );

		//std::cout << " instantiated rotamer set operation, passed it " << lig_it->second.size() << " rb confs for position " << lig_it->first << std::endl;

		task.nonconst_residue_task( lig_it->first ).append_rotamerset_operation( rb_rotsetop );
	}
}

void
AddLigandMotifRotamers::register_options()
{
	basic::options::option.add_relevant( basic::options::OptionKeys::enzdes::run_ligand_motifs );
}

AddLigandMotifRotamers::AddLigandMotifRotamers(){}

AddLigandMotifRotamers::~AddLigandMotifRotamers()= default;

core::pack::task::operation::TaskOperationOP
AddLigandMotifRotamers::clone() const
{
	return core::pack::task::operation::TaskOperationOP( new AddLigandMotifRotamers() );
}

void
AddLigandMotifRotamers::parse_tag( TagCOP /*tag*/ , DataMap & )
{
	/*
	-dtest 2
	-r2 0.4
	-r1 0.6
	-z1 1.0
	-z2 0.97
	-dump_motifs true
	-output_file output.file
	-data_file data.file
	-rotlevel 8
	-motif_filename Pruned_NoCCC.motifs
	*/
	//Those are some options I could add, but it's not really necessary--most users will pass those options on the command line.  I will eventually make them parseable tags. -mdsmith
}

void AddLigandMotifRotamers::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	task_op_schema_empty( xsd, keyname() );
}

void
AddLigandMotifRotamers::apply(
	Pose const & pose,
	PackerTask & task
) const
{
	protocols::motifs::LigandMotifSearchOP motif_search( new protocols::motifs::LigandMotifSearch );
	motif_search->run( pose, task );
}

}//namespace enzdes
}//namespace protocols
