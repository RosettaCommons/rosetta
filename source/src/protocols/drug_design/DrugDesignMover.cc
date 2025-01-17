// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/DrugDesignMover.cc
/// @brief MonteCarlo protocol to design drugs in a protein context
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit Headers
#include <protocols/drug_design/DrugDesignMover.hh>
#include <protocols/drug_design/DrugDesignMoverCreator.hh>

// Package Headers
#include <protocols/moves/Mover.hh>
#include <protocols/chemistries/Chemistry.hh>

// Project Headers
#include <protocols/drug_design/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/chemical/PoseResidueTypeSet.hh>
#include <core/chemical/sdf/mol_writer.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/icoor_support.hh>

#include <protocols/filters/Filter.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/numbers.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Debugging output
#include <protocols/jd2/JobDistributor.hh>
#include <utility/io/ozstream.hh>

// C/C++ headers
#include <string>

namespace protocols {
namespace drug_design {

static basic::Tracer TR("protocols.drug_design.DrugDesignMover");

using namespace core;
using namespace protocols::moves;
using namespace protocols::filters;

std::string
DrugDesignMoverCreator::keyname() const
{
	return DrugDesignMover::mover_name();
}

protocols::moves::MoverOP
DrugDesignMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DrugDesignMover );
}

void DrugDesignMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DrugDesignMover::provide_xml_schema( xsd );
}

/// @brief default constructor
DrugDesignMover::DrugDesignMover():
	Mover(mover_name()),
	redocker_( nullptr ),
	scorer_( nullptr ),
	temperature_( 1.0 ),
	maxtrials_( 10 ),
	chain_( 'X' ),
	restypeset_( new core::chemical::PoseResidueTypeSet(core::chemical::FULL_ATOM_t) )
{}

/// @brief destructor
DrugDesignMover::~DrugDesignMover(){}

/// @brief clone this object
MoverOP
DrugDesignMover::clone() const
{
	return MoverOP( new DrugDesignMover( *this ) );
}

/// @brief create this type of object
MoverOP
DrugDesignMover::fresh_instance() const
{
	return MoverOP( new DrugDesignMover() );
}

/// @brief Return the redocker being used
MoverCOP
DrugDesignMover::redocker() const {
	return redocker_;
}

/// @brief Return the scoring filter
FilterCOP
DrugDesignMover::scorer() const {
	return scorer_;
}

/// @brief Add a chemistry to use
void
DrugDesignMover::add_chemistry( protocols::chemistries::ChemistryOP chemistry, core::Real weight /*= 1.0*/ )
{
	assert( chemistries_.size() == weights_.size() );
	chemistries_.push_back( chemistry );
	weights_.push_back( weight );
}

void
DrugDesignMover::add_before_chemistry( protocols::chemistries::ChemistryOP chemistry )
{
	before_chemistries_.push_back( chemistry );
}

void
DrugDesignMover::add_after_chemistry( protocols::chemistries::ChemistryOP chemistry )
{
	after_chemistries_.push_back( chemistry );
}

/// @brief set redocker
void
DrugDesignMover::redocker( MoverOP const & mover )
{
	redocker_ = mover;
}

/// @brief set scoring filter
void
DrugDesignMover::scorer( FilterOP const & scorer )
{
	scorer_ = scorer;
}

/// @brief Replace residues using only those atoms which are true in mask.
void
replace_residue_type_with_mask( core::pose::Pose & pose, core::Size seqpos, core::chemical::ResidueType const & new_rsd_type, utility::vector1< bool > const & mask ) {
	core::conformation::Residue const & source_rsd( pose.residue( seqpos ) );
	core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( new_rsd_type ) );
	core::conformation::Residue & target_rsd( *new_rsd );
	core::conformation::Conformation & conformation( pose.conformation() );

	Size const natoms( target_rsd.natoms() );

	utility::vector1< bool > missing( natoms, false );
	bool any_missing( false );

	for ( Size i=1; i<= natoms; ++i ) {
		if ( i > mask.size() || mask[i] != true ) {
			any_missing = true;
			missing[i] = true;
			continue;
		}
		std::string const & atom_name( target_rsd.atom_name(i) );
		if ( source_rsd.has( atom_name ) ) {
			target_rsd.atom( i ).xyz( source_rsd.atom( atom_name ).xyz() );
		} else {
			TR.Debug << "replace_residue_type_with_mask: missing atom " << target_rsd.name() << ' ' <<
				atom_name << std::endl;
			any_missing = true;
			missing[i] = true;
		}
	}

	if ( any_missing ) {
		target_rsd.seqpos( source_rsd.seqpos() ); // in case fill_missing_atoms needs context info
		target_rsd.chain ( source_rsd.chain () );
		target_rsd.fill_missing_atoms( missing, conformation );
	}

	pose.replace_residue( seqpos, *new_rsd, false );
}

/// @brief Set the filter to use before redocking
void
DrugDesignMover::prefilter( protocols::filters::FilterOP const & setting ) {
	prefilter_ = setting;
}

/// @brief Set the filter to use after redocking
void
DrugDesignMover::postfilter( protocols::filters::FilterOP const & setting ) {
	postfilter_ = setting;
}

void
DrugDesignMover::apply( Pose & pose )
{
	using namespace core::chemical;
	using namespace core::pose;

	if ( ! chemistries_.size() ) {
		TR.Warning << "No chemistries specified! " << std::endl;
		return;
	}
	if ( chemistries_.size() != weights_.size() ) {
		TR.Warning << "The number of weights aren't equal to the number of chemistries - resizing" << std::endl;
		weights_.resize( chemistries_.size(), 1.0 );
	}
	if ( !redocker_ ) {
		TR.Warning << "Redocking mover is empty! " << std::endl;
		return;
	}
	if ( !scorer_ ) {
		TR.Warning << "No scoring Filter specified! " << std::endl;
		return;
	}

	// Statistics for attempts (MC won't show earlier aborts.)
	std::map< std::string, core::Size > stat_attempt, stat_failed, stat_tested, stat_accepted;

	Real initial_score = scorer_->report_sm( pose );
	assert( temperature() ); // Can't work with a temperature of zero
	core::Size res_pos( find_design_position( pose ) );
	TR << "For round " << 0 << " ligand " << pose.residue( res_pos ).name() << " score " << initial_score << std::endl;
	protocols::moves::MonteCarlo mc( pose, initial_score, temperature() );

	//Setup data for the iteration
	numeric::random::WeightedSampler sampler(weights_);

	// TODO: can I pull the trial number out of the MC object instead?
	for ( Size ii=1; ii<=maxtrials_; ii++ ) {
		TR << "Trial number: " << ii <<std::endl;

		// Make a working copy of the restype of interest. (Restypes in Poses should be treated as const.)
		res_pos = find_design_position( pose );
		MutableResidueTypeOP restype( new MutableResidueType( pose.residue_type( res_pos ) ) );
		IndexVDMapping index_vd_mapping( combine( IndexNameMapping( pose.residue_type( res_pos ) ), NameVDMapping(*restype) ) ); // Starting index to current vds.

		std::string original_name( restype->name() ); // Need this as chemistries aren't necessarily all that good with naming.

		if ( pre_process_restype(restype, index_vd_mapping, pose ) ) {
			// ResidueType is unusable
			continue;
		}

		// Pick chemistry to apply
		Size chem_i( sampler.random_sample(numeric::random::rg()) );
		assert( chem_i != 0 && chem_i <= chemistries_.size() );
		assert( chemistries_[chem_i] );
		protocols::chemistries::Chemistry & chemistry( *chemistries_[chem_i] );

		++stat_attempt[ chemistry.name() ];

		// Apply the selected chemistry
		dump_molecule( *restype, "before_mod" );
		chemistry.apply( *restype, pose ); // Allow context-sensitive chemistries
		dump_molecule( *restype, "after_mod" );
		if ( chemistry.get_last_status() != core::chemical::modifications::SUCCESS ) {
			TR.Warning << "Chemistry " << chemistry.name() << " reported failure on " << original_name << std::endl;
			++stat_failed[ chemistry.name() ];
			continue;
		}

		index_vd_mapping = combine( index_vd_mapping, chemistry.get_mapping() );

		std::string new_name( find_new_res_name( original_name, ii ) );
		//Tracer output for post-processing grep, etc.
		TR << "MODIFICATION: " << original_name << " into " << new_name << " via " << chemistry.name() << std::endl;

		if ( post_process_restype(restype, index_vd_mapping, new_name, pose ) ) {
			// ResidueType is unusable
			++stat_failed[ chemistry.name() ];
			continue;
		}

		dump_molecule( *restype, "after_postprocess" );

		if ( emplace_residue_type(pose, restype, index_vd_mapping ) ) {
			// Emplacing the residue type on the pose failed - reset and continue
			pose = mc.last_accepted_pose();
			++stat_failed[ chemistry.name() ];
			continue;
		}

		// Because the position may change during docking
		res_pos = find_design_position( pose );
		dump_molecule( pose.residue( res_pos ), "tested" );
		++stat_tested[ chemistry.name() ];

		// MonteCarlo
		Real score = scorer_->report_sm( pose );
		TR << "For round " << ii << " ligand " << pose.residue( res_pos ).name() << " score " << score << std::endl;
		if ( mc.boltzmann( pose, score, chemistry.name() ) ) {
			++stat_accepted[ chemistry.name() ];
		}
		mc.show_scores();

	} // i<=maxtrials_

	pose = mc.lowest_score_pose();

	TR << "DockDesignMover finished." << std::endl;

	// Output the final docked ligand as sdf, for convenience
	res_pos = find_design_position( pose );
	dump_molecule( pose.residue( res_pos ), "final" );

	// Stats -- should probably turn this into a full stats sub-object.
	for ( core::Size cc(1); cc <= chemistries_.size(); ++ cc ) {
		std::string const & name( chemistries_[cc]->name() );
		TR << "Stats for " << name << ": Attempts: " << stat_attempt[name] << " (" << stat_attempt[name]/float(maxtrials_) << ")"
			<< " Hard Fails: " << stat_failed[name] << " (" << stat_failed[name]/float(stat_attempt[name]) << ")"
			<< " Tested: " << stat_tested[name] << " (" << stat_tested[name]/float(stat_attempt[name]) << ")"
			<< " Accepted: " << stat_accepted[name] << " (" << stat_accepted[name]/float(stat_attempt[name]) << ")" << std::endl;
	}

	mc.show_scores();
	mc.show_state();
	mc.show_counters();


}// apply

/// @brief Return a unique string representing which job this is
std::string
DrugDesignMover::get_jobname() const {
	// TODO: This needs to be de-JD2-ized at some point.
	return protocols::jd2::JobDistributor::get_instance()->current_output_name();
}


bool // Returns true if we can't use this residuetype
DrugDesignMover::pre_process_restype(
	core::chemical::MutableResidueTypeOP restype, // Can/will be modified!
	core::chemical::IndexVDMapping & index_vd_mapping, // Can/will be modified!
	core::pose::Pose const & pose // Const reference pose for context
) {

	for ( core::Size bc(1); bc <= before_chemistries_.size(); ++bc ) {
		before_chemistries_[bc]->apply( *restype, pose ); // Allow context sensitive chemistries
		if ( before_chemistries_[bc]->get_last_status() != core::chemical::modifications::SUCCESS ) {
			TR.Warning << "Before chemistry " << before_chemistries_[bc]->name() << " reported failure on " << restype->name() << std::endl;
			return true;
		}
		index_vd_mapping = combine( index_vd_mapping, before_chemistries_[bc]->get_mapping() );
	}

	return false;
}

bool // Returns true if we can't use this residuetype
DrugDesignMover::post_process_restype(
	core::chemical::MutableResidueTypeOP restype, // Can/will be modified!
	core::chemical::IndexVDMapping & index_vd_mapping, // Can/will be modified!
	std::string const & new_name,
	core::pose::Pose const & pose // Const reference pose for context
) {
	using namespace core::chemical;

	for ( core::Size ac(1); ac <= after_chemistries_.size(); ++ac ) {
		after_chemistries_[ac]->apply( *restype, pose );
		if ( after_chemistries_[ac]->get_last_status() != core::chemical::modifications::SUCCESS ) {
			TR.Warning << "After chemistry " << after_chemistries_[ac]->name() << " reported failure on " << restype->name() << std::endl;
			return true;
		}
		index_vd_mapping = combine( index_vd_mapping, after_chemistries_[ac]->get_mapping() );
	}

	// Do now, because the after chemistries can possibly change the name
	if ( new_name.size() != 0 ) {
		restype->name( new_name );
	}

	// For placement into the pose, it's best to have the root atom be in the original residue type.
	if ( index_vd_mapping.reverse_lookup( restype->root_atom() ) == index_vd_mapping.invalid_key() ) {
		TR << "Icoord root placed on new atom - rerooting." <<std::endl;
		VD possible_root( index_vd_mapping.invalid_entry() );
		for ( core::chemical::IndexVDMapping::const_iterator itr( index_vd_mapping.begin() ), itr_end( index_vd_mapping.end() );
				itr != itr_end; ++ itr ) {
			// Possibly should have some sort of better central-tendancy type selection - right now just pick an arbitrary heavyatom
			if ( restype->has( itr->second ) && restype->atom( itr->second ).element() != core::chemical::element::H ) {
				possible_root = itr->second;
				break;
			}
		}
		if ( possible_root == index_vd_mapping.invalid_entry() ) {
			TR << "Unable to find acceptable root atom. Skipping rerooting." << std::endl; // Should have reasonable root already
		} else {
			TR << "Re-rooting on atom '" << restype->atom_name( possible_root ) << "'" << std::endl;
			core::chemical::reroot_restype( *restype, restype->graph(), possible_root );
		}
	}

	// Residues with less than three heavy atoms have issues with being properly placed in the pose
	// Register this as a fail
	if ( restype->nheavyatoms() < 3 ) {
		TR << "Chemistry resulted in ligand with less than three heavy atoms. Skipping." <<std::endl;
		dump_molecule( *restype, "hard_reject" );
		return true;
	}

	return false;
}


bool // Returns true if we can't use this residuetype
DrugDesignMover::emplace_residue_type(
	core::pose::Pose & pose, // Can/will be modified!
	core::chemical::MutableResidueTypeOP mutrestype, // Can/will be modified!
	core::chemical::IndexVDMapping const & index_vd_mapping // From original residue type in pose to restype
) {
	using namespace core::chemical;

	core::Size res_pos( find_design_position(pose) );

	// TODO: HACK: need to find a better way to store these.
	restypeset_->add_unpatchable_residue_type( mutrestype );
	core::chemical::ResidueTypeCOP restype( restypeset_->name_mapOP( mutrestype->name() ) );
	restypes_.push_back( restype );

	// Mapping from the indicies in the original residue type to the indicies in the new residue type
	IndexIndexMapping index_index_mapping( index_vd_mapping.downstream_combine( VDNameMapping( *mutrestype ) ).downstream_combine( NameIndexMapping( *restype ) ) );

	dump_molecule( pose.residue( res_pos ), "before_pose" );
	place_new_restype( pose, res_pos, *restype, index_index_mapping );
	dump_molecule( pose.residue( res_pos ), "after_pose" );

	if ( prefilter_ && ! prefilter_->apply( pose ) ) {
		TR << "Prefilter reported failure - rejecting." << std::endl;
		dump_molecule( pose.residue( res_pos ), "prefilter_fail" );
		return true;
	}

	// Apply the redocker
	dump_molecule( pose.residue( res_pos ), "before_dock" );
	redocker_->apply( pose );
	dump_molecule( pose.residue( res_pos ), "after_dock" );

	// Debugging - check to make sure that the we're not leaking constraints
	//pose.constraint_set()->show_definition(TR, pose);

	if ( postfilter_ && ! postfilter_->apply( pose ) ) {
		TR << "Postfilter reported failure - rejecting." << std::endl;
		dump_molecule( pose.residue( res_pos ), "postfilter_fail" );
		return true;
	}

	return false;
}

std::string
DrugDesignMover::get_name() const {
	return DrugDesignMover::mover_name();
}

/// @brief Is position n a valid designable residue according to the settings
bool
DrugDesignMover::check_design_position( Pose const & pose, Size n ) {
	using namespace core::pose;
	// Right now we're doing chain letters only - could be extended later
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) return false; // No chain info - not valid position.
	return pdb_info->chain(n) == chain_;
}

/// @brief Determine which residue on this pose should be designed
Size
DrugDesignMover::find_design_position( Pose const & pose ) {
	using namespace core::pose;
	// Right now we're doing chain letters only - could be extended later
	PDBInfoCOP pdb_info( pose.pdb_info() );
	if ( ! pdb_info ) {
		utility_exit_with_message("Pose does not have PDBInfo - can't determine chain to use.");
	}
	Size res_pos(0);
	for ( Size ii(1); ii <= pdb_info->nres(); ++ii ) {
		if ( pdb_info->chain(ii) == chain_ ) {
			if ( res_pos != 0 ) {
				TR.Warning << "Multiple residues are in designable chain "  << chain_ << " not only " << res_pos << " but also " << ii << std::endl;
			}
			res_pos = ii;
		}
	}
	if ( res_pos == 0 ) {
		utility_exit_with_message( std::string("No designable residue found. Missing chain ") + chain_ );
	}
	return res_pos;
}

/// @brief Find a new residue name, if the name already exists in the residue type set.
std::string
DrugDesignMover::find_new_res_name( std::string original_name, core::Size iteration, core::Size subiteration /*=0*/ ) const {
	// Default is just to append the iteration to the ligand name
	static const std::string letters( "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ" );
	std::string newresname( original_name + "-" + std::to_string(iteration) );
	if ( subiteration != 0 ) {
		newresname += letters[ (subiteration-1)%letters.size() ];
	}
	std::string jobname( get_jobname() );
	// First thing - we prepend the job name if it doesn't have it on the ligand name already.
	if ( ! utility::startswith( newresname, jobname ) ) {
		newresname = jobname + "-" + newresname;
	}
	if ( ! restypeset_->has_name( newresname ) ) {
		return newresname;
	}
	// Okay. Let's just increment a number until we get success
	core::Size mult(2);
	while ( restypeset_->has_name( newresname ) ) {
		newresname = newresname + "x" + std::to_string(mult);
		++mult;
	}
	return newresname;
}

/// @brief Dump the ligand to the given file (appending)
void
DrugDesignMover::dump_molecule( core::chemical::MutableResidueType const & restype, std::string const & stage) const {
	if ( ! debug_prefix_.empty() ) {
		std::string filename( debug_prefix_ + "_" + stage + "_" + get_jobname() + ".sdf" );
		core::chemical::sdf::MolWriter writer;
		utility::io::ozstream output( filename, std::ios_base::out | std::ios_base::app );
		writer.output_residue( output, restype );
		output.close();
	}
}

/// @brief Dump the ligand to the given file (appending)
void
DrugDesignMover::dump_molecule( core::conformation::Residue const & residue, std::string const & stage) const {
	if ( ! debug_prefix_.empty() ) {
		std::string filename( debug_prefix_ + "_" + stage + "_" + get_jobname() + ".sdf" );
		core::chemical::sdf::MolWriter writer;
		utility::io::ozstream output( filename, std::ios_base::out | std::ios_base::app );
		writer.output_residue( output, residue );
		output.close();
	}
}

std::string
DrugDesignMover::mover_name() {
	return "DrugDesignMover";
}

void
DrugDesignMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "trials", xsct_non_negative_integer, "The number of Monte Carlo trials", "10")
		+ XMLSchemaAttribute::attribute_w_default( "temperature", xsct_real, "The Monte Carlo temperature", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "chain", xsct_char, "The chain of the ligand to use", "X" )
		+ XMLSchemaAttribute( "prefilter", xs_string, "If set, apply this filter to the pose prior to the redocker mover." )
		+ XMLSchemaAttribute::attribute_w_default( "redocker", xs_string, "The name of the mover to use prior to scoring", "null_mover" )
		+ XMLSchemaAttribute( "postfilter", xs_string, "If set, apply this filter to the pose after the redocker mover has been applied." )
		+ XMLSchemaAttribute::required_attribute( "scorer", xs_string, "The name of the filter to use for scoring" )
		+ XMLSchemaAttribute( "debug_prefix", xs_string, "If set, use this prefix when dumping intermediate structures" );

	AttributeList add_subelement_attlist;
	add_subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "chemistry", xs_string, "The name of the chemistry to use" )
		+ XMLSchemaAttribute::attribute_w_default( "weight", xsct_real, "Weight for this chemistry when randomly picking.", "1.0" );

	AttributeList ba_subelement_attlist;
	ba_subelement_attlist
		+ XMLSchemaAttribute::required_attribute( "chemistry", xs_string, "The name of the chemistry to use" );

	XMLSchemaSimpleSubelementList subelements;
	subelements
		.add_simple_subelement( "Add", add_subelement_attlist, "Add a chemistry to the list of possibilities to randomly choose from" )
		.add_simple_subelement( "Before", ba_subelement_attlist, "Chemistry to always use before the randomly chosen one" )
		.add_simple_subelement( "After", ba_subelement_attlist, "Chemistry to always use after the randomly chosen one" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements(
		xsd, mover_name(),
		"Run a Monte Carlo drug design protocol on a given ligand.",
		attlist, subelements );

}

// @brief utility function to get chemistry from datamap
protocols::chemistries::ChemistryOP
DrugDesignMover::chemistry_from_subtag( utility::tag::TagCOP const subtag, basic::datacache::DataMap & data  ) const {
	std::string const chemistry_name( subtag->getOption< std::string >( "chemistry" ) );
	if ( ! data.has( "chemistry", chemistry_name ) ) {
		utility_exit_with_message( "Cannot find chemistry named " + chemistry_name + ". Was it defined in the <CHEMISTRY> tag?" );
	}
	return data.get_ptr< protocols::chemistries::Chemistry >( "chemistry", chemistry_name );
}

/// @brief parse xml file
void
DrugDesignMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data )
{
	maxtrials_ = tag->getOption< core::Size >( "trials", 10 );
	temperature_ = tag->getOption< Real >( "temperature", 1.0 );
	chain_ = tag->getOption< char >( "chain", 'X' );
	debug_prefix_ = tag->getOption< std::string >( "debug_prefix", "" );

	if ( ! tag->hasOption( "scorer" ) ) {
		utility_exit_with_message("You must provide a scorer option to the DrugDesignMover for scoring.");
	}
	std::string const filter_name( tag->getOption< std::string >( "scorer" ) );
	scorer( protocols::rosetta_scripts::parse_filter(filter_name, data) );
	std::string const mover_name( tag->getOption< std::string >( "redocker", "null_mover" ) );
	redocker( protocols::rosetta_scripts::parse_mover(mover_name, data) );

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

	for ( TagCOP const & subtag : tag->getTags() ) {
		if ( subtag->getName() == "Add" ) {
			core::Real const weight( subtag->getOption< core::Real >( "weight", 1.0 ) );
			add_chemistry( chemistry_from_subtag( subtag, data ) , weight ) ;
		} else if ( subtag->getName() == "Before" ) {
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
