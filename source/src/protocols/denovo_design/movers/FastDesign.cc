// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/FastDesign.cc
/// @brief Similar to FastRelax, but with design and a few additional options.
/// @details Carries out alternating rounds of design (packing) and minimization, ramping up the repulsive term from a very low value from one iteration to the next.
/// @author Tom Linsky
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added support for D-amino acids.

//Unit Headers
#include <protocols/denovo_design/movers/FastDesign.hh>
#include <protocols/denovo_design/movers/FastDesignCreator.hh>

//Protocol Headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/task_operations/LimitAromaChi2Operation.hh>

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>

//Basic Headers
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>

//Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// ObjexxFCL Headers
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#include <basic/datacache/DataMap.hh> // AUTO IWYU For DataMap
#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask

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

static basic::Tracer TR( "protocols.denovo_design.movers.FastDesign" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace movers {
////////////////////////////////////////////////////////////////////////////////////////////////////



///  ---------------------------------------------------------------------------------
///  FastDesign main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
FastDesign::FastDesign():
	FastRelax(),
	clear_designable_residues_( false ),
	run_count_( 0 ),
	cgs_()
{
	set_enable_design( true );
	//read_script_file( "", default_repeats_ );
}

/// @brief Constructor with some options
///
FastDesign::FastDesign(
	core::scoring::ScoreFunctionOP scorefxn_in,
	core::Size standard_repeats
) :
	FastRelax( scorefxn_in, standard_repeats ),
	clear_designable_residues_(false),
	run_count_(0),
	cgs_()
{
	set_enable_design( true );
}

/// @brief Constructor with some options
///
FastDesign::FastDesign(
	core::scoring::ScoreFunctionOP scorefxn_in,
	std::string const & script_file
) :
	FastRelax( scorefxn_in, script_file ),
	clear_designable_residues_(false),
	run_count_(0),
	cgs_()
{
	set_enable_design( true );
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FastDesign::~FastDesign() = default;



/// @brief return an blank version of ourselves
protocols::moves::MoverOP
FastDesign::fresh_instance() const
{
	return utility::pointer::make_shared< FastDesign >();
}

/// Return a copy of ourselves
protocols::moves::MoverOP
FastDesign::clone() const
{
	return utility::pointer::make_shared< FastDesign >( *this );
}

/// @brief Create the default task factory.  Must be called before design can occur.
void
FastDesign::set_up_default_task_factory() {
	set_task_factory( create_default_task_factory() );
}


void
FastDesign::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	// make sure we create a task factory before parsing FastRelax::parse_my_tag
	// otherwise no design will occur
	set_up_default_task_factory();
	FastRelax::parse_my_tag( tag, data );

	set_clear_designable_residues( tag->getOption< bool >( "clear_designable_residues", clear_designable_residues_ ) );

	// parse constraint generators
	utility::vector1< std::string > const cgs = utility::string_split( tag->getOption< std::string >( "cgs", "" ), ',' );
	for ( auto const & cg : cgs ) {
		if ( cg.empty() ) continue;
		protocols::constraint_generator::ConstraintGeneratorCOP new_cg =
			data.get_ptr< protocols::constraint_generator::ConstraintGenerator const >( "CONSTRAINT_GENERATORS", cg );
		if ( !new_cg ) {
			std::stringstream msg;
			msg << "FastDesign: Could not find a constraint generator named " << cg << " in the data map.  Ensure it has been defined in an AddConstraints mover before being referenced by FastDesign."
				<< std::endl;
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
		}
		add_constraint_generator( new_cg );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
FastDesign::apply( core::pose::Pose & pose )
{
	// increment run counter
	++run_count_;

	TR.Debug   << "========================== FastDesign Run: " << run_count_ << " ===================" << std::endl;
	reset_status();

	// create and print task
	core::pack::task::PackerTaskOP clear_task = get_task_factory()->create_task_and_apply_taskoperations( pose );
	clear_task->show( TR );

	// print movemap
	if ( get_movemap() && TR.Debug.visible() ) {
		TR.Debug << "Movemap: ";
		core::kinematics::MoveMapOP temp_movemap = get_movemap()->clone();
		initialize_movemap( pose, *temp_movemap );
		temp_movemap->show( TR.Debug );
		TR.Debug << std::endl;
	}

	// if requested, clear all residues marked as designable
	if ( clear_designable_residues_ ) {
		TR << "Clearing designable residues...";
		for ( core::Size resid=1; resid<=pose.size(); ++resid ) {
			if ( pose.residue( resid ).is_protein() && clear_task->being_designed( resid ) ) {
				if ( pose.residue( resid ).type().is_alpha_aa() ) {
					if ( pose.residue( resid ).type().is_l_aa() ) { //Note that glycine is skipped, since it's not an L-amino acid.
						// skip proline
						if ( pose.residue( resid ).type().aa() != core::chemical::aa_pro ) {
							TR << resid << "...";
							// mutate to alanine
							protocols::simple_moves::MutateResidue mut_res( resid, core::chemical::aa_ala );
							mut_res.apply( pose );
						} else {
							TR << "<" << resid << ">...";
						}
					} else if ( pose.residue( resid ).type().is_d_aa() ) {
						// skip D-proline
						if ( pose.residue( resid ).type().aa() != core::chemical::aa_dpr ) {
							TR << resid << "...";
							// mutate to D-alanine
							protocols::simple_moves::MutateResidue mut_res( resid, core::chemical::aa_dal );
							mut_res.apply( pose );
						} else {
							TR << "<" << resid << ">...";
						}
					} else { //Achiral case
						TR << "<" << resid << ">...";
					}
				}
			}
		}
		TR << std::endl;
	}

	FastRelax::apply( pose );

	// show final scores
	core::scoring::ScoreFunctionOP local_sfxn( get_scorefxn()->clone() );
	local_sfxn->show( TR, pose );
	TR.flush();
}

/// @brief sets constraint weights -- used with constraint ramping
void
FastDesign::set_constraint_weight(
	core::scoring::ScoreFunctionOP local_scorefxn,
	core::scoring::EnergyMap const & full_weights,
	core::Real const weight,
	core::pose::Pose & pose ) const
{
	runtime_assert( local_scorefxn != nullptr );
	if ( cgs_.size() ) {
		protocols::constraint_generator::RemoveConstraints rm_csts( cgs_, true );

		try {
			rm_csts.apply( pose );
		} catch ( protocols::constraint_generator::EXCN_RemoveCstsFailed const & e ) {
			// if removing constraints fails, we don't really care
			// they are only being removed to clean up before re-adding them below
		}

		if ( weight > 1e-6 ) {
			protocols::constraint_generator::AddConstraints( cgs_ ).apply( pose );
		}

	} else {
		local_scorefxn->set_weight( core::scoring::coordinate_constraint, full_weights[ core::scoring::coordinate_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::atom_pair_constraint, full_weights[ core::scoring::atom_pair_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::angle_constraint, full_weights[ core::scoring::angle_constraint ] * weight );
		local_scorefxn->set_weight( core::scoring::dihedral_constraint, full_weights[ core::scoring::dihedral_constraint ] * weight );
		TR << "[coordinate:atom_pair:angle:dihedral] = "
			<< local_scorefxn->get_weight( core::scoring::coordinate_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::atom_pair_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::angle_constraint )
			<< " : " << local_scorefxn->get_weight( core::scoring::dihedral_constraint ) << std::endl;
	}
}

core::pack::task::TaskFactoryOP
FastDesign::create_default_task_factory() const
{
	using core::pack::task::operation::TaskOperationCOP;

	core::pack::task::TaskFactoryOP local_tf( new core::pack::task::TaskFactory() );
	if ( get_task_factory() ) {
		local_tf = get_task_factory()->clone();
	} else {
		// add command line things
		local_tf->push_back(utility::pointer::make_shared< core::pack::task::operation::InitializeFromCommandline >());

		// add resfile if requested by flags
		if ( basic::options::option[ basic::options::OptionKeys::relax::respect_resfile]() &&
				basic::options::option[ basic::options::OptionKeys::packing::resfile].user() ) {
			local_tf->push_back(utility::pointer::make_shared< core::pack::task::operation::ReadResfile >());
			TR << "Using Resfile for packing step." << std::endl;
		}

		// turn off packing on residues that can't move
		// (Note that this function *doesn't* see any MoveMap info that's set via parse_my_tag().)
		if ( get_movemap() ) {
			core::pack::task::operation::PreventRepackingOP turn_off_packing( new core::pack::task::operation::PreventRepacking() );
			core::kinematics::MoveMap const & mm = *get_movemap();
			for ( auto mm_torsion=mm.movemap_torsion_id_begin();
					mm_torsion!=mm.movemap_torsion_id_end(); ++mm_torsion ) {
				// skip if this residue and torsion are allowed to move
				if ( mm_torsion->second ) continue;
				core::kinematics::MoveMap::MoveMapTorsionID const & torsion_id = mm_torsion->first;

				// skip if this is not chi
				if ( torsion_id.second != core::id::CHI ) continue;

				turn_off_packing->include_residue( torsion_id.first );
			}
			local_tf->push_back( turn_off_packing );
		}
	}

	//Include current rotamer by default - as before.
	local_tf->push_back( utility::pointer::make_shared< core::pack::task::operation::IncludeCurrent >() );

	if ( limit_aroma_chi2() ) {
		local_tf->push_back(utility::pointer::make_shared< protocols::task_operations::LimitAromaChi2Operation >());
	}

	return local_tf;
}

void
FastDesign::add_constraint_generator( protocols::constraint_generator::ConstraintGeneratorCOP generator )
{
	cgs_.push_back( generator );
}

void
FastDesign::clear_constraint_generators()
{
	cgs_.clear();
}

std::string FastDesign::get_name() const {
	return mover_name();
}

std::string FastDesign::mover_name() {
	return "FastDesign";
}

void FastDesign::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaComplexTypeGeneratorOP ct_gen = FastRelax::complex_type_generator_for_fast_relax( xsd );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "clear_designable_residues", xsct_rosetta_bool, "Clear the set of designable residues?")
		+ XMLSchemaAttribute( "cgs", xs_string, "Names of previously specified constraint generators to apply during this run");

	ct_gen->add_attributes( attlist )
		.element_name( mover_name() )
		.description( "FastRelax mover used for design that can take constraint generators" )
		.write_complex_type_to_schema( xsd );


}

/// @brief Provide the citation.
void
FastDesign::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	basic::citation_manager::CitationManager * cm( basic::citation_manager::CitationManager::get_instance() );
	basic::citation_manager::CitationCollectionOP collection( utility::pointer::make_shared< basic::citation_manager::CitationCollection >( get_name(), basic::citation_manager::CitedModuleType::Mover ) );
	collection->add_citation( cm->get_citation_by_doi( "10.1038/nature19791" ) );
	collection->add_citation( cm->get_citation_by_doi( "10.1002/prot.26030" ) );

	citations.add( collection );
	RelaxProtocolBase::provide_citation_info( citations );
}

std::string FastDesignCreator::keyname() const {
	return FastDesign::mover_name();
}

protocols::moves::MoverOP
FastDesignCreator::create_mover() const {
	return utility::pointer::make_shared< FastDesign >();
}

void FastDesignCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FastDesign::provide_xml_schema( xsd );
}


} // namespace movers
} // namespace denovo_design
} // namespace protocols
