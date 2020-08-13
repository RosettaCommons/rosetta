// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/ERRASER2Protocol.cc
/// @brief Run a single-threaded, checkpoint free, RosettaScripts accessible ERRASER2 job
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/rna/movers/ERRASER2Protocol.hh>
#include <protocols/rna/movers/ERRASER2ProtocolCreator.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

#include <protocols/viewer/viewers.hh>

#include <protocols/rna/movers/ErraserMinimizerMover.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.rna.movers.ERRASER2Protocol" );

namespace protocols {
namespace rna {
namespace movers {

using namespace core;
using namespace core::scoring;
using namespace core::pose;
using namespace core::pose::rna;
using namespace protocols;
using namespace protocols::stepwise::monte_carlo;
using namespace basic::options;
using namespace basic::options::OptionKeys;


/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
ERRASER2Protocol::ERRASER2Protocol():
	protocols::moves::Mover( ERRASER2Protocol::mover_name() )
{
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
ERRASER2Protocol::ERRASER2Protocol( ERRASER2Protocol const & src ):
	protocols::moves::Mover( src ),
	minimize_protein_( src.minimize_protein_ ),
	n_rounds_( src.n_rounds_ )
{
	if ( src.scorefxn_ ) scorefxn_ = src.scorefxn_->clone();
	if ( src.rebuild_residue_selector_ ) rebuild_residue_selector_ = src.rebuild_residue_selector_->clone();
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
ERRASER2Protocol::~ERRASER2Protocol(){}


// AMW TODO:
// Goal methods:
// 1. specify a range of residues to rebuild
// 2. specify a range/window of PROTEIN residues to minimize/repack (using residue selector?)



//OPT_KEY( Integer, rounds )
//OPT_KEY( Boolean, minimize_protein )

static basic::Tracer TR( "apps.public.rna.erraser.erraser2" );


utility::vector1< core::Size >
all_pose_residues( core::pose::Pose const & pose ) {
	utility::vector1< core::Size > foo;
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) foo.emplace_back( ii );
	return foo;
}

void show_accuracy_report( pose::Pose const & start_pose, std::string const & tag, core::Size const round ) {
	RNA_SuiteName suite_namer;

	utility::vector1< core::Size > bad_suites;
	for ( core::Size ii = 1; ii <= start_pose.size(); ++ii ) {
		// This should pick up bad suites
		auto suite_assignment = suite_namer.assign( start_pose, ii );
		if ( suite_assignment.name == "!!" ) bad_suites.push_back( ii );
	}

	// TODO: don't output round if # is zero
	// TODO: don't output bad suites if empty vector
	// TODO: output other things too.
	TR << tag << " " << round << ": " << bad_suites << std::endl;
}



bool atoms_have_bond_to_bonded_atoms( pose::Pose const & pose, core::Size const ai, core::Size const ii, core::Size const aj, core::Size const jj ) {
	// Loop through atoms in both residues
	// Accumulate set of atom IDs bonded to each atom.
	// Check if there are bonds between any pairs.

	utility::vector1< core::id::AtomID > bonded_to_ai, bonded_to_aj;

	for ( core::Size ak = 1; ak <= pose.residue( ii ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, ii ) ) ) {
			bonded_to_ai.emplace_back( ak, ii );
		}
		if (  pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, ii ) ) ) {
			bonded_to_aj.emplace_back( ak, ii );
		}
	}
	for ( core::Size ak = 1; ak <= pose.residue( jj ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, jj ) ) ) {
			bonded_to_ai.emplace_back( ak, jj );
		}
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, jj ) ) ) {
			bonded_to_aj.emplace_back( ak, jj );
		}
	}

	for ( auto const & b_to_ai : bonded_to_ai ) {
		for ( auto const & b_to_aj : bonded_to_aj ) {
			if ( pose.conformation().atoms_are_bonded( b_to_ai, b_to_aj ) ) return true;
		}
	}

	return false;
}







bool atoms_have_mutual_bond_to_atom( pose::Pose const & pose, core::Size const ai, core::Size const ii, core::Size const aj, core::Size const jj ) {
	// Loop through atoms in both residues
	// Assumes no 1-atom polymeric residues. That's ok.
	for ( core::Size ak = 1; ak <= pose.residue( ii ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, ii ) )
				&& pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, ii ) ) ) return true;
	}
	for ( core::Size ak = 1; ak <= pose.residue( jj ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, jj ) )
				&& pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, jj ) ) ) return true;
	}

	return false;
}

bool
bump_check( pose::Pose const & pose, core::Size const ii ) {
	for ( core::Size jj = 1; jj <= pose.size(); ++jj ) {
		if ( jj == ii ) continue;
		if ( jj == ii + 1 ) continue;
		if ( jj == ii - 1 ) continue;

		// If first atoms are more than 10A apart no clashes
		if ( pose.residue( ii ).xyz( 1 ).distance_squared( pose.residue( jj ).xyz( 1 ) ) > 100 ) continue;

		for ( core::Size ai = 1; ai <= pose.residue( ii ).natoms(); ++ai ) {
			for ( core::Size aj = 1; aj <= pose.residue( jj ).natoms(); ++aj ) {
				// If the two atoms are bonded ignore.
				if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( aj, jj ) ) ) continue;
				if ( atoms_have_mutual_bond_to_atom( pose, ai, ii, aj, jj ) ) continue;
				if ( atoms_have_bond_to_bonded_atoms( pose, ai, ii, aj, jj ) ) continue;

				Real ai_r = pose.residue_type( ii ).atom_type( ai ).lj_radius();
				Real aj_r = pose.residue_type( jj ).atom_type( aj ).lj_radius();

				if ( pose.residue( ii ).xyz( ai ).distance_squared( pose.residue( jj ).xyz( aj ) ) <
						( aj_r + ai_r - 0.40 ) * ( aj_r + ai_r - 0.40 ) ) {
					//( aj_r + ai_r - 0.35 ) * ( aj_r + ai_r - 0.35 ) ) {
					TR.Trace << "Bump check failed because residue " << pose.pdb_info()->chain( ii ) << ":" << pose.pdb_info()->number( ii ) << " atom " << pose.residue( ii ).atom_name( ai ) << " (lj_radius " << ai_r << ") is within "
						<< ( aj_r + ai_r - 0.40 ) * ( aj_r + ai_r - 0.40 ) << " of residue " << pose.pdb_info()->chain( jj ) << ":" << pose.pdb_info()->number( jj ) << " atom " << pose.residue( jj ).atom_name( aj ) << " (lj_radius " << aj_r << ") -- " << pose.residue( ii ).xyz( ai ).distance_squared( pose.residue( jj ).xyz( aj ) ) << std::endl;
					return true;
				}
			}
		}
	}
	return false;
}


utility::vector1< core::Size >
determine_residues_to_rebuild(
	utility::vector1< core::Size > const & consider_rebuild_res, pose::Pose const & start_pose
) {
	/*if ( utility::file::file_exists( residue_rebuild_log_namer( resample_round, nstruct ) ) ) {
	// process it
	return residues_to_rebuild_log( residue_rebuild_log_namer( resample_round, nstruct ) );
	}*/

	RNA_SuiteName suite_namer;
	utility::vector1< core::Size > residues_to_rebuild;

	using namespace basic::options;
	using namespace core::pose::full_model_info;

	for ( core::Size ii = 1; ii <= start_pose.size(); ++ii ) {
		// Skip any residues that aren't even under consideration
		if ( ! consider_rebuild_res.contains( ii ) ) continue;

		//if ( !option[ minimize_protein ].value() && start_pose.residue_type( ii ).is_protein() ) continue;
		// can minimize protein but can't rebuild it
		if ( start_pose.residue_type( ii ).is_protein() ) continue;

		// This should pick up bad torsions
		Real const torsion_score = start_pose.energies().residue_total_energies( ii )[ rna_torsion ];
		if ( torsion_score > 1.2 ) {
			TR.Debug << "Rebuilding residue " << ii <<" for a bad torsion score of " << torsion_score << std::endl;
			residues_to_rebuild.push_back( ii );
		}

		if ( start_pose.residue_type( ii ).is_RNA() && start_pose.residue_type( ii ).is_canonical_nucleic() ) {
			residues_to_rebuild.push_back( ii );
		}

		// Bad suites, if it hasn't been added already
		// TODO: chain aware
		auto suite_assignment = suite_namer.assign( start_pose, ii );
		if ( suite_assignment.name == "!!" ) {
			TR.Debug << "Rebuilding residue " << ii << " with its neighbors as possible for a bad suite" << std::endl;
			if ( ii > 1 ) {
				residues_to_rebuild.push_back( ii - 1 );
			}
			residues_to_rebuild.push_back( ii );
			if ( ii < start_pose.size() ) {
				residues_to_rebuild.push_back( ii + 1 );
			}
		}

		// This should pick up bad bond lengths and angles.
		// (If cart_bonded is on...)
		Real const cart = start_pose.energies().residue_total_energies( ii )[ cart_bonded ];
		if ( cart > 2 ) {
			TR.Debug << "Rebuilding residue " << ii <<" for a bad cart_bonded score of " << cart << std::endl;
			residues_to_rebuild.push_back( ii );
		}


		// Clashes
		if ( bump_check( start_pose, ii ) ) {
			TR.Debug << "Rebuilding residue " << ii <<" for a bad bump check" << std::endl;
			residues_to_rebuild.push_back( ii );
		}

		// Make sure bad glycosylation sites get picked up here. Need samplers.

		// Make sure branched nucleotides are always resampled.

		// Make sure ions get refined... have to think about how to handle them, water coordination, etc.

		// (at first just Mg?)

		// Make sure cyclic dinucleotide ligands are redocked. Also JUMP_DOCK.
		if ( !start_pose.residue( ii ).is_polymer() ) residues_to_rebuild.push_back( ii );
	}

	// Shuffle for this nstruct. This is your only chance!
	auto tmpset = std::set< core::Size >( residues_to_rebuild.begin(), residues_to_rebuild.end() );
	residues_to_rebuild = utility::vector1< core::Size >( tmpset.begin(), tmpset.end() );
	numeric::random::random_permutation( residues_to_rebuild, numeric::random::rg() );
	TR << "Total to resample this round: " << residues_to_rebuild << std::endl;

	return residues_to_rebuild;
}

mover::StepWiseMasterMover
ERRASER2Protocol::configure_master_mover( Pose const & start_pose ) {
	// initialize options
	options::StepWiseMonteCarloOptionsOP options( new options::StepWiseMonteCarloOptions );
	options->initialize_from_command_line();
	options->set_erraser( true );
	options->set_enumerate( true );
	options->set_skip_deletions( true );
	options->set_output_minimized_pose_list( false );
	options->set_force_moving_res_for_erraser( true );
	options->set_boltzmann_choice_post_enumerated_minimize( true );

	// run StepWiseMasterMover::resample_full_model
	mover::StepWiseMasterMover master_mover( scorefxn_, options );
	// This is essential for rmsd_screen
	master_mover.set_native_pose( utility::pointer::make_shared< pose::Pose >( start_pose ) );
	return master_mover;
}

void
ERRASER2Protocol::resample_full_model( core::Size const resample_round, pose::Pose & start_pose,
	utility::vector1< core::Size > const & consider_rebuild_res, core::Size const nstruct
) {

	// Determine residues to rebuild -- shuffled order, but who knows the checkpointing situation.
	// TODO: if this is slow, figure out a way to skip it if there's a checkpoint file.
	utility::vector1< core::Size > residues_to_rebuild = determine_residues_to_rebuild( consider_rebuild_res, start_pose );

	// OK, what residues might these be? What might be connected?
	for ( auto const elem : residues_to_rebuild ) {
		TR.Debug << elem << ": " << " " << start_pose.pdb_info()->chain( elem ) << start_pose.pdb_info()->number( elem ) << " " << start_pose.residue_type( elem ).name() << std::endl;
	}
	//return;

	pose::Pose seq_rebuild_pose = start_pose;
	protocols::viewer::add_conformation_viewer(seq_rebuild_pose.conformation(), "current", 500, 500);

	mover::StepWiseMasterMover master_mover = configure_master_mover( start_pose );

	master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/,
		residues_to_rebuild, resample_round, nstruct );

	// score seq_rebuild pose
	( *scorefxn_ )( seq_rebuild_pose );

	// show scores of start_pose and full_model_pose
	//if ( option[ show_scores ]() ) {
	TR << "\n Score before seq_rebuild:" << std::endl;
	scorefxn_->show( TR, start_pose );
	TR << "\n Score after seq_rebuild:" << std::endl;
	scorefxn_->show( TR, seq_rebuild_pose );
	//}

	// dump pdb -- should figure out file name based on inputs
	//seq_rebuild_pose.dump_pdb( "SEQ_REBUILD.pdb" );

	// write out to silent file
	//std::string tag = "S_0";
	//stepwise::monte_carlo::output_to_silent_file( tag, silent_file_out, full_model_pose );
	start_pose = seq_rebuild_pose;
	protocols::viewer::clear_conformation_viewers();
}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////


///////////////////////////////////////////////////////////////////////////////
/// @brief Perform the erraser2 protocol in a working directory
/// @details This function also does pose readin (because it could be resuming
/// from checkpoints). In MPI contexts, this function is run on every core. (We
/// use a work-partition job distribution framework by default, where each core
/// is just assigned a fixed set of nstruct.) As a result, it needs to take two
/// ints, which represent its rank and the total number of processors. Since MPI
/// libraries use signed ints, we follow their convention rather than contorting
/// ourselves to make them Sizes.

/// @brief Apply the mover
void
ERRASER2Protocol::apply( core::pose::Pose & start_pose ) {

	// Outline:
	// 1. Read in a pose.
	// 2. Minimize it (using the ErraserMinimizerMover)
	// 3. Rebuild residues that moved a lot plus geometric outliers
	//   a. implement via resample_full_model with a custom sample res
	// 4. [wash again if desired]
	// 5. Minimize it

	// Get rebuild res by applying the residue selector to the pose.
	auto selection = rebuild_residue_selector_ ? rebuild_residue_selector_->apply( start_pose ) : utility::vector1< bool >( start_pose.size(), true );

	utility::vector1< core::Size > consider_rebuild_res;
	for ( core::Size ii = 1; ii <= selection.size(); ++ii ) {
		if ( selection[ ii ] ) consider_rebuild_res.push_back( ii );
	}

	( *scorefxn_ )( start_pose );
	show_accuracy_report( start_pose, "Start", 0/*, TR*/ );

	ErraserMinimizerMover erraser_minimizer;
	// if true: solely responsible for all chunks
	// if false: node 0 coordinates, nodes 1-n min chunks.
	erraser_minimizer.work_partition( true );
	erraser_minimizer.scorefxn( scorefxn_ );
	erraser_minimizer.edens_scorefxn( scorefxn_ );
	erraser_minimizer.minimize_protein( minimize_protein_ /* AMW TODO: change structure */ );

	//using namespace core::scoring::constraints;
	//using namespace core::scoring::func;

	TR << "First pass preminimization..." << std::endl;
	erraser_minimizer.constrain_phosphate( true );
	erraser_minimizer.apply( start_pose );

	for ( core::Size ii = 1; ii <= n_rounds_; ++ii ) {
		erraser_minimizer.constrain_phosphate( false );
		erraser_minimizer.apply( start_pose );
		show_accuracy_report( start_pose, "Minimized", ii/*, TR*/ );

		resample_full_model( ii, start_pose, consider_rebuild_res, 1 );
		show_accuracy_report( start_pose, "Resampled", ii/*, TR*/ );
	}

	erraser_minimizer.apply( start_pose );
	show_accuracy_report( start_pose, "FINAL", 0/*, TR*/ );
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
ERRASER2Protocol::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
ERRASER2Protocol::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	if ( !scorefxn_ && basic::options::option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_function( core::scoring::get_score_function() );
	}
	if ( !scorefxn_ ) {
		// scorefunction unconfigured!
		utility_exit_with_message( "You did not configure the scoring function in tag or -score:weights!" );
	}

	debug_assert( scorefxn_ );

	// AMW TODO: customize the whole structure of the protocol -- let the user
	// 'min, resample, min, min, resample, resample' if they want
	n_rounds( tag->getOption< core::Size >( "n_rounds", true ) );

	// To be rebuilt, residues must be malformed AND in rebuild_res_selector
	if ( tag->hasOption("rebuild_res_selector") ) {
		std::string const selector_name ( tag->getOption< std::string >( "rebuild_res_selector" ) );
		if ( TR.visible() ) TR << "Set rebuild res selector name to " << selector_name << "." << std::endl;
		core::select::residue_selector::ResidueSelectorCOP residue_selector;
		try {
			residue_selector = core::select::residue_selector::get_residue_selector(selector_name, data);
		} catch ( utility::excn::Exception & e ) {
			std::string error_message = "Failed to find ResidueSelector named '" + selector_name + "' from the DataMap from ERRASER2Protocol::parse_tag()\n" + e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_message );
		}
		runtime_assert( residue_selector );
		rebuild_residue_selector( residue_selector );
	}
}

void ERRASER2Protocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "selector" );
	rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute("n_rounds", xsct_non_negative_integer, "Number of ERRASER2 minimize-rebuild rounds." );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Run the ERRASER2 protocol", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ERRASER2Protocol::fresh_instance() const
{
	return utility::pointer::make_shared< ERRASER2Protocol >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
ERRASER2Protocol::clone() const
{
	return utility::pointer::make_shared< ERRASER2Protocol >( *this );
}

std::string ERRASER2Protocol::get_name() const {
	return mover_name();
}

std::string ERRASER2Protocol::mover_name() {
	return "ERRASER2Protocol";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
ERRASER2ProtocolCreator::create_mover() const
{
	return utility::pointer::make_shared< ERRASER2Protocol >();
}

std::string
ERRASER2ProtocolCreator::keyname() const
{
	return ERRASER2Protocol::mover_name();
}

void
ERRASER2ProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ERRASER2Protocol::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, ERRASER2Protocol const & mover )
{
	mover.show(os);
	return os;
}

} //movers
} //rna
} //protocols
