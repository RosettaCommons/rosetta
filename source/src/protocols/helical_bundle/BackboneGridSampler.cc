// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BackboneGridSampler.cc
/// @brief  This mover samples conformations of a repeating chain of a residue type by grid-sampling
/// mainchain torsion values and setting all residues in a range to have the same mainchain torsion values.
/// @details  Note that this mover throws away the input pose and generates new geometry.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/helical_bundle/BackboneGridSampler.hh>
#include <protocols/helical_bundle/BackboneGridSamplerCreator.hh>
#include <protocols/helical_bundle/BackboneGridSamplerHelper.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/Energies.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

#include <numeric/constants.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <numeric/random/random.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/NamedAtomID.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/VariantType.hh>

//JD2:
#include <protocols/jd2/util.hh>

// Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>

// C output headers:
#include <cstdio>
#include <sstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

using basic::Error;
using basic::Warning;

namespace protocols {
namespace helical_bundle {

static basic::Tracer TR("protocols.helical_bundle.BackboneGridSampler");
static basic::Tracer TR_Results("protocols.helical_bundle.BackboneGridSampler.Results");

// XRW TEMP std::string
// XRW TEMP BackboneGridSamplerCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return BackboneGridSampler::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP BackboneGridSamplerCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new BackboneGridSampler );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BackboneGridSampler::mover_name()
// XRW TEMP {
// XRW TEMP  return "BackboneGridSampler";
// XRW TEMP }

///
/// @brief Creator for BackboneGridSampler mover.
BackboneGridSampler::BackboneGridSampler():
	Mover("BackboneGridSampler"),
	nstruct_mode_(false),
	nstruct_mode_repeats_(1),
	select_low_(true),
	max_samples_( 10000 ),
	pre_scoring_mover_(),
	pre_scoring_mover_exists_(false),
	pre_scoring_filter_(),
	pre_scoring_filter_exists_(false),
	dump_pdbs_(false),
	pdb_prefix_("bbs_out"),
	sfxn_set_(false),
	sfxn_(),
	residues_per_repeat_(1),
	torsions_to_sample_(),
	torsions_to_fix_(),
	nres_(20),
	resname_(),
	cap_ends_(false),
	peptide_stub_mover_( new PeptideStubMover ),
	peptide_stub_mover_initialized_(false)
{
	resname_.push_back("ALA");
	torsions_to_fix_.resize(1);
	torsions_to_sample_.resize(1);
}

///
/// @brief Copy constructor for BackboneGridSampler mover.
BackboneGridSampler::BackboneGridSampler( BackboneGridSampler const & src ):
	protocols::moves::Mover( src ),
	nstruct_mode_(src.nstruct_mode_),
	nstruct_mode_repeats_(src.nstruct_mode_repeats_),
	select_low_(src.select_low_),
	max_samples_(src.max_samples_),
	pre_scoring_mover_( src.pre_scoring_mover_ ), //NOTE that we're not cloning this mover, but using it straight
	pre_scoring_mover_exists_(src.pre_scoring_mover_exists_),
	pre_scoring_filter_( src.pre_scoring_filter_ ), //NOTE that we're not cloning this filter, but using it straight
	pre_scoring_filter_exists_(src.pre_scoring_filter_exists_),
	dump_pdbs_(src.dump_pdbs_),
	pdb_prefix_(src.pdb_prefix_),
	sfxn_set_(src.sfxn_set_),
	sfxn_(src.sfxn_), //NOTE that this is also copied without cloning
	residues_per_repeat_( src.residues_per_repeat_ ),
	torsions_to_sample_( src.torsions_to_sample_ ),
	torsions_to_fix_(src.torsions_to_fix_),
	nres_(src.nres_),
	resname_(src.resname_),
	cap_ends_(src.cap_ends_),
	peptide_stub_mover_( utility::pointer::static_pointer_cast< PeptideStubMover >(src.peptide_stub_mover_->clone()) ), //CLONE this mover
	peptide_stub_mover_initialized_( src.peptide_stub_mover_initialized_ )
{
}

///
/// @brief Destructor for BackboneGridSampler mover.
BackboneGridSampler::~BackboneGridSampler() = default;

///
/// @brief Clone operator to create a pointer to a fresh BackboneGridSampler object that copies this one.
protocols::moves::MoverOP BackboneGridSampler::clone() const {
	return protocols::moves::MoverOP( new BackboneGridSampler( *this ) );
}

///
/// @brief Fresh_instance operator to create a pointer to a fresh BackboneGridSampler object that does NOT copy this one.
protocols::moves::MoverOP BackboneGridSampler::fresh_instance() const {
	return protocols::moves::MoverOP( new BackboneGridSampler );
}

////////////////////////////////////////////////////////////////////////////////
//          APPLY FUNCTION                                                    //
////////////////////////////////////////////////////////////////////////////////

///
/// @brief Actually apply the mover to the pose.
void BackboneGridSampler::apply (core::pose::Pose & pose)
{
	using namespace core::chemical;
	using namespace core::id;

	//Calculate the total number of samples:
	core::Size const total_samples( calculate_total_samples() );
	if ( TR.visible() ) { TR << "Starting BackboneGridSampler.  A total of " << total_samples << " mainchain torsion conformations will be sampled." << std::endl; TR.flush(); }

	//Check that the total samples is not greater than the maximum.
	if ( total_samples > max_samples() ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::apply():  The total number of samples exceeds the maximum allowed.  (Note that you can increase the maximum allowed to suppress this error.)\n" );
	}

	//Check that a scorefunction has been set:
	if ( !sfxn_ ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::apply():  No scorefunction has been set for this mover!\n" );
	}

	//Create and initialize the helper object that will keep track of all of the sampling that we want to do:
	BackboneGridSamplerHelperOP helper( new BackboneGridSamplerHelper );
	helper->initialize_data( torsions_to_sample_ );

	//Build the pose
	if ( TR.visible() ) {
		TR << "Building " << nres() << "-repeat pose where each repeating unit consists of ";
		for ( core::Size i=1, imax=residues_per_repeat(); i<=imax; ++i ) {
			TR << resname(i);
			if ( i < imax ) TR << ", ";
		}
		TR << "." << std::endl;
		TR.flush();
	}
	core::pose::PoseOP newpose(new core::pose::Pose);
	if ( !peptide_stub_mover_initialized() ) set_up_peptide_stub_mover();
	peptide_stub_mover_->apply( *newpose );

	//Add termini:
	if ( cap_ends() ) {
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_LOWER, newpose->size() );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_UPPER, newpose->size() );
		core::pose::add_variant_type_to_pose_residue(*newpose, core::chemical::ACETYLATED_NTERMINUS_VARIANT, 1);
		core::pose::add_variant_type_to_pose_residue(*newpose, core::chemical::METHYLATED_CTERMINUS_VARIANT, newpose->size());
	} else {
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_UPPER, 1 );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_LOWER, 1 );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_LOWER, newpose->size() );
		core::pose::remove_variant_type_from_pose_residue(*newpose, core::chemical::CUTPOINT_UPPER, newpose->size() );
		core::pose::add_variant_type_to_pose_residue(*newpose, core::chemical::LOWER_TERMINUS_VARIANT, 1);
		core::pose::add_variant_type_to_pose_residue(*newpose, core::chemical::UPPER_TERMINUS_VARIANT, newpose->size());
	}

	//Do initial checks:
	//Are all of the torsions to fix in the current residue type?
	for ( core::Size ires=1, iresmax=residues_per_repeat(); ires<=iresmax; ++ires ) { //Loop through all residues in the repeating unit.
		for ( core::Size i=1, imax=torsions_to_fix_[ires].size(); i<=imax; ++i ) {
			if ( torsions_to_fix_[ires][i].first > newpose->residue(ires).mainchain_torsions().size() ) {
				utility_exit_with_message(
					"In protocols::helical_bundle::BackboneGridSampler::apply():  Mainchain torsion angles to fix were specified that are not in the " + resname(ires) + " residue type.\n"
				);
			}
		}
		//Are all of the torsions to sample in the current residue type?
		for ( core::Size i=1, imax=torsions_to_sample_[ires].size(); i<=imax; ++i ) {
			if ( torsions_to_sample_[ires][i].first > newpose->residue(ires).mainchain_torsions().size() ) {
				utility_exit_with_message(
					"In protocols::helical_bundle::BackboneGridSampler::apply():  Mainchain torsion angles to sample were specified that are not in the " + resname(ires) + " residue type.\n"
				);
			}
		}

		//Set fixed torsions:
		core::Size rescounter = 0;
		for ( core::Size it=1, itmax=torsions_to_fix_[ires].size(); it<=itmax; ++it ) {
			for ( core::Size ir=1, irmax=newpose->size(); ir<=irmax; ++ir ) {
				++rescounter;
				if ( rescounter > iresmax ) rescounter=1;
				if ( rescounter!=ires ) continue; //Do nothing if the current residue is not the residue in the repeating unit that we're setting.
				newpose->conformation().set_torsion( TorsionID(ir, BB, torsions_to_fix_[ires][it].first), torsions_to_fix_[ires][it].second );
			}
			if ( TR.visible() ) TR << "Set mainchain torsion " << torsions_to_fix_[ires][it].first << " for residue " << ires << " in the repeating unit to " << torsions_to_fix_[ires][it].second << " degrees." << std::endl;
		}

	} //Looping through all residues in the repeating unit.

	newpose->update_residue_neighbors();

	//Loop through all grid samples
	bool at_least_one_success(false);
	core::Size loopstart(1);
	core::Size loopend(total_samples);
	if ( nstruct_mode() ) { //Special case: if we're just doing one set of mainchain torsions per job, we need to figure out which set to do.
		if ( !protocols::jd2::jd2_used() ) {
			utility_exit_with_message(
				"In protocols::helical_bundle::BackboneGridSampler::apply() function: The nstruct_mode option was used, but the current application is not using JD2.");
		}
		core::Size curjob( protocols::jd2::current_nstruct_index() );
		core::Size totaljobs( protocols::jd2::max_nstruct_index() );
		if ( curjob==0 || totaljobs==0 ) {
			utility_exit_with_message(
				"In protocols::helical_bundle::BackboneGridSampler::apply() function: The nstruct_mode option was used, but invalid values were obtained for the current job index or the total number of jobs.");
		}
		if ( totaljobs < total_samples*nstruct_repeats() && TR.Warning.visible() ) {
			TR.Warning << "The BackboneGridSampler mover is in nstruct mode, meaning that one set of mainchain torsion values will be sampled per job.  However, the total number of jobs is less than the total number of samples!  Certain sets of mainchain torsion values will be missed!" << std::endl ;
		}
		//The current job might be greater than the total number of samples, in which case we should wrap around:
		loopstart = ( ( (curjob-1) % total_samples) + 1 ) / nstruct_repeats();
		loopend=loopstart;
	}
	for ( core::Size i=loopstart; i<=loopend; ++i ) {
		if ( nstruct_mode() ) { //If this is one-set-of-Crick-params-per-job mode
			if ( i > 1 ) {
				//If i is greater than 1, increment repeatedly until we reach the current permutation.
				//(Note that we start out on the first permutation, so there's no need to increment the first time).
				for ( core::Size j=2; j<=i; ++j ) {
					helper->increment_cur_indices(); //This is a recursive function that increments indices in a non-trivial way, and so must be called repeatedly to get the right indices.
				}
			}
		} else { //If this is NOT one-set-of-Crick-params-per-job mode
			//If i is greater than 1, increment the current permutation.
			//(Note that we start out on the first permutation, so there's no need to increment the first time).
			if ( i > 1 ) helper->increment_cur_indices();
		}

		//We need a temporary pose copying newpose:
		core::pose::PoseOP temppose( new core::pose::Pose( *newpose ) );

		//Set the variable mainchain torsions of this temporary pose:
		for ( core::Size it=1, itmax=helper->n_torsions_total(); it<=itmax; ++it ) {
			core::Size repeat_index(0);
			for ( core::Size ir=1, irmax=temppose->size(); ir<=irmax; ++ir ) {
				++repeat_index;
				if ( repeat_index > residues_per_repeat() ) repeat_index=1;
				if ( repeat_index != helper->residue_index(it) ) continue; //Only set the torsion if the current residue is the residue in the repeating unit whose torsion we want to set.
				temppose->conformation().set_torsion( TorsionID(ir, BB, helper->torsion_id(it)), helper->torsion_sample_value(it) );
			}
		}
		temppose->update_residue_neighbors();

		if ( TR.visible() ) {
			//Tracer output -- summarize current sample.
			core::Size const refrepeat( nres() > 2 ? 2 : nres() );
			core::Size const refres( (refrepeat-1)*residues_per_repeat() + 1 );
			TR << "Current sample:\t";
			for ( core::Size ir=refres, irmax=refres+residues_per_repeat()-1; ir<=irmax; ++ir ) {
				for ( core::Size it=1, itmax=temppose->residue(ir).mainchain_torsions().size(); it<=itmax; ++it ) {
					TR << "res" << ir-refres+1 << "_tors" << it << "=" << temppose->residue(ir).mainchain_torsions()[it] << "\t";
				}
			}
			TR << std::endl;
		}

		//Apply the preselection mover, if defined:
		if ( prescoring_mover_exists() && pre_scoring_mover_ ) {
			if ( TR.visible() ) {
				TR << "Applying pre-scoring mover." << std::endl;
				TR.flush();
			}
			pre_scoring_mover_->apply( *temppose );
			if ( pre_scoring_mover_->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
				if ( TR.visible() ) {
					TR << "Pre-scoring mover failed.  Discarding current grid sample and moving on to next." << std::endl;
					TR.flush();
				}
				continue;
			}
		}
		if ( prescoring_filter_exists() && pre_scoring_filter_ ) {
			if ( TR.visible() ) {
				TR << "Applying pre-scoring filter." << std::endl;
				TR.flush();
			}
			bool const filterpassed( pre_scoring_filter_->apply( *temppose ) ); //Apply the filter and store the result in filterpassed.
			if ( TR.visible() ) {
				if ( filterpassed ) TR << "Filter passed!" << std::endl;
				else TR << "Filter failed!  Discarding current grid sample and moving on to next." << std::endl;
			}
			if ( !filterpassed ) continue; //Go on to the next grid sample.
		}

		//At this point, at least one sample has had a successful application of the prescoring mover, and has passed the prescoring filter.
		at_least_one_success=true;

		//Score the pose:
		(*sfxn_)(*temppose);

		if ( TR_Results.visible() ) { //Write summary report to separate tracer to make it easy to mute everytihng except the output.
			TR_Results << "Sample " << i << " torsion values:\t";
			core::Size const refrepeat( nres() > 2 ? 2 : nres() );
			core::Size const refres( (refrepeat-1)*residues_per_repeat() + 1 );
			for ( core::Size ir=refres, irmax=refres+residues_per_repeat()-1; ir<=irmax; ++ir ) {
				if ( residues_per_repeat()>1 ) TR_Results << "Res" << ir-refres+1 << ":\t";
				for ( core::Size it=1, itmax=temppose->residue(ir).mainchain_torsions().size(); it<=itmax; ++it ) {
					TR_Results << temppose->residue(ir).mainchain_torsions()[it] << "\t";
				}
			}
			TR_Results << "SCORE:\t" << temppose->energies().total_energy();
			TR_Results << std::endl;
		}

		//Dump a PDB file, if the user has specified that this mover should do so:
		if ( pdb_output() ) {
			char outfile[1024];
			sprintf( outfile, "%s_%05lu.pdb", pdb_prefix().c_str(), static_cast<unsigned long>(i) );
			if ( TR.visible() ) TR << "Writing " << outfile << std::endl;
			temppose->dump_scored_pdb( outfile, *sfxn_ );
		}

		//Store this if it's the lowest-energy pose found so far:
		if (
				i==loopstart ||
				(select_low_ && temppose->energies().total_energy() < newpose->energies().total_energy()) ||
				(!select_low_ && temppose->energies().total_energy() > newpose->energies().total_energy())
				) {
			if ( TR.visible() ) {
				TR << "Current sample is " << (select_low_ ? "lowest" : "highest") << " energy discovered so far." << std::endl;
				TR.flush();
			}
			*newpose = *temppose;
		}
	} //Loop through all states

	//Determine success or failure of this mover:
	if ( at_least_one_success ) {
		if ( TR.visible() ) {
			TR << "Success.  Returning lowest-energy conformation sampled." << std::endl;
			TR.flush();
		}
		pose = *newpose; //DELETE ME
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} else {
		if ( TR.visible() ) {
			TR << "No sampled conformations passed pre-scoring mover and/or filter.  Returning with failed status." << std::endl;
			TR.flush();
		}
		set_last_move_status( protocols::moves::FAIL_DO_NOT_RETRY );
	}

	//Final message:
	if ( TR.visible() ) {
		TR << "Finished BackboneGridSampler apply()." << std::endl;
		TR.flush();
	}

	return;
} //apply()

////////////////////////////////////////////////////////////////////////////////

///
/// @brief Returns the name of this mover ("BackboneGridSampler").
// XRW TEMP std::string BackboneGridSampler::get_name() const{
// XRW TEMP  return "BackboneGridSampler";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////
//          PARSE MY TAG FUNCTION                                            ///
////////////////////////////////////////////////////////////////////////////////

/// @brief parse XML (specifically in the context of the parser/Rosetta_scripting scheme)
///
void
BackboneGridSampler::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &filters,
	protocols::moves::Movers_map const &movers,
	core::pose::Pose const & /*pose*/
) {
	try {

		if ( tag->getName() != "BackboneGridSampler" ) {
			throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible -- the tag name does not match the mover name.");
		}

		if ( TR.visible() ) TR << "Parsing options for BackboneGridSampler (\"" << tag->getOption<std::string>("name" ,"") << "\") mover." << std::endl;

		//Global options for this mover:
		runtime_assert_string_msg( tag->hasOption("scorefxn" ), "In BackboneGridSampler::parse_my_tag(): A \"scorefxn\" option must be specified!");
		set_sfxn(protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data_map )->clone()); // The scorefunction.
		if ( tag->hasOption("max_samples") ) {
			core::Size const val( tag->getOption<core::Size>("max_samples", 10000) );
			if ( TR.visible() ) TR << "Setting maximum number of samples to " << val << "." << std::endl;
			set_max_samples(val);
		}
		if ( tag->hasOption("selection_type") ) {
			std::string const val = tag->getOption<std::string>("selection_type", "");
			runtime_assert_string_msg( val=="high" || val=="low",
				"When parsing options for the BackboneGridSampler mover, could not interpret the selection_type.  This must be set to \"high\" or \"low\"." );
			if ( TR.visible() ) TR << "Setting selection type to " << val << "." << std::endl;
			if ( val=="high" ) set_selection_low(false);
			else set_selection_low(true);
		}
		if ( tag->hasOption("pre_scoring_mover") ) {
			protocols::moves::MoverOP curmover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "pre_scoring_mover" ), movers );
			set_prescoring_mover(curmover);
			if ( TR.visible() ) TR << "BackboneGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" has been assigned mover \"" << tag->getOption< std::string >("pre_scoring_mover") << "\" as a pre-scoring mover that will be applied to all sampled conformations prior to energy evaluation." << std::endl;
		}
		if ( tag->hasOption("pre_scoring_filter") ) {
			protocols::filters::FilterOP curfilter = protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "pre_scoring_filter" ), filters );
			set_prescoring_filter(curfilter);
			if ( TR.visible() ) TR << "BackboneGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" has been assigned filter \"" << tag->getOption< std::string >("pre_scoring_filter") << "\" as a pre-scoring filter that will be applied to all sampled conformations prior to energy evaluation." << std::endl;
		}
		if ( tag->hasOption("dump_pdbs") ) {
			bool const val = tag->getOption<bool>("dump_pdbs", false);
			set_pdb_output(val);
			if ( TR.visible() ) {
				if ( val ) TR << "Setting BackboneGridSampler mover \"" << tag->getOption< std::string >("name", "") << "\" to write out PDB files." << std::endl;
				else TR << "No PDB output will occur from this mover." << std::endl;
			}
		}
		if ( tag->hasOption("pdb_prefix") ) {
			std::string const val = tag->getOption<std::string>("pdb_prefix", "bgs_out");
			set_pdb_prefix(val);
			if ( TR.visible() ) TR << "Setting prefix for PDB output to " << val << "." << std::endl;
		}

		//Set nstruct mode and options:
		bool nstructmode = tag->getOption<bool>("nstruct_mode", false);
		if ( TR.visible() ) TR << "Setting nstruct mode to " << (nstructmode ? "true." : "false.") << "  This means that " << (nstructmode ? "each job will sample a different set of mainchain torsion values." : "every job will sample all sets of mainchain torsion values.") << std::endl;
		set_nstruct_mode(nstructmode);
		core::Size nstructrepeats( tag->getOption<core::Size>( "nstruct_repeats", 1 ) );
		if ( nstructrepeats<1 ) nstructrepeats=1;
		if ( TR.visible() ) TR << "Setting nstruct repeats to " << nstructrepeats << "." << std::endl;
		set_nstruct_repeats(nstructrepeats);

		//Set number of repeats, residues per repeat, and residue type:
		set_residues_per_repeat( tag->getOption<core::Size>( "residues_per_repeat", 1 ));
		if ( TR.visible() ) TR << "Set number of residues per repeating unit to " << residues_per_repeat() << "." << std::endl;

		if ( residues_per_repeat() == 1 && tag->hasOption("residue_count") ) { //For backwards compatibility.
			if ( tag->hasOption("repeat_count") ) {
				utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): Either \"residue_count\" or \"repeat_count\" may be specified if there is one residue per repeating unit, but both options may not be used.  Both were found.\n" );
			}
			set_nres( tag->getOption<core::Size>("residue_count", 12) );
		} else {
			if ( tag->hasOption("residue_count") ) {
				utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): A \"residue_count\" option was found, but more than one residue per repeat was specified.  Use \"repeat_count\" instead.\n" );
			}
			set_nres( tag->getOption<core::Size>("repeat_count", 12) );
		}
		if ( TR.visible() ) TR << "Set number of repeats to " << nres() << "." << std::endl;

		if ( residues_per_repeat() == 1 && tag->hasOption("residue_name") ) { //For backwards compatibility.
			set_resname(1, tag->getOption<std::string>( "residue_name", "ALA" ) );
			if ( TR.visible() ) TR << "Set residue type to " << resname(1) << "." << std::endl;
		} else {
			for ( core::Size i=1, imax=residues_per_repeat(); i<=imax; ++i ) {
				std::stringstream curtagstream;
				curtagstream << "residue_name_" << i;
				set_resname( i, tag->getOption<std::string>(curtagstream.str(), "ALA") );
				if ( TR.visible() ) TR << "Set residue type for residue " << i << " in the repeating unit to " << resname(i) << "." << std::endl;
			}
		}

		//Determine whether the ends should be capped
		set_cap_ends( tag->getOption<bool>( "cap_ends", false ) );
		if ( TR.visible() ) TR << "Set end-capping to " << (cap_ends()?"true":"false") << ".  This means that the termini will " << (cap_ends() ? "have acetylated and carboxyamidated N- and C-termini, respectively." : "have unmodified termini (NH3+/COO-)." ) << std::endl;

		/*** Parse Sub-Tags ***/
		utility::vector1< utility::tag::TagCOP > const branch_tags( tag->getTags() );
		for ( auto const & branch_tag : branch_tags ) { //Loop through all sub-tags
			if ( branch_tag->getName() == "MainchainTorsion" ) { //A mainchain torsion index has been specified.  Parse its options:
				if ( !branch_tag->hasOption("index") ) {
					utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): A <MainchainTorsion> tag must specify a mainchain torsion index with an \"index\" option.\n" );
				}
				core::Size index( branch_tag->getOption<core::Size>( "index", 0 ) );
				core::Size resindex( branch_tag->getOption<core::Size>( "res_index", 1 ) );

				if ( branch_tag->hasOption("value") ) {
					core::Real val( branch_tag->getOption<core::Real>( "value", 0 ) );
					if ( TR.visible() ) {
						TR << "Adding mainchain torsion index " << index << " for residue " << resindex << " in the repeating unit to be fixed to " << val << " degrees." << std::endl;
						TR.flush();
					}
					add_torsion_to_fix(resindex, index, val);
				} else {
					if ( !branch_tag->hasOption("start") ) {
						utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): A <MainchainTorsion> tag must specify a start of the torsion range to be sampled with a \"start\" option.  (Values are in DEGREES.)\n" );
					}
					if ( !branch_tag->hasOption("end") ) {
						utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): A <MainchainTorsion> tag must specify a start of the torsion range to be sampled with a \"end\" option.  (Values are in DEGREES.)\n" );
					}
					if ( !branch_tag->hasOption("samples") ) {
						utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::parse_my_tag(): A <MainchainTorsion> tag must specify the number of samples for a mainchain torsion with a \"samples\" option.\n" );
					}
					core::Real start( branch_tag->getOption<core::Real>( "start", -180.0 ) );
					core::Real end( branch_tag->getOption<core::Real>( "end", 180.0 ) );
					core::Size samples( branch_tag->getOption<core::Real>( "samples", 2 ) );
					if ( TR.visible() ) {
						TR << "Adding mainchain torsion index " << index << " to sample angle range " << start << " to end " << end << " with " << samples << " sample(s)." << std::endl;
						TR.flush();
					}
					add_torsion_to_sample( resindex, index, start, end, samples ); //Actually add the torsion (checking that it hasn't already been added, and that values are reasonable).
				}
			} //if tag name == MainchainTorsion
		} // looping through sub-tags
		/*** End Parse of Sub-Tags ***/

	} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
		TR << "kill me " << e.msg() << std::endl;
	}

	if ( TR.visible() ) TR.flush();
	return;
} //parse_my_tag

////////////////////////////////////////////////////////////////////////////////
//          PUBLIC FUNCTIONS                                                  //
////////////////////////////////////////////////////////////////////////////////

/// @brief Sets the mover that will be applied to all helical bundles generated prior to energy evaluation.
/// @details Note: if this is used, there is no guarantee that the resulting geometry will still lie within the
/// parameter space.  (That is, this mover could move the backbone.)
void BackboneGridSampler::set_prescoring_mover ( protocols::moves::MoverOP mover )
{
	pre_scoring_mover_ = mover;
	pre_scoring_mover_exists_ = true;
	return;
}

/// @brief Sets the filter that will be applied to all helical bundles generated prior to energy evaluation.
/// @details See the pre_scoring_filter_ private member variable for details.
void BackboneGridSampler::set_prescoring_filter ( protocols::filters::FilterOP filter )
{
	pre_scoring_filter_ = filter;
	pre_scoring_filter_exists_ = true;
	return;
}

/// @brief Add a mainchain torsion to sample, the range of values that will be sampled, and the number of samples.
/// @details The residue_index is the index in the repeating unit (1st residue, 2nd residue, etc.).  The torsion_index
/// is the mainchain torsion index in this residue.  Sampled values will go from start_of_range to end_of_range, with the
/// total number of samples given by the samples parameter.  If the range is -180 to 180, the samples are adjusted so that
/// 180 is not sampled twice.
void BackboneGridSampler::add_torsion_to_sample(
	core::Size const residue_index,
	core::Size const torsion_index,
	core::Real const &start_of_range,
	core::Real const &end_of_range,
	core::Size const samples
) {

	//First, let's check some things:
	if ( torsions_to_sample_.size() != residues_per_repeat() ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): Somehow, the number of residues per repeat and the size of the torsions_to_sample_ vector do not match.  This is a program error, so consult a developer or a mortician.\n" );
	}
	if ( residue_index < 1 || residue_index > torsions_to_sample_.size() ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The residue index is out of range.  The value must be greater than zero and less than or equal to the number of residues per repeat.\n" );
	}
	if ( torsion_index<=0 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The mainchain torsion index must be greater than 0.\n" );
	}
	if ( start_of_range < -180 || start_of_range > 180 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The start of the range to sample must lie within (-180,180).\n" );
	}
	if ( end_of_range < -180 || end_of_range > 180 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The end of the range to sample must lie within (-180,180).\n" );
	}
	if ( end_of_range <= start_of_range ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The end of the range to be sampled must be greater than the start of the range to be sampled.\n" );
	}
	if ( samples < 1 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample(): The number of samples must be at least 1.\n" );
	}

	for ( core::Size i=1, imax=torsions_to_sample_[residue_index].size(); i<=imax; ++i ) {
		if ( torsions_to_sample_[residue_index][i].first == torsion_index ) {
			utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample():  Could not add the torsion index to the list of torsions to sample.  It has already been added to this list.\n" );
		}
	}
	for ( core::Size i=1, imax=torsions_to_fix_[residue_index].size(); i<=imax; ++i ) {
		if ( torsions_to_fix_[residue_index][i].first == torsion_index ) {
			utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_sample():  Could not add the torsion index to the list of torsions to sample.  It has already been added to the list of torsions that are fixed.\n" );
		}
	}

	//OK, we're ready to add the torsion value and its range:
	std::pair <core::Real, core::Real> range( start_of_range, end_of_range );
	std::pair < std::pair<core::Real,core::Real>, core::Size > range_and_samples( range, samples );
	std::pair <core::Size, std::pair< std::pair<core::Real,core::Real>, core::Size >  > pair_to_add( torsion_index, range_and_samples );
	torsions_to_sample_[residue_index].push_back( pair_to_add );

	return;
}

/// @brief Add a mainchain torsion to fix, and that torsion's value.
/// @details The residue_index is the index of the residue in the repeating unit (1st, 2nd, 3rd, etc.).
/// The torsion_index is the mainchain torsion index in this residue, and the torsion_value is the
/// value at which to fix this residue.
void BackboneGridSampler::add_torsion_to_fix(
	core::Size const residue_index,
	core::Size const torsion_index,
	core::Real const &torsion_value
) {
	//Checks:
	if ( torsions_to_fix_.size()!=residues_per_repeat() ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix(): Program error!  Somehow, the number of residues per residue does not match the size of the torsions_to_fix_ vector.  Consult a developer or a mortician.\n" );
	}
	if ( residue_index < 1 || residue_index > torsions_to_fix_.size() ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix(): The index of the residue in the repeating unit must be greater than 0 and less than or equal to the number of residues per repeating unit.\n" );
	}
	if ( torsion_index<=0 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix(): The mainchain torsion index must be greater than 0.\n" );
	}
	if ( torsion_value < -180 || torsion_value > 180 ) {
		utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix(): The torsion value must lie within (-180,180).\n" );
	}

	for ( core::Size i=1, imax=torsions_to_fix_[residue_index].size(); i<=imax; ++i ) {
		if ( torsions_to_fix_[residue_index][i].first == torsion_index ) {
			utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix():  Could not add the torsion index to the list of torsions to fix.  It has already been added to this list.\n" );
		}
	}
	for ( core::Size i=1, imax=torsions_to_sample_[residue_index].size(); i<=imax; ++i ) {
		if ( torsions_to_sample_[residue_index][i].first == torsion_index ) {
			utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::add_torsion_to_fix():  Could not add the torsion index to the list of torsions to fix.  It has already been added to the list of torsions to sample.\n" );
		}
	}

	//OK, we're ready to add the torsion.
	std::pair <core::Size, core::Real> pair_to_add(torsion_index, torsion_value);
	torsions_to_fix_[residue_index].push_back(pair_to_add);

	return;
}

/// @brief Set up the PeptideStubMover object that will be used to build geometry.
///
void BackboneGridSampler::set_up_peptide_stub_mover()
{
	peptide_stub_mover_->set_reset_mode(true);
	peptide_stub_mover_->set_update_pdb_numbering_mode(true);

	peptide_stub_mover_->reset_mover_data();
	for ( core::Size i=1, imax=nres(); i<=imax; ++i ) {
		for ( core::Size j=1, jmax=residues_per_repeat(); j<=jmax; ++j ) {
			peptide_stub_mover_->add_residue( "Append", resname(j), 0, false, "", 1, 0, "" );
		}
	}

	peptide_stub_mover_initialized_=true;
	return;
}

/// @brief Reset the torsions_to_sample_, torsions_to_fix_, and resname_ vectors, and
/// initialize them based on the value of residues_per_repeat_.
void BackboneGridSampler::reset_and_initialize_data() {
	resname_.clear();
	torsions_to_fix_.clear();
	torsions_to_fix_.resize(residues_per_repeat());
	torsions_to_sample_.clear();
	torsions_to_sample_.resize(residues_per_repeat());

	for ( core::Size i=1; i<=residues_per_repeat(); ++i ) {
		resname_.push_back("ALA"); //Defaults to a list of alanines.
	}
	return;
}

/// @brief Set the number of residues per repeat.
/// @details This resets the torsions_to_sample_, torsions_to_fix_, and resname_ vectors, and
/// must therefore be called BEFORE setting up the torsions to sample, torsions to fix, or
/// residue types.
void BackboneGridSampler::set_residues_per_repeat( core::Size const val ) {
	if ( val < 1 ) utility_exit_with_message( "In protocols::helical_bundle::BackboneGridSampler::set_residues_per_repeat(): The residues per repeat must be greater than or equal to 1." );
	residues_per_repeat_ = val;
	reset_and_initialize_data();
	return;
}


////////////////////////////////////////////////////////////////////////////////
//          PRIVATE FUNCTIONS                                                 //
////////////////////////////////////////////////////////////////////////////////

/// @brief Is a value in a list?
///
bool BackboneGridSampler::is_in_list( core::Size const val, utility::vector1 < core::Size> const &list ) const {
	core::Size const listsize(list.size());
	if ( listsize==0 ) return false;
	for ( core::Size i=1; i<=listsize; ++i ) {
		if ( list[i]==val ) return true;
	}
	return false;
}

/// @brief Calculate the number of grid points that will be sampled, based on the options set by the user.
///
core::Size BackboneGridSampler::calculate_total_samples() const {

	core::Size total_samples(1);

	for ( core::Size ires=1, iresmax=torsions_to_sample_.size(); ires<=iresmax; ++ires ) {
		for ( core::Size i=1, imax=torsions_to_sample_[ires].size(); i<=imax; ++i ) {
			total_samples *= torsions_to_sample_[ires][i].second.second;
		}
	}

	return total_samples;
}

std::string BackboneGridSampler::get_name() const {
	return mover_name();
}

std::string BackboneGridSampler::mover_name() {
	return "BackboneGridSampler";
}

std::string subtag_for_bbgrid( std::string const & tag_name ) {
	return "stfbbg_" + tag_name;
}

void BackboneGridSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction lohigh;
	lohigh.name( "BackboneGridSampler_lohigh" );
	lohigh.base_type( xs_string );
	lohigh.add_restriction( xsr_enumeration, "low" );
	lohigh.add_restriction( xsr_enumeration, "high" );
	xsd.add_top_level_element( lohigh );

	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("scorefxn", xs_string, "Scorefunction to employ" )
		+ XMLSchemaAttribute::attribute_w_default("max_samples", xsct_non_negative_integer, "Maximum number of total backbone combinations to be sampled.", "10000" )
		+ XMLSchemaAttribute("selection_type", "BackboneGridSampler_lohigh", "Score criterion for selection: \"high\" or \"low\"." )
		+ XMLSchemaAttribute("pre_scoring_mover", xs_string, "A mover to apply after backbone torsions are set but before final scoring and evaluation (like a min mover or something similar)." )
		+ XMLSchemaAttribute("pre_scoring_filter", xs_string, "A filter to apply before scoring, which could help avoid wasteful scoring of bad conformations (like a bump check filter)." )
		+ XMLSchemaAttribute::attribute_w_default("dump_pdbs", xsct_rosetta_bool, "Dump all PDBs, if true; otherwise, there will be no PDB output at all.", "false" )
		+ XMLSchemaAttribute("pdb_prefix", xs_string, "A prefix to apply to all output PDBs." )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_mode", xsct_rosetta_bool, "If true, sample a different set of mainchain torsions for each job; if false, each job consists of the whole mainchain sampling effort.", "false" )
		+ XMLSchemaAttribute::attribute_w_default("nstruct_repeats", xsct_non_negative_integer, "Number of repeats to perform for each nstruct.", "1" )
		+ XMLSchemaAttribute::attribute_w_default("residues_per_repeat", xsct_non_negative_integer, "Number of residues in the minimal repeating unit of this secondary structure.", "1" )
		+ XMLSchemaAttribute::attribute_w_default("residue_count", xsct_non_negative_integer, "Number of residues in the secondary structure exemplar to be sampled.", "12" )
		+ XMLSchemaAttribute::attribute_w_default("repeat_count", xsct_non_negative_integer, "Number of residue-repeats in the secondary structure exemplar to be sampled.", "12" )
		+ XMLSchemaAttribute::attribute_w_default("residue_name", xs_string, "Residue type of which to create the secondary structure, indicated by three-letter code.", "ALA" )
		+ XMLSchemaAttribute::attribute_w_default("residue_name_1", xs_string, "Residue type of which to create the secondary structure, indicated by three-letter code.", "ALA" )
		+ XMLSchemaAttribute::attribute_w_default("residue_name_2", xs_string, "Residue type of which to create the secondary structure, indicated by three-letter code.", "ALA" )
		+ XMLSchemaAttribute::attribute_w_default("residue_name_3", xs_string, "Residue type of which to create the secondary structure, indicated by three-letter code.", "ALA" )
		+ XMLSchemaAttribute::attribute_w_default("residue_name_4", xs_string, "Residue type of which to create the secondary structure, indicated by three-letter code.", "ALA" )
		//More are possible?



		// NOTE residue_name_1 disaster. Have to repair.
		+ XMLSchemaAttribute::attribute_w_default("cap_ends", xsct_rosetta_bool, "If true, adds acetylated and amidated N- and C- termini.", "false" );

	AttributeList subtag_attributes;
	/*** Parse Sub-Tags ***/
	subtag_attributes + XMLSchemaAttribute::required_attribute( "index", xsct_non_negative_integer, "Mainchain torsion index indicated" )
		+ XMLSchemaAttribute::attribute_w_default( "res_index", xsct_non_negative_integer, "Residue whose mainchain torsion is being specified (if there is more than one residue per repeat)", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "value", xsct_real, "A single value in degrees, if this torsion ought to be fixed", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "start", xsct_real, "Starting value of a torsion range in degrees", "-180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "end", xsct_real, "Ending value of a torsion range in degrees", "-180.0" )
		+ XMLSchemaAttribute::attribute_w_default( "samples", xsct_non_negative_integer, "Number of samples to be taken of the dihedral range indicated", "2" );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "MainchainTorsion", subtag_attributes, "Tags describing individual torsions in the helix"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_bbgrid );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Sample mainchain torsions for peptides in a grid, saving good conformations", attlist, ssl );
}

std::string BackboneGridSamplerCreator::keyname() const {
	return BackboneGridSampler::mover_name();
}

protocols::moves::MoverOP
BackboneGridSamplerCreator::create_mover() const {
	return protocols::moves::MoverOP( new BackboneGridSampler );
}

void BackboneGridSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BackboneGridSampler::provide_xml_schema( xsd );
}


} //namespace helical_bundle
} //namespace protocols
