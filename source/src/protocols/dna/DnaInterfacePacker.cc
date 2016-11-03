// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Unit headers
#include <protocols/dna/DnaInterfacePacker.hh>
#include <protocols/dna/DnaInterfacePackerCreator.hh>

#include <protocols/dna/typedefs.hh>
#include <protocols/dna/DnaDesignDef.hh>
#include <protocols/dna/util.hh>
#include <protocols/dna/PDBOutput.hh>
#include <protocols/dna/DnaChains.hh>
#include <protocols/dna/RotamerDNAHBondFilter.hh>
#include <protocols/dna/RestrictDesignToProteinDNAInterface.hh>
#include <protocols/dna/SeparateDnaFromNonDna.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/types.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh> // lead_zero_string_of

// c++ headers
#include <vector> // for rot_to_pack
#include <iostream>
#include <sstream>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace protocols {
namespace dna {

bool ResTypeSequence_lt::operator() ( ResTypeSequence const & a, ResTypeSequence const & b ) const
{
	core::Size const asize( a.size() ), bsize( b.size() );
	if ( asize < bsize ) return true;
	if ( asize > bsize ) return false;
	for ( auto ait( a.begin() ), bit( b.begin() ), aend( a.end() ),
			bend( b.end() ); ait != aend && bit != bend; ++ait, ++bit ) {
		if ( ait->first < bit->first ) return true;
		if ( ait->first > bit->first ) return false;
		if ( ait->second->name1() < bit->second->name1() ) return true;
		if ( ait->second->name1() > bit->second->name1() ) return false;
	}
	return false;
}


int const PRECISION(3); // number of decimal places for floating point numbers

using utility::vector1;
using utility::string_split;
using namespace core;
using namespace chemical;
using namespace conformation;
using namespace optimization;
using namespace basic::options;
using namespace pack;
using namespace task;
using namespace operation;
using namespace pose;
using namespace scoring;

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

using basic::t_warning;
using basic::t_info;
using basic::t_debug;
static THREAD_LOCAL basic::Tracer TR( "protocols.dna.DnaInterfacePacker" );
static THREAD_LOCAL basic::Tracer TR_spec( "protocols.dna.Specificity" );

// for comparing ResTypeSequence, which contain ResidueTypeCOP pointers (must dereference for sorting purposes)

std::string
DnaInterfacePackerCreator::keyname() const
{
	return DnaInterfacePackerCreator::mover_name();
}

protocols::moves::MoverOP
DnaInterfacePackerCreator::create_mover() const {
	return protocols::moves::MoverOP( new DnaInterfacePacker );
}

std::string
DnaInterfacePackerCreator::mover_name()
{
	return "DnaInterfacePacker";
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief lightweight default constructor
DnaInterfacePacker::DnaInterfacePacker()
: protocols::simple_moves::PackRotamersMover( std::string("DnaInterfacePacker") ),
	reference_pose_(/* 0 */),
	dna_chains_(/* 0 */),
	minimize_(false),
	filename_root_( option[ OptionKeys::out::prefix ]() ),
	binding_E_(false),
	probe_specificity_(false),
	reversion_scan_(false),
	protein_scan_(false),
	base_only_(false),
	include_dna_potentials_in_specificity_calculations_(false),
	num_repacks_(1),
	specificity_repacks_(1),
	minimize_options_(/* 0 */),
	min_movemap_(/* 0 */),
	pdboutput_(/* 0 */),
	initialization_state_(false),
	pdbname_("")
{}

/// @brief functional constructor
DnaInterfacePacker::DnaInterfacePacker(
	ScoreFunctionOP scorefxn_in,
	bool minimize,
	std::string filename_root
) : protocols::simple_moves::PackRotamersMover( std::string("DnaInterfacePacker") ),
	reference_pose_(/* 0 */),
	dna_chains_(/* 0 */),
	minimize_( minimize ),
	filename_root_(std::move( filename_root )),
	binding_E_(false),
	probe_specificity_(false),
	reversion_scan_(false),
	protein_scan_(false),
	base_only_(false),
	include_dna_potentials_in_specificity_calculations_(false),
	num_repacks_(1),
	specificity_repacks_(1),
	minimize_options_(/* 0 */),
	min_movemap_(/* 0 */),
	pdboutput_(/* 0 */),
	initialization_state_(false),
	pdbname_("")
{
	score_function( scorefxn_in ); // calls PackRotamersMover setter
	binding_E_ = option[ OptionKeys::dna::design::binding ]();
	num_repacks_ = option[ OptionKeys::packing::ndruns ]();
	if ( option[ OptionKeys::dna::design::probe_specificity ].user() ) {
		probe_specificity_ = true;
		specificity_repacks_ = option[ OptionKeys::dna::design::probe_specificity ]();
	}
	reversion_scan_ = option[ OptionKeys::dna::design::reversion_scan ]();
	base_only_ = option[ OptionKeys::dna::design::base_contacts_only ]();
	include_dna_potentials_in_specificity_calculations_ =
		option[ OptionKeys::dna::design::specificity::include_dna_potentials ]();
	if ( option[ OptionKeys::dna::design::protein_scan ].user() ) {
		protein_scan_ = true;
		allowed_types_ = option[ OptionKeys::dna::design::protein_scan ]();
	}
}

/// @brief destructor
DnaInterfacePacker::~DnaInterfacePacker()= default;

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DnaInterfacePacker::fresh_instance() const
{
	return moves::MoverOP( new DnaInterfacePacker );
}

/// @brief required in the context of the parser/scripting scheme
moves::MoverOP
DnaInterfacePacker::clone() const
{
	return moves::MoverOP( new DnaInterfacePacker( *this ) );
}

std::string
DnaInterfacePacker::get_name() const {
	return DnaInterfacePackerCreator::mover_name();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief runs the packer, with support for efficiently looping over multiple explicit DNA sequences (provided that they can be represented by the RotamerSets/InteractionGraph)
/// @author ashworth
void
DnaInterfacePacker::apply( Pose & pose )
{
	if ( ! initialized() ) init_standard( pose );
	// discard info strings from previous calls to apply
	info().clear();

	if ( protein_scan_ ) protein_scan( pose );
	else standard_packing( pose );
}

void
DnaInterfacePacker::standard_packing( Pose & pose )
{
	// local copy of starting pose
	Pose starting_pose( pose );

	// loop over DNA sequences
	for ( ResTypeSequences::const_iterator dnaseq( dna_sequences_.begin() ),
			end( dna_sequences_.end() ); dnaseq != end; ++dnaseq ) {
		TR << "Packing ";
		print_sequence_pdb_nums( *dnaseq, pose, TR );

		// pack_rotamers trials
		for ( Size trial(0); trial < num_repacks_; ++trial ) {
			// restore starting structure
			if ( dnaseq != dna_sequences_.begin() ) pose = starting_pose;
			// cull over-representative rotamers/energies by masking out unspecified rotameric options
			// (including DNA)
			utility::vector0<int> rot_to_pack;
			if ( !dnaseq->empty() ) restrict_dna_rotamers( rotamer_sets(), *dnaseq, rot_to_pack );
			// run packer (calls PackRotamersMover virtual method)
			run( pose, rot_to_pack );
			// run post-packer routines, output, etc.
			post_packing( pose, *dnaseq, trial );
		}
		TR.flush();
	}
	TR.flush();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
DnaInterfacePacker::post_packing( Pose & pose, ResTypeSequence const & dnaseq, Size trial )
{
	// for collecting information from the various post-packing routines
	Strings info_lines;

	// add descriptors to name of (potential) output pdb file
	std::string seqtag;
	if ( !dnaseq.empty() ) seqtag = "_" + dna_seq_tag( pose, dnaseq );
	std::string pdbroot( filename_root_ + "_pack" + seqtag );
	pdbname_ = pdbroot + "_" + lead_zero_string_of(trial,4);

	// calculate DNA sequence specificity (if enabled)
	// first element: bound specificities, second: binding specificities (if enabled)
	std::pair< SequenceScores, SequenceScores > specificities;

	if ( probe_specificity_ ) {
		specificities = measure_bp_specificities( pose );
		// format specificity results for output
		for ( SequenceScores::const_iterator bound( specificities.first.begin() ),
				bound_end( specificities.first.end() ); bound != bound_end; ++bound ) {
			std::ostringstream os;
			os << std::showpoint << std::fixed << std::setprecision(PRECISION);
			os << "Specificities(bound): ";
			os << seq_pdb_str( bound->first, pose ) << " = " << bound->second;
			TR << os.str() << '\n';
			info_lines.push_back( "REMARK " + os.str() );
		}
		if ( binding_E_ ) {
			for ( SequenceScores::const_iterator binding( specificities.second.begin() ),
					binding_end( specificities.second.end() ); binding != binding_end; ++binding ) {
				std::ostringstream os;
				os << std::showpoint << std::fixed << std::setprecision(PRECISION);
				os << "Specificities(binding): ";
				os << seq_pdb_str( binding->first, pose ) << " = " << binding->second;
				TR << os.str() << '\n';
				info_lines.push_back( "REMARK " + os.str() );
			}
		}
	}

	// minimize bound structure (if minimization enabled)
	if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
		AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
	}

	// score bound structure
	Real const bound_score( (*score_function())(pose) );

	// optional binding energy calculation
	Real binding_score(0.);
	if ( binding_E_ ) {
		// if output_unbound_pdb option is true, output first unbound structure
		bool const output_pdb(
			trial == 0 && option[ OptionKeys::dna::design::output_unbound_pdb ]() );
		binding_score = bound_score - unbound_score( pose, output_pdb, pdbname_ );
		// add binding information to extra pdb info
		std::ostringstream bindingstream;
		bindingstream << std::showpoint << std::fixed << std::setprecision(PRECISION)
			<< "Binding energy: " << binding_score;
		TR << bindingstream.str() << std::endl;
		info_lines.push_back( "REMARK " + bindingstream.str() );
	}

	if ( reversion_scan_ ) {
		// try to revert 'insignificant' amino acid mutations made during design
		std::pair< Real, Real > overall_specificities( 0., 0. );
		ResTypeSequence currseq( current_working_sequence( pose ) );
		// look up overall specificity scores, pass to reversion routine as baseline
		if ( specificities.first.find( currseq ) != specificities.first.end() ) {
			overall_specificities.first = specificities.first[ currseq ];
		}
		if ( specificities.second.find( currseq ) != specificities.second.end() ) {
			overall_specificities.second = specificities.second[ currseq ];
		}
		reversion_scan( pose, bound_score, binding_score, overall_specificities );
	}

	if ( pdboutput_ ) {
		// write out pdb if local pdb outputting is enabled
		// (if using the job distributor, the PDBOutput class should interface with it instead)
		bool const overwrite_old_info(true);
		pdboutput_->add_info( "REMARK DnaInterfacePacker " + pdbroot + ":", info_lines, !overwrite_old_info );
		pdboutput_->note_designed_residues( task() );
		pdboutput_->score_function( *score_function() );
		(*pdboutput_)( pose, pdbname_ + ".pdb" );
	}

	// store collected information for future reference
	info().insert( info().end(), info_lines.begin(), info_lines.end() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief standard initialization of the necessary member data
/// @author ashworth
/// @details pose is nonconst because it is so in pack_rotamers::setup()
void
DnaInterfacePacker::init_standard( Pose & pose )
{

	TR << std::showpoint << std::fixed << std::setprecision(PRECISION);

	reference_residue_types_.clear();
	if ( reference_pose_ ) {
		for ( Size index(1), end( reference_pose_->size() ); index <= end; ++index ) {
			reference_residue_types_.push_back( reference_pose_->residue_type(index).get_self_ptr() );
		}
	} else {
		for ( Size index(1), end( pose.size() ); index <= end; ++index ) {
			reference_residue_types_.push_back( pose.residue_type(index).get_self_ptr() );
		}
	}

	dna_chains_ = DnaChainsOP( new DnaChains );
	find_basepairs( pose, *dna_chains_ );

	TaskFactoryOP my_tf;
	// if there is no initialized TaskFactory, create default one here
	if ( ! task_factory() ) { // PackRotamersMover accessor
		my_tf = TaskFactoryOP( new TaskFactory );
		my_tf->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
		if ( option[ OptionKeys::packing::resfile ].user() ) my_tf->push_back( TaskOperationCOP( new ReadResfile ) );
		RestrictDesignToProteinDNAInterfaceOP rdtpdi( new RestrictDesignToProteinDNAInterface );
		rdtpdi->set_reference_pose( reference_pose_ );
		if ( ! targeted_dna_.empty() ) rdtpdi->copy_targeted_dna( targeted_dna_ );
		rdtpdi->copy_dna_chains( dna_chains_ );
		rdtpdi->set_base_only( base_only_ );
		my_tf->push_back( rdtpdi );
	} else { // TaskFactory already exists, make copy to tamper with
		my_tf = TaskFactoryOP( new TaskFactory( *task_factory() ) );
	}
	// a protein-DNA hbonding filter for ex rotamers that the PackerTask makes available to the rotamer set during rotamer building (formerly known as 'rotamer explosion')
	RotamerDNAHBondFilterOP rot_dna_hb_filter( new RotamerDNAHBondFilter( -0.5, base_only_ ) );
	my_tf->push_back( TaskOperationCOP( new AppendRotamer( rot_dna_hb_filter ) ) );

	task_factory( my_tf ); // PackRotamersMover setter

	setup( pose ); // PackRotamersMover method--fills RotamerSets and IG

	rot_dna_hb_filter->report(); // filtered/accepted statistic for filtered rotamers

	// set up for looping over multiple DNA sequences if task() has DNA positions set to 'scan'
	// makes all canonical combinations that can be described by the rotamer set
	make_dna_sequence_combinations();

	if ( dna_sequences_.empty() ) {
		// no list of DNA sequences to loop over:
		// see if the PackerTask is configured to target any particular DNA sequence, and if so store it
		ResTypeSequence target_seq( get_targeted_sequence( pose ) );
		if ( !target_seq.empty() ) dna_sequences_.push_back( target_seq );
	}
	if ( !dna_sequences_.empty() ) {
		TR << "DNA sequences to be considered:" << '\n';
		print_sequences_pdb_nums( dna_sequences_, pose, TR );
		TR << std::endl;
	} else dna_sequences_.push_back( ResTypeSequence() );
	// no DNA sequences specified by user. add a dummy sequence (placeholder allowing loop)
	// if 'dna_sequences_' is still empty at this point, the SimAnnealer will perform DNA rotamer substitutions using any available DNA rotamers (including mutation). To constrain designed DNA sequences to Watson-Crick pairing, see WatsonCrickRotamerCouplings

	if ( minimize_ ) {
		// set up minimizer
		std::string const min_type( option[ OptionKeys::run::min_type ]() );
		Real const tolerance( option[ OptionKeys::run::min_tolerance ]() );
		minimize_options_ = core::optimization::MinimizerOptionsOP( new MinimizerOptions( min_type, tolerance, true ) );
		min_movemap_ = core::kinematics::MoveMapOP( new kinematics::MoveMap );
		TR << "Chi minimization will be allowed for the following residues:";
		for ( Size index(1); index <= task()->total_residue(); ++index ) {
			if ( !pose.residue_type(index).is_protein() ) continue;
			if ( task()->pack_residue( index ) ) {
				min_movemap_->set_chi( index, true );
				if ( pose.pdb_info() ) {
					TR << " " << pose.pdb_info()->chain(index) << "." << pose.pdb_info()->number(index);
				} else {
					TR << " " << pose.chain(index) << "." << index;
				}
			}
		}
		TR << std::endl;
	}

	if ( pdboutput_ ) pdboutput_->reference_pose( pose );

	initialization_state_ = false;
	runtime_assert( initialized() );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
bool
DnaInterfacePacker::initialized() const
{
	bool is_initialized(
		score_function()
		&& ( task() || task_factory() )
		&& ig()
		&& rotamer_sets()
		&& dna_chains_
		&& !reference_residue_types_.empty()
		&& !initialization_state_
	);
	if ( minimize_ ) is_initialized &= ( minimize_options_ && min_movemap_ );
	return is_initialized;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void
DnaInterfacePacker::clear_initialization()
{
	initialization_state_ = true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void DnaInterfacePacker::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & datamap,
	protocols::filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose
)
{
	if ( tag->hasOption("binding") ) binding_E_ = tag->getOption<bool>("binding");
	if ( tag->hasOption("base_only") ) base_only_ = tag->getOption<bool>("base_only");
	if ( tag->hasOption("minimize") ) minimize_ = tag->getOption<bool>("minimize");
	if ( tag->hasOption("reversion_scan") ) reversion_scan_ = tag->getOption<bool>("reversion_scan");
	if ( tag->hasOption("probe_specificity") ) {
		probe_specificity_ = true;
		specificity_repacks_ = tag->getOption< Size >("probe_specificity");
	}
	if ( tag->hasOption("pdb_output") ) {
		if ( tag->getOption<bool>("pdb_output") ) {
			if ( !pdboutput_ ) pdboutput_ = PDBOutputOP( new PDBOutput );
		}
	}
	if ( tag->hasOption("protein_scan") ) protein_scan_ = tag->getOption<bool>("protein_scan");
	if ( tag->hasOption("allowed_types") ) allowed_types_ = tag->getOption<std::string>("allowed_types");
	if ( protein_scan_ ) runtime_assert( !allowed_types_.empty() );

	// the following are calls to PackRotamersMover methods
	parse_score_function( tag, datamap, filters, movers, pose );
	parse_task_operations( tag, datamap, filters, movers, pose );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details looks for rotable DNA positions in the RotamerSets, generates a list of all canonical sequence combinations for them
/// @author ashworth
void
DnaInterfacePacker::make_dna_sequence_combinations()
{
	runtime_assert( task() != nullptr );
	runtime_assert( dna_chains_ != nullptr );
	dna_sequences_.clear();
	vector1< Size > seq_indices;
	// search the rotamer set for rotable DNA residues
	for ( Size moltenres(1); moltenres <= rotamer_sets()->nmoltenres(); ++moltenres ) {
		Size resid( rotamer_sets()->moltenres_2_resid( moltenres ) );
		if ( !task()->has_behavior( "SCAN", resid ) ) continue;
		if ( !dna_chains_->is_top( resid ) ) continue;
		// just make the top-strand combinations, handle complements later
		seq_indices.push_back( resid );
	}

	if ( seq_indices.empty() ) return; // no positions found at which to make combinations

	ResTypeSequence sequence; // empty starter sequence for the following recursive function
	make_sequence_combinations( seq_indices.begin(), seq_indices, task(), sequence, dna_sequences_ );

	// TR << std::flush << "Single-stranded sequences:" << '\n';
	// print_sequences_pdb_nums( dna_sequences_, pose, TR );
	// TR << std::endl;

	for ( auto & dna_sequence : dna_sequences_ ) {
		add_complementary_sequence( dna_sequence );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details turns single-stranded DNA sequences into double-stranded ones
/// @author ashworth
void
DnaInterfacePacker::add_complementary_sequence( ResTypeSequence & sequence )
{
	runtime_assert( dna_chains_ != nullptr );
	// fill a temporary ResTypeSequence, in order to avoid iterator-related issues with std::map
	ResTypeSequence complement;
	for ( auto & postype : sequence ) {
		Size const index( postype.first );
		DnaPosition const & dnatop( (*dna_chains_)[ index ] );
		runtime_assert( dnatop.top() == index );
		if ( !dnatop.paired() ) continue;
		Size const comppos( dnatop.bottom() );
		// find the complement type in the current typeset
		ResidueTypeCOP type( postype.second );
		ResidueTypeSetCOP typeset( type->residue_type_set() );
		ResidueTypeCOP comptype( typeset->get_representative_type_aa( dna_base_partner( type->aa() ) ) );
		complement[ comppos ] = comptype;
	}
	// append this temporary bottom-stranded sequence to the original top-stranded sequence
	for ( auto & postype : complement ) {
		sequence.insert( postype );
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief
/// @author ashworth
Real
DnaInterfacePacker::unbound_score(
	Pose const & pose,
	bool output_pdb, // = false
	std::string pdbname // = ""
)
{
	// TR << "Calculating unbound score..." << std::endl;
	Pose unbound_pose( pose );
	SeparateDnaFromNonDna unbind;
	unbind.apply( unbound_pose );
	// clone ScoreFunctionCOP so that we can score without distance constraints
	runtime_assert( score_function() != nullptr );
	ScoreFunctionOP nonconst_scorefxn( score_function()->clone() );
	// set distance constraints weights to zero before scoring unbound structure
	nonconst_scorefxn->set_weight( atom_pair_constraint, 0.0 );
	Real const unbound_score( (*nonconst_scorefxn)(unbound_pose) );
	// optional: write out unbound pose to verify proper interface separation
	if ( output_pdb ) {
		PDBOutputOP unbound_outputter;
		if ( pdboutput_ ) unbound_outputter = pdboutput_;
		else unbound_outputter = PDBOutputOP( new PDBOutput );
		unbound_outputter->score_function( *nonconst_scorefxn );
		(*unbound_outputter)( unbound_pose, pdbname + "_unbound.pdb" );
	}
	return unbound_score;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief This method returns the overall bolztmann probabilities for target bound and binding(if enabled) states vs. all single-basepair variant competitors. (For this + specifities for individual basepair positions, use the measure_bp_specificities method instead.)
/// @author ashworth
std::pair< Real, Real > DnaInterfacePacker::measure_specificity( Pose & pose )
{
	if ( ! probe_specificity_ ) return std::make_pair( 0., 0. );
	TR_spec << "\nMeasuring specificity by repacking against other possible DNA states:" << std::endl;

	// figure out current TOP-STRANDED DNA sequence
	ResTypeSequence current_sequence( current_working_sequence( pose ) );

	if ( current_sequence.empty() ) {
		TR << "No double-stranded DNA positions found!" << std::endl;
		return std::make_pair( 0., 0. );
	}
	ResTypeSequences specificity_sequences;
	// add current_sequence
	specificity_sequences.push_back( current_sequence );
	// add competitor DNA sequences
	make_single_mutants( current_sequence, task(), specificity_sequences );
	// add complementary DNA positions to the these top-stranded sequences
	for ( auto & specificity_sequence : specificity_sequences ) {
		add_complementary_sequence( specificity_sequence );
	}

	// for ( ResTypeSequences::iterator sequence( specificity_sequences.begin() ),
	//   end( specificity_sequences.end() ); sequence != end; ++sequence ) {
	//  print_sequence_pdb_nums( *sequence, pose, TR_spec );
	// }
	// TR_spec << std::endl;

	// first element is bound, second binding if binding_ flag is true
	std::pair< SequenceScores, SequenceScores > sequence_scores(
		measure_specificities( pose, specificity_sequences ) );

	Real bound_specificity(0.), binding_specificity(0.);

	// calculate_specificity expects complemented target sequence
	add_complementary_sequence( current_sequence );
	TR_spec << "\nCalculating bound specificity:";
	bound_specificity = calculate_specificity( pose, current_sequence, sequence_scores.first );
	if ( binding_E_ ) {
		TR_spec << "\nCalculating binding specificity:";
		binding_specificity = calculate_specificity( pose, current_sequence, sequence_scores.second );
	}
	TR_spec << std::endl;
	return std::make_pair( bound_specificity, binding_specificity );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Measures and calculates bound and binding specificities of the current protein sequence for its target DNA sequence, vs. single-basepair variant competitors. Also returns specificity scores for each individual basepair in a multiple-basepair target region.
/// @author ashworth
// return values: first maps bound specificities, second maps binding specificities
std::pair< SequenceScores, SequenceScores >
DnaInterfacePacker::measure_bp_specificities( Pose & pose )
{
	if ( ! probe_specificity_ ) return std::pair< SequenceScores, SequenceScores >();
	TR_spec << "\nMeasuring individual basepair specificity by explicitly modeling alternative "
		<< "DNA states:" << std::endl;

	// figure out current TOP-STRANDED DNA sequence
	ResTypeSequence const current_sequence( current_working_sequence( pose ) );

	if ( current_sequence.empty() ) {
		TR << "No targeted double-stranded DNA positions found!" << std::endl;
		return std::pair< SequenceScores, SequenceScores >();
	}

	SequenceScores bound_scores, binding_scores;
	SequenceScores binding_specificities, bound_specificities;
	// separate specificity measurements into individual basepair positions
	for ( auto bppos( current_sequence.begin() ),
			end( current_sequence.end() ); bppos != end; ++bppos ) {
		ResTypeSequences single_bp_variants;
		single_bp_variants.push_back( current_sequence );
		Size index( bppos->first );
		ResidueLevelTask const & rtask( task()->residue_task(index) );
		// add competitors for this position
		for ( auto type( rtask.allowed_residue_types_begin() ),
				end_type( rtask.allowed_residue_types_end() ); type != end_type; ++type ) {
			if ( (*type)->aa() == bppos->second->aa() ) continue; // avoid duplicating input sequence
			// ignore adduct variant types
			if ( (*type)->has_variant_type( chemical::ADDUCT_VARIANT ) ) continue;
			// new copy of current sequence
			ResTypeSequence single_bp_mutant( current_sequence );
			// make change
			single_bp_mutant[ index ] = *type;
			single_bp_variants.push_back( single_bp_mutant );
		}
		// add complements (top-stranded sequences -> double-stranded)
		for ( auto & single_bp_variant : single_bp_variants ) {
			add_complementary_sequence( single_bp_variant ); // uses dna_chains_ information
		}
		// model/score states, calculate specificities
		// first element is bound, second binding if binding_ flag is true
		std::pair< SequenceScores, SequenceScores > sequence_scores(
			measure_specificities( pose, single_bp_variants ) );

		if ( current_sequence.size() > 1 ) {
			// calculate and store specificities for individual basepairs if multiple are involved
			// use single-bp sequence for annotation,
			ResTypeSequence current_single_bp;
			current_single_bp[ bppos->first ] = bppos->second;
			// but need full complemented target sequence for actual calculation
			ResTypeSequence complemented_current_sequence( current_sequence );
			add_complementary_sequence( complemented_current_sequence );
			TR_spec << "\nCalculating bound specificity for " << seq_pdb_str( current_single_bp, pose );
			bound_specificities[ current_single_bp ] =
				calculate_specificity( pose, complemented_current_sequence, sequence_scores.first );
			if ( binding_E_ ) {
				TR_spec << "\nCalculating binding specificity for "
					<< seq_pdb_str( current_single_bp, pose );
				binding_specificities[ current_single_bp ] =
					calculate_specificity( pose, complemented_current_sequence, sequence_scores.second );
			}
		}
		// accumulate sequence scores for overall calculation
		for ( SequenceScores::const_iterator ss_it( sequence_scores.first.begin() ),
				ss_end( sequence_scores.first.end() ); ss_it != ss_end; ++ss_it ) {
			bound_scores[ ss_it->first ] = ss_it->second;
		}
		for ( SequenceScores::const_iterator ss_it( sequence_scores.second.begin() ),
				ss_end( sequence_scores.second.end() ); ss_it != ss_end; ++ss_it ) {
			binding_scores[ ss_it->first ] = ss_it->second;
		}
	}

	// compute overall specificity for target sequence vs. all single-bp variants in design region
	TR_spec << "\nCalculating bound specificity for " << seq_pdb_str( current_sequence, pose );
	// need full complemented target sequence for calculation
	ResTypeSequence complemented_current_sequence( current_sequence );
	add_complementary_sequence( complemented_current_sequence );
	bound_specificities[ current_sequence ] =
		calculate_specificity( pose, complemented_current_sequence, bound_scores );
	if ( binding_E_ ) {
		TR_spec << "\nCalculating binding specificity for " << seq_pdb_str( current_sequence, pose );
		binding_specificities[ current_sequence ] =
			calculate_specificity( pose, complemented_current_sequence, binding_scores );
	}

	// store all of the individual scores as info
	for ( SequenceScores::const_iterator it( bound_scores.begin() ), end( bound_scores.end() );
			it != end; ++it ) {
		std::ostringstream os;
		os << std::showpoint << std::fixed << std::setprecision(PRECISION);
		os << "REMARK SeqScore(bound): " << seq_pdb_str( it->first, pose ) << " = " << it->second;
		info().push_back( os.str() );
	}
	for ( SequenceScores::const_iterator it( binding_scores.begin() ), end( binding_scores.end() );
			it != end; ++it ) {
		std::ostringstream os;
		os << std::showpoint << std::fixed << std::setprecision(PRECISION);
		os << "REMARK SeqScore(binding): " << seq_pdb_str( it->first, pose ) << " = " << it->second;
		info().push_back( os.str() );
	}

	return std::make_pair( bound_specificities, binding_specificities );
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief This requires that all DNA states to be searched are already represented in the rotamer set and interaction graph
/// @author ashworth
std::pair< SequenceScores, SequenceScores >
DnaInterfacePacker::measure_specificities( Pose & pose, ResTypeSequences const & dna_sequences )
{
	// save starting pose to avoid permanently altering it
	Pose starting_pose( pose );

	SequenceScores sequence_scores, sequence_binding_scores;

	for ( auto const & dna_sequence : dna_sequences ) {
		Real best_trial_E(0), best_trial_binding_E(0);
		// restrict packer to current protein sequence and this DNA sequence
		utility::vector0< int > rot_to_pack;
		vector1< ResidueTypeCOP > single_sequence;
		for ( Size index(1), end( pose.size() ); index <= end; ++index ) {
			single_sequence.push_back( pose.residue_type(index).get_self_ptr() );
		}
		// alter dna types according to this dna sequence
		for (const auto & it : dna_sequence) {
			single_sequence[ it.first ] = it.second;
		}
		// populate rot_to_pack with only the rotamers that reflect single_sequence
		restrict_to_single_sequence( rotamer_sets(), single_sequence, rot_to_pack );

		for ( Size trial(0); trial < specificity_repacks_; ++trial ) {
			run( pose, rot_to_pack ); // calls PackRotamersMover method
			if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
				AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
			}
			ScoreFunctionOP nonconst_scorefxn( score_function()->clone() );
			if ( ! include_dna_potentials_in_specificity_calculations_ ) {
				// temporarily disable dna conformational potentials for assessing specificity
				nonconst_scorefxn->set_weight( dna_bp, 0. );
				nonconst_scorefxn->set_weight( dna_bs, 0. );
			}
			Real trial_E( (*nonconst_scorefxn)( pose ) );
			if ( trial == 0 || ( trial_E < best_trial_E ) ) {
				best_trial_E = trial_E;
				if ( binding_E_ ) best_trial_binding_E = trial_E - unbound_score( pose );
			}
			if ( pdboutput_ && option[ OptionKeys::dna::design::specificity::output_structures ]() ) {
				std::string pdbname(
					filename_root_ + "_" + dna_seq_tag( pose, current_working_sequence(pose) ) +
					"_spec_" + dna_seq_tag( pose, dna_sequence ) + ".pdb"
				);
				pdboutput_->score_function( *nonconst_scorefxn );
				(*pdboutput_)( pose, pdbname );
			}
		}
		sequence_scores[ dna_sequence ] = best_trial_E;
		sequence_binding_scores[ dna_sequence ] = best_trial_binding_E;
	}
	pose = starting_pose; // restore starting, unaltered pose
	return std::make_pair( sequence_scores, sequence_binding_scores );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief calculates specificity as a Boltzmann probability of the target sequence in the presence of competitors
/// @author ashworth
Real
DnaInterfacePacker::calculate_specificity(
	Pose const & pose,
	ResTypeSequence const & target_sequence,
	SequenceScores const & sequence_scores
)
{
	TR_spec << std::showpoint << std::fixed << std::setprecision(PRECISION) << '\n';
	// Boltzmann temperature
	Real const temp( option[ OptionKeys::dna::design::Boltz_temp ]() );

	// find low energy
	Real low(0);
	for ( auto iter( sequence_scores.begin() );
			iter != sequence_scores.end(); ++iter ) {
		Real score( iter->second );
		if ( iter == sequence_scores.begin() || ( score < low ) ) low = score;
	}

	Real const inv_temp( 1.0 / temp );
	Real num(0), denom(0);
	for ( auto const & sequence_score : sequence_scores ) {
		ResTypeSequence const & sequence( sequence_score.first );
		Real score( sequence_score.second );
		TR_spec << "\t";
		for ( auto pos( sequence.begin() ); pos != sequence.end(); ++pos ) {
			if ( pos != sequence.begin() ) TR_spec << ", ";
			if ( pose.pdb_info() ) {
				TR_spec << pose.pdb_info()->chain( pos->first ) << "."
					<< pose.pdb_info()->number( pos->first ) << "." << dna_full_name3( pos->second->name3() );
			} else {
				TR_spec << pose.chain( pos->first ) << "."
					<< pos->first << "." << dna_full_name3( pos->second->name3() );
			}
		}
		TR_spec << ": " << score << '\n';
		Real term( std::exp( ( low - score ) * inv_temp ) );
		if ( sequence == target_sequence ) num += term;
		denom += term;
	}
	if ( denom == 0. ) return 0.;
	Real const spec( num / denom );
	TR_spec << "\tspecificity: " << spec << std::endl;
	return spec;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

class Reversion {
public:
	Reversion( Size i = 0, ResidueTypeCOP t = nullptr )
	: index(i), type(std::move(t)), dscore_bound(0.), dspec_bound(0.), dscore_binding(0.), dspec_binding(0.) {}
	~Reversion()= default;
	// assign a number to the effect of this reversion
	Real reversion_score() const { return -1 * dspec_binding; }
	// for sorting when looking for the 'most acceptable' reversion in a set of reversions
	bool operator < ( Reversion const & other ) const
	{ return reversion_score() < other.reversion_score(); }
	Size index;
	ResidueTypeCOP type;
	Real dscore_bound, dspec_bound, dscore_binding, dspec_binding;
};
typedef vector1< Reversion > Reversions;

void
DnaInterfacePacker::reversion_scan(
	Pose & pose,
	Real starting_bound_score, // = 0.
	Real starting_binding_score, // = 0.
	std::pair< Real, Real > starting_specificities // = std::make_pair(0.,0.)
)
{
	if ( ! initialized() ) init_standard( pose );

	TR << std::flush << "Starting reversion scan: using starting scores: " << "bound = "
		<< starting_bound_score << ", binding = " << starting_binding_score
		<< ", specificity.bound = " << starting_specificities.first << ", specificity.binding = "
		<< starting_specificities.second << std::endl << '\n';

	Real current_bound_score( starting_bound_score ), current_binding_score( starting_binding_score );
	std::pair< Real, Real > current_specificities( starting_specificities );

	vector1< ResidueTypeCOP > fixed_residue_types;

	for ( Size index(1), end( pose.size() ); index <= end; ++index ) {
		fixed_residue_types.push_back( pose.residue_type(index).get_self_ptr() );
	}

	// find mutations (positions to revert) based upon comparison to a starting sequence
	Reversions reversions;
	runtime_assert( fixed_residue_types.size() == reference_residue_types_.size() );
	for ( Size index(1), end( fixed_residue_types.size() ); index != end; ++index ) {
		ResidueTypeCOP reference_type( reference_residue_types_[index] );
		if ( reference_type->is_protein() &&
				fixed_residue_types[index]->name3() != reference_type->name3() ) {
			reversions.push_back( Reversion( index, reference_type ) );
		}
	}

	Real const dscore_cutoff( option[ OptionKeys::dna::design::reversion::dscore_cutoff ]() ),
		dspec_cutoff( option[ OptionKeys::dna::design::reversion::dspec_cutoff ]() );

	// do single reversions to wildtype while they do not harm energy/specificity
	Size round(0);
	while ( true ) {
		// assess changes in energy and specificity for each single revertant in parallel
		for ( auto & reversion : reversions ) {
			Size const index( reversion.index );
			ResidueTypeCOP starting_type( fixed_residue_types[ index ] ),
				reference_type( reference_residue_types_[ index ] ); // 'reference' == 'native'
			fixed_residue_types[ index ] = reference_type;

			Real best_score(0.), best_binding_score(0.);
			std::pair< Real, Real > best_specificities;

			for ( Size trial(1); trial <= num_repacks_; ++trial ) {
				// repack bound pose with this reversion
				utility::vector0<int> rot_to_pack;
				restrict_to_single_sequence( rotamer_sets(), fixed_residue_types, rot_to_pack );
				// calls PackRotamersMover method
				run( pose, rot_to_pack );
				if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
					AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
				}
				Real const score( ( *score_function() )( pose ) );
				if ( trial == 1 || score < best_score ) {
					best_score = score;
					if ( binding_E_ ) best_binding_score = score - unbound_score( pose );
					best_specificities = measure_specificity( pose );
				}
			}
			// undo the reversion so that it does not affect the others in the same round
			// (fixed_residue_types will restore the starting type in the next call to the packer)
			fixed_residue_types[ reversion.index ] = starting_type;

			Real const dscore_bound( best_score - current_bound_score ),
				dscore_binding( best_binding_score - current_binding_score ),
				dspec_bound( best_specificities.first - current_specificities.first ),
				dspec_binding( best_specificities.second - current_specificities.second );

			TR << "Scores for reversion from " << starting_type->name3() << " to "
				<< reference_type->name3() << " at ";
			if ( pose.pdb_info() ) {
				TR << pose.pdb_info()->chain( index ) << "." << pose.pdb_info()->number( index ) << ":";
			} else {
				TR << pose.chain( index ) << "." << index << ":";
			}
			TR << " bound = " << best_score << " (" << dscore_bound << ")"
				<< ", binding = " << best_binding_score << " (" << dscore_binding << ")"
				<< ", specificity.bound = " << best_specificities.first << " (" << dspec_bound << ")"
				<< ", specificity.binding = " << best_specificities.second << " (" << dspec_binding
				<< ")\n";

			reversion.dscore_bound = dscore_bound;
			reversion.dscore_binding = dscore_binding;
			reversion.dspec_bound = dspec_bound;
			reversion.dspec_binding = dspec_binding;
		}
		// sort reversions, and then keep the first acceptable reversion from the current round
		// (there could have been multiple 'acceptable' reversions)
		std::sort( reversions.begin(), reversions.end() );
		auto rev( reversions.begin() );
		for ( Reversions::const_iterator end( reversions.end() ); rev != end; ++rev ) {
			// ignore reversion if it results in the loss of too much binding energy or specificity
			if ( rev->dscore_binding > dscore_cutoff || rev->dspec_binding < dspec_cutoff ) continue;
			// make 'best' reversion 'permanent'
			Size const index( rev->index );
			ResidueTypeCOP starting_type( fixed_residue_types[ index ] ),
				reference_type( reference_residue_types_[ index ] ); // 'reference' == 'native'
			fixed_residue_types[ index ] = reference_type;
			TR << "(round " << round << ") Reversion from " << starting_type->name3()
				<< " to " << reference_type->name3() << " at ";
			if ( pose.pdb_info() ) {
				TR << pose.pdb_info()->chain( index ) << "." << pose.pdb_info()->number( index );
			} else {
				TR << pose.chain( index ) << "." << index;
			}
			TR << " is acceptable and is now fixed.\n";
			break;
		}

		// remove from list of remaining reversions to try
		if ( rev != reversions.end() ) reversions.erase( rev );
		// no acceptable reversion remains, stop scan
		else {
			TR << "No (more) acceptable reversions found." << std::endl;
			break;
		}
		// repeat whole process: mst re-assess the rest of the reversions in the new context
		++round;
	}
	// repack one last time in case the last reversion tried was not acceptable
	utility::vector0<int> rot_to_pack;
	restrict_to_single_sequence( rotamer_sets(), fixed_residue_types, rot_to_pack );
	// calls PackRotamersMover method
	run( pose, rot_to_pack );
	if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
		AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
	}
}

std::string
DnaInterfacePacker::allowed_types() const
{
	std::string protein("ACDEFGHIKLMNPQRSTVWY");
	if ( allowed_types_ == "protein" )  return protein;
	if ( allowed_types_ == "PROTEIN" )  return protein;
	if ( allowed_types_ == "std" )      return protein;
	if ( allowed_types_ == "standard" ) return protein;
	if ( allowed_types_ == "" )         return protein;
	return allowed_types_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief runs a single-residue scan of user-defined amino acid possibilities to estimate affinity and specificity of single mutants w/ respect to relevant DNA
/// @author ashworth
void
DnaInterfacePacker::protein_scan( Pose & pose )
{
	std::string const typestring( allowed_types() );
	TR << "Starting protein_scan with allowed types " << typestring << "." << std::endl;
	// parse allowed_types string into residue types
	ResidueTypeCOPs allowed_type_caps;
	ResidueTypeSetCOP rts( pose.residue(1).residue_type_set() );
	for ( char typechar : typestring ) {
		ResidueTypeCOP aa_type( rts->get_representative_type_aa( aa_from_oneletter_code( typechar ) ) );
		if ( ! aa_type ) {
			TR(t_warning) << "no ResidueType found in ResidueTypeSet for " << typechar << std::endl;
			runtime_assert(false);
		}
		allowed_type_caps.push_back( aa_type );
	}
	runtime_assert( !allowed_type_caps.empty() );

	// get list of positions to scan from PackerTask
	std::list< Size > scan_positions;
	for ( Size index(1); index <= task()->total_residue(); ++index ) {
		if ( !task()->design_residue( index ) ) continue;
		if ( pose.residue_type( index ).is_DNA() ) continue; // skip DNA
		scan_positions.push_back( index );
	}

	if ( option[ OptionKeys::dna::design::checkpoint ].user() ) {
		// skip previously completed positions found in log file
		utility::io::izstream file;
		std::string filename( filename_root_ + ".protein_scan" );
		file.open( filename.c_str() );
		if ( file ) {
			TR << "Reading existing (incomplete?) protein scan results file\n";
			// remove completed positions as they are found in file
			std::string line;
			while ( file.getline( line ) ) {
				utility::vector1< std::string > words( string_split( line ) );
				if ( words.front() != "Done" ) continue; // skip to next line
				// expected line format: "Done scanning at index #"
				if ( words.size() < 5 ) continue;
				std::istringstream ss_index( words.back() );
				Size index;
				ss_index >> index;
				TR << "skipping previously completed scan position " << index << '\n';
				scan_positions.remove( index );
			}
			file.close();
		}
	}

	// open file for output
	std::string outfilename( filename_root_ + ".protein_scan" );
	// APPPEND to file (note: must be careful when parsing results)
	utility::io::ozstream outfile( outfilename.c_str(), std::ios::app );
	if ( !outfile ) {
		std::cerr << "trouble opening file " << outfilename << " for writing" << std::endl;
		assert( false ); // die here in debug mode
		return;
	}

#ifndef WIN32
	// this statement causes a build error on Windows.
	outfile << std::showpoint << std::fixed << std::setprecision(PRECISION);
#endif

	Pose const input_pose( pose );

	vector1< ResidueTypeCOP > pose_residue_types;
	for ( Size index(1), end( pose.size() ); index <= end; ++index ) {
		pose_residue_types.push_back( pose.residue_type(index).get_self_ptr() );
	}

	// get native energies
	Real best_native_score(0.), best_native_dG(0.);
	std::pair< Real, Real > best_native_specificities;

	// repack [and minimize] native interface in the relevant region, measure specificity
	{ // scope
		Pose best_pose( pose );
		for ( Size trial(1); trial <= num_repacks_; ++trial ) {
			// repack bound pose
			utility::vector0<int> native_rot_to_pack;
			restrict_to_single_sequence( rotamer_sets(), pose_residue_types, native_rot_to_pack );
			// calls PackRotamersMover method
			run( pose, native_rot_to_pack );

			if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
				AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
			}
			// native bound score
			Real const native_score( ( *score_function() )( pose ) );

			if ( trial == 1 || native_score < best_native_score ) {
				best_native_score = native_score;
				// native binding score (bound - unbound)
				if ( binding_E_ ) best_native_dG = native_score - unbound_score( pose );
				best_native_specificities  = measure_specificity( pose );
				best_pose = pose;
			}
		}
		// the interface is now repacked [and minimized] in the relevant region for subsequent calculations
		pose = best_pose;

		outfile << "Scanning protein positions that interface with DNA position(s) "
			<< dna_seq_tag( pose, current_working_sequence( pose ) ) << '\n';
		outfile << "Using native scores from best trial: " << "bound = " << best_native_score
			<< ", binding = " << best_native_dG << ", specificity.bound = "
			<< best_native_specificities.first << ", specificity.binding = "
			<< best_native_specificities.second << '\n';

	} // end scope

	// for each designing protein residue (in the interface)
	for ( std::list< Size >::const_iterator index( scan_positions.begin() ),
			end( scan_positions.end() ); index != end; ++index ) {

		outfile << "current designable residues are";
		for ( Size i(1); i <= task()->total_residue(); ++i ) {
			if ( !task()->design_residue(i) ) continue;
			if ( pose.pdb_info() ) {
				outfile << " " << pose.pdb_info()->chain(i) << "." << pose.pdb_info()->number(i);
			} else {
				outfile << " " << pose.chain(i) << "." << i;
			}
			outfile << " " << dna_full_name3( pose.residue_type(i).name3() );
		}
		outfile << '\n';
		// save the native type
		ResidueTypeCOP native_type( pose_residue_types[ *index ] );

		for ( ResidueTypeCOPs::const_iterator scan_type( allowed_type_caps.begin() );
				scan_type != allowed_type_caps.end(); ++scan_type ) {

			// ensure that this type was allowed by the PackerTask/RotamerSets/I.G. before proceeding
			ResidueLevelTask const & rtask( task()->residue_task(*index) );
			ResidueLevelTask::ResidueTypeCOPList const & art( rtask.allowed_residue_types() );
			// maybe this should be a ResidueLevelTask method, and maybe the ResidueLevelTask should clear
			// its allowed_residue_types if !being_packed
			if ( !rtask.being_packed() ) {
				TR << "packing was disabled at " << pose.pdb_info()->chain(*index) << "."
					<< pose.pdb_info()->number(*index) << std::endl;
				runtime_assert(false);
			}
			if ( std::find( art.begin(), art.end(), *scan_type ) == art.end() ) {
				TR << (*scan_type)->name() << " not allowed at " << pose.pdb_info()->chain(*index) << "."
					<< pose.pdb_info()->number(*index) << std::endl;
				//runtime_assert(false);
				continue;
			}
			// set scan type at this residue
			pose_residue_types[ *index ] = *scan_type;

			Real best_score(0.), best_dG(0.);
			std::pair< Real, Real > best_specificities;

			for ( Size trial(1); trial <= num_repacks_; ++trial ) {
				utility::vector0<int> rot_to_pack;
				restrict_to_single_sequence( rotamer_sets(), pose_residue_types, rot_to_pack );
				// calls PackRotamersMover method
				run( pose, rot_to_pack );
				if ( min_movemap_ != nullptr && minimize_options_ != nullptr ) {
					AtomTreeMinimizer().run( pose, *min_movemap_, *score_function(), *minimize_options_ );
				}
				Real const score( ( *score_function() )( pose ) );
				if ( trial == 1 || score < best_score ) {
					best_score = score;
					if ( binding_E_ ) best_dG = score - unbound_score( pose );
					best_specificities = measure_specificity( pose );
				}
			}
			outfile << "Scores for mutation to " << pose_residue_types[ *index ]->name3() << " at ";
			if ( pose.pdb_info() ) {
				outfile << pose.pdb_info()->chain( *index ) << "." << pose.pdb_info()->number( *index );
			} else {
				outfile << pose.chain( *index ) << "." << *index;
			}
			outfile << "." <<  native_type->name() << ":" << " bound = " << best_score << " ("
				<< best_score - best_native_score << ")" << ", binding = " << best_dG <<  " ("
				<< best_dG - best_native_dG << ")" << ", specificity.bound = "
				<< best_specificities.first << " ("
				<< best_specificities.first - best_native_specificities.first << ")"
				<< ", specificity.binding = " << best_specificities.second
				<< " (" << best_specificities.second - best_native_specificities.second << ")\n";
		}
		// restore pose_residue_types to native state at this position
		pose_residue_types[ *index ] = native_type;
		outfile << "Done scanning at index " << *index << '\n';
	}
	outfile.close();
	// rename output file so that future runs do not read/append to it accidentally
	std::string newname( outfilename + ".done" );
	std::rename( outfilename.c_str(), newname.c_str() );

	pose = input_pose;
	TR.flush();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief makes hard copy to guarantee that the reference pose isn't changed from elsewhere
void DnaInterfacePacker::reference_pose( Pose const & pose ) { reference_pose_ = PoseCOP( PoseOP( new Pose( pose ) ) ); }
PoseCOP DnaInterfacePacker::reference_pose() const { return reference_pose_; }

void DnaInterfacePacker::targeted_dna( DnaDesignDefOPs const & defs ) { targeted_dna_ = defs; }
DnaDesignDefOPs const & DnaInterfacePacker::targeted_dna() const { return targeted_dna_; }

void DnaInterfacePacker::pdboutput( PDBOutputOP pdboutput ) { pdboutput_ = pdboutput; }

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details similar to basic::dna_seq_str, but returns only top stranded sequence, delimited by "_". (safe for filenames)
/// @author ashworth
std::string
DnaInterfacePacker::dna_seq_tag( Pose const & pose, ResTypeSequence const & sequence ) const
{
	std::ostringstream ss;
	bool sep(false);
	for ( auto const & pos : sequence ) {
		Size const seqpos( pos.first );
		if ( !dna_chains_->is_top( seqpos ) ) continue;
		if ( sep ) ss << "_";
		if ( pose.pdb_info() ) {
			ss << pose.pdb_info()->chain( seqpos ) << "." << pose.pdb_info()->number( seqpos );
		} else {
			ss << pose.chain( seqpos ) << "." << seqpos;
		}
		ss << "." << dna_full_name3( pos.second->name3() );
		sep = true;
	}
	return ss.str();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details returns DNA sequence that the PackerTask is configured to target (if any)
/// @author ashworth
ResTypeSequence
DnaInterfacePacker::get_targeted_sequence( Pose const & pose ) const
{
	runtime_assert( task() != nullptr );
	runtime_assert( dna_chains_ != nullptr );
	ResTypeSequence sequence;
	for ( DnaPositions::const_iterator it( dna_chains_->begin() ), end( dna_chains_->end() );
			it != end; ++it ) {
		DnaPosition const & pos( it->second );
		Size const resid( pos.top() );
		ResidueLevelTask const & rtask( task()->residue_task(resid) );
		if ( rtask.has_behavior("TARGET") ) {
			// dna position whose PackerTask behavior is 'targeted'
			// if PackerTask indicates nucleotide type to target, add this to local targeted sequence
			if ( rtask.target_type() != nullptr ) sequence[ resid ] = rtask.target_type();
			// otherwise use the residue type at this position in the current pose
			else sequence[ resid ] = pose.residue_type( resid ).get_self_ptr();

			if ( !pos.paired() ) continue;
			// similar treatment for paired 'lower-strand' nucleotides
			Size const comp_resid( pos.bottom() );
			ResidueLevelTask const & comp_rtask( task()->residue_task(comp_resid) );
			if ( comp_rtask.target_type() != nullptr ) sequence[ comp_resid ] = comp_rtask.target_type();
			else sequence[ comp_resid ] = pose.residue_type( comp_resid ).get_self_ptr();
		}
	}
	return sequence;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief current TOP-STRANDED DNA sequence of the pose, at PackerTask's 'targeted' or 'scan' positions
/// @author ashworth
ResTypeSequence
DnaInterfacePacker::current_working_sequence( Pose const & pose ) const
{
	ResTypeSequence current_sequence;
	for ( DnaPositions::const_iterator it( dna_chains_->begin() ); it != dna_chains_->end(); ++it ) {
		if ( ! it->second.paired() ) continue; // skip unpaired DNA positions
		Size const resid( it->first );
		ResidueLevelTask const & rtask( task()->residue_task(resid) );
		if ( !rtask.has_behavior("TARGET") && !rtask.has_behavior("SCAN") ) continue;
		current_sequence[ resid ] = pose.residue( resid ).type().get_self_ptr();
	}
	return current_sequence;
}

std::string
DnaInterfacePacker::current_dna_design_string( Pose const & pose ) const
{
	return dna_seq_tag( pose, current_working_sequence( pose ) );
}

} // namespace dna
} // namespace protocols
