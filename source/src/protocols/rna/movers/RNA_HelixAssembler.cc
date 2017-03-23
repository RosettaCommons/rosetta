// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_HelixAssembler.cc
/// @brief make a helix that is low energy in Rosetta
/// @details Adds ideal base pairs one at a time and minimizes.
/// @author Rhiju Das

#include <protocols/rna/movers/RNA_HelixAssembler.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/sequence/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/rna/util.hh>
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/types.hh>

// External library headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>


//C++ headers
#include <vector>
#include <string>
#include <sstream>

#include <core/pose/annotated_sequence.hh>

///////////////////////////////////////////////////////////////////////////////////////
//
// Purpose: create an RNA helix with A-form torsions, but low in Rosetta energy.
//
//  * Uses a hacky weights file with high rna_torsion, otherwise short helices gets wonky.
//  * Updated to allow modeling of bulges & dangling ends.
//  * Typically called through rna_helix.py available in Rosetta/tools/rna_tools/
//
// Note: does *not* currently guarantee symmetry, i.e. 5'-cc-3'/5'-gg-3' is not
//       the exact same as 5'-gg-3'/5'-cc-3'. This is due to buildup strategy,
//       which guarantees low energy but not much else.
//
///////////////////////////////////////////////////////////////////////////////////////

using namespace core;
using ObjexxFCL::string_of;
using basic::T;

static THREAD_LOCAL basic::Tracer TR( "protocols.rna.movers.RNA_HelixAssembler" );

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace rna {
namespace movers {


// silly helper function
bool
is_blank_seq( char const & c ){

	static utility::vector1< char > blank_chars;
	static bool init( false );

	if ( !init ) {
		blank_chars.push_back( '-' );
		blank_chars.push_back( '0' );
		blank_chars.push_back( '.' );
		init = true;
	}

	return blank_chars.has_value( c );
}

RNA_HelixAssembler::RNA_HelixAssembler():
	dump_( false ),
	random_perturbation_( false ),
	minimize_all_( false ),
	minimize_jump_( false ),
	use_phenix_geo_( false ),
	rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) ),
	torsion_info_(),
	perturb_amplitude_( 10.0 ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "stepwise/rna/rna_helix" ) ),
	model_and_remove_capping_residues_( true ),
	capping_residues_( "gc" ),
	full_sequence_( "" )
{
	Mover::type("RNA_HelixAssembler");
	initialize_minimizer();
}

RNA_HelixAssembler::~RNA_HelixAssembler(){}

protocols::moves::MoverOP RNA_HelixAssembler::clone() const
{
	return protocols::moves::MoverOP( new RNA_HelixAssembler(*this) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details  Apply the RNA HelixAssembler
///
void RNA_HelixAssembler::apply( core::pose::Pose & pose )
{
	apply( pose, pose.sequence() );
}

std::string
RNA_HelixAssembler::get_name() const {
	return "RNA_HelixAssembler";
}

void RNA_HelixAssembler::set_scorefxn( core::scoring::ScoreFunctionOP setting )
{
	scorefxn_ = setting;
}


void RNA_HelixAssembler::set_finish_scorefxn( core::scoring::ScoreFunctionOP setting )
{
	finish_scorefxn_ = setting;
}


void RNA_HelixAssembler::use_phenix_geo( bool const setting )
{
	use_phenix_geo_ = setting;
	if ( setting ) {
		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	} else {
		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	}
}

////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::apply( core::pose::Pose & pose, std::string const & full_sequence_in )
{
	// handle weird non-naturals, including g[IGU] for isoguanosine, etc.
	full_sequence_ = full_sequence_in;
	non_standard_residues_ = core::sequence::parse_out_non_standard_residues( full_sequence_ );

	// figure out if there are dangling ends -- just extract helix portion.
	figure_out_and_remove_dangling_ends();

	TR << "Will build helix sequence: " << full_sequence_ << std::endl;

	// build helix portion only.
	build_helix( pose );

	// build on dangling ends.
	build_dangling_ends( pose );

	fill_chain_info( pose );
}


void
RNA_HelixAssembler::add_capping_base_pairs_to_full_sequence()
{
	Size L( full_sequence_.size() );
	std::string sequence1( full_sequence_.substr( 0, L/2 ) );
	std::string sequence2( full_sequence_.substr( L/2, L/2 ) );

	TR << "adding capping residues (g-c base pair)" << std::endl;
	sequence1 = capping_residues_[0] + sequence1 + capping_residues_[1];
	sequence2 = capping_residues_[0] + sequence2 + capping_residues_[1]; // will remove these base pairs later.
	full_sequence_ = sequence1 + sequence2;

	// full name info for isoguanosine, etc.
	std::map< Size, std::string > non_standard_residues_new;
	for ( auto const & elem : non_standard_residues_ ) {
		Size res( elem.first );
		if ( res <= L/2 ) {
			non_standard_residues_new[ res + 1 ] = elem.second;
		} else {
			non_standard_residues_new[ res + 3 ] = elem.second;
		}
	}
	non_standard_residues_ = non_standard_residues_new;
}

////////////////////////////////////////////////////////////////////////////////
std::string
RNA_HelixAssembler::get_sequence( core::Size const n ){
	std::string rsd( "" );
	rsd += full_sequence_[ n-1 ];
	if ( non_standard_residues_.find( n ) == non_standard_residues_.end() ) return rsd;
	return rsd + "[" + non_standard_residues_[n] + "]";
}

////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_helix( core::pose::Pose & pose ){

	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::id;

	if  ( model_and_remove_capping_residues_ ) add_capping_base_pairs_to_full_sequence();

	Size const numres = full_sequence_.size() / 2;

	pose = *build_init_pose( get_sequence( 1 ), get_sequence( 2 * numres ) );

	for ( Size n = 2; n <= numres; n++ ) {
		TR << "Building on base pair: " << n << std::endl;
		build_on_base_pair( pose, n, get_sequence( n ), get_sequence( 2 * numres - n + 1 ) );
		if ( dump_ ) pose.dump_pdb( "helix_extend" + string_of(n) + ".pdb" );

		put_constraints_on_base_step( pose, n );
		//TR << "Finished adding constraints..." << std::endl;
		minimize_base_step( pose, n, scorefxn_ );

		if ( dump_ ) pose.dump_pdb( "helix_min" + string_of(n) + ".pdb" );
	}

	if ( finish_scorefxn_ )  minimize_base_step( pose, 0, finish_scorefxn_ );

	if  ( model_and_remove_capping_residues_ ) {
		TR << "removing capping residues" << std::endl;
		get_rid_of_capping_base_pairs( pose );
	}

	if ( finish_scorefxn_ ) ( *finish_scorefxn_)( pose );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// @brief build pose with correct fold-tree for modeling helix, possibly incl. dangling ends.
//
// @details
//
// in use by RECCES -- note that dangles are always at 3'-end of seq1, or 5'-end of seq2,
//  e.g if you specify  -seq1 gu -seq2 aac, you get:
//
//      5'-gu-3'
//         ||
//      3'-caa-5'
//
//   or if you specify  -seq1 gua -seq2 ac, you get:
//
//      5'-gua-3'
//         ||
//      3'-ca-5'
//
// Only used to make initial base pair by rna_helix.cc.
//
pose::PoseOP
RNA_HelixAssembler::build_init_pose( std::string const & seq1, std::string const & seq2 ){
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::chemical::rna;

	PoseOP pose( new Pose() );
	core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ );
	make_pose_from_sequence( *pose, seq1 + seq2, *rsd_set );
	for ( Size i = 1; i <= pose->size(); ++i ) set_Aform_torsions( *pose, i );
	if ( seq1 == "" || seq2 == "" ) return pose;
	Size const len1( get_sequence_len( seq1 ) );

	FoldTree f( pose->size() );
	f.new_jump( 1, pose->size(), len1 );
	f.set_jump_atoms(
		1,
		chi1_torsion_atom( pose->residue_type( 1 ) ),
		chi1_torsion_atom( pose->residue_type( pose->size() ) )
	);
	pose->fold_tree( f );

	// AMW: whoa buddy would this be a problem for NC basepairs?
	// Get reference jump from database
	std::string path( basic::database::full_name("sampling/rna/") );
	std::string const & pose_seq( pose->sequence() );
	std::string first_bp( pose_seq.substr(0, 1) );
	first_bp += pose_seq[pose_seq.size()];
	// Use gc jump for non-WC pairs
	if (
			first_bp != "au" &&
			first_bp != "ua" &&
			first_bp != "gc" &&
			first_bp != "cg"
			) first_bp = "gc";

	path += first_bp + "_std.pdb";
	Pose ref_pose;
	core::io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, path );
	FoldTree f1( 2 );
	f1.new_jump( 1, 2, 1 );
	f1.set_jump_atoms(
		1,
		chi1_torsion_atom( ref_pose.residue_type( 1 ) ),
		chi1_torsion_atom( ref_pose.residue_type( 2 ) )
	);
	ref_pose.fold_tree( f1 );
	Jump const & jump( ref_pose.jump( 1 ) );

	// Now set the jump
	pose->set_jump( 1, jump );

	return pose;
}
//////////////////////////////////////////////////////
void
RNA_HelixAssembler::set_Aform_torsions( pose::Pose & pose, Size const & n ) const
{
	using namespace core::id;
	using namespace core::pose::rna;
	using namespace core::chemical::rna;
	pose.set_torsion( TorsionID( n, BB, 1), torsion_info_.alpha_aform() );
	pose.set_torsion( TorsionID( n, BB, 2), torsion_info_.beta_aform() );
	pose.set_torsion( TorsionID( n, BB, 3), torsion_info_.gamma_aform() );
	pose.set_torsion( TorsionID( n, BB, 4), torsion_info_.delta_north() );
	pose.set_torsion( TorsionID( n, BB, 5), torsion_info_.epsilon_aform() );
	pose.set_torsion( TorsionID( n, BB, 6), torsion_info_.zeta_aform() );

	apply_pucker(pose, n, NORTH, false /*skip_same_state*/, use_phenix_geo_);
}

/////////////////////////////////////////////////
void
RNA_HelixAssembler::build_on_base_pair( pose::Pose & pose, Size const & n, std::string const & seq1, std::string const & seq2 ) const {
	append_Aform_residue(  pose, n - 1, seq1 );
	prepend_Aform_residue( pose, n + 1, seq2 );
	core::pose::full_model_info::update_full_model_info_from_pose( pose );
}


/////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::initialize_minimizer(){

	using namespace core::optimization;

	minimizer_ = core::optimization::AtomTreeMinimizerOP( new AtomTreeMinimizer );

	float const dummy_tol( 0.0000025);
	bool const use_nblist( false );
	minimizer_options_ = core::optimization::MinimizerOptionsOP( new MinimizerOptions( "lbfgs_armijo_nonmonotone", dummy_tol, use_nblist, false, false ) );
	minimizer_options_->nblist_auto_update( true );
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_base_step( pose::Pose & pose, Size const n, core::scoring::ScoreFunctionOP scorefxn ) const {

	using namespace core::scoring;
	using namespace core::optimization;

	if ( !scorefxn_->has_nonzero_weight( core::scoring::base_pair_constraint ) ) {
		TR << TR.Magenta << "Scorefunction does not have base_pair_constraint on -- setting weight to 5.0." << TR.Reset << std::endl;
		scorefxn_->set_weight( core::scoring::base_pair_constraint, 5.0 );
	}

	Real const cst_weight = scorefxn->get_weight( base_pair_constraint );
	//runtime_assert( cst_weight != 0 );

	kinematics::MoveMap mm;

	if ( minimize_all_ || n == 0 ) {
		mm.set_bb( true );
		mm.set_chi( true );
		if ( minimize_jump_ ) mm.set_jump( true );
	} else {
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );
		for ( Size i = n-1; i <= n+2; i++ ) {
			mm.set_bb( i, true );
			mm.set_chi( i, true );
		}
	}

	minimizer_->run( pose, mm, *scorefxn, *minimizer_options_ );

	scorefxn->set_weight( base_pair_constraint, 0.0 );
	minimizer_->run( pose, mm, *scorefxn, *minimizer_options_ );

	scorefxn->set_weight( base_pair_constraint, cst_weight );
	(*scorefxn)( pose );
}

//////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::put_constraints_on_base_step( pose::Pose & pose, Size const & n ) const {
	utility::vector1< std::pair< Size, Size > > pairings;
	pairings.push_back( std::make_pair( n, n+1)  );
	pairings.push_back( std::make_pair( n-1, n+2) );

	scoring::constraints::ConstraintSetOP new_cst_set;
	pose.constraint_set( new_cst_set ); //blank out cst set.

	protocols::rna::setup_base_pair_constraints( pose, pairings );
}

//////////////////////////////////////////////////////
void
RNA_HelixAssembler::get_rid_of_capping_base_pairs( pose::Pose & pose ) {

	using namespace core::kinematics;

	//Need to fix up fold_tree
	Size const nres = pose.size();
	FoldTree f( nres );
	f.new_jump( nres/2-1, nres/2+2, nres/2);  // put jump across second-to-last base pair
	f.set_jump_atoms( 1,
		core::chemical::rna::chi1_torsion_atom( pose.residue_type(nres/2-1) ),
		core::chemical::rna::chi1_torsion_atom( pose.residue_type(nres/2+2) ) );
	f.reorder( nres/2-1 ); // root cannot be deleted, I think, so place it at this second-to-last base pair.
	pose.fold_tree( f );

	runtime_assert( get_cutpoint( pose ) == nres/2 );

	// get rid of last base pair
	pose.delete_polymer_residue( nres/2 );
	pose.delete_polymer_residue( nres/2 );

	// get rid of first base pair
	pose.delete_polymer_residue( 1 );
	pose.delete_polymer_residue( pose.size() );

	runtime_assert( get_cutpoint( pose ) == pose.size()/2 );

	std::string sequence_helix1( full_sequence_.substr( 0,      nres/2 ) );
	std::string sequence_helix2( full_sequence_.substr( nres/2, nres ) );
	remove_first_base_pair( full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );
	remove_last_base_pair(  full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );

	core::pose::full_model_info::update_full_model_info_from_pose( pose );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::remove_first_base_pair( std::string & full_sequence,
	std::map< Size, std::string > & non_standard_residues,
	std::string & sequence_helix1,
	std::string & sequence_helix2 ) const {
	Size L = full_sequence.size()/2;
	runtime_assert( L >= 2 );
	full_sequence = full_sequence.substr( 1, full_sequence.size()-2 );
	sequence_helix1 = sequence_helix1.substr( 1,  L-1 );
	sequence_helix2 = sequence_helix2.substr( 0,  L-1 );

	std::map< Size, std::string > non_standard_residues_new;
	for ( auto const & elem : non_standard_residues ) {
		Size res( elem.first );
		if ( res > 1 || res < full_sequence_.size() ) non_standard_residues_new[ res-1 ] = elem.second;
	}
	non_standard_residues = non_standard_residues_new;

	Size const seq_length( full_sequence.size() );
	runtime_assert( seq_length == 2*L - 2 );
	runtime_assert( sequence_helix1 == full_sequence_.substr( 0, seq_length/2 ) );
	runtime_assert( sequence_helix2 == full_sequence_.substr( seq_length/2, seq_length ) );
}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::remove_last_base_pair( std::string & full_sequence,
	std::map< Size, std::string > & non_standard_residues,
	std::string & sequence_helix1,
	std::string & sequence_helix2 ) const {
	Size L = full_sequence.size()/2;
	runtime_assert( L >= 2 );
	full_sequence = full_sequence.substr( 0, L-1 ) + full_sequence.substr( L+1, L-1 ) ;
	sequence_helix1 = sequence_helix1.substr( 0,  L-1 );
	sequence_helix2 = sequence_helix2.substr( 1,  L-1 );

	std::map< Size, std::string > non_standard_residues_new;
	for ( auto const & elem : non_standard_residues ) {
		Size res( elem.first );
		if ( res < L  ) {
			non_standard_residues_new[ res ] = elem.second;
		} else if ( res > (L + 1) ) {
			non_standard_residues_new[ res - 2 ] = elem.second;
		}
	}
	non_standard_residues = non_standard_residues_new;

	Size const seq_length( full_sequence.size() );
	runtime_assert( seq_length == 2*L - 2 );
	runtime_assert( sequence_helix1 == full_sequence_.substr( 0, seq_length/2 ) );
	runtime_assert( sequence_helix2 == full_sequence_.substr( seq_length/2, seq_length ) );
}

////////////////////////////////////////////////////////////////////////////////////////
// dangles are defined by nucleotides at either helix end which are "base-paired" to
// gaps. Treat them separately.
//
void
RNA_HelixAssembler::figure_out_and_remove_dangling_ends(){

	Size const seq_length( full_sequence_.size() );
	std::string const sequence1( full_sequence_.substr( 0, seq_length/2 ) );
	std::string const sequence2( full_sequence_.substr( seq_length/2, seq_length ) );

	// sequence_helix1 & sequence_helix2 aren't really necessary anymore, but carry them forward
	//  to run consistency checks...
	std::string sequence_helix1( sequence1 );
	std::string sequence_helix2( sequence2 );

	dangle_seq1_5prime_ = "";
	dangle_seq1_3prime_ = "";
	dangle_seq2_5prime_ = "";
	dangle_seq2_3prime_ = "";

	if ( is_blank_seq( sequence_helix1[ 0 ] ) ) {
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) );
		dangle_seq2_3prime_ = get_sequence( full_sequence_.size() ); //sequence_helix2[ sequence2.size()-1 ];
		remove_first_base_pair( full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );
	}

	if ( is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) ) {
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix2[ 0 ] ) );
		dangle_seq2_5prime_ = get_sequence( full_sequence_.size()/2 + 1 ); // sequence_helix2[ 0 ];
		remove_last_base_pair( full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );
	}

	if ( is_blank_seq( sequence_helix2[ 0 ] ) ) {
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) );
		dangle_seq1_3prime_ = get_sequence( full_sequence_.size()/2 ); //sequence_helix1[ sequence_helix1.size()-1  ];
		remove_last_base_pair( full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );
	}

	if ( is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) ) {
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix1[ 0 ] ) );
		dangle_seq1_5prime_ = get_sequence( 1 ); //sequence_helix1[ 0  ];
		remove_first_base_pair( full_sequence_, non_standard_residues_, sequence_helix1, sequence_helix2 );
	}

	// can't have more than one dangle at the moment.
	runtime_assert( !is_blank_seq( sequence_helix1[ 0 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix2[ 0 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) );

	runtime_assert( ! ( dangle_seq1_5prime_.size()>0 && dangle_seq2_3prime_.size()>0 ) );
	runtime_assert( ! ( dangle_seq1_3prime_.size()>0 && dangle_seq2_5prime_.size()>0 ) );

	runtime_assert( sequence_helix1.size() == sequence_helix2.size() );
}


////////////////////////////////////////////////////////////////////////////////////////
// build on dangling ends. [perhaps this could be a different class...]
void
RNA_HelixAssembler::build_dangling_ends( pose::Pose & pose ) const {
	build_dangle_seq1_5prime( pose, dangle_seq1_5prime_ );
	build_dangle_seq2_5prime( pose, dangle_seq2_5prime_ );
	build_dangle_seq1_3prime( pose, dangle_seq1_3prime_ );
	build_dangle_seq2_3prime( pose, dangle_seq2_3prime_ );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq1_5prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	prepend_Aform_residue( pose, 1, dangle_seq );
	minimize_prepend_res( pose, 1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq2_5prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const n = get_cutpoint( pose ); // boundary between two strands
	prepend_Aform_residue( pose, n+1, dangle_seq );
	minimize_prepend_res( pose, n+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq1_3prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const n = get_cutpoint( pose ); // boundary between two strands
	append_Aform_residue( pose, n, dangle_seq );
	minimize_append_res( pose, n+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq2_3prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const nres = pose.size();
	append_Aform_residue( pose, nres, dangle_seq );
	minimize_append_res( pose, nres+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_HelixAssembler::get_cutpoint( pose::Pose const & pose ) const {

	Size n( 1 );
	for ( n = 1; n < pose.size(); n++ ) if ( pose.fold_tree().is_cutpoint(n) ) break;
	runtime_assert( pose.fold_tree().is_cutpoint( n ) ) ;

	return n;
}

////////////////////////////////////////////////////////////////////////////////////////
core::conformation::ResidueOP
RNA_HelixAssembler::get_residue( std::string const & nt ) const {
	// rsd1 =  ResidueFactory::create_residue( *(rsd_set_->aa_map( aa_from_oneletter_code( nt ) )[1] ) );
	return core::conformation::ResidueFactory::create_residue(
		*core::pose::residue_types_from_sequence( nt, *rsd_set_, false /*auto_termini*/ )[1] );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::append_Aform_residue( pose::Pose & pose, Size const & n, std::string const & nt ) const {

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::id;

	runtime_assert( pose.fold_tree().is_cutpoint(n) || n == pose.size() );

	/////////////////////////////////////
	ResidueOP rsd1 = get_residue( nt );
	runtime_assert( rsd1->is_NA() );
	pose.append_polymer_residue_after_seqpos( *rsd1, n, true /*build_ideal_geometry*/ );
	pose.set_torsion( TorsionID( n, BB, 5), torsion_info_.epsilon_aform());
	pose.set_torsion( TorsionID( n, BB, 6), torsion_info_.zeta_aform());
	set_Aform_torsions( pose, n+1 );

	if ( random_perturbation_ ) {
		utility::vector1<TorsionID> id_list;
		id_list.push_back( TorsionID( n, BB, 5) );
		id_list.push_back( TorsionID( n, BB, 6) );
		id_list.push_back( TorsionID( n+1, BB, 1) );
		id_list.push_back( TorsionID( n+1, BB, 2) );
		id_list.push_back( TorsionID( n+1, BB, 3) );
		id_list.push_back( TorsionID( n+1, BB, 4) );
		id_list.push_back( TorsionID( n+1, BB, 5) );
		id_list.push_back( TorsionID( n+1, BB, 6) );
		id_list.push_back( TorsionID( n+1, id::CHI, 1) );
		perturb_torsion(pose, id_list);
	}

	core::pose::full_model_info::update_full_model_info_from_pose( pose );
}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::prepend_Aform_residue( pose::Pose & pose, Size const & n, std::string const & nt ) const {

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::id;

	runtime_assert( n == 1 || pose.fold_tree().is_cutpoint(n-1) );

	ResidueOP rsd2 = get_residue( nt );
	runtime_assert( rsd2->is_NA() );
	pose.prepend_polymer_residue_before_seqpos( *rsd2, n, true /*build_ideal_geometry*/ );

	set_Aform_torsions( pose, n );
	pose.set_torsion( TorsionID( n+1, BB, 1), torsion_info_.alpha_aform());
	pose.set_torsion( TorsionID( n+1, BB, 2), torsion_info_.beta_aform());
	pose.set_torsion( TorsionID( n+1, BB, 3), torsion_info_.gamma_aform());


	// perturb all torsion angles except those at sugar puckers (DELTA, NU1, NU2 )
	if ( random_perturbation_ ) {
		utility::vector1<TorsionID> id_list;
		id_list.push_back( TorsionID( n, BB, 1) );
		id_list.push_back( TorsionID( n, BB, 2) );
		id_list.push_back( TorsionID( n, BB, 3) );
		id_list.push_back( TorsionID( n, BB, 4) );
		id_list.push_back( TorsionID( n, BB, 5) );
		id_list.push_back( TorsionID( n, BB, 6) );
		id_list.push_back( TorsionID( n, id::CHI, 1) );
		id_list.push_back( TorsionID( n+1, BB, 1) );
		id_list.push_back( TorsionID( n+1, BB, 2) );
		id_list.push_back( TorsionID( n+1, BB, 3) );
		perturb_torsion(pose, id_list);
	}

	core::pose::full_model_info::update_full_model_info_from_pose( pose );
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::perturb_torsion(
	pose::Pose & pose,
	utility::vector1<id::TorsionID> const & id_list
) const {
	for ( auto const & id : id_list ) {
		pose.set_torsion( id, pose.torsion(id) + perturb_amplitude_ * numeric::random::rg().gaussian() );
	}
}
//////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_append_res( pose::Pose & pose, Size const n ) const {

	kinematics::MoveMap mm;

	runtime_assert( n > 1 );
	if ( minimize_all_  ) {
		mm.set_bb( true );
		mm.set_chi( true );
		mm.set_jump( true );
	} else {
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		mm.set_bb( n, true );
		mm.set_chi( n, true );

		mm.set_bb( n-1, true ); //perhaps should whittle down to suite residues...
	}

	core::scoring::ScoreFunctionOP scorefxn_minimize = scorefxn_;
	if ( finish_scorefxn_ )  scorefxn_minimize = finish_scorefxn_;
	minimizer_->run( pose, mm, *scorefxn_minimize, *minimizer_options_ );
	(*scorefxn_minimize)( pose );
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_prepend_res( pose::Pose & pose, Size const n ) const {

	kinematics::MoveMap mm;

	runtime_assert( n < pose.size() );
	if ( minimize_all_  ) {
		mm.set_bb( true );
		mm.set_chi( true );
		mm.set_jump( true );
	} else {
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );

		mm.set_bb( n, true );
		mm.set_chi( n, true );

		mm.set_bb( n+1, true ); //perhaps should whittle down to suite residues...
	}

	core::scoring::ScoreFunctionOP scorefxn_minimize = scorefxn_;
	if ( finish_scorefxn_ )  scorefxn_minimize = finish_scorefxn_;
	minimizer_->run( pose, mm, *scorefxn_minimize, *minimizer_options_ );
	(*scorefxn_minimize)( pose );
}

/////////////////////////////////////////////////////////////////////////////////////
// Could make chains distinct ('A' and 'B'), but this was causing
// more confusion than it was worth. Better to have 'default' values (blank chains).
/////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::fill_chain_info( pose::Pose & pose ){

	using namespace core::pose;
	using namespace core::pose::full_model_info;

	utility::vector1< char > chains( pose.size(), ' ' );
	PDBInfoOP pdb_info( new PDBInfo( pose ) );
	pdb_info->set_chains( chains );
	pdb_info->obsolete( false ); // this is silly.

	pose.pdb_info( pdb_info );

	FullModelInfoOP full_model_info( new FullModelInfo( pose ) );
	set_full_model_info( pose, full_model_info );
}


} //movers
} //rna
} //protocols
