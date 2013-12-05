// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_HelixAssembler
/// @detailed
/// @author Rhiju Das

#include <protocols/rna/RNA_HelixAssembler.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/vector1.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/io/pdb/file_data.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <basic/database/open.hh>

//Minimizer stuff
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <numeric/xyzVector.hh>

#include <core/types.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

// External library headers
#include <numeric/random/random.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>

#include <core/pose/annotated_sequence.hh>

using namespace core;
using ObjexxFCL::string_of;
using basic::T;

static basic::Tracer TR( "protocols.rna.RNA_HelixAssembler" ) ;

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace rna {

static numeric::random::RandomGenerator RG(17720);  // <- Magic number, do not change it!

// silly helper function
bool
is_blank_seq( char const & c ){

	static utility::vector1< char > blank_chars;
	static bool init( false );

	if ( !init ){
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
	rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA ) ),
	torsion_info_(),
	perturb_amplitude_( 10.0 ),
	scorefxn_( core::scoring::ScoreFunctionFactory::create_score_function( "rna/rna_helix" ) ),
	model_and_remove_capping_residues_( true ),
	capping_residues_( "gc" )
{
	Mover::type("RNA_HelixAssembler");
	//	scorefxn_->set_weight( core::scoring::atom_pair_constraint, 5.0 );
	//	scorefxn_->set_weight( core::scoring::rna_torsion, 20.0 );
	// why are we getting clashes?
	//	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
	//	scorefxn_->set_weight( core::scoring::fa_intra_rep, 0.1 );
	initialize_minimizer();
}

RNA_HelixAssembler::~RNA_HelixAssembler(){}

protocols::moves::MoverOP RNA_HelixAssembler::clone() const
{
	return new RNA_HelixAssembler(*this);
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
	if (setting) {
		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna_phenix" );
	} else {
		rsd_set_ = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );
	}
}

////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::apply( core::pose::Pose & pose, std::string const & full_sequence )
{

	// figure out if there are dangling ends -- just extract helix portion.
	std::string const sequence_helix = figure_out_and_remove_dangling_ends( full_sequence );

	TR << "Will build helix sequence: " << sequence_helix << std::endl;

	// build helix portion only.
	build_helix( pose, sequence_helix );

	// build on dangling ends.
	build_dangling_ends( pose );

	fill_chain_info( pose, full_sequence );

}

////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_helix( core::pose::Pose & pose, std::string const & full_sequence ){

	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::id;

	Size seq_length( full_sequence.size() );
	std::string sequence1( full_sequence.substr( 0, seq_length/2 ) );
	std::string sequence2( full_sequence.substr( seq_length/2, seq_length ) );

	if 	( model_and_remove_capping_residues_ ){
		TR << "adding capping residues (g-c base pair)" << std::endl;
		sequence1 = capping_residues_[0] + sequence1 + capping_residues_[1];
		sequence2 = capping_residues_[0] + sequence2 + capping_residues_[1]; // will remove these base pairs later.
	}


  TR << "Sequence of first strand:  " << sequence1 << std::endl;
	TR << "Sequence of second strand: " << sequence2 << std::endl;

	Size const numres = sequence1.size();
	assert( sequence2.size() == numres );

	pose = *build_init_pose( sequence1.substr( 0, 1 ), sequence2.substr( numres-1, numres ) );

	//	pose.dump_pdb( "helix_init.pdb" );

	for ( Size n = 2; n <= numres; n++ ) {
		TR << "Building on base pair: " << n << std::endl;
		build_on_base_pair( pose, n, sequence1[n-1], sequence2[numres-n]);
		if ( dump_ ) pose.dump_pdb( "helix_extend" + string_of(n) + ".pdb" );

		put_constraints_on_base_step( pose, n );
        //TR << "Finished adding constraints..." << std::endl;

		minimize_base_step( pose, n, scorefxn_ );

		if ( dump_ ) pose.dump_pdb( "helix_min" + string_of(n) + ".pdb" );
	}

	if ( finish_scorefxn_ > 0 )  minimize_base_step( pose, 0, finish_scorefxn_ );

	if 	( model_and_remove_capping_residues_ ){
		TR << "removing capping residues" << std::endl;
		get_rid_of_capping_base_pairs( pose );
	}

	if ( finish_scorefxn_ ) (*finish_scorefxn_)( pose );

}


//////////////////////////////////////////////////////
pose::PoseOP
RNA_HelixAssembler::build_init_pose( std::string const & seq1, std::string const & seq2 ){
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::chemical::rna;

	PoseOP pose( new Pose() );
	make_pose_from_sequence( *pose, seq1 + seq2,	*rsd_set_ );
	for ( Size i = 1; i <= pose->n_residue(); ++i ) set_Aform_torsions( *pose, i );
	if ( seq1 == "" || seq2 == "" ) return pose;
	Size const len1( get_sequence_len( seq1 ) );

	FoldTree f( pose->n_residue() );
	f.new_jump( 1, pose->n_residue(), len1 );
	f.set_jump_atoms(
		1,
		chi1_torsion_atom( pose->residue( 1 ) ),
		chi1_torsion_atom( pose->residue( pose->n_residue() ) )
	);
	pose->fold_tree( f );

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
	core::io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set_, path );
	FoldTree f1( 2 );
	f1.new_jump( 1, 2, 1 );
	f1.set_jump_atoms(
		1,
		chi1_torsion_atom( ref_pose.residue( 1 ) ),
		chi1_torsion_atom( ref_pose.residue( 2 ) )
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
RNA_HelixAssembler::build_on_base_pair( pose::Pose & pose, Size const & n, char const & seq1, char const & seq2 ) const {

 	using namespace core::conformation;
 	using namespace core::chemical;
 	using namespace core::id;

	append_Aform_residue(  pose, n - 1, seq1 );
	prepend_Aform_residue( pose, n + 1, seq2 );
}


/////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::initialize_minimizer(){

 	using namespace core::optimization;

	minimizer_ = new AtomTreeMinimizer;

	float const dummy_tol( 0.0000025);
	bool const use_nblist( false );
	minimizer_options_ = new MinimizerOptions( "dfpmin", dummy_tol, use_nblist, false, false );
	minimizer_options_->nblist_auto_update( true );
}

/////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_base_step( pose::Pose & pose, Size const n, core::scoring::ScoreFunctionOP scorefxn ) const {

	using namespace core::scoring;
	using namespace core::optimization;

	Real const cst_weight = scorefxn->get_weight( atom_pair_constraint );
	//runtime_assert( cst_weight != 0 );

	kinematics::MoveMap mm;

	if ( minimize_all_ || n == 0 ){
		mm.set_bb( true );
		mm.set_chi( true );
		if ( minimize_jump_ )	mm.set_jump( true );
	} else {
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );
		for (Size i = n-1; i <= n+2; i++ ) {
			mm.set_bb( i, true );
			mm.set_chi( i, true );
		}
	}

	minimizer_->run( pose, mm, *scorefxn, *minimizer_options_ );

	scorefxn->set_weight( atom_pair_constraint, 0.0 );
	minimizer_->run( pose, mm, *scorefxn, *minimizer_options_ );

	scorefxn->set_weight( atom_pair_constraint, cst_weight );
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
RNA_HelixAssembler::get_rid_of_capping_base_pairs( pose::Pose & pose ) const {

	using namespace core::kinematics;

	//Need to fix up fold_tree
	Size const nres = pose.total_residue();
	FoldTree f( nres );
	f.new_jump( nres/2-1, nres/2+2, nres/2);  // put jump across second-to-last base pair
	f.set_jump_atoms( 1,
										core::chemical::rna::chi1_torsion_atom( pose.residue(nres/2-1) ),
										core::chemical::rna::chi1_torsion_atom( pose.residue(nres/2+2) ) );
	f.reorder( nres/2-1 ); // root cannot be deleted, I think, so place it at this second-to-last base pair.
	pose.fold_tree( f );

	runtime_assert( get_cutpoint( pose ) == nres/2 );

 // get rid of last base pair
	pose.delete_polymer_residue( nres/2 );
	pose.delete_polymer_residue( nres/2 );

 // get rid of first base pair
	pose.delete_polymer_residue( 1 );
	pose.delete_polymer_residue( pose.total_residue() );

	runtime_assert( get_cutpoint( pose ) == pose.total_residue()/2 );

}


////////////////////////////////////////////////////////////////////////////////////////
// dangles are defined by nucleotides at either helix end which are "base-paired" to
// gaps. Treat them separately.
//
std::string
RNA_HelixAssembler::figure_out_and_remove_dangling_ends( std::string const & full_sequence ){

	Size const seq_length( full_sequence.size() );
	std::string const sequence1( full_sequence.substr( 0, seq_length/2 ) );
	std::string const sequence2( full_sequence.substr( seq_length/2, seq_length ) );

	std::string sequence_helix1( sequence1 );
	std::string sequence_helix2( sequence2 );

	dangle_seq1_5prime_ = "";
	dangle_seq1_3prime_ = "";
	dangle_seq2_5prime_ = "";
	dangle_seq2_3prime_ = "";

	if ( is_blank_seq( sequence_helix1[ 0 ] ) ){
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) );
		dangle_seq2_3prime_ = sequence_helix2[ sequence2.size()-1 ];
		sequence_helix1 = sequence_helix1.substr( 1,  sequence_helix1.size()-1 );
		sequence_helix2 = sequence_helix2.substr( 0,  sequence_helix2.size()-1 );
	}

	if ( is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) ){
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix2[ 0 ] ) );
		dangle_seq2_5prime_ = sequence_helix2[ 0 ];
		sequence_helix1 = sequence_helix1.substr( 0,  sequence_helix1.size()-1 );
		sequence_helix2 = sequence_helix2.substr( 1,  sequence_helix2.size()-1 );
	}

	if ( is_blank_seq( sequence_helix2[ 0 ] ) ){
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) );
		dangle_seq1_3prime_ = sequence_helix1[ sequence_helix1.size()-1  ];
		sequence_helix1 = sequence_helix1.substr( 0,  sequence_helix1.size()-1 );
		sequence_helix2 = sequence_helix2.substr( 1,  sequence_helix2.size()-1 );
	}

	if ( is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) ){
		// can't have double gaps...
		runtime_assert( !is_blank_seq( sequence_helix1[ 0 ] ) );
		dangle_seq1_5prime_ = sequence_helix1[ 0  ];
		sequence_helix1 = sequence_helix1.substr( 1,  sequence_helix1.size()-1 );
		sequence_helix2 = sequence_helix2.substr( 0,  sequence_helix2.size()-1 );
	}

	// can't have more than one dangle at the moment.
	runtime_assert( !is_blank_seq( sequence_helix1[ 0 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix2[ 0 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix1[ sequence_helix1.size()-1 ] ) );
	runtime_assert( !is_blank_seq( sequence_helix2[ sequence_helix2.size()-1 ] ) );

	runtime_assert( ! ( dangle_seq1_5prime_.size()>0 && dangle_seq2_3prime_.size()>0 ) );
	runtime_assert( ! ( dangle_seq1_3prime_.size()>0 && dangle_seq2_5prime_.size()>0 ) );

	std::string full_sequence_helix = sequence_helix1  + sequence_helix2;
	return full_sequence_helix;

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
	prepend_Aform_residue( pose, 1, dangle_seq[ 0 ] );
	minimize_prepend_res( pose, 1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq2_5prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const n = get_cutpoint( pose ); // boundary between two strands
	prepend_Aform_residue( pose, n+1, dangle_seq[ 0 ] );
	minimize_prepend_res( pose, n+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq1_3prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const n = get_cutpoint( pose ); // boundary between two strands
	append_Aform_residue( pose, n, dangle_seq[ 0 ] );
	minimize_append_res( pose, n+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::build_dangle_seq2_3prime( pose::Pose & pose, std::string const & dangle_seq ) const {
	if ( dangle_seq.size() == 0 ) return;
	Size const nres = pose.total_residue();
	append_Aform_residue( pose, nres, dangle_seq[ 0 ] );
	minimize_append_res( pose, nres+1 );
}

////////////////////////////////////////////////////////////////////////////////////////
Size
RNA_HelixAssembler::get_cutpoint( pose::Pose const & pose ) const {

	Size n( 1 );
	for ( n = 1; n < pose.total_residue(); n++ ) if ( pose.fold_tree().is_cutpoint(n) ) break;
	runtime_assert( pose.fold_tree().is_cutpoint( n ) ) ;

	return n;
}

////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::append_Aform_residue( pose::Pose & pose, Size const & n, char const & nt ) const {

 	using namespace core::conformation;
 	using namespace core::chemical;
 	using namespace core::id;

	runtime_assert( pose.fold_tree().is_cutpoint(n) || n == pose.total_residue() );

	/////////////////////////////////////
	ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set_->aa_map( aa_from_oneletter_code( nt ) )[1] ) ) );
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
		id_list.push_back( TorsionID( n+1, CHI, 1) );
		perturb_torsion(pose, id_list);
	}
}


////////////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::prepend_Aform_residue( pose::Pose & pose, Size const & n, char const & nt ) const {

 	using namespace core::conformation;
 	using namespace core::chemical;
 	using namespace core::id;

	runtime_assert( n == 1 || pose.fold_tree().is_cutpoint(n-1) );

	ResidueOP rsd2( ResidueFactory::create_residue( *(rsd_set_->aa_map( aa_from_oneletter_code( nt ) )[1] ) ) );
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
		id_list.push_back( TorsionID( n, CHI, 1) );
		id_list.push_back( TorsionID( n+1, BB, 1) );
		id_list.push_back( TorsionID( n+1, BB, 2) );
		id_list.push_back( TorsionID( n+1, BB, 3) );
		perturb_torsion(pose, id_list);
	}
}

//////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::perturb_torsion(
	pose::Pose & pose,
	utility::vector1<id::TorsionID> const & id_list
) const {
	for (Size i = 1; i <= id_list.size(); ++i) {
		pose.set_torsion( id_list[i], pose.torsion(id_list[i]) + perturb_amplitude_ * RG.gaussian() );
	}
}
//////////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_append_res( pose::Pose & pose, Size const n ) const {

	kinematics::MoveMap mm;

	runtime_assert( n > 1 );
	if ( minimize_all_  ){
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

	runtime_assert( n < pose.total_residue() );
	if ( minimize_all_  ){
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

///////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::fill_chain_info( pose::Pose & pose, std::string const & full_sequence ){

	using namespace core::pose;

	utility::vector1< char > chains;

	Size const nres = full_sequence.size();

	for ( Size i = 1; i <= (nres/2); i++ ){
		if ( !is_blank_seq( full_sequence[i-1] ) ) chains.push_back( 'A' );
	}

	for ( Size i = (nres/2 + 1); i <= nres; i++ ){
		if ( !is_blank_seq( full_sequence[i-1] ) ) chains.push_back( 'B' );
	}

	PDBInfoOP pdb_info = new PDBInfo( pose );
	pdb_info->set_chains( chains );
	pdb_info->obsolete( false ); // this is silly.

	pose.pdb_info( pdb_info );

}




} // namespace rna
} // namespace protocols
