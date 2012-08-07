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
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <core/sequence/Sequence.fwd.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_StructureParameters.fwd.hh>
// AUTO-REMOVED #include <protocols/rna/RNA_StructureParameters.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>

//Minimizer stuff
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>


// External library headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <numeric/random/random.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>

#include <core/pose/annotated_sequence.hh>


using namespace core;
using basic::T;

static basic::Tracer TR( "protocols.rna.HelixAssembler" ) ;

typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace rna {

static numeric::random::RandomGenerator RG(17720);  // <- Magic number, do not change it!

RNA_HelixAssembler::RNA_HelixAssembler():
	verbose_( true ),
	random_perturbation_( false ),
	minimize_all_( false ),
	use_phenix_geo_( false ),
	ideal_jump( "RT -0.994805 -0.0315594 0.0967856 -0.0422993 0.992919 -0.111004 -0.092597 -0.114522 -0.989096 6.34696 -0.449942 0.334582 " ),
	rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA ) ),
	ALPHA_A_FORM( -64.11),
	BETA_A_FORM( 176.33),
	GAMMA_A_FORM( 53.08),
	DELTA_A_FORM( 82.90),
	EPSILON_A_FORM( -150.17),
	ZETA_A_FORM( -71.45),
	CHI_A_FORM( 79.43),
	NU2_A_FORM( 38.82),
	NU1_A_FORM( 95.34),
	perturb_amplitude_( 10.0 ),
	scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::RNA_HIRES_WTS ) )
{
	Mover::type("RNA_HelixAssembler");
	// scorefxn->set_weight( core::scoring::atom_pair_constraint, 0.01 );
	// scorefxn->set_weight( core::scoring::rna_torsion, 5.0 );
	scorefxn->set_weight( core::scoring::atom_pair_constraint, 5.0 );
	scorefxn->set_weight( core::scoring::rna_torsion, 20.0 );
	//scorefxn->set_weight( core::scoring::fa_stack, 0.125 );

	// why are we getting clashes?
	//	scorefxn->set_weight( core::scoring::fa_rep, 1.0 );
	//	scorefxn->set_weight( core::scoring::fa_intra_rep, 0.1 );

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

void RNA_HelixAssembler::use_phenix_geo( bool const setting )
{
	use_phenix_geo_ = setting;
	if (setting) {
		core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna_phenix" );
	} else {
		core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );
	}
}
////////////////////////////////////////////////////////////////////////////////
void RNA_HelixAssembler::apply( core::pose::Pose & pose, std::string const & full_sequence )
{

	using namespace core::pose;
	using namespace core::kinematics;

	Size const seq_length( full_sequence.size() );

	std::string const sequence1( full_sequence.substr( 0, seq_length/2 ) );
	std::string const sequence2( full_sequence.substr( seq_length/2, seq_length ) );
	if ( verbose_ ) {
		std::cout << "SEQ1 " << sequence1 << std::endl;
		std::cout << "SEQ2 " << sequence2 << std::endl;
	}
	Size const numres = sequence1.size();
	assert( sequence2.size() == numres );

	Pose pose_scratch;
	core::pose::make_pose_from_sequence( pose_scratch, sequence1.substr(0,1)+sequence2.substr( numres-1,numres) ,	*rsd_set );
	pose = pose_scratch;

	/////////////////////////////////////
	FoldTree f( 2 );
	f.new_jump( 1, 2, 1);
	f.set_jump_atoms( 1,
										core::scoring::rna::chi1_torsion_atom( pose.residue(1) ),
										core::scoring::rna::chi1_torsion_atom( pose.residue(2) ) );
	pose.fold_tree( f );

	/////////////////////////////////////
	//Need sample jumps for a-u, g-c, g-u.
	Jump j;
	std::stringstream jump_stream( ideal_jump );
	jump_stream >> j;
	if (verbose_) std::cout << j << std::endl;

	pose.set_jump( 1, j );

	using namespace core::id;
	set_Aform_torsions( pose, 1 );
	set_Aform_torsions( pose, 2 );
	//	pose.dump_pdb( "helix_init.pdb" );

	for ( Size n = 2; n <= numres; n++ ) {
		std::cout << "Building on base pair: " << n << std::endl;
		build_on_base_pair( pose, n, sequence1[n-1], sequence2[numres-n]);
		//		pose.dump_pdb( "helix_extend"+string_of(n)+".pdb" );

		put_constraints_on_base_step( pose, n );

		minimize_base_step( pose, n );
		//		pose.dump_pdb( "helix_min"+string_of(n)+".pdb" );
	}

}


//////////////////////////////////////////////////////
void
RNA_HelixAssembler::set_Aform_torsions( pose::Pose & pose, Size const & n )
{
	using namespace core::id;
	pose.set_torsion( TorsionID( n, BB, 1),  ALPHA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 2),   BETA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 3),  GAMMA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 4),  DELTA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 5),  EPSILON_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 6),  ZETA_A_FORM);

	if (use_phenix_geo_) {
		ideal_coord_.apply(pose, n);
	} else {
		pose.set_torsion( TorsionID( n, CHI, 1),  CHI_A_FORM);
		pose.set_torsion( TorsionID( n, CHI, 2),  NU2_A_FORM);
		pose.set_torsion( TorsionID( n, CHI, 3),  NU1_A_FORM);
	}
}

/////////////////////////////////////////////////
void
RNA_HelixAssembler::build_on_base_pair( pose::Pose & pose, Size const & n, char const & seq1, char const & seq2 ) {

 	using namespace core::conformation;
 	using namespace core::chemical;
 	using namespace core::id;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	/////////////////////////////////////
	ResidueOP rsd1( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( seq1 ) )[1] ) ) );
	pose.append_polymer_residue_after_seqpos(   *rsd1, n - 1, true /*build_ideal_geometry*/ );


	ResidueOP rsd2( ResidueFactory::create_residue( *(rsd_set->aa_map( aa_from_oneletter_code( seq2 ) )[1] ) ) );
	pose.prepend_polymer_residue_before_seqpos( *rsd2, n + 1, true /*build_ideal_geometry*/ );

	pose.set_torsion( TorsionID( n-1, BB, 5),  EPSILON_A_FORM);
	pose.set_torsion( TorsionID( n-1, BB, 6),  ZETA_A_FORM);

	pose.set_torsion( TorsionID( n, BB, 1),  ALPHA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 2),   BETA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 3),  GAMMA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 4),  DELTA_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 5),  EPSILON_A_FORM);
	pose.set_torsion( TorsionID( n, BB, 6),  ZETA_A_FORM);

	pose.set_torsion( TorsionID( n+1, BB, 1),  ALPHA_A_FORM);
	pose.set_torsion( TorsionID( n+1, BB, 2),   BETA_A_FORM);
	pose.set_torsion( TorsionID( n+1, BB, 3),  GAMMA_A_FORM);
	pose.set_torsion( TorsionID( n+1, BB, 4),  DELTA_A_FORM);
	pose.set_torsion( TorsionID( n+1, BB, 5),  EPSILON_A_FORM);
	pose.set_torsion( TorsionID( n+1, BB, 6),  ZETA_A_FORM);

	pose.set_torsion( TorsionID( n+2, BB, 1),  ALPHA_A_FORM);

	if (use_phenix_geo_) {
		ideal_coord_.apply(pose, n);
		ideal_coord_.apply(pose, n+1);
	} else {
		pose.set_torsion( TorsionID( n, CHI, 1),  CHI_A_FORM);
		pose.set_torsion( TorsionID( n, CHI, 2),  NU2_A_FORM);
		pose.set_torsion( TorsionID( n, CHI, 3),  NU1_A_FORM);
		pose.set_torsion( TorsionID( n+1, CHI, 1),  CHI_A_FORM);
		pose.set_torsion( TorsionID( n+1, CHI, 2),  NU2_A_FORM);
		pose.set_torsion( TorsionID( n+1, CHI, 3),  NU1_A_FORM);
	}

	// perturb all torsion angles except those at sugar puckers (DELTA, NU1, NU2 )
	if ( random_perturbation_ ) {
		pose.set_torsion( TorsionID( n-1, BB, 5),  EPSILON_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n-1, BB, 6),  ZETA_A_FORM + perturb_amplitude_ * RG.gaussian() );

		pose.set_torsion( TorsionID( n, BB, 1),  ALPHA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n, BB, 2),   BETA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n, BB, 3),  GAMMA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n, BB, 5),  EPSILON_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n, BB, 6),  ZETA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n, CHI, 1),  CHI_A_FORM + perturb_amplitude_ * RG.gaussian() );

		pose.set_torsion( TorsionID( n+1, BB, 1),  ALPHA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n+1, BB, 2),   BETA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n+1, BB, 3),  GAMMA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n+1, BB, 5),  EPSILON_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n+1, BB, 6),  ZETA_A_FORM + perturb_amplitude_ * RG.gaussian() );
		pose.set_torsion( TorsionID( n+1, CHI, 1),  CHI_A_FORM + perturb_amplitude_ * RG.gaussian() );

		pose.set_torsion( TorsionID( n+2, BB, 1),  ALPHA_A_FORM + perturb_amplitude_ * RG.gaussian() );
	}

}

/////////////////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::minimize_base_step( pose::Pose & pose, Size const n ){

	using namespace core::scoring;
	using namespace core::optimization;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	kinematics::MoveMap mm;

	if (minimize_all_){
		mm.set_bb( true );
		mm.set_chi( true );
		mm.set_jump( true );
	} else {
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );
		for (Size i = n-1; i <= n+2; i++ ) {
			mm.set_bb( i, true );
			mm.set_chi( i, true );
		}
		//mm.set_jump( true );
	}

	minimizer.run( pose, mm, *scorefxn, options );
	(*scorefxn)( pose );

}

//////////////////////////////////////////////////////////////////
void
RNA_HelixAssembler::put_constraints_on_base_step( pose::Pose & pose, Size const & n ){
	utility::vector1< std::pair< Size, Size > > pairings;
	pairings.push_back( std::make_pair( n, n+1)  );
	pairings.push_back( std::make_pair( n-1, n+2) );

	scoring::constraints::ConstraintSetOP new_cst_set;
	pose.constraint_set( new_cst_set ); //blank out cst set.

	protocols::rna::setup_base_pair_constraints( pose, pairings );
}


} // namespace rna
} // namespace protocols
