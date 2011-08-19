// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to RNA_FragmentMover
/// @detailed
/// @author Rhiju Das


#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>


#include <core/pose/Pose.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>

// External library headers
// AUTO-REMOVED #include <numeric/random/random.hh>
#include <utility/io/ozstream.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
// AUTO-REMOVED #include <ctime>

//Auto Headers
#include <core/pose/annotated_sequence.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>


using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;
using basic::T;

static basic::Tracer TR( "protocols.rna.rna_fragment_mover" ) ;

namespace protocols {
namespace rna {

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// A bunch of helper functions used in rna apps...
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void
figure_out_reasonable_rna_fold_tree( pose::Pose & pose )
{
	using namespace core::conformation;

	//Look for chainbreaks in PDB.
	Size const nres = pose.total_residue();
	kinematics::FoldTree f( nres );

	Size m( 0 );

	for (Size i=1; i < nres; ++i) {

		if ( scoring::rna::is_rna_chainbreak( pose, i ) ){

			f.new_jump( i, i+1, i );
			m++;

			Residue const & current_rsd( pose.residue( i   ) ) ;
			Residue const &    next_rsd( pose.residue( i+1 ) ) ;
			//			Size dummy( 0 ), jump_atom1( 0 ), jump_atom2( 0 );
			//rna_basepair_jump_atoms( current_rsd.aa(), jump_atom1, dummy, dummy );
			//rna_basepair_jump_atoms( next_rsd.aa(), jump_atom2, dummy, dummy );
			//f.set_jump_atoms( m, current_rsd.atom_name( jump_atom1 ), next_rsd.atom_name( jump_atom2 ) );

			f.set_jump_atoms( m,
												core::scoring::rna::chi1_torsion_atom( current_rsd ),
												core::scoring::rna::chi1_torsion_atom( next_rsd )   );

		}

	}

	pose.fold_tree( f );
}

///////////////////////////////////////////////////////////////////////////////
void
get_base_pairing_info( pose::Pose const & pose,
											 Size const & seqpos,
											 char & secstruct,
											 FArray1D <bool> & edge_is_base_pairing ){

	using namespace core::scoring::rna;
	using namespace core::chemical;
	using namespace core::conformation;

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	edge_is_base_pairing.dimension( 3 );
	edge_is_base_pairing = false;

	bool forms_canonical_base_pair( false ), forms_base_pair( false );

	Size k( 0 ), m( 0 );
	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		if ( i == seqpos ){
			k = base_pair.edge1;
			m = base_pair.edge2;
		} else if ( j == seqpos ){
			k = base_pair.edge2;
			m = base_pair.edge1;
		} else {
			continue;
		}

		edge_is_base_pairing( k ) = true;
		forms_base_pair = true;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) ) {
			std::string atom1, atom2;
			get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
			if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
				forms_canonical_base_pair = true;
			}
		}
	}

	secstruct = 'N';
	if (forms_canonical_base_pair ) {
		secstruct = 'H';
	} else if (forms_base_pair ){
		secstruct = 'P';
	}
}



///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose, utility::io::ozstream & torsions_out)
{

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace protocols::rna;

	Size const total_residue = pose.total_residue();

	figure_out_reasonable_rna_fold_tree( pose );

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	bool const idealize_frag( false );
	protocols::idealize::IdealizeMover idealizer;
	idealizer.fast( false /* option[ fast ] */ );

	//	EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	initialize_non_main_chain_sugar_atoms();

	for (Size i=1; i <= total_residue; ++i) {

		torsions_out << pose.residue( i ).name1() << " " ;

		if (idealize_frag ) {

			////////////////////////////////////////////////
			// NEW: can we idealize this thing?
			pose::Pose mini_pose;
			Size offset = 1;
			if ( i > 1 ) {
				mini_pose.append_residue_by_bond( pose.residue( i-1 ) );
				offset = 2;
			}
			mini_pose.append_residue_by_bond( pose.residue( i ) );
			if ( i < total_residue ){
				mini_pose.append_residue_by_bond( pose.residue( i+1 ) );
			}

			figure_out_reasonable_rna_fold_tree( mini_pose );
			idealizer.apply( mini_pose );
			idealizer.apply( mini_pose );

			for (Size j=1; j <= core::scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j) {
				id::TorsionID my_ID( offset, id::BB, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) );
			}

			for (Size j=1; j <= core::scoring::rna::NUM_RNA_CHI_TORSIONS; ++j) {
				id::TorsionID my_ID( offset, id::CHI, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) ) << " ";
			}

		} else {
			for (Size j=1; j <= core::scoring::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j) {
				id::TorsionID my_ID( i, id::BB, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) );
			}

			for (Size j=1; j <= core::scoring::rna::NUM_RNA_CHI_TORSIONS; ++j) {
				id::TorsionID my_ID( i, id::CHI, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) ) << " ";
			}

			//New (Feb., 2009) ...
			// x-y-z of coordinates of C2*, C1*, and O4*, in a local coordiante system defined
			// by C3*, C4*, and C5* (as "stub" atoms).
			conformation::Residue rsd = pose.residue( i );
			kinematics::Stub const input_stub( rsd.xyz( " C3*" ), rsd.xyz( " C3*" ), rsd.xyz( " C4*" ), rsd.xyz( " C5*" ) );

			torsions_out << " S  " ;
			for (Size n = 1; n <= non_main_chain_sugar_atoms.size(); n++  ) {
				Vector v = input_stub.global2local( rsd.xyz( non_main_chain_sugar_atoms[ n ] ) );
				torsions_out << F( 12, 6, v.x() ) << " " ;
				torsions_out << F( 12, 6, v.y() ) << " " ;
				torsions_out << F( 12, 6, v.z() ) << "  " ;
			}

		}

		// "Secondary structure" ==>
		//    H = forming canonical Watson/Crick base pair
		//    P = forming some kind of base pair
		//    N = no base pairs.
		FArray1D < bool > edge_is_base_pairing( 3, false );
		char secstruct( 'X' );
		get_base_pairing_info( pose, i, secstruct, edge_is_base_pairing );
		torsions_out << secstruct << " " <<
			edge_is_base_pairing( 1 ) << " " <<
			edge_is_base_pairing( 2 ) << " " <<
			edge_is_base_pairing( 3 ) << " ";

		torsions_out << pose.fold_tree().is_cutpoint( i ) << I(6, i)  << std::endl;
	}


}

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose, std::string const outfile )
{
	utility::io::ozstream torsions_out ( outfile );
	create_rna_vall_torsions( pose, torsions_out );

}



//////////////////////////////////////////////////////////////////////////////////////
Real
get_o1p_o2p_sign( pose::Pose const & pose ) {
	Real sign( 0 );
	for (Size i = 2; i <= pose.total_residue(); i++ ) {

		conformation::Residue const & rsd( pose.residue(i)  );
		if (!rsd.is_RNA() ) continue;

		sign = dot( rsd.xyz( " O5*" ) - rsd.xyz( " P  " ),
				 cross( rsd.xyz( " O2P" ) - rsd.xyz( " P  " ),
								rsd.xyz( " O1P" ) - rsd.xyz( " P  " ) ) );
		break;
	}
	return sign;
}

////////////////////////////////////////////////////////////////
void
ensure_phosphate_nomenclature_matches_mini( pose::Pose & pose )
{
	runtime_assert( pose.total_residue() > 1 );
	Real sign1 = get_o1p_o2p_sign( pose );

	pose::Pose mini_pose;
	core::pose::make_pose_from_sequence( mini_pose, "aa", pose.residue(1).residue_type_set() );
	Real sign2 = get_o1p_o2p_sign( mini_pose );

	if ( sign1 * sign2 > 0 ) return;

	std::cout << "*************************************************************" << std::endl;
	std::cout << " Warning ... flipping O1P <--> O2P to match mini convention  " << std::endl;
	std::cout << "*************************************************************" << std::endl;

	for (Size i = 1; i <= pose.total_residue(); i++ ) {

		conformation::Residue const & rsd( pose.residue(i) );
		if (!rsd.is_RNA() ) continue;

		Vector const temp1 = rsd.xyz( " O1P" );
		Vector const temp2 = rsd.xyz( " O2P" );
		pose.set_xyz( id::AtomID( rsd.atom_index( " O1P" ), i ), temp2 );
		pose.set_xyz( id::AtomID( rsd.atom_index( " O2P" ), i ), temp1 );
	}

}

//////////////////////////////////////////////////////////////////////////////////////
Real
get_o1p_o2p_sign_parin( pose::Pose const & pose , Size res_num) {
	Real sign( 0 );

		conformation::Residue const & rsd( pose.residue(res_num)  );

//		std::cout << "O3*= " << rsd.xyz( " O3*" );
//		std::cout << "  O2P= " << rsd.xyz( " O2P" );
//		std::cout << "  O1P= " << rsd.xyz( " O1P" );
//		std::cout << "  P= " << rsd.xyz( " P  " ) << std::endl;

		sign = dot( rsd.xyz( " O5*" ) - rsd.xyz( " P  " ),
				 cross( rsd.xyz( " O2P" ) - rsd.xyz( " P  " ),
								rsd.xyz( " O1P" ) - rsd.xyz( " P  " ) ) );


	return sign;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
void
assert_phosphate_nomenclature_matches_mini( pose::Pose const & pose){

	runtime_assert( pose.total_residue() > 1 );


	for(Size res_num=1; res_num<=pose.total_residue(); res_num++){


		Real sign1 = get_o1p_o2p_sign_parin( pose,  res_num);

		pose::Pose mini_pose;
		core::pose::make_pose_from_sequence( mini_pose, "aa", pose.residue(res_num).residue_type_set() );
		Real sign2 = get_o1p_o2p_sign( mini_pose);

		if ( sign1 * sign2 < 0 ) {

			std::cout << " In the assert_phosphate_nomenclature_matches_mini function, phosphate_nomenclature_matches does not match mini! " << std::endl;
			exit(1);

			conformation::Residue const & rsd( pose.residue(res_num) );
			if (!rsd.is_RNA() ) {
				std::cout << " residue is not a RNA" << std::endl;
				exit (1);
			};

		}
	}
}


void
ensure_phosphate_nomenclature_matches_mini_parin( pose::Pose & pose)
{


	//std::cout << "Enter ensure_phosphate_nomenclature_matches_mini_parin function" << std::endl;

	runtime_assert( pose.total_residue() > 1 );


	for(Size res_num=1; res_num<=pose.total_residue(); res_num++){


		Real sign1 = get_o1p_o2p_sign_parin( pose,  res_num);

		pose::Pose mini_pose;
		core::pose::make_pose_from_sequence( mini_pose, "aa", pose.residue(res_num).residue_type_set() );
		Real sign2 = get_o1p_o2p_sign( mini_pose);


		if ( sign1 * sign2 < 0 ) {

			std::cout << "*************************************************************" << std::endl;
			std::cout << " Warning ... flipping O1P <--> O2P to match mini convention  " << std::endl;
			std::cout << "*************************************************************" << std::endl;
			std::cout << "res_num " << res_num << std::endl;
			std::cout << "sign1: " << sign1 << std::endl;
			std::cout << "sign2: " << sign2 << std::endl;


			conformation::Residue const & rsd( pose.residue(res_num) );
			if (!rsd.is_RNA() ) {
				std::cout << " residue is not a RNA" << std::endl;
				exit (1);
			};

			Vector const temp1 = rsd.xyz( " O1P" );
			Vector const temp2 = rsd.xyz( " O2P" );
			pose.set_xyz( id::AtomID( rsd.atom_index( " O1P" ), res_num ), temp2 );
			pose.set_xyz( id::AtomID( rsd.atom_index( " O2P" ), res_num ), temp1 );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
void
export_packer_results( 	utility::vector1< std::pair< Real, std::string > > & results,
												utility::vector1< pose::PoseOP > pose_list,
												scoring::ScoreFunctionOP & scorefxn,
												std::string const & outfile,
												bool const dump )
{

	utility::io::ozstream out( outfile );
	for (Size n = 1; n <= results.size() ; n++ ){
		out << F(8,3,results[n].first) << " " << results[n].second << std::endl;
	}
	out.close();

	using namespace core::io::silent;
	core::io::silent::SilentFileData silent_file_data;

	std::string silent_file( outfile );
	Size pos( silent_file.find( ".txt" ) );
	silent_file.replace( pos, 4, ".out" );

	for (Size n = 1; n <= results.size() ; n++ ){
		pose::Pose & pose( *pose_list[n] );
		(*scorefxn)( pose );
		std::string const tag( "S_"+lead_zero_string_of( n, 4 ) );
		RNA_SilentStruct s( pose, tag );
		if ( dump ) pose.dump_pdb( tag+".pdb");
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
check_base_pair( pose::Pose & pose, FArray1D_int & struct_type )
{

	using namespace core::scoring::rna;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
 	(*scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	Size const nres( pose.total_residue() );
	FArray1D_bool forms_noncanonical_base_pair( nres, false );
	FArray1D_bool forms_canonical_base_pair( nres, false );
	FArray1D_int  WC_base_pair_partner( nres, 0 );

	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		bool WC_base_pair( false );
		if ( ( base_pair.edge1 == WATSON_CRICK && base_pair.edge2 == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) ) 			{
			std::string atom1, atom2;
			//		get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
			//			if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
			WC_base_pair = true;
				//			}
		}

		std::cout << rsd_i.name1() << I(3,i) << "--" << rsd_j.name1() << I(3,j) << "  WC:" << WC_base_pair << "  score: " << F(10,6,it->first) << std::endl;

		if ( WC_base_pair ) {
			forms_canonical_base_pair( i ) = true;
			forms_canonical_base_pair( j ) = true;
			WC_base_pair_partner( i ) = j;
			WC_base_pair_partner( j ) = i;
		} else {
			forms_noncanonical_base_pair( i ) = true;
			forms_noncanonical_base_pair( j ) = true;
		}

	}


	struct_type.dimension( nres );
	for (Size i = 1; i <= nres; i++ ) {
		if ( !forms_noncanonical_base_pair(i) && !forms_canonical_base_pair(i) ) {
			struct_type( i ) = 0;
		} else if (forms_canonical_base_pair(i) && !forms_noncanonical_base_pair( WC_base_pair_partner(i) )  ) {
			struct_type( i ) = 1;
		} else {
			struct_type( i ) = 2;
		}
	}

}


/////////////////////////////////////////////////////////
void
setup_base_pair_constraints(
														pose::Pose & pose,
														utility::vector1< std::pair< Size, Size > > const &  pairings )
{

	using namespace core::scoring::constraints;
	using namespace core::scoring::rna;

 	Real const WC_distance( 1.9 );
 	Real const distance_stddev( 0.25 ); //Hmm. Maybe try linear instead?
	FuncOP const distance_func( new HarmonicFunc( WC_distance, distance_stddev ) );

	// Need to force base pairing -- not base stacking!
 	Real const C1star_distance( 10.5 );
 	Real const C1star_distance_stddev( 1.0 ); //Hmm. Maybe try linear instead?
	FuncOP const C1star_distance_func( new HarmonicFunc( C1star_distance, C1star_distance_stddev ) );

	for ( Size n = 1; n <= pairings.size(); n++ ) {

		Size const & i = pairings[n].first;
		Size const & j = pairings[n].second;

		Size const atom1 = pose.residue(i).type().atom_index( " C1*" ) ;
		Size const atom2 = pose.residue(j).type().atom_index( " C1*" ) ;
		pose.add_constraint( new AtomPairConstraint(
																								id::AtomID(atom1,i),
																								id::AtomID(atom2,j),
																								C1star_distance_func ) );

		utility::vector1< std::string > atom_ids1, atom_ids2;
		get_watson_crick_base_pair_atoms( pose.residue(i).aa(), pose.residue(j).aa(), atom_ids1, atom_ids2 );

		for (Size p = 1; p <= atom_ids1.size(); p++ ){

			Size const atom1 = pose.residue(i).type().atom_index( atom_ids1[p] ) ;
			Size const atom2 = pose.residue(j).type().atom_index( atom_ids2[p] ) ;

			TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue(i).name1() << I(3,i) << " <-->  " <<
				pose.residue(j).name1() << I(3,j) << "   " <<
				atom_ids1[p] << " <--> " <<
				atom_ids2[p] << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;

			pose.add_constraint( new AtomPairConstraint(
																									id::AtomID(atom1,i),
																									id::AtomID(atom2,j),
																									distance_func ) );

		}
	}


}



} // namespace rna
} // namespace protocols
