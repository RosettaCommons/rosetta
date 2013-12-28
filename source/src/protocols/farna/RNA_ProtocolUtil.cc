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


#include <protocols/farna/RNA_ProtocolUtil.hh>
#include <protocols/farna/RNA_SecStructInfo.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/forge/methods/fold_tree_functions.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>


#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>



// External library headers
#include <utility/io/ozstream.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

#include <protocols/farna/RNA_MatchType.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>


using namespace core;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
using basic::T;

static basic::Tracer TR( "protocols.rna.RNA_ProtocolUtil" ) ;

namespace protocols {
namespace farna {

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

		if ( !pose.residue(i).is_RNA() && !pose.residue(i+1).is_RNA() )  continue;

		if ( (  pose.residue(i).is_RNA() && !pose.residue(i+1).is_RNA() ) ||
				 ( !pose.residue(i).is_RNA() &&  pose.residue(i+1).is_RNA() ) ) {
			f.new_jump( i, i+1, i );
			m++;
			continue;
		}

		if ( pose::rna::is_rna_chainbreak( pose, i ) ){

			//std::cout << "CHAINBREAK between " << i << " and " << i+1 << std::endl;

			f.new_jump( i, i+1, i );
			m++;

			Residue const & current_rsd( pose.residue( i   ) ) ;
			Residue const &    next_rsd( pose.residue( i+1 ) ) ;
			//			Size dummy( 0 ), jump_atom1( 0 ), jump_atom2( 0 );
			//rna_basepair_jump_atoms( current_rsd.aa(), jump_atom1, dummy, dummy );
			//rna_basepair_jump_atoms( next_rsd.aa(), jump_atom2, dummy, dummy );
			//f.set_jump_atoms( m, current_rsd.atom_name( jump_atom1 ), next_rsd.atom_name( jump_atom2 ) );

			f.set_jump_atoms( m,
												chemical::rna::chi1_torsion_atom( current_rsd ),
												chemical::rna::chi1_torsion_atom( next_rsd )   );

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

			if ( !rsd_i.is_coarse() ) { // doesn't work for coarse-grained RNA
				get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
				if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() > 3.5 ) continue;
			}

			forms_canonical_base_pair = true;

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
get_base_pairing_list( pose::Pose & pose,
											 utility::vector1< std::pair<Size, Size> > & base_pairing_list ){

	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::rna;

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list const & scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	//	bool forms_canonical_base_pair( false );

	Size k( 0 ), m( 0 );

	std::list < std::pair< Size,Size > > base_pair_list0;

	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		if ( i > j ) continue;

		k = base_pair.edge1;
		m = base_pair.edge2;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) ) {
			std::string atom1, atom2;

			if ( !rsd_i.is_coarse() ) { // doesn't work for coarse-grained RNA
				get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
				if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() > 3.5 ) continue;
			}

			base_pair_list0.push_back( std::make_pair( i, j ) );

		}
	}


	base_pair_list0.sort();
	base_pairing_list.clear();
	for ( std::list< std::pair<Size,Size > >::const_iterator it = base_pair_list0.begin();
				it != base_pair_list0.end(); ++it ){
		base_pairing_list.push_back( *it );
	}

}


///////////////////////////////////////////////////////////////////////////////
void
figure_out_secstruct( pose::Pose & pose ){
	using namespace core::scoring;
	using namespace ObjexxFCL;

	// Need some stuff to figure out which residues are base paired. First score.
	ScoreFunctionOP scorefxn( new ScoreFunction );
	scorefxn->set_weight( rna_base_pair, 1.0 );
	(*scorefxn)( pose );

	std::string secstruct = "";
	FArray1D < bool > edge_is_base_pairing( 3, false );
	char secstruct1( 'X' );
	for (Size i=1; i <= pose.total_residue() ; ++i) {
		get_base_pairing_info( pose, i, secstruct1, edge_is_base_pairing );
		secstruct += secstruct1;
	}

	std::cout << "SECSTRUCT: " << secstruct << std::endl;

	set_rna_secstruct( pose, secstruct );

}


///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose,
												  utility::io::ozstream & torsions_out,
													utility::vector1 <Size> const & exclude_res_list )
{

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace protocols::farna;

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

	for (Size i=1; i <= total_residue; ++i) {
		if ( is_num_in_list(i, exclude_res_list) ) continue;

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

			for (Size j=1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j) {
				id::TorsionID my_ID( offset, id::BB, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) );
			}

			for (Size j=1; j <= chemical::rna::NUM_RNA_CHI_TORSIONS; ++j) {
				id::TorsionID my_ID( offset, id::CHI, j );
				torsions_out << F( 12, 6, mini_pose.torsion( my_ID ) ) << " ";
			}

		} else {
			for (Size j=1; j <= chemical::rna::NUM_RNA_MAINCHAIN_TORSIONS; ++j) {
				id::TorsionID my_ID( i, id::BB, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) );
			}

			for (Size j=1; j <= chemical::rna::NUM_RNA_CHI_TORSIONS; ++j) {
				id::TorsionID my_ID( i, id::CHI, j );
				torsions_out << F( 12, 6, pose.torsion( my_ID ) ) << " ";
			}

			//New (Feb., 2009) ...
			// x-y-z of coordinates of C2', C1', and O4', in a local coordiante system defined
			// by C3', C4', and C5' (as "stub" atoms).
			conformation::Residue rsd = pose.residue( i );
			kinematics::Stub const input_stub( rsd.xyz( " C3'" ), rsd.xyz( " C3'" ), rsd.xyz( " C4'" ), rsd.xyz( " C5'" ) );

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

		bool is_cutpoint = false;
		if ( pose.fold_tree().is_cutpoint( i ) ||
				 is_num_in_list(i + 1, exclude_res_list) ) {
			is_cutpoint = true;
		}

		torsions_out << is_cutpoint << I(6, i)  << std::endl;
	}


}

///////////////////////////////////////////////////////////////////////////////
void
create_rna_vall_torsions( pose::Pose & pose,
													std::string const outfile,
													utility::vector1 <Size> const & exclude_res_list )
{
	utility::io::ozstream torsions_out ( outfile );
	create_rna_vall_torsions( pose, torsions_out, exclude_res_list );

}



//////////////////////////////////////////////////////////////////////////////////////
Real
get_o1p_o2p_sign( pose::Pose const & pose ) {

	Real sign= 0;
	bool found_valid_sign=false;

	for (Size i = 2; i <= pose.total_residue(); i++ ) {

		conformation::Residue const & rsd( pose.residue(i)  );
		if (!rsd.is_RNA() ) continue;

		sign = dot( rsd.xyz( " O5'" ) - rsd.xyz( " P  " ), cross( rsd.xyz( " OP1" ) - rsd.xyz( " P  " ), rsd.xyz( " OP2" ) - rsd.xyz( " P  " ) ) );

		found_valid_sign=true;

		break;
	}

	if(found_valid_sign==false) utility_exit_with_message("found_valid_sign==false");

	return sign;
}

//////////////////////////////////////////////////////////////////////////////////////
//This version used to be called get_o1p_o2p_sign_parin()
Real
get_o1p_o2p_sign( pose::Pose const & pose , Size res_num) {

	if(res_num > pose.total_residue()) utility_exit_with_message("res_num > pose.total_residue()");

	conformation::Residue const & rsd( pose.residue(res_num)  );

	//SML PHENIX conference cleanup
	if (basic::options::option[basic::options::OptionKeys::rna::rna_prot_erraser].value()){
		if ( !rsd.is_RNA() ) return 0.0;
	} else {
		if(rsd.is_RNA()==false) utility_exit_with_message("rsd.is_RNA()==false!");
	}

	Real const sign = dot( rsd.xyz( " O5'" ) - rsd.xyz( " P  " ), cross( rsd.xyz( " OP1" ) - rsd.xyz( " P  " ), rsd.xyz( " OP2" ) - rsd.xyz( " P  " ) ) );

	return sign;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
void
assert_phosphate_nomenclature_matches_mini( pose::Pose const & pose){

	runtime_assert( pose.total_residue() > 1 );

	for(Size res_num=1; res_num<=pose.total_residue(); res_num++){

		Real sign1 = get_o1p_o2p_sign( pose,  res_num);

		pose::Pose mini_pose; //Could move this part outside the for loop
		make_pose_from_sequence( mini_pose, "aa", pose.residue(res_num).residue_type_set() );
		Real const sign2 = get_o1p_o2p_sign( mini_pose);

		if ( sign1 * sign2 < 0 ) {

			std::cout << "In the assert_phosphate_nomenclature_matches_mini function: phosphate_nomenclature_matches does not match mini! " << std::endl;
			utility_exit_with_message("In the assert_phosphate_nomenclature_matches_mini function: phosphate_nomenclature_matches does not match mini!");

			conformation::Residue const & rsd( pose.residue(res_num) );

			if(rsd.is_RNA()==false){ //Consistency check!
				std::cout << "residue # " << res_num << " should be a RNA nucleotide" << std::endl;
				utility_exit_with_message("residue # " + string_of(res_num)+ " should be a RNA nucleotide!");
			};


		}
	}
}


////////////////////////////////////////////////////////////////
//Jan 08, 2012 Parin S:
//To Rhiju, I think this version is buggy and should be deprecated (use make_phosphate_nomenclature_matches_mini instead!)
void
ensure_phosphate_nomenclature_matches_mini( pose::Pose & pose )
{
	runtime_assert( pose.total_residue() > 1 );
	Real sign1 = get_o1p_o2p_sign( pose );

	pose::Pose mini_pose;
	make_pose_from_sequence( mini_pose, "aa", pose.residue(1).residue_type_set() );
	Real sign2 = get_o1p_o2p_sign( mini_pose );

	if ( sign1 * sign2 > 0 ) return;

	std::cout << "*************************************************************" << std::endl;
	std::cout << " Warning ... flipping OP2 <--> OP1 to match mini convention  " << std::endl;
	std::cout << "*************************************************************" << std::endl;

	for (Size i = 1; i <= pose.total_residue(); i++ ) {

		conformation::Residue const & rsd( pose.residue(i) );
		if (!rsd.is_RNA() ) continue;

		if (!rsd.type().has( " OP2")) continue;
		if (!rsd.type().has( " OP1")) continue;

		Vector const temp1 = rsd.xyz( " OP2" );
		Vector const temp2 = rsd.xyz( " OP1" );
		pose.set_xyz( id::AtomID( rsd.atom_index( " OP2" ), i ), temp2 );
		pose.set_xyz( id::AtomID( rsd.atom_index( " OP1" ), i ), temp1 );
	}

}

////////////////////////////////////////////////////////////////
void
make_phosphate_nomenclature_matches_mini( pose::Pose & pose)
{


	for(Size res_num=1; res_num<=pose.total_residue(); res_num++){

		if ( !pose.residue( res_num ).is_RNA() ) continue;

		pose::Pose mini_pose; //Could move this part outside of the for loop
		make_pose_from_sequence( mini_pose, "aa", pose.residue( res_num ).residue_type_set());
		Real const sign2 = get_o1p_o2p_sign( mini_pose);

		Real sign1 = get_o1p_o2p_sign( pose,  res_num);

		if ( sign1 * sign2 < 0 ) {

			//std::cout << " Flipping OP2 <--> OP1 " << "res_num " << res_num << " | sign1: " << sign1 << " | sign2: " << sign2 << std::endl;

			conformation::Residue const & rsd( pose.residue(res_num) );

			if(rsd.is_RNA()==false){ //Consistency check!
				std::cout << "residue # " << res_num << " should be a RNA nucleotide!" << std::endl;
				utility_exit_with_message("residue # " + string_of(res_num)+ " should be a RNA nucleotide!");
			};

			Vector const temp1 = rsd.xyz( " OP2" );
			Vector const temp2 = rsd.xyz( " OP1" );
			pose.set_xyz( id::AtomID( rsd.atom_index( " OP2" ), res_num ), temp2 );
			pose.set_xyz( id::AtomID( rsd.atom_index( " OP1" ), res_num ), temp1 );
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
														utility::vector1< std::pair< Size, Size > > const &  pairings,
														Real const suppress_factor /* = 1.0 */ )
{

	using namespace core::scoring::constraints;
	using namespace core::scoring::rna;

 	Real const WC_distance( 1.9 );
 	Real const distance_stddev( 0.25 / suppress_factor ); //Hmm. Maybe try linear instead?
	core::scoring::func::FuncOP const distance_func( new core::scoring::func::HarmonicFunc( WC_distance, distance_stddev ) );

	// Need to force base pairing -- not base stacking!
 	Real const C1prime_distance( 10.5 );
 	Real const C1prime_distance_stddev( 1.0 / suppress_factor ); //Hmm. Maybe try linear instead?
	core::scoring::func::FuncOP const C1prime_distance_func( new core::scoring::func::HarmonicFunc( C1prime_distance, C1prime_distance_stddev ) );

	for ( Size n = 1; n <= pairings.size(); n++ ) {

		Size const & i = pairings[n].first;
		Size const & j = pairings[n].second;

		if ( !pose.residue(i).is_RNA() ) continue;
		if ( !pose.residue(j).is_RNA() ) continue;


		if ( !pose.residue(i).is_coarse() ) { //fullatom
			Size const atom1 = pose.residue(i).type().atom_index( " C1'" ) ;
			Size const atom2 = pose.residue(j).type().atom_index( " C1'" ) ;
			pose.add_constraint( new AtomPairConstraint(
																									id::AtomID(atom1,i),
																									id::AtomID(atom2,j),
																									C1prime_distance_func ) );

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

		} else { //coarse-grained RNA

			static Real const coarse_WC_SUG_distance_min( 12.3 );
			static Real const coarse_WC_SUG_distance_max( 14.3 );
			static Real const coarse_WC_SUG_distance_stddev( 2.0 );
			static Real const coarse_WC_SUG_bonus( -2.5 );
			static core::scoring::func::FuncOP const coarse_SUG_distance_func( new core::scoring::func::FadeFunc( coarse_WC_SUG_distance_min - coarse_WC_SUG_distance_stddev,
																																  coarse_WC_SUG_distance_max + coarse_WC_SUG_distance_stddev,
																																	coarse_WC_SUG_distance_stddev,
																																	coarse_WC_SUG_bonus ) );

			if ( true ){
				Size const & atom1 = pose.residue(i).atom_index( " S  " );
				Size const & atom2 = pose.residue(j).atom_index( " S  " );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue(i).name1() << I(3,i) << " <-->  " <<
					pose.residue(j).name1() << I(3,j) << "   " <<
					" S  " << " <--> " <<
					" S  " << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( new AtomPairConstraint(
																										id::AtomID(atom1,i),
																										id::AtomID(atom2,j),
																										coarse_SUG_distance_func ) );
			}

			static Real const coarse_WC_CEN_distance( 5.5 );
			static Real const coarse_WC_CEN_distance_stddev( 3.0 );
			static Real const coarse_WC_CEN_bonus( -5.0 );
			static core::scoring::func::FuncOP const coarse_CEN_distance_func( new core::scoring::func::FadeFunc( coarse_WC_CEN_distance - coarse_WC_CEN_distance_stddev,
																																	coarse_WC_CEN_distance + coarse_WC_CEN_distance_stddev,
																																	coarse_WC_CEN_distance_stddev,
																																	coarse_WC_CEN_bonus ) );

			if ( true ){
				Size const & atom1 = pose.residue(i).atom_index( " CEN" );
				Size const & atom2 = pose.residue(j).atom_index( " CEN" );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue(i).name1() << I(3,i) << " <-->  " <<
					pose.residue(j).name1() << I(3,j) << "   " <<
					" CEN" << " <--> " <<
					" CEN" << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( new AtomPairConstraint(
																										id::AtomID(atom1,i),
																										id::AtomID(atom2,j),
																										coarse_CEN_distance_func ) );
			}

			static Real const coarse_WC_X_distance( 3.5 );
			static Real const coarse_WC_X_distance_stddev( 2.0 );
			static Real const coarse_WC_X_bonus( -5.0 );
			static core::scoring::func::FuncOP const coarse_X_distance_func( new core::scoring::func::FadeFunc( coarse_WC_X_distance - coarse_WC_X_distance_stddev,
																																coarse_WC_X_distance + coarse_WC_X_distance_stddev,
																																coarse_WC_X_distance_stddev, coarse_WC_X_bonus ) );

			{
				Size const & atom1 = pose.residue(i).atom_index( " X  " );
				Size const & atom2 = pose.residue(j).atom_index( " X  " );
				TR << "BASEPAIR: Adding rna_force_atom_pair constraint: " << pose.residue(i).name1() << I(3,i) << " <-->  " <<
					pose.residue(j).name1() << I(3,j) << "   " <<
					" X  " << " <--> " <<
					" X  " << ".  [ " << atom1 << "-" << atom2 << "]" << std::endl;
				pose.add_constraint( new AtomPairConstraint(
																										id::AtomID(atom1,i),
																										id::AtomID(atom2,j),
																										coarse_X_distance_func ) );
			}
		}
	}
}


/////////////////////////////////////////////////////////
// Following only works for coarse-grained RNA poses --
//  perhaps should set up something similar for regular RNA
//  as an alternative to linear_chainbreak.
void
setup_coarse_chainbreak_constraints( pose::Pose & pose, Size const & n )
{

	using namespace core::scoring::constraints;
	using namespace core::scoring::rna;

	if ( !pose.residue_type( n ).has( " S  " ) ) return;
	if ( !pose.residue_type( n ).is_RNA() ) return;

	std::cout << "Adding CHAINBREAK constraints to " << n << " " << n+1 << std::endl;

 	static Real const S_P_distance_min( 3 );
 	static Real const S_P_distance_max( 4 );
 	static Real const S_P_distance_fade( 1.0 );
 	static Real const S_P_distance_bonus( -10.0 );
	static core::scoring::func::FuncOP const S_P_distance_func( new core::scoring::func::FadeFunc( S_P_distance_min - S_P_distance_fade,
																								S_P_distance_max + S_P_distance_fade,
																								S_P_distance_fade,
																								S_P_distance_bonus ) );

	// loose long-range constraint too
 	static Real const S_P_distance( 4.0 );
 	static Real const S_P_distance_stdev( 10.0 );
	static core::scoring::func::FuncOP const S_P_harmonic_func( new core::scoring::func::HarmonicFunc( S_P_distance,
																													 S_P_distance_stdev ) );


 	static Real const S_S_distance_min( 4.5 );
 	static Real const S_S_distance_max( 7.5 );
 	static Real const S_S_distance_fade( 2.0 );
 	static Real const S_S_distance_bonus( -5.0 );
	static core::scoring::func::FuncOP const S_S_distance_func( new core::scoring::func::FadeFunc( S_S_distance_min - S_S_distance_fade,
																								S_S_distance_max + S_S_distance_fade,
																								S_S_distance_fade,
																								S_S_distance_bonus ) );

 	static Real const P_P_distance_min( 4.5 );
 	static Real const P_P_distance_max( 7.5 );
 	static Real const P_P_distance_fade( 2.0 );
 	static Real const P_P_distance_bonus( -5.0 );
	static core::scoring::func::FuncOP const P_P_distance_func( new core::scoring::func::FadeFunc( P_P_distance_min - P_P_distance_fade,
																								P_P_distance_max + P_P_distance_fade,
																								P_P_distance_fade,
																								P_P_distance_bonus ) );

	Size const & atom_S1 = pose.residue( n ).atom_index( " S  " );
	Size const & atom_P1 = pose.residue( n ).atom_index( " P  " );
	Size const & atom_S2 = pose.residue( n+1 ).atom_index( " S  " );
	Size const & atom_P2 = pose.residue( n+1 ).atom_index( " P  " );

	pose.add_constraint( new AtomPairConstraint(
																							id::AtomID(atom_S1, n),
																							id::AtomID(atom_P2, n+1),
																									S_P_distance_func ) );

	pose.add_constraint( new AtomPairConstraint(
																							id::AtomID(atom_S1, n),
																							id::AtomID(atom_P2, n+1),
																									S_P_harmonic_func ) );

	pose.add_constraint( new AtomPairConstraint(
																							id::AtomID(atom_P1, n),
																							id::AtomID(atom_P2, n+1),
																							P_P_distance_func ) );

	pose.add_constraint( new AtomPairConstraint(
																							id::AtomID(atom_S1, n),
																							id::AtomID(atom_S2, n+1),
																							S_S_distance_func ) );

}

/////////////////////////////////////////////////////////////////////////////////////////////////
std::string const
convert_based_on_match_type( std::string const RNA_string, Size const type ){

		std::string RNA_string_local = RNA_string;

		Size const size = RNA_string.length();

		//Obey orders to match exactly, match pyrimidine/purine, or match all.
		if (type == MATCH_ALL){

			for (Size i = 0; i < size; i++) 	RNA_string_local[ i ] = 'n';

		} else if ( type == MATCH_YR ) {

			for (Size i = 0; i < size; i++) {
				if (RNA_string[ i ] == 'g' || RNA_string[ i ] == 'a' ){
					RNA_string_local[ i ] = 'r';
				} else {
					runtime_assert( RNA_string[ i ] == 'u' || RNA_string[ i ] == 'c' );
					RNA_string_local[ i ] = 'y';
				}
			}

		}

		return RNA_string_local;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	compare_RNA_char( char const char1, char const char2 ) {
		//Man this is silly, there must be a more elegant way to do this.
		if (char1 == char2) return true;
		if (char1 == 'n' || char2 == 'n') return true;
		if (char1 == 'r' && (char2 == 'a' || char2 == 'g')) return true;
		if (char1 == 'y' && (char2 == 'c' || char2 == 'u')) return true;
		if (char2 == 'r' && (char1 == 'a' || char1 == 'g')) return true;
		if (char2 == 'y' && (char1 == 'c' || char1 == 'u')) return true;
		return false;
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	compare_RNA_secstruct( char const char1, char const char2 ) {
		if (char1 == char2) return true;
		if (char1 == 'X' || char2 == 'X' ) return true;
		if (char1 == 'L' && ( char2 == 'N' || char2 == 'P') ) return true;
		if (char2 == 'L' && ( char1 == 'N' || char1 == 'P') ) return true;
		return false;
	}

///////////////////////////////////////////////////////////////////////////
Vector
get_sugar_centroid( core::conformation::Residue const & rsd ){
	Vector cen( 0.0 );
	cen += rsd.xyz(  rsd.atom_index( " C1'" ) );
	cen += rsd.xyz(  rsd.atom_index( " C2'" ) );
	cen += rsd.xyz(  rsd.atom_index( " C3'" ) );
	cen += rsd.xyz(  rsd.atom_index( " C4'" ) );
	cen += rsd.xyz(  rsd.atom_index( " O4'" ) );
	cen /= 5.0;
	return cen;
}

/////////////////////////////////////////////////
void
make_extended_coarse_pose( pose::Pose & coarse_pose, std::string const & full_sequence ){

 	using namespace core::chemical;
 	using namespace core::id;

	ResidueTypeSetCAP rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );

	make_pose_from_sequence( coarse_pose, full_sequence, *rsd_set_coarse );

	for ( Size n = 1; n <= coarse_pose.total_residue(); n++ ) {
		coarse_pose.set_torsion( TorsionID( n, BB, 1 ), -150.0 );
		coarse_pose.set_torsion( TorsionID( n, BB, 2 ),  180.0 );
	}

}


///////////////////////////////////////////////////////////////////////////
void
make_coarse_pose( pose::Pose const & pose, pose::Pose & coarse_pose ){

 	using namespace core::chemical;
 	using namespace core::id;
 	using namespace core::scoring::rna;

	ResidueTypeSetCAP rsd_set, rsd_set_coarse;
	RNA_CentroidInfo rna_centroid_info;

	rsd_set_coarse = ChemicalManager::get_instance()->residue_type_set( COARSE_RNA );
	make_extended_coarse_pose( coarse_pose, pose.sequence() );

	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		coarse_pose.set_xyz(  NamedAtomID( " P  ", n ),  pose.xyz( NamedAtomID( " P  ", n )) );

		//coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  pose.xyz( NamedAtomID( " C4'", n )) );
		Vector sugar_centroid = get_sugar_centroid( pose.residue( n ) );
		coarse_pose.set_xyz(  NamedAtomID( " S  ", n ),  sugar_centroid );

		Vector base_centroid = rna_centroid_info.get_base_centroid( pose.residue( n ) );
		kinematics::Stub stub = rna_centroid_info.get_base_coordinate_system( pose.residue( n ), base_centroid );

		coarse_pose.set_xyz(  NamedAtomID( " CEN", n ),  base_centroid );
		coarse_pose.set_xyz(  NamedAtomID( " X  ", n ),  base_centroid + stub.M.col_x() );
		coarse_pose.set_xyz(  NamedAtomID( " Y  ", n ),  base_centroid + stub.M.col_y() );

	}

	core::kinematics::FoldTree f( pose.fold_tree() );
	for ( Size n = 1; n <= pose.num_jump(); n++ ) {
		f.set_jump_atoms( n, " Y  ", " Y  " );
	}
	coarse_pose.fold_tree( f );

}

////////////////////////////////////////////////////////
void
remove_cutpoint_closed( pose::Pose & pose, Size const i ){

	using namespace core::chemical;
	using namespace core::kinematics;

	remove_variant_type_from_pose_residue( pose, CUTPOINT_LOWER, i );
	remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i+1 );

	//	using namespace protocols::forge::methods;
	// 	FoldTree f( pose.fold_tree() );
	// 	remove_cutpoint( i, f );
	// 	pose.fold_tree( f );

	utility::vector1< int > const & cutpoints = pose.fold_tree().cutpoints();

	Size const num_jump = 	pose.fold_tree().num_jump();
	utility::vector1< Size > upstream_pos, downstream_pos;
	for ( Size n = 1; n <= num_jump; n++ ) {
		upstream_pos.push_back( pose.fold_tree().upstream_jump_residue( n ) );
		downstream_pos.push_back( pose.fold_tree().downstream_jump_residue( n ) );
	}

	ObjexxFCL::FArray1D<int> cuts( num_jump-1 );
	Size count( 0 );
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( cutpoints[n] == static_cast<int>( i ) ) continue;
		count++;
		cuts( count ) = cutpoints[ n ];
	}

	Size const nres( pose.total_residue() );

	// Just brute-force iterate through to find a jump we can remove.
	Size jump_to_remove( 0 );
	for ( Size k = 1; k <= num_jump; k++ ) {
		FArray1D< bool > partition_definition( nres, false );
		pose.fold_tree().partition_by_jump( k, partition_definition );
		if ( partition_definition( i ) != partition_definition( i+1 ) ){
			jump_to_remove = k; break;
		}
	}

	bool success( false );

	ObjexxFCL::FArray2D<int> jump_point( 2, num_jump-1 );

	count = 0;
	for ( Size n = 1; n <= num_jump; n++ ) {
		if ( n == jump_to_remove ) continue;
		count++;
		if ( upstream_pos[ n ] < downstream_pos[ n ] ){
			jump_point( 1, count ) = upstream_pos[ n ];
			jump_point( 2, count ) = downstream_pos[ n ];
		} else {
			jump_point( 1, count ) = downstream_pos[ n ];
			jump_point( 2, count ) = upstream_pos[ n ];
		}
	}

	FoldTree f( nres );

	success = f.tree_from_jumps_and_cuts( nres, num_jump-1, jump_point, cuts, 1, false /*verbose*/ );



	if ( !success ) utility_exit_with_message( "FAIL to remove cutpoint "+string_of( i ) );

	pose.fold_tree( f );

}

////////////////////////////////////////////////////////
void
remove_cutpoints_closed( pose::Pose & pose ){

	// Make a list of each cutpoint_closed.
	for ( Size i = 1; i < pose.total_residue(); i++ ) {

		if ( pose.fold_tree().is_cutpoint( i ) &&
				 pose.residue( i   ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
				 pose.residue( i+1 ).has_variant_type( chemical::CUTPOINT_UPPER ) ){
			remove_cutpoint_closed( pose, i ); // this will cycle through to find a jump that is removable.
		}
	}

}

////////////////////////////////////////////////////////
void
virtualize_5prime_phosphates( pose::Pose & pose ){

	for ( Size i = 1; i <= pose.total_residue(); i++ ) {

		if ( i==1 || (( pose.fold_tree().is_cutpoint( i-1 ) &&
									 !pose.residue( i-1 ).has_variant_type( chemical::CUTPOINT_LOWER ) &&
									 !pose.residue( i   ).has_variant_type( chemical::CUTPOINT_UPPER ) ) &&
				 pose.residue(i).is_RNA()) ){
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", i );
		}

	}

}

///////////////////////////////////////////////////////////////////////
// One of these days I'll put in residue 1 as well.
void
print_internal_coords( core::pose::Pose const & pose ) {

	using namespace core::id;
	using numeric::conversions::degrees;

	for (Size i = 2;  i <= pose.total_residue(); i++ ){

		conformation::Residue const & rsd( pose.residue( i ) ) ;

		std::cout << "----------------------------------------------------------------------" << std::endl;
		std::cout << "RESIDUE: " << rsd.name3() << " " << rsd.seqpos() << std::endl;

		for (Size j = 1; j <= rsd.natoms(); j++ ){
			core::kinematics::tree::AtomCOP current_atom ( & pose.atom_tree().atom( AtomID(j,i) ) );
			core::kinematics::tree::AtomCOP input_stub_atom1( current_atom->input_stub_atom1() );
			core::kinematics::tree::AtomCOP input_stub_atom2( current_atom->input_stub_atom2() );
			core::kinematics::tree::AtomCOP input_stub_atom3( current_atom->input_stub_atom3() );

			if ( !(current_atom && input_stub_atom1 && input_stub_atom2 && input_stub_atom3) )  continue;

			if ( current_atom->is_jump() ) {

				core::kinematics::tree::AtomCOP jump_stub_atom1( current_atom->stub_atom1() );
				core::kinematics::tree::AtomCOP jump_stub_atom2( current_atom->stub_atom2() );
				core::kinematics::tree::AtomCOP jump_stub_atom3( current_atom->stub_atom3() );

				// Now try to reproduce this jump based on coordinates of atoms.
				std::cout << 	 "JUMP! " <<
					//					A( 5, rsd.atom_name( j )) << " " <<
					"STUB ==> " <<
					"FROM " << input_stub_atom1->id().rsd() << " " <<
					pose.residue( (input_stub_atom1->id()).rsd() ).atom_name( (input_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue( (input_stub_atom2->id()).rsd() ).atom_name( (input_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue( (input_stub_atom3->id()).rsd() ).atom_name( (input_stub_atom3->id()).atomno() );
				std::cout << "TO " << current_atom->id().rsd() << " " <<
					pose.residue( (jump_stub_atom1->id()).rsd() ).atom_name( (jump_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue( (jump_stub_atom2->id()).rsd() ).atom_name( (jump_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue( (jump_stub_atom3->id()).rsd() ).atom_name( (jump_stub_atom3->id()).atomno() ) << std::endl;
				std::cout << " MY JUMP: " << current_atom->jump() << std::endl;

				kinematics::Stub const input_stub( input_stub_atom1->xyz(), input_stub_atom1->xyz(), input_stub_atom2->xyz(), input_stub_atom3->xyz());
				kinematics::Stub const jump_stub ( jump_stub_atom1->xyz(), jump_stub_atom1->xyz(), jump_stub_atom2->xyz(), jump_stub_atom3->xyz());

				std::cout << " MY JUMP: " << 	kinematics::Jump( input_stub, jump_stub ) << std::endl;


			} else {
				std::cout << "ICOOR_INTERNAL  " <<
					A( 5, rsd.atom_name( j )) << " " <<
					F(11,6, degrees(	pose.atom_tree().dof( DOF_ID( current_atom->id(), id::PHI ) ) ) )  << " " <<
					F(11,6, degrees(  pose.atom_tree().dof( DOF_ID( current_atom->id(), id::THETA ) ) ) ) << " " <<
					F(11,6, pose.atom_tree().dof( DOF_ID( current_atom->id(), id::D ) ) )    << "  " <<
					pose.residue( (input_stub_atom1->id()).rsd() ).atom_name( (input_stub_atom1->id()).atomno() ) << "  " <<
					pose.residue( (input_stub_atom2->id()).rsd() ).atom_name( (input_stub_atom2->id()).atomno() ) << "  " <<
					pose.residue( (input_stub_atom3->id()).rsd() ).atom_name( (input_stub_atom3->id()).atomno() ) << "  " <<
					" [" <<
					" " << (input_stub_atom1->id()).rsd() <<
					" " << (input_stub_atom2->id()).rsd() <<
					" " << (input_stub_atom3->id()).rsd() <<
					"]" <<
					std::endl;
			}

		}

// 		/////////////////////////////////
// 		if ( option[ rsd_type_set]() == "rna" ){
// 			std::cout << "Give me LOWER-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( pose.residue(i-1).xyz( "O3'" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;
// 			std::cout << "Give me OP2-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( rsd.xyz( "OP2" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;
// 			std::cout << "Give me OP1-P-O5' " << 180.0 - numeric::conversions::degrees( numeric::angle_radians( rsd.xyz( "OP1" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ) ) ) << std::endl;

// 			Real main_torsion = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i-1).xyz( "O3'" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );
// 			Real main_torsion1 = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i).xyz( "OP2" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );
// 			Real main_torsion2 = numeric::conversions::degrees(  numeric::dihedral_radians( pose.residue(i).xyz( "OP1" ), rsd.xyz( "P" ), rsd.xyz( "O5'" ), rsd.xyz( "C5'" ) ) );

// 			std::cout << "OFFSETS: " << main_torsion1 - main_torsion << " " << main_torsion2 - main_torsion1 << std::endl;
// 		}

	}
}


////////////////////////////////////////////////////////////////////////////////////
// Apparently can only reroot the tree at a "vertex", i.e. beginning or end of an "edge".
bool
possible_root( core::kinematics::FoldTree const & f, core::Size const & n ){
	return f.possible_root( n );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
get_rigid_body_jumps( core::pose::Pose const & pose ) {

	utility::vector1< Size > rigid_body_jumps;

	TR.Debug << "Initialize RigidBodyMover: Is last residue virtual? " <<  pose.residue( pose.total_residue() ).name3()  << std::endl;

	if ( pose.residue( pose.total_residue() ).name3() != "XXX" ) return rigid_body_jumps; // must have a virtual anchor residue.

	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ){
		TR.Debug << "checking jump: " <<  pose.fold_tree().upstream_jump_residue( n ) << " to " <<  pose.fold_tree().downstream_jump_residue( n ) << std::endl;
		if ( pose.fold_tree().upstream_jump_residue( n ) == (int)pose.total_residue()  ||
				 pose.fold_tree().downstream_jump_residue( n ) == (int)pose.total_residue()  ){
			TR.Debug << "found jump to virtual anchor at: " << n << std::endl;

			rigid_body_jumps.push_back( n );
			if ( rigid_body_jumps.size() > 1 ) {
				TR.Debug << "found moveable jump!" << std::endl;
			}
		}
	}

	return rigid_body_jumps;
}


////////////////////////////////////////////////////////////////////////////////////
bool
let_rigid_body_jumps_move( core::kinematics::MoveMap & movemap,
													 pose::Pose const & pose,
													 bool const move_first_rigid_body /* = false */ ){

	utility::vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
	Size const found_jumps = rigid_body_jumps.size();
	if ( found_jumps <= 1 )	 return false; // nothing to rotate/translate relative to another object.

	Size start = ( move_first_rigid_body ) ? 1 : 2;
	for ( Size n = start; n <= rigid_body_jumps.size(); n++ ) movemap.set_jump( rigid_body_jumps[n], true );

	return true;
}

/////////////////////////////////////////////////////////////////////////////////////
void
translate_virtual_anchor_to_first_rigid_body( pose::Pose & pose ){

	utility::vector1< Size > const rigid_body_jumps = get_rigid_body_jumps( pose );
	if ( rigid_body_jumps.size() == 0 ) return;

	Size const nres = pose.total_residue(); // This better be a virtual residue -- checked in get_rigid_body_jumps() above.

	Size anchor_rsd = pose.fold_tree().downstream_jump_residue( rigid_body_jumps[1] );
	if ( anchor_rsd == nres ) anchor_rsd = pose.fold_tree().upstream_jump_residue( rigid_body_jumps[1] );

	Vector anchor1 = pose.xyz( id::AtomID( 1, anchor_rsd ) );
	Vector root1   = pose.xyz( id::AtomID( 1, nres ) );
	Vector const offset = anchor1 - root1;

	for ( Size j = 1; j <= pose.residue( nres ).natoms(); j++ ){
		id::AtomID atom_id( j, nres );
		pose.set_xyz( atom_id, pose.xyz( atom_id ) + offset );
	}

}

////////////////////////////////////////////////////////////////////////////////////////
bool
involved_in_phosphate_torsion( std::string atomname )
{
	utility::vector1< std::string > const & atoms_involved = core::chemical::rna::atoms_involved_in_phosphate_torsion;

	for ( Size n = 1; n <= atoms_involved.size(); n++ ){
		if (  atomname == atoms_involved[ n ] ) return true;
	}
	return false;
}


////////////////////////////////////////////////////////////////////////////////////////
bool
mutate_position( pose::Pose & pose, Size const i, char const & new_seq ){

	using namespace core::conformation;
	using namespace core::chemical;

	if ( new_seq == pose.sequence()[i-1] ) return false;

	ResidueTypeSet const & rsd_set = pose.residue( i ).residue_type_set();

	ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( new_seq ).match_variants( pose.residue(i).type() ).select( rsd_set )[1] );
	ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );

	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );

	return true;
}

///////////////////////////////////////////////////////////////////////////////
void
set_output_res_num( pose::Pose & extended_pose,
										utility::vector1< Size > const & output_res_num ){
	using namespace pose;
	if ( output_res_num.size() == 0 ) return;
	runtime_assert( output_res_num.size() == extended_pose.total_residue() );

	PDBInfoOP pdb_info = new PDBInfo( extended_pose );
	pdb_info->set_numbering( output_res_num );
	extended_pose.pdb_info( pdb_info );
}



//////////////////////////////////////////////////////////////////////////////////////
void
figure_out_base_pair_partner( pose::Pose & pose, std::map< Size, Size > & partner,
															bool const strict )
{

	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::chemical::rna;
	using namespace core::chemical;
	using namespace core::conformation;

	partner.clear();

	ScoreFunctionOP lores_scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );
	(*lores_scorefxn)( pose );
	lores_scorefxn->show( std::cout, pose );

	RNA_ScoringInfo const & rna_scoring_info( rna_scoring_info_from_pose( pose ) );
	RNA_FilteredBaseBaseInfo const & rna_filtered_base_base_info( rna_scoring_info.rna_filtered_base_base_info() );
	Energy_base_pair_list scored_base_pair_list( rna_filtered_base_base_info.scored_base_pair_list() );

	Size k( 0 ), m( 0 );
	for ( Energy_base_pair_list::const_iterator it = scored_base_pair_list.begin();
				it != scored_base_pair_list.end(); ++it ){

		Base_pair const base_pair = it->second;

		Size const i = base_pair.res1;
		Size const j = base_pair.res2;

		k = base_pair.edge1;
		m = base_pair.edge2;

		Residue const & rsd_i( pose.residue( i ) );
		Residue const & rsd_j( pose.residue( j ) );

		if ( ( k == WATSON_CRICK && m == WATSON_CRICK
					 && base_pair.orientation == 1 )  &&
				 possibly_canonical( rsd_i.aa(), rsd_j.aa() ) &&
				 pose.torsion( id::TorsionID( i, id::CHI, 1 ) ) > 0  && //Need to check syn/anti
				 pose.torsion( id::TorsionID( j, id::CHI, 1 ) ) > 0     //Need to check syn/anti
				 )
			{

				if (strict && !possibly_canonical_strict( rsd_i.aa(), rsd_j.aa() ) ) continue;

				std::string atom1, atom2;
				get_watson_crick_base_pair_atoms( rsd_i.aa(), rsd_j.aa(), atom1, atom2 );
				if ( ( rsd_i.xyz( atom1 ) - rsd_j.xyz( atom2 ) ).length() < 3.5 ) {
					partner[ i ] = j;
					partner[ j ] = i;
				}

			}

	}

}

//////////////////////////////////////////////////////////////////////////////////////////////
void
process_input_file( std::string const & input_file,
										utility::vector1< pose::PoseOP > & pose_list,
										bool is_pdb /*= false*/,
										bool coarse_rna /* = false */)
{
	using namespace core::io::silent;
	using namespace protocols::farna;

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA );

	if ( is_pdb ){

		pose::PoseOP pose_op( new pose::Pose );
		core::import_pose::pose_from_pdb( *pose_op, *rsd_set, input_file );
		//			ensure_phosphate_nomenclature_matches_mini( *pose_op );
		figure_out_reasonable_rna_fold_tree( *pose_op );
		pose_list.push_back( pose_op );

	} else { //its a silent file.

			SilentFileData silent_file_data;
			silent_file_data.read_file( input_file );
			for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
							end = silent_file_data.end(); iter != end; ++iter ) {
				pose::PoseOP pose_op( new pose::Pose );
				iter->fill_pose( *pose_op );
				pose_list.push_back( pose_op );
			}

	}

	// further cleanup.
	for (Size n = 1; n <= pose_list.size(); n++ ){

		pose::PoseOP pose_op = pose_list[ n ];

		remove_cutpoints_closed( *pose_op );

		if ( coarse_rna && !pose_op->residue(1).is_coarse() ){
			pose::Pose coarse_pose;
			make_coarse_pose( *pose_op, coarse_pose );
			*pose_op = coarse_pose;
		}

		virtualize_5prime_phosphates( *pose_op );
	}

	if ( pose_list.size() < 1)  {
		utility_exit_with_message(  "No structure found in input file  " + input_file );
	}

}


} //farna
} //protocols
