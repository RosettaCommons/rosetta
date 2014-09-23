// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody/H3CterInsert.cc
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#include <protocols/antibody/H3CterInsert.hh>

#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <basic/database/open.hh>




static thread_local basic::Tracer TR( "protocols.antibody.H3CterInsert" );
using namespace core;



namespace protocols {
namespace antibody {





H3CterInsert::H3CterInsert() : Mover() {
	user_defined_ = false;
}



H3CterInsert::H3CterInsert(antibody::AntibodyInfoOP  antibody_info, bool camelid ) : Mover() {
	user_defined_ = true;
	init(antibody_info, camelid, false);
}




// H3CterInsert default destructor
H3CterInsert::~H3CterInsert() {}




void H3CterInsert::init(AntibodyInfoOP antibody_info, bool camelid, bool benchmark) {
	Mover::type( "H3CterInsert" );

	set_default();

	if ( user_defined_ ) {
		is_camelid_ = camelid;
		benchmark_  = benchmark;
		ab_info_= antibody_info;
	}

	//    setup_objects();

	read_H3_cter_fragment( ) ;

}





void H3CterInsert::set_default() {
	is_camelid_ = false;
	benchmark_ = false;
}


//void H3CterInsert::setup_objects(){
//}



//void H3CterInsert::finalize_setup(core::pose::Pose & pose ){

//}



void H3CterInsert::apply(pose::Pose & pose) {

//    finalize_setup(pose);




	TR <<  "Inserting CDR H3 C-ter Fragments" << std::endl;

	Size loop_begin(0), loop_end(0), cutpoint(0);



	loop_begin = ab_info_->get_CDR_loop(h3).start()-1;  //JQX: to match R2_Antibody
	cutpoint   = ab_info_->get_CDR_loop(h3).cut(); // keep the cutpoint unchanged
	loop_end   = ab_info_->get_CDR_loop(h3).stop()+2; //JQX: to match R2_Antibody

	pose.fold_tree( * ab_info_->setup_simple_fold_tree(loop_begin, cutpoint, loop_end, pose) );
	TR<<pose.fold_tree()<<std::endl;

	bool success = false;
	while(success == false) {

		fragment::FragData f;
		Size random_H3_ter(0);
		random_H3_ter = numeric::random::rg().random_range( 1, H3_base_library_.size() );
		f = H3_base_library_[ random_H3_ter ];
		//    TR<<H3_base_library_.size()<<std::endl;
		// f = H3_base_library_[ 1 ]; //JQX: a test case to match R2_Antibody
		//JQX: realize R2_antibody the fragment data was vector, not vector1
		//JQX: in another word, R2->20 match R3->21
		//TODO:
		//JQX: for the M18 case, the H3_base_library_[ 30 ] has a problem in both R2 and R3
		//     OK, this is due to the database used omega-phi-psi, but Rosetta use phi-psi-omega
		//     the omega is the the omega of last residue
		//TR<<f<<std::endl;

		//inserting base dihedrals
		Size cter_insertion_pos( is_camelid_ ? 4 : 2 ); // not sure why 4 for camelid


		if( (ab_info_->get_CDR_loop(h3).stop()-cter_insertion_pos) <= ab_info_->get_CDR_loop(h3).start() ) {
			TR << "H3 LOOP IS TOO SHORT: CAN NOT USE C-TERM INFORMATION"<< std::endl;
		} else {
			f.apply( pose, ab_info_->get_CDR_loop(h3).stop() - cter_insertion_pos, ab_info_->get_CDR_loop(h3).stop() + 1 );

		}

		//JQX: it seems movemap is not required..... kind of weird, but it works! maybe there's a default movemap


		success = CDR_H3_cter_filter(pose, ab_info_);
		//JQX: this seems ridiculous, right? You are inserting a kinked or extended structure, why you want
		//     to check again? The reason is that ... the H3_CTERM files is wrong! It is supposed to be
		//     phi-psi-omega, but the info there is omega-phi-psi, in another word, it is using a wrong omega.
		//     If the H3_CTERM file is fixed, this filter is certainly not necessary for sure.
	}


	TR <<  "Finished Inserting CDR H3 C-ter Fragments" << std::endl;

	return;




}




std::string H3CterInsert::get_name() const {
	return "H3CterInsert";
}







/// FIXME: JQX the get_CDR_loop is so weird here
void H3CterInsert::read_H3_cter_fragment( ) {
	using namespace fragment;

	TR <<  "Reading CDR H3 C-ter Fragments" << std::endl;

	bool is_kinked = false;
	bool is_extended  = false;
	if (ab_info_->get_H3_kink_type() == Kinked) is_kinked = true;
	if (ab_info_->get_H3_kink_type() == Extended) is_extended = true;


	// extract single letter aa codes for the chopped loop residues
	TR<< ab_info_->get_CDR_loop(h3)<<std::endl;

	Size cdr_h3_size = (     ab_info_->get_CDR_loop(h3).stop() - ab_info_->get_CDR_loop(h3).start()  ) + 1;

	TR<<"cdr_h3_size="<<cdr_h3_size<<std::endl;
	std::string aa_1name = ab_info_->get_CDR_sequence_with_stem(h3,2,2);

	TR<<"aa_1name.size()="<<aa_1name.length()<<std::endl;

	TR<<aa_1name<<std::endl;





	// used only when no length & kink match are found
	utility::vector1< FragData > H3_base_library_seq_kink;

	// used only when no (length & kink) or (length & seq) are found
	utility::vector1< FragData > H3_base_library_kink;


	// file is read in from where other contraints are supposed to exist
	if( ab_info_->is_camelid() ) {
		H3_ter_library_filename_ = basic::database::full_name( "sampling/antibodies/Fragments/camelid_H3_CTERM") ;
	} else {
		H3_ter_library_filename_ = basic::database::full_name( "sampling/antibodies/Fragments/H3_CTERM") ;
	}

	// Read the file defined by command line option
	utility::io::izstream H3_ter_library_stream( H3_ter_library_filename_ );

	// Check to see if file exists
	if ( !H3_ter_library_stream ) {
		TR << "[Error]: Could not open H3 base library file: "
		   << H3_ter_library_filename_ << std::endl
		   << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
		std::exit( EXIT_FAILURE );
	}

	std::string pdb_name;
	std::string res_no;
	char res_name;
	Real phi(0.0);
	Real psi(0.0);
	Real omega(0.0);
	Size H3_length(0);
	Real resolution(0.0);
	std::string base_type;

	Size pdb_H3_length = cdr_h3_size;
	Size h3_base_frag_size( ab_info_->is_camelid() ? 6 : 4 );

	bool end_not_reached(true);
	while(end_not_reached) {
		bool seq_match( true );
		bool base_match( false );

		FragData f;
		f.set_valid( true );

		for ( Size i = 1; i <= h3_base_frag_size; ++i ) {
			H3_ter_library_stream >> pdb_name  >> res_no     >> res_name
			                      >> omega     >> phi        >> psi
			                      >> H3_length >> resolution >> base_type
			                      >> std::skipws;

			if(H3_ter_library_stream.eof()) {
				end_not_reached = false;
				break;
			}

			if( res_name != aa_1name[aa_1name.size() - 5 + i] ) {
				seq_match = false;
			}
			//TR<<res_name<<"   .vs.    "<<aa_1name[aa_1name.size() - 5 + i]<<"    "<<seq_match<<std::endl;

			BBTorsionSRFDOP res_torsions(
			    new BBTorsionSRFD( 3, 'L', res_name ) ); // 3 protein torsions

			// ugly numbers 1-3, but pose.set_phi also uses explicit numbers
			res_torsions->set_torsion   ( 1, phi   );
			res_torsions->set_torsion   ( 2, psi   );
			res_torsions->set_torsion   ( 3, omega );
			res_torsions->set_secstruct ( 'L' );
			f.add_residue( res_torsions );
		}


		//regular antibody
		if( ab_info_->is_camelid() == false) {
			if( is_kinked && base_type == "KINK" ) {
				base_match = true;
			} else if( is_extended && base_type == "EXTENDED" ) {
				base_match = true;
			} else if( !is_kinked && !is_extended && base_type == "NEUTRAL" ) {
				base_match = true;
			}
		}
		//camelid antibody
		else {
			if( is_extended && base_type == "EXTENDED" ) {
				base_match = true;
			} else if( ( is_kinked && base_type == "KINK" ) || ( is_kinked && base_type == "EXTENDED" ) ) {
				base_match = true;
			} else if( !is_kinked && !is_extended ) {
				base_match = true;
			}
		}


		//camelid antibody
		if( ab_info_->is_camelid() && end_not_reached && base_match ) {

			H3_base_library_.push_back( f );
		}

		// regular antibody
		else if( end_not_reached && (H3_length==pdb_H3_length) && base_match ) {
			//TR<< "111111111111      H3_length = "<<H3_length<<"     pdb_h3_length="<<pdb_H3_length<<"    base_match="<<base_match<<std::endl;
			//TR<<pdb_name<<"   "<<omega<<"    "<<phi<<"    "<<psi<<std::endl;

			H3_base_library_.push_back( f );
		}


		if( end_not_reached && seq_match && base_match ) {

			H3_base_library_seq_kink.push_back( f );
		}
		if( end_not_reached && base_match  ) {

			H3_base_library_kink.push_back( f );
		}
	}




	//for (int i=1;i<=H3_base_library_.size();i++) {
	//    TR<<H3_base_library_[i]<<std::endl;
	//}

	H3_ter_library_stream.close();
	H3_ter_library_stream.clear();

	// if no match found based on sequence and kink match criterion
	// then choose based on size and kink match criterion
	// if still no match, then choose based only on kink
	if( H3_base_library_.size() == 0 ) {
		H3_base_library_ = H3_base_library_seq_kink;
	}
	if( H3_base_library_.size() == 0 ) {
		H3_base_library_ = H3_base_library_kink;
	}

	TR <<  "Finished reading CDR H3 C-ter Fragments" << std::endl;
	return;


}

















}//antibody
}//protocols





