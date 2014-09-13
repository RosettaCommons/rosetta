// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/GraftOneCDRLoop.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @detailed
/// @author Jianqing Xu (xubest@gmail.com)

#include <protocols/antibody/GraftOneCDRLoop.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/antibody/Ab_TemplateInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/util.tmpl.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>



static basic::Tracer TRG("protocols.antibody.GraftOneCDRLoop");

namespace protocols {
namespace antibody {
using namespace core;



GraftOneCDRLoop::GraftOneCDRLoop() {}



GraftOneCDRLoop::GraftOneCDRLoop( CDRNameEnum const & cdr_name,
                                  AntibodyInfoOP antibody_info,
                                  Ab_TemplateInfoOP ab_t_info) : Mover( "GraftOneCDRLoop" ) {

	set_default();

	cdr_name_  = cdr_name;
	ab_info_   = antibody_info;
	ab_t_info_ = ab_t_info;


	init();

}


GraftOneCDRLoop::~GraftOneCDRLoop() {}


void GraftOneCDRLoop::set_default() {
	ab_info_=NULL;
	ab_t_info_=NULL;
	stem_copy_size_ = 2;
	flank_size_ = 4 ;
	// JQX:  the default value of flank_size_ is equle to 4, meaning there are 4 residues
	//       on the C-ter and N-ter of the actual loop
	//       However, based on the old R2 antibody code, only 3 residues on each stem

	stem_not_graft_ = false;
	benchmark_ = false;
	preprocessing_script_version_ = "R3_Python"; //JQX: R2_Perl or R3_Python

}



void GraftOneCDRLoop::init() {

	//idealize the loop
	idealize::IdealizeMover idealizer;
	idealizer.fast( false );

}




void GraftOneCDRLoop::finalize_setup() {

	if(scorefxn_) {
		scorefxn_=core::scoring::get_score_function();
	}

	if(stem_not_graft_) {
		stem_copy_size_ = 0 ;
	}

}




void GraftOneCDRLoop::apply( pose::Pose & pose_in ) {


	TRG<<"Start to Graft CDRs   "<< ab_info_->get_CDR_name(cdr_name_) <<" ............"<<std::endl;
	TRG<<"flank_size: "<<flank_size_<<std::endl;

	finalize_setup();



	core::Size query_start = ab_info_->get_CDR_loop(cdr_name_).start();
	core::Size query_end   = ab_info_->get_CDR_loop(cdr_name_).stop();
	Size query_size = ( query_end - query_start )+1;


	// create a sub pose with  4 flanking residues on either side of CDR loop
	// TRG<<"query_start="<<query_start<<std::endl;
	// TRG<<"query_end="<<query_end<<std::endl;
	// TRG<<"flank_size="<<flank_size_<<std::endl;
	// TRG<<"query_size="<<query_size<<std::endl;
	// TRG<<"truncated_pose will be from "<<query_start-flank_size_<<" to "<<query_end+flank_size_<<std::endl;
	pose::Pose truncated_pose( pose_in, query_start-flank_size_, query_end+flank_size_ );
	pose::Pose template_pose = ab_t_info_->get_one_template_pose(ab_info_->get_CDR_name(cdr_name_)) ;
	TRG<< "template_pose: "<<std::endl;
	TRG<< template_pose<<std::endl;

	// Just want to make life a litter easier
	//TRG<<"**************&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*************"<<std::endl;
	//TRG<<pose_in.sequence()<<std::endl;
	//TRG<<"  the template sequence:  "<<template_pose.sequence()<<std::endl;
	//TRG<<"trucated_query sequence:  "<<truncated_pose.sequence()<<std::endl;
	//TRG<<"**************&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*************"<<std::endl;


	// create atom map for superimposing 2 flanking resiudes
	id::AtomID_Map< id::AtomID > atom_map;
	pose::initialize_atomid_map( atom_map, template_pose, id::BOGUS_ATOM_ID );



	//   ****AAAAAAAAAAAAAAAAAAA****  the template pose should have 4 residues each side
	//    @@@AAAAAAAAAAAAAAAAAAA@@@   the real alignment only based on 3 residues each side due to error in R2_Perl
	//   @@@@AAAAAAAAAAAAAAAAAAA@@@@  the corrected alignment based on 4 residues each side in R3_Python

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ OptionKeys::antibody::preprocessing_script_version ].user() ) {
		preprocessing_script_version_ = option[ OptionKeys::antibody::preprocessing_script_version ]() ;
	}

	Size correction = 0;
	if(preprocessing_script_version_ == "R2_Perl") correction = 1;
	if(preprocessing_script_version_ == "R3_Python") correction = 0;

	for( Size start_stem = 1+correction; start_stem <= flank_size_; ++start_stem ) {
		/// starting from the 2nd residue in the l1-3.pdb, H1-3.pdb
		Size const ref_stem ( start_stem  );
		for( Size j=1; j <= 4; j++ ) {    /// four backbone heavy atoms
//          		TRG<<"j="<<j<<"  start_stem_in_template_pose="<<start_stem<<"  ref_stem_in_truncated_pose_="<<ref_stem<<std::endl;
			id::AtomID const id1( j, start_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}

	// start at the end of the actual loop
	for( Size end_stem = query_size+flank_size_+1; end_stem <= query_size+flank_size_+flank_size_-correction; ++end_stem ) {
		Size const ref_stem ( end_stem);
		//   if(template_name_ == "h3") Size const ref_stem(end_stem+1);
		for( Size j=1; j <= 4; j++ ) {    /// four backbone heavy atoms
//          	TRG<<"j="<<j<<"  end_stem_in_template_pose="<<end_stem<<"  ref_stem_in_truncated_pose_="<<ref_stem<<std::endl;
			id::AtomID const id1( j, end_stem );
			id::AtomID const id2( j, ref_stem );
			atom_map[ id1 ] = id2;
		}
	}


	scoring::superimpose_pose( template_pose, truncated_pose, atom_map );


	// TODO:
	// JQX:
	// when copying the backbone heavy atom coordinates, you copy two more residues on each side
	// pay attention to this, it is "2", not "3" here. "3" was used to superimpose like the code did above.
	// No reason for that, just to match R2_antibody.
	// This may be the reason why you need hfr to reasign the phi, psi, or omega angle of h3_start-1
	// when doing H3 loop modeling. But this is still weird, because you need to take care
	// of 2 residues on each side then, not just one more on h3_start-1 or +1 position using hfr.pdb.
	// Maybe I should remove this number 2, just copy the loop itself !!!!!

	for ( Size i=query_start-stem_copy_size_; i <= query_end+stem_copy_size_; ++i ) {
		Size template_num ( i - (query_start-5) );  // this "5" is the default in the L1.pdb, you have 4 residues before the 1st one
//        TRG<<" i="<<i<<"    template_num="<<template_num<<std::endl;
		conformation::Residue const & source_rsd( pose_in.residue( i ) );
		conformation::Residue const & target_rsd( template_pose.residue( template_num) );
//        TRG<<source_rsd<<std::endl;
//        TRG<<target_rsd<<std::endl;

		Size const natoms( source_rsd.natoms() );

		bool any_missing( false );
		id::AtomID_Mask missing( false );
		// dimension the missing-atom mask
		pose::initialize_atomid_map( missing, pose_in );

		if( source_rsd.name() != target_rsd.name() ) pose_in.set_secstruct( i, 'X' );
		for ( Size j=1; j<= natoms; ++j ) {
			std::string const & atom_name( source_rsd.atom_name(j) );
			if ( target_rsd.has( atom_name ) ) {
				pose_in.set_xyz( id::AtomID( source_rsd.atom_index( atom_name),i ), target_rsd.xyz( atom_name ) );
				//make sure the source_rsd in pose_in only has 4 backbone atoms
				//my script should take care of it
			} else {
				TRG<<"watch out !!! missing "<<atom_name << "   !!!!!" <<std::endl;
				any_missing = true;
				missing[ id::AtomID( pose_in.residue_type(i).atom_index(source_rsd.atom_name(j)), i ) ] = true;
			}
		}

		if ( any_missing ) {
			pose_in.conformation().fill_missing_atoms( missing );
		}
	}

//    pose_in.dump_pdb(template_name_+"_graft");


}









std::string GraftOneCDRLoop::get_name() const {
	return "GraftOneCDRLoop";
}

// copy ctor
GraftOneCDRLoop::GraftOneCDRLoop( GraftOneCDRLoop const & rhs ) : Mover(rhs) {
	initForEqualOperatorAndCopyConstructor(*this, rhs);
}


///@brief assignment operator
GraftOneCDRLoop & GraftOneCDRLoop::operator=( GraftOneCDRLoop const & rhs ) {
	//abort self-assignment
	if (this == &rhs) return *this;
	Mover::operator=(rhs);
	initForEqualOperatorAndCopyConstructor(*this, rhs);
	return *this;
}

void GraftOneCDRLoop::initForEqualOperatorAndCopyConstructor(GraftOneCDRLoop & lhs, GraftOneCDRLoop const & rhs) {
	lhs.flank_size_ = rhs.flank_size_;
	lhs.stem_copy_size_=rhs.stem_copy_size_;
	lhs.cdr_name_=rhs.cdr_name_;
	lhs.ab_info_=rhs.ab_info_;
	lhs.ab_t_info_=rhs.ab_t_info_;
	lhs.benchmark_=rhs.benchmark_;
	lhs.stem_not_graft_=rhs.stem_not_graft_;
}





}  // namespace antibody
}  // namespace protocols
