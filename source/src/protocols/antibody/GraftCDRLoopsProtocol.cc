// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/GraftCDRLoopsProtocol.cc
/// @brief Build a homology model of an antibody
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )



#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <protocols/loops/loops_main.hh>

#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/antibody/GraftOneCDRLoop.hh>
#include <protocols/antibody/CloseOneCDRLoop.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/Ab_TemplateInfo.hh>
#include <protocols/antibody/GraftCDRLoopsProtocol.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/CDRsMinPackMin.hh>
#include <protocols/antibody/GraftedStemOptimizer.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::format;


using basic::Error;
using basic::Warning;

static basic::Tracer TR( "protocols.antibody.GraftCDRLoopsProtocol" );
using namespace core;

namespace protocols {
namespace antibody {

// default constructor
GraftCDRLoopsProtocol::GraftCDRLoopsProtocol() : Mover() {
	init();
}

// default destructor
GraftCDRLoopsProtocol::~GraftCDRLoopsProtocol() = default;

//clone
protocols::moves::MoverOP GraftCDRLoopsProtocol::clone() const {
	return( protocols::moves::MoverOP( new GraftCDRLoopsProtocol() ) );
}


void GraftCDRLoopsProtocol::init() {
	Mover::type( "GraftCDRLoopsProtocol" );

	set_default();
	register_options();
	init_from_options();
	setup_objects();

}


void GraftCDRLoopsProtocol::set_default() {
	TR <<  "Setting default settings" << std::endl;
	graft_l1_  = true;
	graft_l2_  = true;
	graft_l3_  = true;
	graft_h1_  = true;
	graft_h2_  = true;
	graft_h3_  = true;
	h3_no_stem_graft_=false;
	packonly_after_graft_=false;
	benchmark_ = false;
	camelid_   = false;
	camelid_constraints_ = false;
	sc_min_ = false;
	rt_min_ = false;
	stem_optimize_ = true;
	cst_weight_ = 0.0;

}


void GraftCDRLoopsProtocol::register_options() {
	using namespace basic::options;


	option.add_relevant( OptionKeys::antibody::graft_l1 );
	option.add_relevant( OptionKeys::antibody::graft_l2 );
	option.add_relevant( OptionKeys::antibody::graft_l3 );
	option.add_relevant( OptionKeys::antibody::graft_h1 );
	option.add_relevant( OptionKeys::antibody::graft_h2 );
	option.add_relevant( OptionKeys::antibody::graft_h3 );
	option.add_relevant( OptionKeys::antibody::h3_no_stem_graft );
	option.add_relevant( OptionKeys::antibody::camelid );
	option.add_relevant( OptionKeys::antibody::camelid_constraints );
	option.add_relevant( OptionKeys::antibody::packonly_after_graft);
	option.add_relevant( OptionKeys::antibody::stem_optimize);
	option.add_relevant( OptionKeys::antibody::stem_optimize_size);
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::run::benchmark );
	option.add_relevant( OptionKeys::in::file::native );
	//option.add_relevant( OptionKeys::antibody::sc_min);
	//option.add_relevant( OptionKeys::antibody::rt_min);
}


void GraftCDRLoopsProtocol::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	TR <<  "Reading Options" << std::endl;

	if ( option[ OptionKeys::antibody::graft_l1 ].user() ) {
		set_graft_l1( option[ OptionKeys::antibody::graft_l1 ]() );
	}
	if ( option[ OptionKeys::antibody::graft_l2 ].user() ) {
		set_graft_l2( option[ OptionKeys::antibody::graft_l2 ]() );
	}
	if ( option[ OptionKeys::antibody::graft_l3 ].user() ) {
		set_graft_l3( option[ OptionKeys::antibody::graft_l3 ]() );
	}
	if ( option[ OptionKeys::antibody::graft_h1 ].user() ) {
		set_graft_h1( option[ OptionKeys::antibody::graft_h1 ]() );
	}
	if ( option[ OptionKeys::antibody::graft_h2 ].user() ) {
		set_graft_h2( option[ OptionKeys::antibody::graft_h2 ]() );
	}
	if ( option[ OptionKeys::antibody::graft_h3 ].user() ) {
		set_graft_h3( option[ OptionKeys::antibody::graft_h3 ]() );
	}
	if ( option[ OptionKeys::antibody::h3_no_stem_graft ].user()  ) {
		set_h3_stem_graft( option[ OptionKeys::antibody::h3_no_stem_graft ]()  );
	}
	if ( option[ OptionKeys::antibody::packonly_after_graft ].user()  ) {
		set_packonly_after_graft( option[ OptionKeys::antibody::packonly_after_graft ]()  );
	}
	if ( option[ OptionKeys::antibody::stem_optimize ].user()  ) {
		set_stem_optimize( option[ OptionKeys::antibody::stem_optimize ]()  );
	}
	if ( option[ OptionKeys::antibody::camelid ].user() ) {
		set_camelid( option[ OptionKeys::antibody::camelid ]() );
	}
	if ( option[ OptionKeys::antibody::camelid_constraints ].user() ) {
		set_camelid_constraints( option[ OptionKeys::antibody::camelid_constraints ]() );
	}
	if ( option[ OptionKeys::run::benchmark ].user() ) {
		set_benchmark( option[ OptionKeys::run::benchmark ]() );
	}
	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		set_cst_weight( option[ OptionKeys::constraints::cst_weight ]() );
	}
	//if ( option[ OptionKeys::antibody::sc_min_ ].user() ) {
	// set_sc_min( option[ OptionKeys::antibody::sc_min_ ]() );
	//}
	//if ( option[ OptionKeys::antibody::rt_min_ ].user() ) {
	// set_rt_min( option[ OptionKeys::antibody::rt_min_ ]() );
	//}


	//set native pose if asked for
	if ( option[ OptionKeys::in::file::native ].user() ) {
		core::pose::PoseOP native_pose( new core::pose::Pose() );
		core::import_pose::pose_from_file( *native_pose, option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		set_native_pose( native_pose );
	} else {
		set_native_pose(nullptr);
	}


	if ( camelid_ ) {
		graft_l1_ = false;
		graft_l2_ = false;
		graft_l3_ = false;
	}

}


void GraftCDRLoopsProtocol::setup_objects() {

	ab_info_ = nullptr;
	ab_t_info_ = nullptr;

	graft_sequence_ = nullptr;

	tf_ = nullptr;

	scorefxn_pack_ = core::scoring::get_score_function();

	scorefxn_pack_->set_weight( core::scoring::chainbreak, 1.0 );
	scorefxn_pack_->set_weight( core::scoring::overlap_chainbreak, 10./3. );
	scorefxn_pack_->set_weight( core::scoring::atom_pair_constraint, 1.00 );

}


void GraftCDRLoopsProtocol::finalize_setup( pose::Pose & frame_pose ) {
	TR<<"AAAAAAAA     cst_weight: "<<cst_weight_<<std::endl;

	// check for native and input pose
	if ( !get_input_pose() ) {
		pose::PoseOP input_pose( new pose::Pose(frame_pose) );  //JQX: QUESTION: why owning pointer here
		set_input_pose( input_pose );   // JQX: pass the input_pose to the mover.input_pose_
	}


	pose::PoseOP native_pose;
	if ( !get_native_pose() ) {
		TR << "Danger Will Robinson! Native is an impostor!" << std::endl;
		TR << "   'native_pose' is just a copy of the 'input_pose'    " << std::endl;
		TR << "    since you didn't sepcifiy the native pdb name"<<std::endl;
		native_pose = pose::PoseOP( new pose::Pose(frame_pose) );
	} else {
		native_pose = pose::PoseOP( new pose::Pose( *get_native_pose() ) );
	}

	// JQX: this is the secondary structure from the native pose
	pose::set_ss_from_phipsi( *native_pose );

	set_native_pose( native_pose ); // pass the native pose to the mover.native_pose_


	ab_info_ = AntibodyInfoOP( new AntibodyInfo(frame_pose) );
	ab_t_info_ = Ab_TemplateInfoOP( new Ab_TemplateInfo(graft_l1_, graft_l2_, graft_l3_,
		graft_h1_, graft_h2_, graft_h3_, camelid_) );

	graft_sequence_ = protocols::moves::SequenceMoverOP( new moves::SequenceMover() );


	TR<<" Checking AntibodyInfo object: "<<std::endl<<*ab_info_<<std::endl<<std::endl;
	TR<<" Checking Ab_TemplateInfo object: "<<std::endl<<*ab_t_info_<<std::endl<<std::endl;

	for ( CDRNameEnum it = start_cdr_loop; it <= ab_info_->get_total_num_CDRs(); it=CDRNameEnum(it+1) ) {
		TR << "Creating movers for " << ab_info_->get_CDR_name(it) << std::endl;
		TR << "                  start (chothia): "<<ab_info_->get_CDR_loop(it).start()<<std::endl;
		TR << "                   stop (chothia): "<<ab_info_->get_CDR_loop(it).stop()<<std::endl;

		GraftOneCDRLoopOP graft_one_cdr( new GraftOneCDRLoop( it, ab_info_, ab_t_info_) ) ;
		if ( it == h3 && h3_no_stem_graft_ ) {
			graft_one_cdr->set_no_stem_copy(true);
		}


		graft_one_cdr->enable_benchmark_mode( benchmark_ );
		graft_sequence_->add_mover( graft_one_cdr);


		/*
		CloseOneCDRLoopOP closeone( new Clos eOneCDRLoop( ab_info.get_loop(it->first)->start(),
		ab_info.get_loop(it->first)->stop()   )     );
		closeone->enable_benchmark_mode( benchmark_ );
		graft_sequence_->add_mover( closeone );
		*/
	}

	// TR<<" frame_pose.size()  =  "<< frame_pose.size()<<std::endl;


	// Exact match Aroop's old code in Rosetta 2:
	// graft all CDRs by superimpose stems, then pack the whole new pose

	// When do packing, pack the whole pose, but minimize the CDRs
	tf_ = setup_packer_task(frame_pose);

	cdrs_min_pack_min_ = CDRsMinPackMinOP( new CDRsMinPackMin(ab_info_) );
	cdrs_min_pack_min_ -> set_task_factory(tf_);
	// the tf_ include all the residues, the movemap is to use the deafult one in CDRsMinPackMin, which is the CDRs
	cdrs_min_pack_min_->set_sc_min(sc_min_);
	cdrs_min_pack_min_->set_rt_min(rt_min_);
	cdrs_min_pack_min_->set_turnoff_minimization(packonly_after_graft_);


	optimize_sequence_ = protocols::moves::SequenceMoverOP( new moves::SequenceMover() );

	for ( CDRNameEnum it = start_cdr_loop; it <= ab_info_->get_total_num_CDRs(); it=CDRNameEnum(it+1) ) {
		if (  !(it == h3 && h3_no_stem_graft_)   ) {
			GraftedStemOptimizerOP optimize_oneCDR_stems( new GraftedStemOptimizer(it, ab_info_) );
			optimize_oneCDR_stems->enable_benchmark_mode( benchmark_ );
			optimize_sequence_->add_mover(optimize_oneCDR_stems);
		}
	}


}


//APPLY
void GraftCDRLoopsProtocol::apply( pose::Pose & frame_pose ) {

	using namespace chemical;
	using namespace id;
	using namespace scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::moves;


	TR<<" in the apply function "<<std::endl;


	finalize_setup(frame_pose);


	basic::prof_reset();
	// utility::exit( EXIT_FAILURE, __FILE__, __LINE__);

	pose::set_ss_from_phipsi( frame_pose );

	// display constraints and return
	if ( camelid_constraints_ ) {
		display_constraint_residues( frame_pose );
		return;
	}


	Size nres = frame_pose.size();

	// Storing secondary structure
	utility::vector1<char> secondary_struct_storage;
	for ( Size i = 1; i <= nres; i++ ) {
		secondary_struct_storage.push_back( frame_pose.secstruct( i ) );
		//        TR<<"JQX:   residue: "<<i<<"       secstruct: "<<frame_pose.secstruct(i)<<std::endl;
	}

	graft_sequence_->apply( frame_pose );

	frame_pose.dump_pdb("grafted.pdb");

	cdrs_min_pack_min_->apply(frame_pose);
	frame_pose.dump_pdb("grafted_minimized.pdb");

	if ( stem_optimize_ ) {
		optimize_sequence_->apply(frame_pose);
		frame_pose.dump_pdb("grafted_minimized_stem-optimized.pdb");

		cdrs_min_pack_min_->apply(frame_pose);
		frame_pose.dump_pdb("grafted_minimized_stem-optimized_minimized.pdb");
	}


	// Recover secondary structures
	for ( Size i = 1; i <= nres; i++ ) frame_pose.set_secstruct( i, secondary_struct_storage[ i ] );

	// align pose to native pose
	pose::Pose native_pose;
	if ( get_native_pose() ) native_pose = *get_native_pose();
	else                    native_pose = frame_pose;

	AntibodyInfoOP native_ab_info( new AntibodyInfo( native_pose ) );


	align_to_native( frame_pose, native_pose, ab_info_, native_ab_info );


	basic::prof_show();

}// end apply


std::string GraftCDRLoopsProtocol::get_name() const {
	return "GraftCDRLoopsProtocol";
}


void GraftCDRLoopsProtocol::display_constraint_residues( core::pose::Pose & pose ) {
	// Detecting di-sulfide bond

	Size H1_Cys(0), H3_Cys(0);

	if (      pose.residue( pose.pdb_info()->pdb2pose('H',32 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( 'H', 32 );
	} else if ( pose.residue( pose.pdb_info()->pdb2pose('H',33 ) ).name3() == "CYS" ) {
		H1_Cys = pose.pdb_info()->pdb2pose( 'H', 33 );
	}

	for ( Size ii = ab_info_->get_CDR_loop(h3).start(); ii <= ab_info_->get_CDR_loop(h3).stop(); ii++ ) {

		if ( pose.residue(ii).name3() == "CYS" ) H3_Cys = ii;
	}

	if ( ( H1_Cys != 0 ) && ( H3_Cys != 0 ) ) {
		TR << "CONSTRAINTS: "<< "AtomPair CA " << H1_Cys << " CA " << H3_Cys
			<< " BOUNDED 4.0 6.1 0.6 BOND; mean 5.6 sd 0.6" << std::endl;
	}

	// Specifying extended kink

	Size hfr_46(0), h3_closest(0);
	hfr_46 = pose.pdb_info()->pdb2pose( 'H', 46 );
	if ( ab_info_->get_H3_kink_type() == Extended ) h3_closest = ab_info_->get_CDR_loop(h3).stop() - 5;
	if ( h3_closest != 0 ) {
		TR << "CONSTRAINTS: " << "AtomPair CA " << hfr_46 << " CA " << h3_closest
			<< " BOUNDED 6.5 9.1 0.7 DISTANCE; mean 8.0 sd 0.7" << std::endl;
	}

	return;
} // display_constraint_residues


/// @details  Show the complete setup of the docking protocol
void GraftCDRLoopsProtocol::show( std::ostream & out ) const {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const GraftCDRLoopsProtocol & ab_m_2 ) {
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta 3 Antibody Modeler" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << A( 47, "     CDR Assembling       " ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;
	out << line_marker << " Graft_l1:  " << ab_m_2.graft_l1_<<std::endl;
	out << line_marker << " Graft_l2:  " << ab_m_2.graft_l2_<<std::endl;
	out << line_marker << " Graft_l3:  " << ab_m_2.graft_l3_<<std::endl;
	out << line_marker << " Graft_h1:  " << ab_m_2.graft_h1_<<std::endl;
	out << line_marker << " Graft_h2:  " << ab_m_2.graft_h2_<<std::endl;
	out << line_marker << " Graft_h3:  " << ab_m_2.graft_h3_<<std::endl;
	out << line_marker << "  camelid:  " << ab_m_2.camelid_ <<std::endl;
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}

} // end antibody
} // end protocols
