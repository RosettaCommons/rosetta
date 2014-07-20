// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/smlewis/InterfaceAnalyzer.cc
/// @brief Q&D protocol to run InterfaceAnalyzerMover as protocol
/// @author Steven Lewis, Bryan Der, Ben Stranges

// Unit Headers
#include <devel/anchored_design/InterfaceAnalyzerMover.hh>

// Project Headers
#include <ObjexxFCL/format.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>

#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("apps.pilot.smlewis.InterfaceAnalyzer");

//define local options
basic::options::IntegerOptionKey const jumpnum("jumpnum");
//basic::options::BooleanOptionKey const is_compute_hbond_unsat("is_compute_hbond_unsat");
basic::options::BooleanOptionKey const compute_packstat("compute_packstat");
basic::options::BooleanOptionKey const write_pymol_selections_to_file("write_pymol_selections_to_file");
basic::options::BooleanOptionKey const write_pymol_selections_to_tracer("write_pymol_selections_to_tracer");
basic::options::BooleanOptionKey const tracer_data_print("tracer_data_print");
basic::options::BooleanOptionKey const force_output("force_output");
basic::options::StringVectorOptionKey const fixedchains( "fixedchains" );
basic::options::StringOptionKey const ia_stats_filename( "ia_stats_filename" );
basic::options::BooleanOptionKey const pack_input("pack_input");


// mover deffinition
class IAMover : public protocols::moves::Mover {
public:

  IAMover();

  virtual void apply( core::pose::Pose& pose );

	virtual
	std::string
	get_name() const {
		return "IAMover";
	}

  virtual void assign_IA_mover( devel::anchored_design::InterfaceAnalyzerMoverOP & moverOP, core::pose::Pose & pose);

private:
  //bool firstline_printed_;
	std::string output_name_;
  //std::ofstream output_statsFile_;
	utility::io::ozstream output_statsFile_;
  devel::anchored_design::InterfaceAnalyzerMoverOP IAM_;
	core::scoring::ScoreFunctionOP scorefxn_;
	const char* name_star_;
};

IAMover::IAMover() {
  // variable definitions
  //firstline_printed_=false;
	output_name_=basic::options::option[ia_stats_filename].value();
	name_star_ = output_name_.c_str() ;
	scorefxn_ = core::scoring::get_score_function();
	//add to the output file
	//output_statsFile_.open( name_star );

}
//assign the correct constructor for the mover
void IAMover::assign_IA_mover(devel::anchored_design::InterfaceAnalyzerMoverOP & moverOP, core::pose::Pose & pose){

  //
  if(pose.conformation().num_chains() <= 2){
    moverOP = new  devel::anchored_design::InterfaceAnalyzerMover(
								 basic::options::option[ jumpnum ].value(),
								 basic::options::option[ tracer_data_print ].value(),
								 core::scoring::get_score_function(),
								 basic::options::option[ compute_packstat ].value(),
								 basic::options::option[ pack_input ].value(),
								 true, //pack the separated poses
								 false //use the input pose name instead of jobname
								 );
    return;
  }

  else if(pose.conformation().num_chains() >= 2 && !basic::options::option[fixedchains].active()){
    TR << "WARNING more than two chains present but no -fixedchains declared.  Interface calculations unreliable!" << std::endl;
    moverOP = new  devel::anchored_design::InterfaceAnalyzerMover(
								 basic::options::option[ jumpnum ].value(),
								 basic::options::option[ tracer_data_print ].value(),
								 core::scoring::get_score_function(),
								 basic::options::option[ compute_packstat ].value(),
								 basic::options::option[ pack_input ].value(),
								 true, //pack the separated poses
								 false //use the input pose name instead of jobname
								 );
    return;
  }

  else if(pose.conformation().num_chains() >= 2 && basic::options::option[fixedchains].active()){
    utility::vector1<std::string> fixed_chains_string (basic::options::option[fixedchains].value());
    //parse the fixed chains to figure out pose chain nums
    std::set< int > fixed_chains;
    TR << "Fixed chains are: " ;
    for(core::Size j = 1; j <= fixed_chains_string.size(); ++j){
      char this_chain (fixed_chains_string[ j ][0]);
      for (core::Size i = 1; i<=pose.total_residue(); ++i){
	if (pose.pdb_info()->chain( i ) == this_chain){
	  fixed_chains.insert( pose.chain(i) );
	  break;
	}
      }
      TR << this_chain << ", ";
    }
    TR << "these will be moved together." << std::endl;
    moverOP = new devel::anchored_design::InterfaceAnalyzerMover(
								fixed_chains,
								basic::options::option[ tracer_data_print ].value(),
								core::scoring::get_score_function(),
								basic::options::option[ compute_packstat ].value(),
								basic::options::option[ pack_input ].value()
								);
    return;
  }
  else
    utility_exit_with_message_status( "Something weird with number of chains/options...exiting", 1 );

  return;
} //end assign_IA_mover

///begin apply
void IAMover::apply( core::pose::Pose & pose ) {
	using namespace std;
	using namespace ObjexxFCL::format;
  //get filename and some info
  utility::file::FileName filename = pose.pdb_info()->name();
	string pdbname( filename.base() );


  //check to make sure there are enough chains
  if(pose.conformation().num_chains() < 2){
    string exitname( pdbname + " has only one chain, exiting...");
    utility_exit_with_message_status( exitname, 1 );
  }

  //fill the interface analyzer mover
  assign_IA_mover( IAM_ , pose );
	(*scorefxn_)(pose);

	//now apply and get cool data and stuff
	IAM_->apply(pose);
	//core::Real total_sasa( IAM_->get_total_complexed_sasa() );
	core::Real whole_protein_energy( IAM_->get_complex_energy() );
	core::Real interface_delta_sasa( IAM_->get_interface_delta_sasa() );
	core::Real separated_interface_energy( IAM_->get_separated_interface_energy() );
	core::Real separated_interface_energy_ratio( IAM_->get_separated_interface_energy_ratio() );
	core::Real interface_packstat( IAM_->get_interface_packstat() );
	core::Size interface_delta_hbond_unsat( IAM_->get_interface_delta_hbond_unsat() );
	string pymol_sel_interface( IAM_->get_pymol_sel_interface() );
	string pymol_sel_hbond_unsat( IAM_->get_pymol_sel_hbond_unsat() );
	string pymol_sel_packing( IAM_->get_pymol_sel_packing() );
	core::Real per_res_E( IAM_->get_per_residue_energy() );
	core::Size n_residues( IAM_->get_num_interface_residues() ) ;
	core::Real gly_dG( IAM_->get_gly_interface_energy() );
	core::Real centroid_dG (IAM_->get_centroid_dG() );
	core::Real hbond_sasa (IAM_->get_interface_Hbond_sasa());
	core::Real hb_exposure(IAM_->get_Hbond_exposure_ratio());
	core::Real hbond_E(IAM_->get_total_Hbond_E());

	//add the first line to the output file
	//if( !firstline_printed_ ){
	if ( !utility::file::file_exists( output_name_ ) ){
		output_statsFile_.open( name_star_ );
		output_statsFile_ << A(45, "InputFile" )<< " "
											<< A(14, "TotalScore" )<< " "
											<< A(10, "Energy/Res" )<< " "
											<< A(10, "N_residues" )<< " "
											<< A(10, "dGbind" ) << " "
											<< A(10, "dSASA") << " "
											<< A(10, "dG/dSASA" ) << " "
											<< A(10, "HBond_E" ) << " "
											<< A(15, "HbondE/dG" ) << " "
											<< A(15, "hb_SASA/dSASA") << " "
											<< A(15, "HB_exposure") << " "
											<< A(15, "gly_dG/dSASA" ) << " "
											<< A(20, "centorid_dG/dSASA" ) << " "
											<< A(10, "packstat" )<< " "
											<< A(12, "unsatHbonds" ) << " "
											<< A(12, "dSASA/unsat" ) << " "
											<< A(12, "unsat/dSASA" )<< std::endl;
		//firstline_printed_=true;
	}
	else
		output_statsFile_.open_append( name_star_ );

	//other values and the like
	core::Real sasa_per_unsat( interface_delta_sasa );
	if (interface_delta_hbond_unsat != 0)
		sasa_per_unsat = interface_delta_sasa / interface_delta_hbond_unsat;

	//Critical output included all in one line (packstat and hbond_unsat will equal 0 if not included in the options)
	int const precision(3);
	output_statsFile_ << A(45, pdbname ) << " "
										<< F(14, precision,  whole_protein_energy ) << " "
										<< F(10, precision, per_res_E) << " "
										<< I(10, 1        , n_residues) << " "
										<< F(10, precision,  separated_interface_energy )<< " "
										<< F(10, 1        ,   interface_delta_sasa )<< " "
										<< F(10, precision,  separated_interface_energy_ratio )<< " "
										<< F(10, precision, hbond_E ) << " "
										<< F(15, precision, hbond_E/ separated_interface_energy )
										<< F(15, 5        , hbond_sasa / interface_delta_sasa )<< " "
										<< F(15, 5        , hb_exposure ) << " "
										<< F(15, 5        , gly_dG / interface_delta_sasa ) << " "
										<< F(20, 5        , centroid_dG / interface_delta_sasa ) << " "
										<< F(10, precision,   interface_packstat )<< " "
										<< I(12, interface_delta_hbond_unsat )<< " "
										<< F(12, precision,  sasa_per_unsat )<< " "
										<< F(12, 4,  interface_delta_hbond_unsat / interface_delta_sasa )<< std::endl;


	//Tracer output of pymol selections
	//if you prefer not to create a lot of little files
	if( basic::options::option[write_pymol_selections_to_tracer].value() ) {
		std::string output_filename = "pymol_sel_" + pdbname + ".pml";

		TR << "cmd.hide(\"everything\",\"" << pdbname << "\")" << std::endl;
		TR << "cmd.show(\"cartoon\"   ,\"" << pdbname << "\")" << std::endl;
		TR << "cmd.show(\"lines\"     ,\"" << pdbname << "\")" << std::endl;
		TR << "cmd.hide(\"(" << pdbname << " and hydro)\")" << std::endl; // one can choose to hide hydrogens
		TR << "cmd.remove(\"(" << pdbname << ") and hydro\")" << std::endl; // I like to remove hydrogens
		TR << "util.color_chains(\"(" << pdbname << " and elem c)\")" << std::endl;
		TR << "cmd.center(\"" << pdbname << "\",animate=-1)" << std::endl;

		TR << pymol_sel_packing << std::endl;
		TR << pymol_sel_interface << std::endl;
		TR << pymol_sel_hbond_unsat << std::endl;
	}

	if( basic::options::option[write_pymol_selections_to_file].value() ) { //if you prefer the convenience of being able to immediately run the pymol script
		std::string output_filename = "pymol_sel_" + pdbname + ".pml";
		utility::io::ozstream OutputFile( output_filename );

		OutputFile << "cmd.hide(\"everything\",\"" << pdbname << "\")" << std::endl;
		OutputFile << "cmd.show(\"cartoon\"   ,\"" << pdbname << "\")" << std::endl;
		OutputFile << "cmd.show(\"lines\"     ,\"" << pdbname << "\")" << std::endl;
		OutputFile << "cmd.hide(\"(" << pdbname << " and hydro)\")" << std::endl; // one can choose to hide hydrogens
		OutputFile << "cmd.remove(\"(" << pdbname << ") and hydro\")" << std::endl; // I like to remove hydrogens
		OutputFile << "util.color_chains(\"(" << pdbname << " and elem c)\")" << std::endl;
		OutputFile << "cmd.center(\"" << pdbname << "\",animate=-1)" << std::endl;

		OutputFile << pymol_sel_packing << std::endl;
		OutputFile << pymol_sel_interface << std::endl;
		OutputFile << pymol_sel_hbond_unsat << std::endl;

	}
	//if you want all the usual output from the job distributor
 	if( basic::options::option[ force_output ].value() )
		TR << "Forcing output of PDBs" << std::endl;
	else
		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY); //hacky way to prevent output

	return;

}//end apply



int
main( int argc, char* argv[] )
{

	try {

	using basic::options::option;
	option.add( jumpnum, "jump between chains of interface" ).def(1);
	option.add( compute_packstat, "packstat of interface residues only" ).def(false);
	option.add( write_pymol_selections_to_file, "write pymol selections to file" ).def(false);
	option.add( write_pymol_selections_to_tracer, "write pymol selections to tracer" ).def(false);
	option.add( tracer_data_print, "print pymol selections" ).def(true);
	option.add( force_output, "forces output of scorefile and pdbs" ).def(false);
	option.add( fixedchains, "Which chain(s) is/are moved away from the others." );
	option.add( ia_stats_filename,"Name of the file output with interface stats." ).def("interface_stats.txt");
	option.add( pack_input,"Run pack rots on the input pose?").def(false);

	devel::init(argc, argv);

	protocols::jd2::JobDistributor::get_instance()->go(new IAMover);

	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
