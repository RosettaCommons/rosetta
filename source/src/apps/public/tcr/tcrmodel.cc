// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/tcr/tcrmodel.cc
///
/// @brief application for Tcr (T cell receptor) modelling
/// @author Ragul Gowthaman


#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <string>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRmodel.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

using std::string;
using namespace protocols::tcr;
static basic::Tracer TR("apps.public.tcr.tcrmodel");

OPT_KEY( String, alpha )
OPT_KEY( String, beta )

int main(int argc, char * argv [])
{
  try {

		NEW_OPT( alpha, "Sequence of Tcr alpha chain", "" );
		NEW_OPT( beta, "Sequence of Tcr beta chain", "" );

		devel::init(argc, argv);



		string alpha_chain_sequence, beta_chain_sequence;

		if ( basic::options::option[ basic::options::OptionKeys::alpha ].user() ) {
			alpha_chain_sequence = basic::options::option[ basic::options::OptionKeys::alpha ]();
		}

		if ( basic::options::option[ basic::options::OptionKeys::beta ].user() ) {
			beta_chain_sequence = basic::options::option[ basic::options::OptionKeys::beta ]();
		}
		if( !alpha_chain_sequence.size() and !beta_chain_sequence.size() ) {
			utility_exit_with_message("Error reading input sequences.");
		}

		TCRseqInfoOP seqinfo( new TCRseqInfo(alpha_chain_sequence, beta_chain_sequence) );
		TCRmodel tcrinfo(seqinfo);

		//Alpha
		TR <<"Tcr Alpha truncated Domain sequence: "<<seqinfo->atcr().truncdomain<<std::endl;
		TR << "Tcr Alpha germline sequence : "<<seqinfo->atcr().gm<<std::endl;
		TR << "Tcr Alpha Framework sequence : "<<seqinfo->atcr().fr<<std::endl;
		TR << "Tcr Alpha Domain CDR1 sequence: "<<seqinfo->atcr().cdr1<<std::endl;
		TR << "Tcr Alpha Domain CDR2 sequence: "<<seqinfo->atcr().cdr2<<std::endl;
		TR << "Tcr Alpha Domain CDR3 sequence: "<<seqinfo->atcr().cdr3<<std::endl;

		if (tcrinfo.use_gma_templates()) {
			TR << "Alpha germline template : "<< tcrinfo.atmplt().gm.tid<< std::endl;
		} else {
			TR << "Alpha Framework template : "<< tcrinfo.atmplt().fr.tid<< std::endl;
			TR << "Alpha CDR1 template : "<<tcrinfo.atmplt().cdr1.tid<<std::endl;
			TR << "Alpha CDR2 template : "<<tcrinfo.atmplt().cdr2hv4.tid<<std::endl;
		}
		TR << "Alpha CDR3 template : "<<tcrinfo.atmplt().cdr3.tid<<std::endl;

		//Beta
		TR << "Tcr Beta truncated Domain : "<<seqinfo->btcr().truncdomain<<std::endl;
		TR << "Tcr Beta germline sequence : "<<seqinfo->btcr().gm<<std::endl;
		TR << "Tcr Beta Framework sequence : "<<seqinfo->btcr().fr<<std::endl;
		TR << "Tcr Beta Domain CDR1 sequence: "<<seqinfo->btcr().cdr1<<std::endl;
		TR << "Tcr Beta Domain CDR2 sequence: "<<seqinfo->btcr().cdr2<<std::endl;
		TR << "Tcr Beta Domain CDR3 sequence: "<<seqinfo->btcr().cdr3<<std::endl;

		if (tcrinfo.use_gmb_templates()) {
			TR << "Beta germline template : "<< tcrinfo.btmplt().gm.tid<< std::endl;
		} else {
			TR << "Beta Framework template : "<< tcrinfo.btmplt().fr.tid<< std::endl;
			TR << "Beta CDR1 template : "<<tcrinfo.btmplt().cdr1.tid<<std::endl;
			TR << "Beta CDR2 template : "<<tcrinfo.btmplt().cdr2hv4.tid<<std::endl;
		}
		TR << "Beta CDR3 template : "<<tcrinfo.atmplt().cdr3.tid<<std::endl;
		//orientation template
		TR << "Alpha Beta orientation template : "<< tcrinfo.atmplt().ori.tid << " " << tcrinfo.btmplt().ori.tid << std::endl;
		std::string const prefix = basic::options::option[basic::options::OptionKeys::out::prefix];
		if ( !(tcrinfo.tcr_graft_model()->empty()) ) {
			std::string out_tag_tcrgraftmodel = prefix + "tcrgraftmodel.pdb";
			TR <<"Output grafted model:"<<"\t"<<out_tag_tcrgraftmodel<<std::endl;
			tcrinfo.tcr_graft_model()->dump_pdb(out_tag_tcrgraftmodel);
		}
		if ( !(tcrinfo.tcr_loop_model()->empty()) ) {
			std::string out_tag_tcrloopmodel = prefix + "tcrloopmodel.pdb";
			TR <<"Output loop model:"<<"\t"<<out_tag_tcrloopmodel<<std::endl;
			tcrinfo.tcr_loop_model()->dump_pdb(out_tag_tcrloopmodel);
		}
		if ( !(tcrinfo.tcr_model()->empty()) ) {
			std::string out_tag_tcrmodel = prefix + "tcrmodel.pdb";
			TR <<"Output final model:"<<"\t"<<out_tag_tcrmodel<<std::endl;
			tcrinfo.tcr_model()->dump_pdb(out_tag_tcrmodel);
		}
	} catch (utility::excn::Exception const & e ) {
    e.display();
    return 1;
  }
}

#else // __ANTIBODY_GRAFTING__

#include <iostream>

int main( int /* argc */, char * /* argv */ [] )
{
	std::cerr << "Antibody/TCRmodel protocol require to be build with full C++11 support! Please use GCC-4.9+ or Clang-3.6+ and add extras=cxx11 to your scons build command line." <<std::endl;
	return 1;
}

#endif // __ANTIBODY_GRAFTING__
