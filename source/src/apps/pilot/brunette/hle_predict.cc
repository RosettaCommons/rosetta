// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/brunette/minimalCstRelax.cc
///
/// @brief  For the first structure in the alignment compare the deviation between native and the relaxed native structure. Then constraints are modified and coordinate constraints are added. If an alignment is given only the residues from the beginning to end of the alignment are used.

/// @usage: -in:file:s <pdb files> [options: -minimalCstRelax:coordinate_cst_gap in first round of relax gap between CA coordinate constraints] -in::file::alignmen [alignment.filt file]
/// @author TJ Brunette

#include <devel/init.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <protocols/ss_prediction/SS_predictor.hh>

#include <ctime>
#include <iostream>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/format.hh>

using namespace ObjexxFCL::format;
using std::string;
using utility::vector1;
using core::Real;
using core::Size;
using core::SSize;

static thread_local basic::Tracer tr( "hle_predict" );

char get_label(vector1 <Real> ss_pred_pos){
	char label = 'X';
	if((ss_pred_pos[1] >= ss_pred_pos[2]) && (ss_pred_pos[1] >= ss_pred_pos[3]))
		label = 'H';
	else
		if((ss_pred_pos[2] >= ss_pred_pos[1]) && (ss_pred_pos[2] >= ss_pred_pos[3]))
			label = 'C';
		else
			if((ss_pred_pos[3] >= ss_pred_pos[1]) && (ss_pred_pos[3] >= ss_pred_pos[2]))
				label = 'E';
	return(label);
}

void predict_ss(vector1< vector1 <Real> > ss_pred, string fasta, string id){
	std::string out_file_name_str( id + ".che_ss2");
	utility::io::ozstream output(out_file_name_str);
	output << "# PSIPRED VFORMATED but done with Rosetta v1.0\n" << std::endl;
	tr << "NOTE HLE produces CHE output identical to psipred" << std::endl;
	for(Size ii = 1; ii<= ss_pred.size(); ++ii){
		char label = get_label(ss_pred[ii]);
		output << I(4,ii) << " " << fasta.at(ii-1) << " " << label <<" " << F(7,3,ss_pred[ii][2]) << F(7,3,ss_pred[ii][1])  << F(7,3,ss_pred[ii][3]) << std::endl;
 	}
	output.close();
}


int main( int argc, char * argv [] ) {
	try {

	using namespace protocols::ss_prediction;
	using namespace basic;
	using namespace core::sequence;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	devel::init(argc, argv);
    typedef vector1< SequenceOP > seqlist;
    SS_predictor *hle_predictor = new SS_predictor("HLE");
	for(Size ii = 1; ii<= option[in::file::fasta]().size(); ++ii){
		seqlist query_seq = read_fasta_file(option[in::file::fasta ]()[ii]);
		string fasta = query_seq[1]->sequence();
		string id = option[in::file::fasta ]()[ii];
		vector1< vector1 <Real> > ss_pred = hle_predictor->predict_ss(fasta);
		predict_ss(ss_pred,fasta,id);
    }
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

