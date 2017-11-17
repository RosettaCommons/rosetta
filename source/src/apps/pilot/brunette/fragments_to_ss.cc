// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/fragments_to_ss
/// @brief  outputs the fragments as if they were a psipred secondary structure prediction. This will be used to calibrate
//          rosetta abinitio
/// @author TJ Brunette

#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FragData.hh>

#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <ObjexxFCL/format.hh>

#include <utility/io/ozstream.hh>

class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}

OPT_KEY( File, fragment_file )
OPT_KEY( Integer, fragment_depth )
OPT_1GRP_KEY( File, out, ss )

using namespace core;
using utility::vector1;

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( fragment_file, "fragment file", "" );
    NEW_OPT( fragment_depth, "number of top fragments", 25);
    NEW_OPT( out::ss, "secondary structure out", "");
}

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

int main( int argc, char * argv [] ) {
	try {
	using namespace basic;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
    using namespace fragment;
    using namespace ObjexxFCL::format;
    ThisApplication::register_options();
	devel::init(argc, argv);
    utility::io::ozstream out( option[ out::ss ] );
    out << "# PSIPRED VFORMATED but done with Rosetta fragments\n" << std::endl;
    FragSetOP orig_frags;
    orig_frags = FragmentIO().read_data( option[ fragment_file ]() );
    vector1<Size> per_pos_count;
    vector1<Size> per_pos_count_coil;
    vector1<Size> per_pos_count_helix;
    vector1<Size> per_pos_count_sheet;
    for(Size ii =1; ii <= orig_frags->max_pos(); ++ii){
        per_pos_count.push_back(0);
        per_pos_count_coil.push_back(0);
        per_pos_count_helix.push_back(0);
        per_pos_count_sheet.push_back(0);
    }
    for ( ConstFrameIterator frame = orig_frags->begin(), eframe=orig_frags->end(); frame != eframe; ++frame ) {
        assert( frame->length() == 9);
        for ( Size ii=1; ii<=frame->nr_frags()&& ii<=(Size)option[ fragment_depth ]() ; ++ii ) {
            for( Size kk=0; kk< frame->fragment_ptr(ii)->secstruct().size(); ++kk){
                per_pos_count[frame->start()+kk]+=1;
                if(frame->fragment_ptr(ii)->secstruct().at(kk) == 'L'){
                    per_pos_count_coil[frame->start()+kk]+=1;
                }
                if(frame->fragment_ptr(ii)->secstruct().at(kk) == 'H'){
                    per_pos_count_helix[frame->start()+kk]+=1;
                }
                if(frame->fragment_ptr(ii)->secstruct().at(kk) == 'E'){
                    per_pos_count_sheet[frame->start()+kk]+=1;
                }

            }
        }
    }
    for(Size ii=1; ii<=per_pos_count.size(); ++ii){
        vector1 <Real> pred;
        pred.push_back(((Real)per_pos_count_helix[ii])/((Real)per_pos_count[ii]));
        pred.push_back(((Real)per_pos_count_coil[ii])/((Real)per_pos_count[ii]));
        pred.push_back(((Real)per_pos_count_sheet[ii])/((Real)per_pos_count[ii]));
        char label = get_label(pred);
        out << I(4,ii) << " " << "A" << " " << label <<" " << F(7,3,pred[2]) << F(7,3,pred[1])  << F(7,3,pred[3]) << std::endl;
    }
    out.close();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

