// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/fragments_to_ss
/// @brief  outputs the fragments as if they were a psipred secondary structure prediction. This will be used to calibrate
//          rosetta abinitio
/// @author TJ Brunette
// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <protocols/jumping/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <ObjexxFCL/format.hh>
#include <basic/options/option_macros.hh>

#include <utility/io/ozstream.hh> 


using namespace core;
using utility::vector1;

class ThisApplication  {
public:
	ThisApplication();
	static void register_options();
private:
};

ThisApplication::ThisApplication()
{}
OPT_1GRP_KEY( File, out, ss )

void ThisApplication::register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
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

vector1<vector1 <Real> > gather_ss_pct(){
		using namespace core::chemical;
		using namespace core::import_pose::pose_stream;
		using core::import_pose::pose_from_pdb;
		MetaPoseInputStream input = streams_from_cmd_line();
    ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );
		vector1<vector1 <Size> > counts; //H=0,L=1,E=2
		vector1<vector1 <Real> > pcts;
		bool first_pose=true;
		while(input.has_another_pose()){
				core::pose::PoseOP poseOP;
				poseOP = core::pose::PoseOP( new core::pose::Pose() );
				input.fill_pose(*poseOP,*rsd_set);
				if(first_pose){
						for(Size ii=1; ii<=poseOP->total_residue(); ++ii){
								vector1<Size> tmp_pos;
								for(Size kk=1; kk<=3; ++kk){
										tmp_pos.push_back(0);
								}
								counts.push_back(tmp_pos);
						}
						first_pose=false;
				}
				protocols::jumping::assign_ss_dssp( *poseOP );
				for (Size ii=1; ii <= poseOP->total_residue(); ++ii) {	
						if(poseOP->secstruct(ii) == 'H')
								counts[ii][1]+=1;
						if(poseOP->secstruct(ii) == 'L')
								counts[ii][2]+=1;
						if(poseOP->secstruct(ii) == 'E')
								counts[ii][3]+=1;
				}
		}


		for(Size ii=1; ii<=counts.size();++ii){
				Size total_count = 0;
				for(Size kk=1; kk<=3; ++kk){
						total_count += counts[ii][kk];
				}
				vector1<Real> tmp_pct;
				for(Size kk=1; kk<=3; ++kk){
						tmp_pct.push_back((Real)counts[ii][kk]/(Real)total_count);
				}
				pcts.push_back(tmp_pct);
		}
		return(pcts);
}

int main( int argc, char * argv [] ) {
	try {
		using namespace basic;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
    using namespace ObjexxFCL::format;
		ThisApplication::register_options();
		devel::init(argc, argv);
    utility::io::ozstream out( option[ out::ss ] );
    out << "# PSIPRED VFORMATED but done with Rosetta fragments\n" << std::endl;
		vector1<vector1<Real> > pct = gather_ss_pct();
    for(Size ii=1; ii<=pct.size(); ++ii){
        char label = get_label(pct[ii]);
        out << I(4,ii) << " " << "A" << " " << label <<" " << F(7,3,pct[ii][2]) << F(7,3,pct[ii][1])  << F(7,3,pct[ii][3]) << std::endl;
    }
    out.close();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

