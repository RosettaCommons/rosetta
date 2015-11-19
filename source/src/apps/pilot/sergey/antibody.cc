// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file antibody.cc
///
/// @brief C++ port of Python antibody script
/// @author Sergey Lyskov

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/antibody/grafting/antibody_sequence.hh>
#include <protocols/antibody/grafting/regex_based_cdr_detection.hh>
#include <protocols/antibody/grafting/chothia_numberer.hh>
#include <protocols/antibody/grafting/scs_blast.hh>
#include <protocols/antibody/grafting/exception.hh>
#include <protocols/antibody/grafting/scs_functor.hh>


// Grafting util function related includes
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/simple_moves/MutateResidue.hh>


#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

//#include <utility/CSI_Sequence.fwd.hh>

#include <devel/init.hh>

#include <cassert>
#include <cstdlib>
#include <iostream>

#include <string>
#include <regex>

//#include <typeinfo>

using std::string;

static THREAD_LOCAL basic::Tracer TR("antibody");

// Command-line options relevant only to this application
OPT_KEY( String, heavy )
OPT_KEY( String, light )

void register_options()
{
	NEW_OPT( heavy, "Sequence of antibody heavy chain", "" );
	NEW_OPT( light, "Sequence of antibody light chain", "" );
}


int antibody_main()
{
	using namespace utility;
	using namespace protocols::antibody::grafting;
	using protocols::antibody::grafting::uint;

	using std::string;

	string heavy_chain_sequence, light_chain_sequence;

	if ( basic::options::option[ basic::options::OptionKeys::heavy ].user() ) {
		heavy_chain_sequence = basic::options::option[ basic::options::OptionKeys::heavy ]();
	}

	if ( basic::options::option[ basic::options::OptionKeys::light ].user() ) {
		light_chain_sequence = basic::options::option[ basic::options::OptionKeys::light ]();
	}

	if( !heavy_chain_sequence.size()  and  !light_chain_sequence.size() ) {
		// 2BRR, default values for testing
		heavy_chain_sequence = "AVQLEQSGPELKKPGETVKISCKASGYTFTNYGMNWVKQAPGKGLKWMGWINTYTGEPTYADDFKERFAFSLETSASAAYLQINNLKNEDTATYFCARDYYGSTYPYYAMDYWGQGTTVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPRGG";
		light_chain_sequence = "ENVLTQSPAIMSASPGEKVTMTCRASSSVSSSYLHWYQQKSGASPKLWIYSTSNLASGVPARFSGSGSGTSYSLTISSVEAEDAATYYCQQYSGYPYTFGGGTKLEIKRADAAPTVSIFPPSSEQLTSGGASVVCFLNNFYPKDINVKWKIDGSERQNGVLNSWTDQDSKDSTYSMSSTLTLTKDEYERHNSYTCEATHKTSTSPIVKSFNRNEC";
	}

	{
		TR << TR.Red << "Antibody grafting protocol" << TR.Reset << " Running with input:" << std::endl;
		TR << TR.Green << "Heavy chain sequence: " << TR.Bold << heavy_chain_sequence << TR.Reset << std::endl;
		TR << TR.Blue  << "Light chain sequence: " << TR.Bold << light_chain_sequence << TR.Reset << std::endl;
	}

	AntibodySequence as(heavy_chain_sequence, light_chain_sequence);

	RegEx_based_CDR_Detector().detect(as);
	std::cout << std::endl;

	AntibodyFramework heavy_fr = as.heavy_framework();
	AntibodyFramework light_fr = as.light_framework();

	std::cout << "Heavy" << heavy_fr << std::endl;
	std::cout << "Light" << light_fr << std::endl;
	std::cout << std::endl;

	try {
		// AntibodyNumbering an( Chothia_Numberer().number(as, heavy_fr, light_fr) );
		// std::cout << "Heavy Numbering: " << an.heavy << std::endl << std::endl;
		// std::cout << "Light Numbering: " << an.heavy << std::endl << std::endl;

		SCS_BlastPlus scs;

		scs.add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_sequence_length) );
		scs.add_filter( SCS_FunctorCOP(new SCS_BlastFilter_by_alignment_length) );

		scs.set_sorter( SCS_FunctorCOP(new SCS_BlastComparator_BitScore_Resolution) );

		SCS_ResultsOP scs_results { scs.select(as) };

		if( scs_results->frh.size()  and  scs_results->frl.size()  and  scs_results->orientation.size() ) {
			TR << "SCS frh best template: " << TR.Bold << TR.Blue  << scs_results->frh[0]->pdb << TR.Reset << std::endl;
			TR << "SCS frl best template: " << TR.Bold << TR.Green << scs_results->frl[0]->pdb << TR.Reset << std::endl;
			TR << "SCS Orientation best template: " << TR.Bold << TR.Yellow << scs_results->orientation[0]->pdb << TR.Reset << std::endl;

			string frh_pdb_name = "./../../tools/antibody/antibody_database/pdb" + scs_results->frh[0]->pdb + "_chothia.pdb";
			string frl_pdb_name = "./../../tools/antibody/antibody_database/pdb" + scs_results->frh[0]->pdb + "_chothia.pdb";
			string orientation_pdb_name = "./../../tools/antibody/antibody_database/pdb" + scs_results->orientation[0]->pdb + "_chothia.pdb";

			core::pose::PoseOP frh = core::import_pose::pose_from_pdb( frh_pdb_name );
			core::pose::PoseOP frl = core::import_pose::pose_from_pdb( frh_pdb_name );

			frh->dump_pdb("../_ab_/frh.pdb");
			frl->dump_pdb("../_ab_/frl.pdb");

			AntibodyFramework trimmed_heavy_fr = as.heavy_framework();
			AntibodyFramework trimmed_light_fr = as.light_framework();

			trim_framework(as, trimmed_heavy_fr, trimmed_light_fr);

			AntibodyNumbering an( Chothia_Numberer().number(as, trimmed_heavy_fr, trimmed_light_fr) );

			//std::map<char, char> opposite_chain { {'L', 'H'},  {'H', 'L'} };

			struct {
				char chain;
				core::pose::PoseOP &pose;
				//AntibodyChainNumbering::NumberingVector numbering;
				AntibodyChainNumbering numbering;
				string trimmed_sequence;
				string color;
			} J[] {
				{'H', frh, an.heavy, trimmed_heavy_fr.fr1 + as.h1_sequence() + trimmed_heavy_fr.fr2 + as.h2_sequence() + trimmed_heavy_fr.fr3 + as.h3_sequence() + trimmed_heavy_fr.fr4, TR.Blue},
				{'L', frl, an.light, trimmed_light_fr.fr1 + as.l1_sequence() + trimmed_light_fr.fr2 + as.l2_sequence() + trimmed_light_fr.fr3 + as.l3_sequence() + trimmed_light_fr.fr4, TR.Green},
			};


			for(auto &j : J) {
				TR << "Adjusting fr" << char(std::tolower(j.chain)) << " template sequence [" << j.pose->pdb_info()->name() << "]..." << std::endl;

				AntibodyChainNumbering::NumberingVector numbering = j.numbering.all();

				TR << "By using numbering: " << utility::join( numbering, " ") << std::endl;

				int chain = -1;
				// Finding right chain in SCS template and delete eveything else
				for(uint i=1; i<=j.pose->total_residue(); ++i) {
					if( j.pose->pdb_info()->chain(i) == j.chain ) { chain=j.pose->chain(i);  break; }
				}

				TR.Debug << "Sequence before (all chains): " << j.pose->sequence() << std::endl;

				if( chain >= 0 ) j.pose = j.pose->split_by_chain(chain);
				else throw _AE_grafting_failed_( string("Could not find chain: ") + j.chain + " in template " + j.pose->pdb_info()->name() );

				// for(uint i=1; i<=j.pose->total_residue(); ++i) {
				// 	TR << "New pose pdb info for res " << i << ":" << j.pose->pdb_info()->pose2pdb(i) << " i:" << j.pose->pdb_info()->icode(i) << std::endl;
				// }

				TR << "Sequence before: " << j.color << j.pose->sequence() << TR.Reset << std::endl;

				for(uint i=j.pose->total_residue()-1; i>=1; --i) {
					string pdb_res_n = std::to_string(j.pose->pdb_info()->number(i)) + ( j.pose->pdb_info()->icode(i) == ' ' ?  "" : string(1, j.pose->pdb_info()->icode(i) ) );
					auto np = std::find(numbering.begin(), numbering.end(), pdb_res_n);

					if( np < numbering.end() ) {
						char aa = j.trimmed_sequence[ np-numbering.begin() ];
						TR.Trace << "Replacing pose residue " << i << " with " << aa << std::endl;
						protocols::simple_moves::MutateResidue(i, aa).apply( *j.pose.get() );

						j.trimmed_sequence.erase( j.trimmed_sequence.begin() + ( np - numbering.begin() ) );
						numbering.erase( np );
					}
				}
				if( numbering.size() ) {
					TR << TR.Red << "WARNING: Was not able to adjust all residue in chain: " << j.chain << "!!!" << TR.Reset << std::endl;
					TR << TR.Red << "Leftovers numbering: " << numbering  << "  AA: " << j.trimmed_sequence << TR.Reset << std::endl;
					TR << TR.Red << TR.Underline << "Original numbering was:" << TR.Reset << ' ' << j.numbering << std::endl;
				}

				TR << "Sequence after:  " << TR.Bold << j.color << j.pose->sequence() << TR.Reset << std::endl;
			}

			frh->dump_pdb("../_ab_/frh_after_grafting.pdb");
			frl->dump_pdb("../_ab_/frl_after_grafting.pdb");

			frh->append_pose_by_jump(*frl, 1);
			frh->pdb_info()->obsolete(false);
			//frh->update_pose_chains_from_pdb_chains();

			frh->dump_pdb("../_ab_/frh_frl.pdb");

			// for(uint i=1; i<=frh->total_residue(); ++i) {
			// 	TR << "New pose pdb info for res " << i << ":" << frh->pdb_info()->pose2pdb(i) << " i:" << frh->pdb_info()->icode(i) << std::endl;
			// }



			core::pose::PoseOP orientation = core::import_pose::pose_from_pdb( orientation_pdb_name );
			orientation->dump_pdb("../_ab_/orientation.pdb");

			//TR << "Pose: " << std::endl;
			//TR <<  *orientation.get() << std::endl;


		}

	}
	catch(Grafting_Base_Exception e) {
		std::cout << e << std::endl;
	}

	// CDR1_Sequence cdr1 = get_cdr<AntibodyLightChain, CDR1_Sequence>()(as);
	// std::cout << "AntibodyLightChain, CDR1_Sequence: (" << cdr1.start << ", " << cdr1.stop << ")" << " valid:" << cdr1.valid(light_chain_sequence) << std::endl;


	// std::cout << CDR_Bounds() << std::endl;

	return 0;
}

int main(int argc, char * argv [])
{
	try {
		register_options();
		devel::init(argc, argv);
		return antibody_main();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return 1;
	}
}


#else // __ANTIBODY_GRAFTING__

#include <iostream>

int main( int /* argc */, char * /* argv */ [] )
{
	std::cerr << "Antibody protocol require to be build with full C++11 support! Please use GCC-4.9+ or Clang-3.6+ and add extras=cxx11 to your scons build command line.";
	return 1;
}

#endif // __ANTIBODY_GRAFTING__
