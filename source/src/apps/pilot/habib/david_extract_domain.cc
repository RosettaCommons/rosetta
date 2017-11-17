// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk + dj

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sstream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <basic/MetricValue.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/rms_util.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <numeric/random/random.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>


// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.david_recompute_score_and_rmsd.main" );


OPT_KEY( Integer, startCut )
OPT_KEY( Integer, endCut )
OPT_KEY( Integer, domainNum )
OPT_KEY( String, chain )
OPT_KEY( String, input_complex )

//set to store pdb info keys
std::vector<std::string> surface;
std::set <std::string> interface;


/// General testing code
int
main( int argc, char * argv [] )
{
	try {

		std::vector<std::string> surface;
		std::set <std::string> interface;
		NEW_OPT ( startCut, "Starting position of domain to be extracted",0);
		NEW_OPT ( endCut, "Ending position of domain to be extracted",0);
		NEW_OPT ( domainNum, "Number to append to output_tag for domain name",0);
		NEW_OPT ( chain, "Chain to extract the domain from","");
		NEW_OPT ( input_complex, "File name specifying the protein-RNA complex (the other PDB only has the specified chain","");

		using namespace core;
		using namespace core::scoring;

		devel::init(argc, argv);
		pose::Pose whole_pose;
		pose::Pose whole_pose_complex;
		pose::Pose chain_pose;
		std::string const input_pdb_name ( basic::options::start_file() );
		core::import_pose::pose_from_file( chain_pose, input_pdb_name );
		std::string const complex_pdb_name ( option[input_complex] );
		core::import_pose::pose_from_file( whole_pose, complex_pdb_name );
		core::import_pose::pose_from_file( whole_pose_complex, complex_pdb_name );

		core::Size const start = (core::Size)option[ startCut ];
		core::Size const end = (core::Size)option[ endCut ];
		core::Size const domain = (core::Size)option[ domainNum ];
		std::string const chain_string = option[ chain ];
		if ( start <= 0 ) {
			std::cerr<< "Error: startCut must be greater than 0"<<std::endl;
			return -100;
		}
		if ( start >= end ) {
			std::cerr<< "Error: endCut must be greater than startCut"<<std::endl;
			return -101;
		}
		if ( chain_string.length() != 1 ) {
			std::cerr<<"Error: invalid chain specified"<<std::endl;
			return -102;
		}
		char chain = chain_string.c_str()[0];
		Size chain_count=0;
		for ( Size i = 1; i <= chain_pose.size(); ++i ) {
			if ( chain_pose.residue(i).is_protein() ) {
				chain_count++;
			}
		}
		if ( start>= chain_count ) {
			std::cerr<< "Error: startCut must be less than the total residues in the chain"<<std::endl;
			return -103;
		}
		if ( end > chain_count ) {
			std::cerr<< "Error: endCut must be less than the total residues in the chain"<<std::endl;
			return -104;
		}

		bool found = false;
		for ( Size i = 1; i <= whole_pose.size(); ++i ) {
			std::cout<<i<<":"<<whole_pose.size()<<" ";
			if ( !found && whole_pose.pdb_info()->chain(i) == chain ) {
				found=true;
				whole_pose_complex.conformation().delete_residue_range_slow(i+start-1, i+end-1);
				std::cout<<i<<"\n";
				for ( Size j=1; j<start; j++ ) {
					std::cout<<j<<"\n";
					whole_pose.conformation().delete_residue_slow(i);
				}
				std::cout<<i<<"\n";
				i+=(end-start);
				std::cout<<i<<"\n";
				continue;
			}
			//if (whole_pose.pdb_info()->chain(i) != chain && !whole_pose.residue(i).is_RNA()){
			if ( !whole_pose.residue(i).is_RNA() ) {
				whole_pose.conformation().delete_residue_slow(i);
				i--;
			}
		}
		std::cout<<"printing stuff\n";
		std::stringstream domfilename;
		domfilename<<option[ OptionKeys::out::output_tag ]()<<"_D"<<domain<<".pdb";
		whole_pose.dump_pdb(domfilename.str());

		std::cout<<"deleting stuff\n";
		if ( end - start+1 == chain_pose.size() ) {
			chain_pose.clear();
		} else {
			chain_pose.conformation().delete_residue_range_slow(start, end);
		}

		std::cout<<"printing stuff\n";
		std::stringstream leftfilename;
		leftfilename<<option[ OptionKeys::out::output_tag ]()<<"_leftover.pdb";
		chain_pose.dump_pdb(leftfilename.str());

		std::stringstream leftcompfilename;
		leftcompfilename<<option[ OptionKeys::out::output_tag ]()<<"_leftover_complex.pdb";
		std::cout<<"printing stuff\n";
		whole_pose_complex.dump_pdb(leftcompfilename.str());

		/*    while (ifs.good()){
		std::string intres;
		ifs >> intres;
		if (intres.length() >0){
		if (intres[0]!= '>'){
		seq.append(intres);
		}
		}
		}
		}

		if (seq.length()==0){
		std::cout<< "Fasta file contains no sequence: "<<ffilename<<std::endl;
		return -101;

		}

		int last = -100;
		core::chemical::ResidueTypeSet const & rsd_set( pose.residue(1).residue_type_set() );

		for( Size i = 1; i <= pose.size(); ++i){
		if (last == -100){
		last = pose.pdb_info()->number(i);
		std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
		continue;
		}
		if (last > pose.pdb_info()->number(i)+1){
		for (int j = last-1; j > pose.pdb_info()->number(i); j--){
		std::cout<<seq[j-1]<<"\n";
		core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueSelector().set_name1( seq[j-1] ).exclude_variants().select( rsd_set )[1] );
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
		pose.append_polymer_residue_after_seqpos(*new_rsd, i, false );
		pose.pdb_info()->number(i+2, pose.pdb_info()->number(i)-1);
		pose.pdb_info()->chain(i+2, pose.pdb_info()->chain(i));
		}

		}else if (last < pose.pdb_info()->number(i)-1){
		for (int j = last+1; j < pose.pdb_info()->number(i); j++){
		std::cout<<seq[j-1]<<"\n";
		core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueSelector().set_name1( seq[j-1] ).exclude_variants().select( rsd_set )[1] );
		core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );
		pose.append_polymer_residue_after_seqpos(*new_rsd, i-1, false );
		pose.pdb_info()->number(i, pose.pdb_info()->number(i-1)+1);
		pose.pdb_info()->chain(i, pose.pdb_info()->chain(i-1));
		i++;
		}


		}
		std::cout<<i<<" "<< pose.pdb_info()->number(i)<<" "<<pose.residue(i).name1()<<"\n" ;
		last = pose.pdb_info()->number(i);
		}

		for( Size i = 1; i <= pose.size(); ++i){
		std::cout<<pose.pdb_info()->number(i)<<"\n";
		}
		for( int i = 1; i <pose.pdb_info()->number(1); ++i){
		std::cout<<" ";
		}
		for( Size i = 1; i <= pose.size(); ++i){
		std::cout<<pose.residue(i).name1();
		}
		std::cout<<"\n"<<seq<<"\n";

		pose.pdb_info()->obsolete(false);
		pose.dump_pdb("test.pdb");
		*/
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}



