// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


#include <protocols/jd2/JobDistributor.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DiagnosticData.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/ScoreMap.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/assembly.OptionKeys.gen.hh>
#include <string>
/// ObjexxFCL headers
#include <numeric/random/random.hh>
#include <ObjexxFCL//string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace ObjexxFCL::format;

using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
using namespace protocols;
using namespace protocols::moves;
using namespace core::scoring;

using utility::vector1;

static THREAD_LOCAL basic::Tracer TR( "ediff" );


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
    try {
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace scoring;

	devel::init(argc, argv);

	std::string pdb1_filename =  option[ assembly::pdb1 ]();
	std::string pdb2_filename =  option[ assembly::pdb2 ]();

	core::chemical::ResidueTypeSetCAP rsd_set_full = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueTypeSetCAP rsd_set_cen = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );

	core::pose::Pose pose1;
	core::pose::Pose pose2;

	core::import_pose::pose_from_file( pose1, *rsd_set_full, pdb1_filename , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( pose2, *rsd_set_full, pdb2_filename , core::import_pose::PDB_file);

	runtime_assert( pose1.size() == pose2.size() );

	ScoreFunctionOP score_function_( get_score_function() );

	(*score_function_)(pose1);
	(*score_function_)(pose2);

	score_function_->show( TR, pose1);
	score_function_->show( TR, pose2);
	EnergyMap totals1 = pose1.energies().total_energies();
	EnergyMap totals2 = pose2.energies().total_energies();

	for ( EnergyMap::const_iterator it=totals1.begin(),
				it_end = totals1.end(); it != it_end; ++it ) {
			ScoreType const scoretype = ScoreType( it - totals1.begin() + 1 ); // hacky
			if( 	score_function_->get_weight( scoretype ) != 0 )
			TR << LJ(24,ScoreType( scoretype ))  << " " << (totals1[ scoretype ] - totals2[ scoretype ])* score_function_->get_weight( scoretype ) << std::endl;
	}

	EnergyMap total;
	core::Real total_energy2 = 0.0;
	for( core::Size ir=1; ir <= pose1.size(); ++ ir ){
		EnergyMap emap1 = pose1.energies().residue_total_energies(ir);
		EnergyMap emap2 = pose2.energies().residue_total_energies(ir);

		TR << "Res: " << ir << "  ";
		core::Real sum=0;
		for ( EnergyMap::const_iterator it=emap1.begin(),
				it_end = emap1.end(); it != it_end; ++it ) {
			ScoreType const scoretype = ScoreType( it - emap1.begin() + 1 ); // hacky

			core::Real diff = (emap1[ scoretype ] - emap2[ scoretype ]) * score_function_->get_weight( scoretype );

			if( score_function_->get_weight( scoretype ) > 0 ){
				TR << F(5,1,diff) << " ";
				total[ scoretype ] += diff;
				sum+=diff;
			}

			if( scoretype == fa_dun ){
				if( ir <= pose1.pdb_info()->nres() ){
					for( core::Size iat = 1; iat <= pose1.pdb_info()->natoms(ir); iat++ ) pose1.pdb_info()->temperature(  ir, iat, 5*diff+50.0 );
				}
				if( ir <= pose2.pdb_info()->nres() ){
					for( core::Size iat = 1; iat <= pose2.pdb_info()->natoms(ir); iat++ ) pose2.pdb_info()->temperature(  ir, iat, -5*diff+50.0 );
				}
			}

		}
		TR << "| " <<  F(7,1,sum);
		total_energy2 += sum;
		TR << std::endl;
	}
	TR << "--------------------------------" << std::endl;

	core::Real total_sum1 = 0;
	for ( EnergyMap::const_iterator it=total.begin(),
	        it_end = total.end(); it != it_end; ++it ) {

		ScoreType const scoretype = ScoreType( it - total.begin() + 1 );
		if( score_function_->get_weight( scoretype ) > 0 ){
		    TR << F(5,1,total[ scoretype ]) << " ";
				total_sum1 += total[ scoretype ];
		}
	}
	TR << "| " << F(7,1,total_sum1) << std::endl;
	TR << "| " << F(7,1,total_energy2) << std::endl;


	TR <<  std::endl;

	pose1.dump_pdb("pdb1.pdb");
	pose2.dump_pdb("pdb2.pdb");

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
