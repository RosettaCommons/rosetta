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


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/chemical/rna/util.hh>
#include <devel/init.hh>
#include <protocols/viewer/viewers.hh>
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/id/TorsionID.hh>
#include <utility/vector1.hh>

#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_NucleosideStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SuiteStepWiseSampler.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_StepWiseSamplerGeneratorWrapper.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>

//////////////////////////////////////////////////////////


// C++ headers
#include <string>

using namespace core;
using namespace core::chemical::rna;
using namespace protocols::stepwise::sampler::rna;
using namespace protocols::stepwise::sampling::rna;

void test() {
	using namespace core::pose;
	using core::id::TorsionID;

	chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = chemical::ChemicalManager::get_instance()->
	          residue_type_set( chemical::RNA );
	Pose pose;
	make_pose_from_sequence( pose, "aaa", *rsd_set );

	//A-form torsion
	chemical::rna::RNA_FittedTorsionInfo const torsion_info;
	for (Size i = 1; i <= pose.size(); ++i) {
		pose.set_torsion( TorsionID( i, id::BB, 1 ), torsion_info.alpha_aform() );
		pose.set_torsion( TorsionID( i, id::BB, 2 ), torsion_info.beta_aform() );
		pose.set_torsion( TorsionID( i, id::BB, 3 ), torsion_info.gamma_aform() );
		pose.set_torsion( TorsionID( i, id::BB, 5 ), torsion_info.epsilon_aform() );
		pose.set_torsion( TorsionID( i, id::BB, 6 ), torsion_info.zeta_aform() );
	}


/*
	RNA_ChiStepWiseSampler chi_rotamer(2, WHATEVER, NORTH);
	RNA_SugarStepWiseSampler sugar_rotamer(2, WHATEVER);

	//chi_rotamer.set_extra_chi(true);
	chi_rotamer.init();

	sugar_rotamer.init();

	while (chi_rotamer.has_more()) {
		++chi_rotamer;
		std::cout << chi_rotamer.value() << std::endl;
	}

	while (sugar_rotamer.has_more()) {
		++sugar_rotamer;
		std::cout << sugar_rotamer.pucker() << std::endl;
	}

	chi_rotamer.apply(pose);
	sugar_rotamer.apply(pose);
*/
	utility::vector1< core::Size > const suite_list( 1, 2 );

	Pose pose_copy = pose;
	StepWiseRNA_StepWiseSamplerGeneratorWrapper rotamer_generator(
			pose_copy, suite_list, false, true );
	rotamer_generator.set_include_syn_chi( true );
	rotamer_generator.set_allow_syn_pyrimidine( true );
	rotamer_generator.initialize_rotamer_generator_list();


	RNA_NucleosideStepWiseSampler n_rotamer( 3, WHATEVER, WHATEVER );
	n_rotamer.set_idealize_coord( false );
	n_rotamer.set_skip_same_pucker( false );
	n_rotamer.init();

	RNA_SuiteStepWiseSampler suite_rotamer( 2, NORTH, WHATEVER, WHATEVER, WHATEVER );
	suite_rotamer.set_sample_nucleoside_lower( false );
	suite_rotamer.set_sample_nucleoside_upper( true );
	suite_rotamer.set_idealize_coord( false );
	suite_rotamer.set_skip_same_pucker( false );
	suite_rotamer.init();

	utility::vector1< Torsion_Info > rotamer;
	utility::vector1< Real > data;
	Size count1( 0 ), count2( 0 );
	Real sum1( 0 ), sum2( 0 ), val( 0 );

	for ( suite_rotamer.reset(); suite_rotamer.not_end(); ++suite_rotamer ) {
		++count1;
		suite_rotamer.apply( pose );

		data.clear();
		data.push_back( pose.torsion( TorsionID( 3, id::BB, 1 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::BB, 2 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::BB, 3 ) ) );
		data.push_back( pose.torsion( TorsionID( 2, id::BB, 5 ) ) );
		data.push_back( pose.torsion( TorsionID( 2, id::BB, 6 ) ) );
		//data.push_back( pose.torsion( TorsionID( 2, id::CHI, 1 ) ) );
		//data.push_back( pose.torsion( TorsionID( 2, id::CHI, 2 ) ) );
		//data.push_back( pose.torsion( TorsionID( 2, id::CHI, 3 ) ) );
		//data.push_back( pose.torsion( TorsionID( 2, id::BB, 4 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::CHI, 1 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::CHI, 2 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::CHI, 3 ) ) );
		data.push_back( pose.torsion( TorsionID( 3, id::BB, 4 ) ) );

		//std::cout << "NEW ";
		for ( Size i = 1; i <= data.size(); ++i ) {
			val = data[i];
			if ( val > 170.0001 ) {
				val -= 360;
			} else if ( val < -190.0001 ) {
				val += 360;
			}
			//std::cout << val << ' ';
			sum1 += val;
		}
		//std::cout << std::endl;

		if ( rotamer_generator.has_another_rotamer() ) {
			utility::vector1< Torsion_Info > const rot = rotamer_generator.get_next_rotamer();
			//std::cout << "OLD ";
			for ( Size i = 1; i <= rot.size(); ++i ) {
				val = rot[i].value;
				if ( val > 170.0001 ) {
					val -= 360;
				} else if ( val < -190.0001 ) {
					val += 360;
				}
				//std::cout << val << ' ';
				sum2 += val;
			}
			//std::cout << std::endl;
		} else {
			std::cout << "BUG" << std::endl;
		}
	}
	if ( rotamer_generator.has_another_rotamer() )
		std::cout << "BUG" << std::endl;

	std::cout << sum1 << ' ' << sum2 << std::endl;
	pose.dump_pdb( "test.pdb" );
}

///////////////////////////////////////////////////////////////
void*
my_main ( void* ) {
	test();
}


///////////////////////////////////////////////////////////////////////////////
int
main ( int argc, char * argv [] ) {
  try {
		devel::init ( argc, argv );
		protocols::viewer::viewer_main ( my_main );
  } catch ( utility::excn::EXCN_Base const & e ) {
    std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }
}
