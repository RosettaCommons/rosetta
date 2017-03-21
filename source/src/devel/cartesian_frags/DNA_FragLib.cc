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


#include <devel/cartesian_frags/DNA_FragLib.hh>

// libRosetta headers


#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


// // // C++ headers

// //silly using/typedef


namespace devel {
namespace cartesian_frags {

/// @details Auto-generated virtual destructor
DNA_FragLib::~DNA_FragLib() = default;

using namespace core;

static THREAD_LOCAL basic::Tracer tt( "devel.cartesian_frags.DNA_FragLib", basic::t_trace );
static THREAD_LOCAL basic::Tracer td( "devel.cartesian_frags.DNA_FragLib", basic::t_debug );
static THREAD_LOCAL basic::Tracer ti( "devel.cartesian_frags.DNA_FragLib", basic::t_info );
static THREAD_LOCAL basic::Tracer tw( "devel.cartesian_frags.DNA_FragLib", basic::t_warning );


core::pose::Pose &
DNA_FragLib::suite_pose() const
{
	if ( suite_pose_.size() == 0 ) {
		setup_suite_pose( suite_pose_ );
	}
	return suite_pose_;
}


///////////////////////////////////////////////////////////////////////////////
void
build_frag_libraries(
	utility::vector1< std::string > const & files,
	DNA_FragLib & lib
)
{
	using namespace pose;
	using namespace kinematics;
	using namespace id;
	using namespace conformation;
	using namespace chemical;
	using namespace scoring::dna;


	for ( Size n=1; n<= files.size(); ++n ) {
		Pose pose;
		td << "Reading " << files[n] << std::endl;

		core::import_pose::pose_from_file( pose, files[n] , core::import_pose::PDB_file);

		set_base_partner( pose ); // fills base partner info

		BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

		Conformation const & conf( pose.conformation() );

		Size const nres( pose.size() );
		for ( Size i=1; i<= nres; ++i ) {
			Residue const & rsd( pose.residue(i) );

			if ( !rsd.is_DNA() ) continue;
			if ( !partner[i] ) continue; // unpaired
			if ( rsd.is_terminus() ) continue;

			// the takeoff and landing stubs for the sugar fragment:

			// new way:
			TorsionStubID const     chi( TorsionID( i, CHI, 1 ), Backward );
			TorsionStubID const   gamma( TorsionID( i,  BB, 3 ), Backward );
			TorsionStubID const epsilon( TorsionID( i,  BB, 5 ), Forward );

			utility::vector1< TorsionStubID > outgoing;
			outgoing.push_back( gamma   ); // first gamma,
			outgoing.push_back( epsilon ); // then epsilon

			lib.sugars.push_back( CartesianFragment( chi, SafeAtomID( "C1*", i ), outgoing, conf ) );

			SafeAtomID const N3( "N3", i ), N3_partner( "N3", partner[i] ), N3_next( "N3", i+1 );

			//if ( partner[i] > i ) {
			{
				// the basepair
				std::string const bp( std::string() + rsd.name1() + pose.residue( partner[i] ).name1() );
				TorsionStubID const partner_chi( TorsionID( partner[i], CHI, 1 ), Backward );
				lib.add_base_pair( bp, CartesianFragment( chi, N3, partner_chi, SafeBondID( N3, N3_partner ), conf ) );
			}


			// include the suite from i to i+1 if:
			//
			// -- i+1 not terminus
			// -- i+1 is paired
			// -- i+1 is paired to partner[i] - 1
			//
			if ( !pose.residue(i+1).is_terminus() && partner[i+1] && ( partner[i+1] == partner[i] - 1 ) ) {
				TorsionStubID const gamma_next( TorsionID( i+1, BB, 3 ), Backward );

				lib.forward_suites.push_back(  CartesianFragment( epsilon   , SafeAtomID( "O3*", i  ), gamma_next, conf));
				lib.backward_suites.push_back( CartesianFragment( gamma_next, SafeAtomID( "C5*", i+1),    epsilon, conf));

				// the basestep
				std::string const bs( std::string() + rsd.name1() + pose.residue(i+1).name1() );
				TorsionStubID const next_chi( TorsionID( i+1, CHI, 1 ), Backward );
				lib.add_base_step( bs, CartesianFragment( chi, N3, next_chi, SafeBondID( N3, N3_next ), conf ) );
			}
		} // i
	} // files

	ti << "FragLib: sugars= " << lib.sugars.size() << " suites= " << lib.forward_suites.size() <<
		" basepairs: at= " << lib.base_pairs("at").size() << " gc= " << lib.base_pairs("gc").size() << std::endl;

}
///////////////////////////////////////////////////////////////////////////////

void
setup_suite_pose(
	pose::Pose & suite_pose
)
{
	using namespace conformation;
	using namespace chemical;

	suite_pose.clear();

	ResidueTypeSet const & rsd_set( *chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	ResidueType const & vrt1_type( rsd_set.name_map("VRT1BB" ) );
	ResidueType const & vrt3_type( rsd_set.name_map("VRT3BB" ) );
	ResidueOP vrt1( ResidueFactory::create_residue( vrt1_type ) );
	ResidueOP vrt3( ResidueFactory::create_residue( vrt3_type ) );
	suite_pose.append_residue_by_bond( *vrt3 );
	vrt1->set_xyz("ORIG", Vector(4,0,0) );
	suite_pose.append_residue_by_bond( *vrt1 );
	vrt1->set_xyz("ORIG", Vector(4,2,0) );
	suite_pose.append_residue_by_bond( *vrt1 );
	vrt1->set_xyz("ORIG", Vector(6,2,0) );
	suite_pose.append_residue_by_bond( *vrt3 );
}


} // ns cartesian_frags
} // ns devel

//   { // some chainbreak diagnostics
//    using namespace scoring;
//    setup_dna_chainbreak_constraints( pose );
//    ScoreFunction scorefxn;
//    scorefxn.set_weight( atom_pair_constraint, 1.0 );
//    scorefxn.set_weight(     angle_constraint, 1.0 );
//    Size nstep(0.0);
//    for ( Size i=1; i< pose.size(); ++i ) {
//     Residue const & rsd1( pose.residue( i   ) );
//     Residue const & rsd2( pose.residue( i+1 ) );
//     if ( rsd1.is_DNA() && !rsd1.is_upper_terminus() && rsd2.is_DNA() && !rsd2.is_lower_terminus() ) {
//      ++nstep;
//      if ( rsd1.xyz("O3*").distance( rsd2.xyz("P") ) > 2.0 ) {
//       td.Warning << "long bond in frag input file: " << files[n] << ' ' << i << ' ' <<
//        rsd1.xyz("O3*").distance( rsd2.xyz("P") ) << std::endl;
//      }
//     }
//    }
//    Real const cstscore( scorefxn( pose ) );
//    tt << "cstscore " << files[n] << ' ' << cstscore << " nstep: " << nstep << " score/step " <<
//     cstscore/nstep << std::endl;
//    scorefxn.show( tt, pose );
//   }
