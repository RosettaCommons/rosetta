// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/jkleman/mp_lipid_acc.cc
/// @brief  compute lipid accessibility from a membrane protein structure
/// @author  JKLeman (julia.koehler.leman@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

//#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

//#include <core/conformation/membrane/MembraneInfo.hh>
//#include <core/conformation/membrane/util.hh>
//#include <core/conformation/membrane/SpanningTopology.hh>
//#include <core/conformation/membrane/Span.hh>
//#include <protocols/membrane/util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/membrane/MPLipidAccessibility.hh>
//#include <protocols/moves/MoverContainer.hh>
//#include <protocols/membrane/AddMembraneMover.hh>
//// I have no idea why this is needed??? but it doesn't compile without it
//#include <protocols/membrane/TranslationRotationMover.hh>
//#include <basic/options/keys/mp.OptionKeys.gen.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
//#include <core/scoring/sasa/SasaCalc.hh>
//#include <core/scoring/sasa/util.hh>
//
//// Utility Headers
//#include <core/types.hh>
//#include <core/pose/util.hh>
#include <utility/pointer/owning_ptr.hh>
//#include <numeric/xyzVector.hh>
//#include <numeric/xyz.functions.hh>
//#include <numeric/conversions.hh>
//#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
//#include <numeric/numeric.functions.hh>
//#include <utility/numbers.hh>
//
//#include <core/membrane/hull.hh>

// C++ headers
#include <iostream>
#include <cstdlib>

static basic::Tracer TR( "apps.pilot.jkleman.mp_lipid_acc" );

using namespace core;
using namespace core::pose;
//using namespace numeric;
//using namespace protocols::moves;
//using namespace core::conformation::membrane;
//using namespace protocols::simple_moves;
using namespace protocols::membrane;
//using namespace basic::options;

/// @brief Load Membrane Mover: Load and Initialize a Membrane Protein
/// using RosettaMP
//class MPLipidAccessibility : public Mover {
//public:
//
// /////////////////////
// /// Constructors  ///
// /////////////////////
//
// /// @brief Default constructor
// MPLipidAccessibility();
//
// /// @brief Copy constructor (not needed unless you need deep copies)
// // MPLipidAccessibility( MPLipidAccessibility const & src );
//
// /// @brief Destructor (important for properly forward-declaring smart-pointer members)
// virtual ~MPLipidAccessibility();
//
// /////////////////////
// /// Mover Methods ///
// /////////////////////
//
//public:
// /// @brief Apply the mover
// virtual void
// apply( core::pose::Pose & pose );
//
// virtual void
// show( std::ostream & output = std::cout ) const;
//
// /// @brief Get the name of the Mover
// virtual std::string
// get_name() const;
//
// ///////////////////////////////
// /// Rosetta Scripts Support ///
// ///////////////////////////////
//
// /// @brief parse XML tag (to use this Mover in Rosetta Scripts)
// virtual void
// parse_my_tag(
//     utility::tag::TagCOP tag,
//     basic::datacache::DataMap & data,
//     protocols::filters::Filters_map const & filters,
//     protocols::moves::Movers_map const & movers,
//     core::pose::Pose const & pose );
//
// //MPLipidAccessibility & operator=( MPLipidAccessibility const & src );
//
// /// @brief required in the context of the parser/scripting scheme
// virtual protocols::moves::MoverOP
// fresh_instance() const;
//
// /// @brief required in the context of the parser/scripting scheme
// virtual protocols::moves::MoverOP
// clone() const;
//
//private: // methods
//
// /// @brief set defaults
// void set_defaults();
//
// /// @brief Register Options from Command Line
// void register_options();
//
// /// @brief Initialize from commandline
// void init_from_cmd();
//
// /// @brief finalize setup
// void finalize_setup( Pose & pose );
//
// /// @brief check whether protein is in membrane
// bool protein_in_membrane( Pose & pose );
//
// /// @brief fill up slice arrays with protein data
// void fill_up_slices( Pose & pose );
//
// /// @brief compute slice COM
// void compute_slice_com();
//
//private: // data
//
// /// @brief original data from first implementation
// // core::Real user_radius_angstrom_;
// // core::Real user_rel_sasa_;
// // bool user_;
// core::Real angle_cutoff_;
// // core::Real relative_radius_;
// core::Real slice_width_;
//
// // define variables, the outer vector goes through the slices
// // inner vector goes through residues in each slice
// utility::vector1< utility::vector1< core::Size > > resi_;
// utility::vector1< utility::vector1< core::Vector > > ca_coord_, cb_coord_;
// utility::vector1< core::Vector > slice_com_;
// utility::vector1< core::Real > slice_zmin_;
// core::Real shell_radius_;
//
//};
//
//std::ostream &
//operator<<( std::ostream & os, MPLipidAccessibility const & mover );
//
///////////////////////////////////////////
//
///////////////////////
///// Constructors  ///
///////////////////////
//
///// @brief Default constructor
//MPLipidAccessibility::MPLipidAccessibility() : protocols::moves::Mover() {
//
// register_options();
// set_defaults();
// init_from_cmd();
//}
//
//////////////////////////////////////////////////////////////////////////////////
///// @brief Copy constructor
////MPLipidAccessibility::MPLipidAccessibility( MPLipidAccessibility const & src ):
//// protocols::moves::Mover( src )
////{
////
////}
//
//////////////////////////////////////////////////////////////////////////////////
///// @brief Destructor (important for properly forward-declaring smart-pointer members)
//MPLipidAccessibility::~MPLipidAccessibility(){}
//
//////////////////////////////////////////////////////////////////////////////////
///// @brief Show the contents of the Mover
//void
//MPLipidAccessibility::show(std::ostream & output) const
//{
// protocols::moves::Mover::show(output);
//}
//
///// @brief Get the name of the Mover
//std::string
//MPLipidAccessibility::get_name() const {
// return "MPLipidAccessibility";
//}
//
//////////////////////////////////////////////////////////////////////////////////
///// Rosetta Scripts Support ///
/////////////////////////////////
//
///// @brief parse XML tag (to use this Mover in Rosetta Scripts)
//void
//MPLipidAccessibility::parse_my_tag(
//           utility::tag::TagCOP ,
//           basic::datacache::DataMap& ,
//           protocols::filters::Filters_map const & ,
//           protocols::moves::Movers_map const & ,
//           core::pose::Pose const & )
//{
//
//}
//
//////////////////////////////////////////////////////////////////////////////////
///// @brief required in the context of the parser/scripting scheme
//protocols::moves::MoverOP
//MPLipidAccessibility::fresh_instance() const
//{
// return protocols::moves::MoverOP( new MPLipidAccessibility );
//}
//
///// @brief required in the context of the parser/scripting scheme
//protocols::moves::MoverOP
//MPLipidAccessibility::clone() const
//{
// return protocols::moves::MoverOP( new MPLipidAccessibility( *this ) );
//}
//
//std::ostream &
//operator<<( std::ostream & os, MPLipidAccessibility const & mover )
//{
// mover.show(os);
// return os;
//}
//
//////////////////////////////////////////////////////////////////////////////////
///// MOVER METHODS ///
///////////////////////
//void MPLipidAccessibility::apply( core::pose::Pose & pose ){
//
// using namespace core::id;
//
// // finalize setup
// finalize_setup( pose );
//
// ////////////////// SET B-FACTOR TO ZERO FOR ALL ATOMS //////////////////
//
// for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {
//  for ( core::Size a = 1; a <= pose.residue( r ).natoms(); ++a ) {
//   pose.pdb_info()->bfactor( r, a, 0 );
//  }
// }
//
// //////////////////////// COMPUTE SASA //////////////////////////////////
//
// // compute relative SASA for the protein,
// // be aware: the SASA vector goes over all residues
// // utility::vector1< core::Real > rel_res_sasa = core::scoring::sasa::rel_per_res_sc_sasa( pose );
//
//
// //////////////////////// GO THROUGH SLICES /////////////////////////////
//
// // compute concave hull with shell boundary
// core::id::AtomID_Map< bool > shell = core::membrane::concave_shell( pose, slice_zmin_[ 1 ], slice_zmin_[ slice_zmin_.size() ] + slice_width_, slice_width_, shell_radius_ );
//
// // get lipid-exposed residues from the AtomID_map
// utility::vector1< core::Size > boundary;
// for ( core::Size r = 1; r <= nres_protein( pose ); ++r ) {
//
//  // create AtomID for CA atom
//  AtomID aid = AtomID( 2, r );
//
//  // if atomID has a true in the map, add to boundary vector
//  if ( shell.get( aid ) == true ) {
//   boundary.push_back( r );
//  }
// }
//
// // go through slices
// for ( core::Size s = 1; s < slice_zmin_.size(); ++s ) {
//
//  // go through residues
//  for ( core::Size r = 1; r <= resi_[ s ].size(); ++r ) {
//
//   bool lipid_exposed( false );
//
//   // if residues are on boundary, then lipid exposed
//   if ( boundary.has_value( resi_[ s ][ r ] ) ) {
//    lipid_exposed = true;
//   }
//
//   // if angle between COM-CA-CB is smaller than cutoff, then not lipid-exposed
//   // compute angle between COM-CA-CB
//   // if angle < angle_cutoff, then sidechain face COM
//   // if angle >= angle_cutoff, then sidechain faces away from COM
//   core::Real angle = numeric::angle_degrees( slice_com_[ s ], ca_coord_[ s ][ r ], cb_coord_[ s ][ r ] );
//
//   // go through atoms in residue and set B-factor
//   for ( core::Size a = 1; a <= pose.residue( resi_[ s ][ r ] ).natoms(); ++a ) {
//
//    // debug
//    //    if ( a == 2 ) {
//    //     TR << "residue " << resi_[ s ][ r ] << " angle " << angle << " accessibility " << rel_res_sasa[ resi_[ s ][ r ] ] << " ca_com " << ca_com.length() << " radius " << radius*relative_radius << std::endl;
//    //    }
//
//    // is this a 2TM protein? if yes, all residues are lipid exposed
//    if ( pose.conformation().membrane_info()->spanning_topology()->nspans() <= 2 ) {
//     pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    }
//
//    // if lipid exposed and COM-CA-CB larger than cutoff, i.e. facing outwards
//    if ( lipid_exposed == true && angle > angle_cutoff_ ) {
//     pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    }
//
//
//    //    // user-defined, i.e. uses SASA and radius in Å instead of relative radius
//    //    if ( user == true ) {
//    //
//    //     // facing outwards and larger than certain radius with user-defined
//    //     // radius and accessibility
//    //     if ( angle > angle_cutoff && ca_com.length() >= user_radius_angstrom && rel_res_sasa[ resi_[ s ][ r ] ] >= user_rel_sasa ) {
//    //      pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    //     } else if ( angle > angle_cutoff && ca_com.length() >= user_radius_angstrom && rel_res_sasa[ resi_[ s ][ r ] ] == 0 && pose.residue( resi_[ s ][ r ] ).name3() == "GLY" ) {
//    //      // same, just for GLY
//    //      pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    //     }
//    //    } else {
//    //     // doesn't use SASA (too finicky) and uses relative radius, default app
//    //
//    //     // facing outwards and larger than certain radius
//    //     if ( angle > angle_cutoff && ca_com.length() >= radius * relative_radius ) {
//    //      pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    //     }
//    //
//    //     //  && rel_res_sasa[ resi_[ s ][ r ] ] > 0.2
//    //     // for single helices and 2 helices, mostly lipid exposed
//    //     if ( pose.conformation().membrane_info()->spanning_topology()->nspans() <= 2 ) {
//    //      pose.pdb_info()->bfactor( resi_[ s ][ r ], a, 50.0 );
//    //     }
//    //    } // user or not
//   } // atoms
//  } // residues
// } // slices
//
//}// apply
//
//////////////////////////////////////////////////////////////////////////////////
///// private methods ///
/////////////////////////
//
///// @brief set defaults
//void MPLipidAccessibility::set_defaults() {
//
// // original data
// // user_radius_angstrom_ = 10;
// // user_rel_sasa_ = 0.2;
// // user_ = false;
// angle_cutoff_ = 65.0;
// // relative_radius_ = 0.65;
//
// // slices from hull data
// slice_width_ = 5.0;
// shell_radius_ = 5.0;
//
//}// set defaults
//
////////////////////////////////////////////
///// @brief Register Options from Command Line
//void MPLipidAccessibility::register_options() {
//
// using namespace basic::options;
// // option.add_relevant( OptionKeys::mp::lipid_acc::radius_cutoff );
// // option.add_relevant( OptionKeys::mp::lipid_acc::rel_sasa_cutoff );
// option.add_relevant( OptionKeys::mp::lipid_acc::angle_cutoff );
//
//}// register options
//
////////////////////////////////////////////
///// @brief Initialize from commandline
//void MPLipidAccessibility::init_from_cmd() {
//
// using namespace basic::options;
//
// // if ( option[ OptionKeys::mp::lipid_acc::radius_cutoff ].user() ) {
// //  user_radius_angstrom = option[ OptionKeys::mp::lipid_acc::radius_cutoff ]();
// //  user = true;
// // }
// // if ( option[ OptionKeys::mp::lipid_acc::rel_sasa_cutoff ].user() ) {
// //  user_rel_sasa = option[ OptionKeys::mp::lipid_acc::rel_sasa_cutoff ]();
// //  user = true;
// // }
// if ( option[ OptionKeys::mp::lipid_acc::angle_cutoff ].user() ) {
//  angle_cutoff_ = option[ OptionKeys::mp::lipid_acc::angle_cutoff ]();
////  user = true;
// }
//
//}// init from cmd
//
////////////////////////////////////////////
///// @brief finalize setup
//void MPLipidAccessibility::finalize_setup( Pose & pose ){
//
// ////////// SWITCH TO FULL-ATOM, REQUIRED FOR VECTOR DEFINITIONS ////////
// using namespace protocols::simple_moves;
// SwitchResidueTypeSetMoverOP full_atom( new SwitchResidueTypeSetMover( "fa_standard" ) );
// full_atom->apply( pose );
//
// // add membrane
// using namespace protocols::membrane;
// AddMembraneMoverOP addmem( new AddMembraneMover() );
// addmem->apply( pose );
//
// // check whether protein is in membrane
// protein_in_membrane( pose );
//
//}// finalize setup
//
////////////////////////////////////////////
///// @brief check whether protein is transformed into membrane
//bool MPLipidAccessibility::protein_in_membrane( Pose & pose ){
//
// bool in_membrane( true );
//
// // get minimum and maximum z-coordinate
// // check whether structure is transformed into membrane, get CA position
// core::Real min_z( 10000 );
// core::Real max_z( -10000 );
//
// for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {
//
//  // skip for membrane residue
//  if ( pose.residue( i ).name3() == "MEM" ) {
//   continue;
//  }
//
//  if ( pose.residue( i ).xyz( "CA" ).z() < min_z ) {
//   min_z = pose.residue( i ).xyz( "CA" ).z();
//  } else if ( pose.residue( i ).xyz( "CA" ).z() > max_z ) {
//   max_z = pose.residue( i ).xyz( "CA" ).z();
//  }
// }
//
// // crude check: if no CA atoms with positive AND negative z-coordinates,
// // then the protein doesn't span the membrane
// if ( !( min_z < 0 && max_z > 0 ) ) {
//  in_membrane = false;
//  TR.Warning << "YOUR PROTEIN DOES NOT SPAN THE MEMBRANE! EITHER YOU KNOW WHAT YOU ARE DOING OR YOUR PROTEIN IS NOT TRANSFORMED INTO MEMBRANE COORDINATES!!!" << std::endl;
//  //  utility_exit_with_message( "Your protein is not transformed into membrane coordinates! Cannot compute lipid accessibility. Quitting." );
// }
// return in_membrane;
//
//}// protein in membrane?
//
////////////////////////////////////////////
///// @brief fill up slice arrays with protein data
//void MPLipidAccessibility::fill_up_slices( Pose & pose ) {
//
// core::Real iter = - pose.conformation().membrane_info()->membrane_thickness();
//
// // go through protein in the membrane and get slices and arrays for them
// while ( iter <= pose.conformation().membrane_info()->membrane_thickness() ) {
//
//  slice_zmin_.push_back( iter );
//
//  // temp utility vectors for each slice
//  utility::vector1< core::Size > slice_res;
//  utility::vector1< core::Vector > slice_ca, slice_cb;
//
//  // go through protein residues
//  for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {
//
//   // skip membrane residue
//   if ( pose.residue( i ).name3() == "MEM" ) {
//    continue;
//   }
//
//   // get CA and CB coordinates
//   core::Vector ca = pose.residue( i ).xyz( "CA" );
//   core::Vector cb;
//
//   // use 2HA instead of CB for glycine
//   if ( pose.residue( i ).name3() == "GLY" ) {
//    cb = pose.residue( i ).xyz( "2HA" );
//   } else {
//    cb = pose.residue( i ).xyz( "CB" );
//   }
//
//   // if CA-coordinate is within slice, add resi, CA and CB to arrays
//   if ( ( ca.z() >= iter && ca.z() < iter + slice_width_ ) or
//    ( cb.z() >= iter && cb.z() < iter + slice_width_ ) ) {
//
//    slice_res.push_back( i );
//    slice_ca.push_back( ca );
//    slice_cb.push_back( cb );
//
//   }
//  }
//
//  // push back into slices
//  resi_.push_back( slice_res );
//  ca_coord_.push_back( slice_ca );
//  cb_coord_.push_back( slice_cb );
//
//  // go to next slice
//  iter += slice_width_;
//
// } // get arrays of slices
//
//}// fill up slice data
//
////////////////////////////////////////////
///// @brief compute slice COM
//void MPLipidAccessibility::compute_slice_com(){
//
// // go through slices and compute COMs
// for ( core::Size s = 1; s < slice_zmin_.size(); ++s ) {
//
//  core::Vector com( 0, 0, 0 );
//  core::Real radius = 0;
//
//  // go through residues, compute COM per 5Å slice
//  for ( core::Size r = 1; r <= resi_[ s ].size(); ++r ) {
//   com += ca_coord_[ s ][ r ];
//  }
//  com /= resi_[ s ].size();
//
//  slice_com_.push_back( com );
// }
//
//}// compute slice COM


/////////////////////////////////////////

using MPLipidAccessibilityOP = utility::pointer::shared_ptr<MPLipidAccessibility>;

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		devel::init(argc, argv);

		using namespace protocols::membrane;

		protocols::jd2::register_options();

		// Create and kick off a new load membrane mover
		MPLipidAccessibilityOP lipid_acc( new MPLipidAccessibility() );
		protocols::jd2::JobDistributor::get_instance()->go( lipid_acc );

		return 0;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

