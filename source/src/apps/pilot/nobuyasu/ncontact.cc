// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief Nobuyasu Koga

// projects headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <devel/init.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/fldsgn/sspot/SS_Info2.hh>

// utility headers
#include <utility/vector1.hh>

// C++ headers
#include <cmath>

//
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::jobdist;
using namespace protocols::moves;


static thread_local basic::Tracer TT( "ncontact" );

namespace ncontact
{
	RealOptionKey dist( "ncontact:dist" );
	IntegerOptionKey seqsep( "ncontact:sep" );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class MyAtom {

public:
	MyAtom(){

		for ( Size j=1; j<= 50 ; ++j ) {
			atom_hydrophobic_.push_back(false);
		}

		atom_hydrophobic_[ 3 ] = true;  // CH1
		atom_hydrophobic_[ 4 ] = true;  // CH2
		atom_hydrophobic_[ 5 ] = true;  // CH3
		atom_hydrophobic_[ 6 ] = true;  // aroC
		atom_hydrophobic_[ 19 ] = true; // CAbb

	}

	bool is_hydrophobic( Size const index ){
		return atom_hydrophobic_[ index ];
	}

	private:
	utility::vector1< bool > atom_hydrophobic_;

}; // MyAtom

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class CountContact : public protocols::moves::Mover {

	typedef core::Real Real;
	typedef core::Size Size;
	typedef ObjexxFCL::FArray2D< Size > FArray2D_Size;

	typedef core::conformation::Residue Residue;
	typedef core::scoring::ScoreFunctionOP ScoreFunctionOP;
	typedef core::scoring::dssp::Dssp Dssp;
	typedef protocols::fldsgn::sspot::SS_Info2 SS_Info2;

public:

	/// @brief default constructor
	CountContact():
		Mover( "CountCotact" ),
		condist_( 36.0 ),
		isep_residue_( 4 )
	{
		scorefxn_ = core::scoring::get_score_function();
		initialize();
		//write_header();
	}

	/// @brief value constructor
	CountContact( Real const cdist,
								Size const sep ):
		Mover( "CountCotact" ),
		condist_( numeric::square(cdist) ),
		isep_residue_( sep )
	{
		scorefxn_ = core::scoring::get_score_function();
		initialize();
		//write_header();
	}

	/// @brief
	virtual void apply ( core::pose::Pose & );

	/// @brief
	virtual ~CountContact(){};

private:

	/// @brief
	void initialize() {
		ncontact_allatm_ = 0;
		ncontact_hpatm_ = 0;
		ncontact_atm_hpres_ = 0;
		ncontact_res_hh_ = 0;
		ncontact_res_eh_ = 0;
	}

	/// @brief
	/*
	void write_header() const {
		TT << "rms" << " "	<< "ncon_res" << " " << "ncon_atm_hpres" << " "
			 << "ncon_atm" << " " << "ncon_hpatm" << std::endl;
	}
	*/

private:

	/// @brief distact**2 used for juding contact pair
	Real condist_;

	/// @brief residue pairs of i < i+isep_residue_ are used for counting #countacts
	Size isep_residue_;

	/// @brief #atom-contacts among all heavy atoms
	Size ncontact_allatm_;

	/// @brief #atom-contacts among hydrophobic heavy atoms
	Size ncontact_hpatm_;

	/// @brief #atom-contacts among heavy atoms of sidechains of hydrophobic residues
	Size ncontact_atm_hpres_;

	/// @brief #residue-contacts among helices, defeined only by Ca
	Size ncontact_res_hh_;

	/// @brief #residue-contacts between helices and strands, defeined only by Ca
	Size ncontact_res_eh_;

	/// @brief score function
	ScoreFunctionOP scorefxn_;

}; //CountContact

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void CountContact::apply ( pose::Pose & pose ){

	// secondary structure setting
	Dssp dssp( pose );
	dssp.insert_ss_into_pose( pose );
	SS_Info2 ssinfo( pose.secstruct() );

	// intialize
	Size max_ssele = ssinfo.ss_element_id( pose.total_residue() );
	utility::vector1< utility::vector1< Size > > ncon_sselements( max_ssele, utility::vector1<core::Size>(max_ssele, 0) );
	for ( Size i=1 ;i<=max_ssele; ++i ) {
		for ( Size j=1 ;j<=max_ssele; ++j ) {
			ncon_sselements[i][j] = 0;
		}
	}

	initialize();

	// clac score
	(*scorefxn_)( pose );

	// calc number of contacts
	MyAtom myatom;
	for ( Size iaa=1 ;iaa<= pose.total_residue()-isep_residue_; ++iaa ) {
		Size iaa_ssele( ssinfo.ss_element_id( iaa ) );

		for ( Size jaa=iaa+isep_residue_ ;jaa<=pose.total_residue(); ++jaa ) {
			Size jaa_ssele( ssinfo.ss_element_id( jaa ) );

			Residue const & ires( pose.residue( iaa ) );
			for ( Size iatm=1, iatm_end=ires.natoms(); iatm<=iatm_end; ++iatm ) {

				if( ires.atom_type( int(iatm) ).is_heavyatom() ){
					Residue const & jres( pose.residue( jaa ) );

					for ( Size jatm=1, jatm_end=jres.natoms(); jatm<=jatm_end; ++jatm ) {

						if( jres.atom_type(int(jatm)).is_heavyatom() ){

							conformation::Atom const & iatom( ires.atom( iatm ) );
							conformation::Atom const & jatom( jres.atom( jatm ) );
							Size const iadex ( ires.atom_type_index( int(iatm) ) );
							Size const jadex ( jres.atom_type_index( int(jatm) ) );

							Real const dsq( distance_squared( iatom.xyz(), jatom.xyz() ));

							if( dsq <= condist_ ){
								ncontact_allatm_ ++;

								if( myatom.is_hydrophobic(iadex) && myatom.is_hydrophobic(jadex) ){
									ncontact_hpatm_ ++;

									if( !ires.atom_is_backbone( int(iatm) ) && !jres.atom_is_backbone( int(jatm) ) &&
        						   ires.type().is_polar() != 1 && jres.type().is_polar() != 1 ){
										ncontact_atm_hpres_ ++;

										ncon_sselements[iaa_ssele][jaa_ssele] ++ ;
										//std::cout << iaa  << " "  << jaa  << " " << iatm << " "  << jatm << " "
										//<< ires.name() << " "
										//<< jres.name() << " "
										//<< ires.atom_type( iatm ).name() << " "
										//<< jres.atom_type( jatm ).name() << " "
										//<< dsq << std::endl;
									}
								}
							} // dsq

						} //fi jatm is_hydrogen ?
					} // loop for jatm
				} // fi iatm is_hydrogen ?
			} // loop for iatm

		} // loop for jaa
	}// loop for iaa

	//
	Size tot( 0 );
	for ( Size i=1 ;i<=max_ssele; ++i ) {
		for ( Size j=1 ;j<=max_ssele; ++j ) {
			tot += ncon_sselements[i][j];
		}
	}

	//

	float ss_entrpy( 0.0 );
	for ( Size i=1 ;i<=max_ssele; ++i ) {
		for ( Size j=1 ;j<=max_ssele; ++j ) {
			if( ncon_sselements[i][j] > 0 ){
				float prob = Real(ncon_sselements[i][j])/Real(tot);
				ss_entrpy += -prob*std::log( prob );
			}
		}
	}


	float ncatm   ( Real(ncontact_allatm_)/Real(pose.total_residue()));
	float nchpatm ( Real(ncontact_hpatm_)/Real(pose.total_residue())) ;
  float nchpres ( Real(ncontact_atm_hpres_)/Real(pose.total_residue()));

	setPoseExtraScore( pose, "nres", float(pose.total_residue()) );
	setPoseExtraScore( pose, "ncon_atm", ncatm  );
	setPoseExtraScore( pose, "ncon_hpatm", nchpatm );
	setPoseExtraScore( pose, "ncon_atm_hpres", nchpres );
	setPoseExtraScore( pose, "ss_entrpy", ss_entrpy );

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	option.add( ncontact::dist, "distance of contact" );
	option.add( ncontact::seqsep,  "sequence separation distance of sequence to be counted" );

	devel::init( argc, argv );

	Size sep( 4 );
	Real dist( 6.0 );
	if( option[ ncontact::dist ].user() ){
		dist = option[ ncontact::dist ];
	}
	if( option[ ncontact::seqsep ].user() ){
		sep = option[ ncontact::seqsep ];
	}

	MoverOP protocol;
  protocol = new CountContact( dist, sep );

	universal_main( *protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
