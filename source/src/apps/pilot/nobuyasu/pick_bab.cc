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
/// @brief
/// @author Nobuyasu Koga

// Package headers
#include <protocols/jd2/parser/BluePrint.hh>
#include <core/util/ABEGOManager.hh>
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <protocols/fldsgn/topology/BetaAlphaBetaMotif.hh>

// Project headers
#include <basic/MetricValue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MakePolyXMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>
#include <utility/excn/Exceptions.hh>

// C++ header
#include <fstream>
#include <map>

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

typedef std::string String;


static basic::Tracer TR("pick_bab");

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, output )
OPT_KEY( File, blue )
OPT_KEY( Integer, abego )
OPT_KEY( Real, pore_radius )
OPT_KEY( Integer, maxread )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( output, "output file name for log files", "bab" );
	NEW_OPT( abego, "output level of abego discription", 3 );
	NEW_OPT( blue, "blueprint", "" );
	NEW_OPT( pore_radius, "pore radius for sasa", 2.0 );
	NEW_OPT( maxread, "max number to read structures", 0 );
}


/// local mover for testing purposes
class PickBAB : public protocols::moves::Mover {
public: // typedef

	typedef protocols::jd2::parser::BluePrint BluePrint;
	typedef protocols::jd2::parser::BluePrintOP BluePrintOP;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;
	typedef protocols::fldsgn::topology::StrandCOP StrandCOP;
	typedef protocols::fldsgn::topology::HelixCOP HelixCOP;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::SheetSet SheetSet;
	typedef protocols::fldsgn::topology::SheetSetOP SheetSetOP;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotif BetaAlphaBetaMotif;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotifs BetaAlphaBetaMotifs;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotifSet BetaAlphaBetaMotifSet;


public: // constructor/deconstructor


	PickBAB()
	{
		// read blueprint
		ssinput_ = false;
		if( option[ blue ].user() ) {
			BluePrintOP blueprint = new BluePrint( option[ blue ]() );
			ssinfo_ = new SS_Info2( blueprint->secstruct() );
			ssinput_ = true;
		} else {
			ssinfo_ = new SS_Info2;
		}

		// output file
		std::ostringstream filename;
		filename <<  option[ output ]() << ".out";
		output_.open( filename.str().c_str() ,std::ios::out );

		// output file for phi psi
		filename.str("");
		filename <<  option[ output ]() << ".phipsi.out";
		output2_.open( filename.str().c_str() ,std::ios::out );

		output_ << "id " << "filename " << "updown1 " << "updown2 " << "positions " << "babmotif " << "cross_over "
				<< "hs_dist " << "hs_angle1 " << "hs_angle2 " << "dsasa " << "broken_hb " << "lp_hb " << "helix_cycle "
				<< "helix_bend " << "rama1 " << "rama2 " << "lp1E " << "lp2E "
				<< "len_s1 " << "len_l1 " << "len_h " << "len_l2 " << "len_s2 "
				<< "abe_s1e " << "abe_l1 " << "abe_hb "
				<< "seq_s1e " << "seq_l1 " << "seq_hb "
				<< "abe_he " << "abe_l2 " << "abe_s2b "
				<< "seq_he " << "seq_l2 " << "seq_s2b "
				<< "spairs " << std::endl;

		// set scorefxn
		scorefxn_centroid_ = core::scoring::get_score_function( false );
		scorefxn_fullatom_ = core::scoring::get_score_function();

		num_bab_ = 0;
		nstruct_ = 0;

	}
	virtual ~PickBAB(){ output_.close(); }

	virtual std::string get_name() const { return ""; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return new PickBAB;
	}


public: // apply


	virtual
	void
	apply( core::pose::Pose & pose )
	{
		using core::util::switch_to_residue_type_set;
		using namespace ObjexxFCL::format;
		using namespace protocols::jd2;
		JobOP job_me( JobDistributor::get_instance()->current_job() );
		String me( JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

		nstruct_++;
		if( nstruct_	> option[ maxread ]() && option[ maxread ].user() ) {
			return;
		}

		core::scoring::dssp::Dssp dssp( pose );
		if( ssinput_ ) {
			ssinfo_->set_SSorient( pose );
		} else {
			ssinfo_->initialize( pose, dssp.get_dssp_secstruct() );
		}

		// calc strand pairing
		StrandPairingSetOP spairset
			= new StrandPairingSet( protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo_ ) );

			// calc sheet
		SheetSetOP sheet_set = new SheetSet( ssinfo_, spairset );

		// calc bab
		BetaAlphaBetaMotifSet babset( ssinfo_, sheet_set );

		// make poly ala
		protocols::simple_moves::MakePolyXMover make_poly_ala( "ALA", false, false, false );
		make_poly_ala.apply( pose );

		// set rama to score function
		scorefxn_fullatom_->set_weight( core::scoring::rama, 1.0 );

		// calc fullatom energy
		( *scorefxn_fullatom_ )( pose );

		// make poly gly
		protocols::simple_moves::MakePolyXMover make_poly_gly( "GLY", false, false, false );
		Pose gly_pose( pose );
		make_poly_gly.apply( gly_pose );
		( *scorefxn_fullatom_ )( gly_pose );

//		read abego setup
		Size level( option[ abego ]() );

		for( BetaAlphaBetaMotifs::const_iterator iter = babset.bab_motifs().begin(),
					 iter_end = babset.bab_motifs().end(); iter != iter_end; ++iter ) {
			BetaAlphaBetaMotif & bab( **iter );

			if ( ! bab.is_lefthanded() ) {
				StrandCOP strand1 = ssinfo_->strand( bab.strand1() );
				StrandCOP strand2 = ssinfo_->strand( bab.strand2() );
				HelixCOP helix = ssinfo_->helix( bab.helix() );

				using protocols::forge::build::Interval;
				utility::vector1< Interval > intervals;
				Interval interval( strand1->end()+1, strand2->begin()-1 );
				intervals.push_back( interval );
				Real dsasa = protocols::fldsgn::topology::calc_delta_sasa( pose, intervals, option[ pore_radius ]() );

				if( helix->length() < 5 ) continue;
				if( strand1->length() < 3 ) continue;
				if( strand2->length() < 3 ) continue;

				Size loop1_length = helix->begin() - strand1->end() - 1;
				Size loop2_length = strand2->begin() - helix->end() - 1;
				if( loop1_length > 10 ) continue;
				if( loop2_length > 10 ) continue;

				Size pdblen = pose.pdb_info()->number( strand2->end() ) - pose.pdb_info()->number( strand1->begin() );
				if( pdblen != ( strand2->end() - strand1->begin() ) ) continue;

				// get inout of pleats
				if( pose.aa( strand1->end() ) == core::chemical::aa_gly &&
						pose.aa( strand1->end()-1 ) == core::chemical::aa_gly ) continue;
				if( pose.aa( strand2->begin() ) == core::chemical::aa_gly &&
						pose.aa( strand2->begin()+1 ) == core::chemical::aa_gly ) continue;

				String inout1("");
				if( pose.aa( strand1->end() ) == core::chemical::aa_gly ) {
					Size io = bab.calc_inout( ssinfo_, strand1->end()-1 );
					if( io == 1 ) {
						inout1 = "DOWN";
					} else {
						inout1 = "UP";
					}
				} else {
					Size io = bab.calc_inout( ssinfo_, strand1->end() );
					if( io == 1 ) {
						inout1 = "UP";
					} else {
						inout1 = "DOWN";
					}
				}

				String inout2("");
				if( pose.aa( strand2->begin() ) == core::chemical::aa_gly ) {
					Size io = bab.calc_inout( ssinfo_, strand2->begin()+1 );
					if( io == 1 ) {
						inout2 = "DOWN";
					} else {
						inout2 = "UP";
					}
				} else {
					Size io = bab.calc_inout( ssinfo_, strand2->begin() );
					if( io == 1 ) {
						inout2 = "UP";
					} else {
						inout2 = "DOWN";
					}
				}

				if( inout2 == "UP" && pose.aa( strand2->begin() - 1 ) == core::chemical::aa_pro ) {
					inout2 = "PDWN";
				}

				// get cross over
				Size cross_over = bab.cross_over();

				// get abego for loop1
				Size loop1_begin = strand1->end() + 1;
				Size loop1_end = helix->begin() - 1;
				utility::vector1< String > abegop1_b = am_.get_symbols( pose, loop1_begin-2, loop1_begin-1, level );
				utility::vector1< String > abegop1 = am_.get_symbols( pose, loop1_begin, loop1_end, level );
				utility::vector1< String > abegop1_a = am_.get_symbols( pose, loop1_end+1, loop1_end+2, level );

				Real rama1( 0.0 );
				for( Size ii=loop1_begin; ii<=loop1_end; ii++ ) {
					Real r1( pose.energies().residue_total_energies( ii )[ core::scoring::rama ] );
					Real r2( gly_pose.energies().residue_total_energies( ii )[ core::scoring::rama ] );
					if( r1 > r2 ) {
						rama1 += r2;
					} else {
						rama1 += r1;
					}
				}

				// get abego for loop2
				Size loop2_begin = helix->end() + 1;
				Size loop2_end = strand2->begin() - 1;
				utility::vector1< String > abegop2_b = am_.get_symbols( pose, loop2_begin-2, loop2_begin-1, level );
				utility::vector1< String > abegop2 = am_.get_symbols( pose, loop2_begin, loop2_end, level );
				utility::vector1< String > abegop2_a = am_.get_symbols( pose, loop2_end+1, loop2_end+2, level);


				Real rama2( 0.0 );
				for( Size ii=loop2_begin; ii<=loop2_end; ii++ ) {
					Real r1( pose.energies().residue_total_energies( ii )[ core::scoring::rama ] );
					Real r2( gly_pose.energies().residue_total_energies( ii )[ core::scoring::rama ] );
					if( r1 > r2 ) {
						rama2 += r2;
					} else {
						rama2 += r1;
					}
				}

				// make loop pose
				// Pose loop1_pose( pose ), loop2_pose( pose );
				// loop1_pose.conformation().delete_residue_range_slow( helix->begin() + 1 , pose.total_residue() );
				// loop1_pose.conformation().delete_residue_range_slow( 1, strand1->end() - 1 );
				// make_poly_gly.apply( loop1_pose );
				// (*scorefxn_fullatom_)( loop1_pose );
				// Real loop1E( loop1_pose.energies().total_energies()[ core::scoring::fa_rep ] );

				Real loop1E( 0.0 );
				for( Size ii=loop1_begin-1; ii<=loop1_end+5; ii++ ) {
					loop1E += gly_pose.energies().residue_total_energies( ii )[ core::scoring::fa_rep ];
				}

				// loop2_pose.conformation().delete_residue_range_slow( strand2->begin() + 1, pose.total_residue() );
				// loop2_pose.conformation().delete_residue_range_slow( 1, helix->end() - 1 );
				// make_poly_gly.apply( loop2_pose );
				// (*scorefxn_fullatom_)( loop2_pose );
				// Real loop2E( loop2_pose.energies().total_energies()[ core::scoring::fa_rep ] );

				Real loop2E( 0.0 );
				for( Size ii=loop2_begin-5; ii<=loop2_end+1; ii++ ) {
					loop2E += gly_pose.energies().residue_total_energies( ii )[ core::scoring::fa_rep ];
				}


				using core::chemical::oneletter_code_from_aa;

				Size broken_hydrogen = protocols::fldsgn::topology::check_kink_helix( pose, helix->begin()-1, helix->end()-5 );
				utility::vector1< core::scoring::hbonds::HBond > lp_hbonds = protocols::fldsgn::topology::check_internal_hbonds( pose, strand1->end(), helix->begin() );

				std::ostringstream name;
				name << pose.pdb_info()->number( strand1->begin() ) << "-" << pose.pdb_info()->number( strand1->end() ) << ";"
					 << pose.pdb_info()->number( helix->begin() ) << "-" << pose.pdb_info()->number( helix->end() ) << ";"
					 << pose.pdb_info()->number( strand2->begin() ) << "-" << pose.pdb_info()->number( strand2->end() );

				String id( "A00000" );
				std::ostringstream num;
				num << ++num_bab_;
				id.replace( 6 - num.str().length(), num.str().length(), num.str() );

				output_ << A( 6, id ) << " " << LJ( 20, me ) << A( 5, inout1 ) << " " << A( 5, inout2 ) << " "
						<< LJ( 27, name.str() ) << " " <<  LJ( 8, bab.name() ) << " "
						<< I( 2, cross_over ) << "  "
						<< F( 5, 2, bab.hsheet_dist() ) << " " << F( 7, 2, bab.hs_angle() ) << " " << F( 7, 2, bab.hsheet_elev_angle() ) << " "
						<< F( 6, 1, dsasa ) << " "
						<< I( 2, broken_hydrogen ) << " " << I( 2, lp_hbonds.size() ) << " " << A( 3, bab.helix_cycle_as_string() ) << " "
						<< F( 5, 2, helix->bend() ) << " " << F( 5, 2, rama1 ) << " " << F( 5, 2, rama2 ) << " "
				        << F( 6, 2, loop1E ) << " " << F( 6, 2, loop2E ) << " "
						<< I( 3, strand1->length() ) << " " << I( 3, loop1_length ) << " "
						<< I( 3, helix->length() )   << " " << I( 3, loop2_length ) << " "
						<< I( 3, strand2->length() ) << " ";

				output2_ << A(6, id ) << " " << A( 20, me ) << " ";
				// write abego for loop1
				name.str("");
				output_ << abegop1_b[ 1 ] << abegop1_b[ 2 ] << " ";
				for( Size i=1; i<=abegop1.size(); i++ ) {
					name << abegop1[ i ];
				}
				if( name.str() == "" ) name << ".";
				output_ << LJ( 10, name.str() ) << " ";
				output_ << abegop1_a[ 1 ] << abegop1_a[ 2 ] << "  ";

				// write phi psi for loop1
				output2_ << LJ( 7, name.str() ) << " ";
				for( Size i=1; i<=loop1_length; i++ ) {
					output2_ <<  F( 7, 2, pose.phi( strand1->end() + i ) ) << " " << F( 7, 2, pose.psi( strand1->end() + i ) ) << " ";
				}
				output2_ << " ";

				// write sequence for loop1
				name.str("");
				output_ << oneletter_code_from_aa( pose.aa( loop1_begin-2 ) ) << oneletter_code_from_aa( pose.aa( loop1_begin-1 ) ) << " ";
				for( Size i=1; i<=abegop1.size(); i++ ) {
					name << oneletter_code_from_aa( pose.aa( loop1_begin+i-1 ) );
				}
				if( name.str() == "" ) name << ".";
				output_ << LJ( 10, name.str() ) << " ";
				output_ << oneletter_code_from_aa( pose.aa( loop1_end+1 ) ) << oneletter_code_from_aa( pose.aa( loop1_end+2 ) ) << "  ";

				// write abego for loop2
				name.str("");
				output_ << abegop2_b[ 1 ] << abegop2_b[ 2 ] << " ";
				for( Size i=1; i<=abegop2.size(); i++ ) {
					name << abegop2[ i ];
				}
				if( name.str() == "" ) name << ".";
				output_ << LJ( 10, name.str() ) << " ";
				output_ << abegop2_a[ 1 ] << abegop2_a[ 2 ] << "  ";

                // write phi psi for loop2
				output2_ << LJ( 7, name.str() ) << " ";
				for( Size i=1; i<=loop2_length; i++ ) {
					output2_ <<  F( 7, 2, pose.phi( helix->end() + i ) ) << " " << F( 7, 2, pose.psi( helix->end() + i ) ) << " ";
                }

				// write phi psi for 1st res of 2nd strand
				output2_ <<  F( 7, 2, pose.phi( strand2->begin() ) ) << " " << F( 7, 2, pose.psi( strand2->begin() ) ) << " ";
				output2_ << std::endl;

				// write sequence for loop2
				name.str("");
				output_ << oneletter_code_from_aa( pose.aa( loop2_begin-2 ) ) << oneletter_code_from_aa( pose.aa( loop2_begin-1 ) )<< " ";
				for( Size i=1; i<=abegop2.size(); i++ ) {
					name << oneletter_code_from_aa( pose.aa( loop2_begin+i-1 ) );
				}
				if( name.str() == "" ) name << ".";
				output_ << LJ( 10, name.str() ) << " ";
				output_ << oneletter_code_from_aa( pose.aa( loop2_end+1 ) ) << oneletter_code_from_aa( pose.aa( loop2_end+2 ) );

				// output spairset name
				output_ << " " << spairset->name() << std::endl;

			}
		}

	}

private: // data

	protocols::fldsgn::topology::SS_Info2_OP ssinfo_;
	std::ofstream output_;
    std::ofstream output2_;
	core::scoring::ScoreFunctionOP scorefxn_centroid_;
	core::scoring::ScoreFunctionOP scorefxn_fullatom_;
	core::util::ABEGOManager am_;
	bool ssinput_;
	Size num_bab_;
	Size nstruct_;
};

typedef utility::pointer::owning_ptr< PickBAB > PickBABOP;


int
main( int argc, char * argv [] )
{
	try{
	ThisApplication::register_options();

	// init
  devel::init(argc, argv);

	// mover
	protocols::moves::MoverOP protocol;
	protocol = new PickBAB();

	// run
	protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

//		using protocols::forge::build::Interval;
//		utility::vector1< Interval > intervals;
//		Interval interval( ssinfo_->strand( babset.bab_motif( 1 )->strand1() )->end()+1, ssinfo_->strand( babset.bab_motif( 1 )->strand2() )->begin()-1 );
//		intervals.push_back( interval );
//		Real dsasa2 = protocols::fldsgn::topology::calc_delta_sasa( pose, intervals, option[ pore_radius ]() );
//		std::cout << "OK" << dsasa << " " << dsasa2 << std::endl;

//		Real twist( 0.0 );
//		utility::vector1< Real > surfaces( 2, 0.0 );
//		if( sheet_set->sheets().size() == 1 ) {
//			surfaces = sheet_set->sheets()[ 1 ]->calc_sasa_bothsides( pose, ssinfo_, option[ pore_radius ]() );
//			twist = surfaces[ 1 ]/surfaces[ 2 ];
//		}
//		output_ << dsasa << " " << twist << " ";



//    String spairs_string = option[ spairs ]();
//		utility::vector1< String > sp( utility::string_split( spairs_string, '.' ) );
//		if( sp.size() == 2 ) {
//			if( spairs_string != "" && spairset->name_wo_rgstr() != spairs_string ) return;
//		} else if ( sp.size() == 3 ) {
//			if( spairs_string != "" && spairset->name() != spairs_string ) return;
//		}

//		    using namespace protocols::toolbox::pose_metric_calculators;
//				BuriedUnsatisfiedPolarsCalculator bu_calc( "default", "default" );
//				bu_calc.special_region( region );
//				basic::MetricValue< Size > val;
//				bu_calc.get( "special_region_bur_unsat_polars", val, pose );
//				std::cout << val.value() << std::endl;

