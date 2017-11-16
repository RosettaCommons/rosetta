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
/// @author Nobuyasu Koga

// Package headers
#include <protocols/parser/BluePrint.hh>
#include <protocols/fldsgn/topology/SheetFoldTypeManager.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/Sheet.hh>
#include <protocols/fldsgn/topology/BetaAlphaBetaMotif.hh>

// Project headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <core/id/AtomID_Map.hh>
#include <core/scoring/sasa.hh>
#include <core/pose/util.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/excn/Exceptions.hh>

#include <ObjexxFCL/format.hh>

// C++ header
#include <fstream>
#include <map>
#include <cmath>


using basic::Error;
using basic::Warning;

using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

typedef std::string String;


static basic::Tracer TR( "foldptn" );

class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, output )
OPT_KEY( File, output2 )
OPT_KEY( File, blue )
OPT_KEY( Boolean, bab )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( output, "output file name", "foldptn.log" );
	NEW_OPT( output2, "output file name", "loop_rama.dat" );
	NEW_OPT( blue, "blueprint", "" );
	NEW_OPT( bab, "if the data is bab motif ", false );
}


/// local mover for testing purposes
class Foldptn : public protocols::moves::Mover {
public: // typedef


	typedef protocols::parser::BluePrint BluePrint;
	typedef protocols::parser::BluePrintOP BluePrintOP;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetOP StrandPairingSetOP;
	typedef protocols::fldsgn::topology::SheetSet SheetSet;
	typedef protocols::fldsgn::topology::SheetSetOP SheetSetOP;
	typedef protocols::fldsgn::topology::SheetFoldType SheetFoldType;
	typedef protocols::fldsgn::topology::SheetFoldTypeManager SheetFoldTypeManager;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotif BetaAlphaBetaMotif;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotifs BetaAlphaBetaMotifs;
	typedef protocols::fldsgn::topology::BetaAlphaBetaMotifSet BetaAlphaBetaMotifSet;


public: // constructor/deconstructor


	Foldptn()
	{
		// read blueprint
		ssinput_ = false;
		if ( option[ blue ].user() ) {
			BluePrintOP blueprint( new BluePrint( option[ blue ]() ) );
			ssinfo_ = protocols::fldsgn::topology::SS_Info2_OP( new SS_Info2( blueprint->secstruct() ) );
			ssinput_ = true;
		} else {
			ssinfo_ = protocols::fldsgn::topology::SS_Info2_OP( new SS_Info2 );
		}

		// set scorefxn
		scorefxn_fullatom_ = get_score_function();
		scorefxn_centroid_ = core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen" );

		// output file
		std::ostringstream filename;
		filename <<  option[ output ]();
		output_.open( filename.str().c_str() ,std::ios::out );

		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( native_, option[ basic::options::OptionKeys::in::file::native ] , core::import_pose::PDB_file);
		}

		std::ostringstream header;
		header << "filename " << "foldname " << "strand_pairs " << "fullE " << "cenE ";

		if ( option[ in::file::native ].user() ) {
			header << "rms";
		}

		output_ << header.str() << std::endl;

	}

	virtual ~Foldptn(){};

	virtual std::string get_name() const { return ""; }

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return protocols::moves::MoverOP( new Foldptn );
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

		core::scoring::dssp::Dssp dssp( pose );
		if ( ssinput_ ) {
			ssinfo_->set_SSorient( pose );
		} else {
			ssinfo_->initialize( pose, dssp.get_dssp_secstruct() );
		}

		StrandPairingSetOP spairset( new StrandPairingSet( protocols::fldsgn::topology::calc_strand_pairing_set( pose, ssinfo_ ) ) );

		SheetFoldType sfold = sm_.foldtype_from_spairs( spairset->name_wo_rgstr() );

		// calc sheet
		SheetSetOP sheet_set( new SheetSet( ssinfo_, spairset ) );
		sheet_set->calc_geometry( ssinfo_ );

		// calc bab
		bool lefty = false;
		BetaAlphaBetaMotifSet babset( ssinfo_, sheet_set );
		for ( BetaAlphaBetaMotifs::const_iterator iter = babset.bab_motifs().begin(),
				iter_end = babset.bab_motifs().end(); iter != iter_end; ++iter ) {
			BetaAlphaBetaMotif bab( **iter );
			bab.calc_geometry( ssinfo_, sheet_set );
			if ( bab.is_lefthanded() ) lefty = true;
		}

		Size bab_inout = 0;
		if ( option[ bab ]() && spairset->num_strands() == 2 && babset.size() >= 1 ) {
			bab_inout = babset.bab_motif( 1 )->calc_inout( ssinfo_, ssinfo_->strand( 1 )->end() );
		}

		// calc score
		( *scorefxn_fullatom_ )( pose );
		Real fullE = pose.energies().total_energies()[ core::scoring::total_score ];

		// switch to centroid model
		Pose copy_pose( pose );
		switch_to_residue_type_set( copy_pose, core::chemical::CENTROID );
		( *scorefxn_centroid_ )( copy_pose );
		Real cenE = copy_pose.energies().total_energies()[ core::scoring::total_score ];

		// calc rmsd
		Real rms( -1.0 );
		if ( option[ in::file::native ].user() ) {
			rms = core::scoring::CA_rmsd( pose, native_ );
		}

		String foldname;
		if ( lefty ) {
			foldname = "LeftHanded";
		} else {
			foldname = sm_.name_from_foldtype( sfold );
		}

		// output
		String spair_string;
		if ( foldname == "NO_STRANDS" ) {
			spair_string = "-";
		} else {
			spair_string = spairset->name();
		}


		String inout("");
		if ( bab_inout != 0 ) {
			if ( bab_inout == 1 ) {
				inout = "UP";
			} else {
				inout = "DOWN";
			}

			output_ << me << " " << foldname << "_" << inout << " " << spair_string << " "
				<< fullE << " " << cenE << " ";

		} else {

			output_ << me << " " << foldname << " " << spair_string << " "
				<< fullE << " " << cenE << " ";
		}

		if ( rms > 0.0 ) {
			output_ << " " << rms << std::endl;
		} else {
			output_ << " " << std::endl;
		}

	}


private: // data


	SheetFoldTypeManager sm_;
	protocols::fldsgn::topology::SS_Info2_OP ssinfo_;
	core::scoring::ScoreFunctionOP scorefxn_fullatom_;
	core::scoring::ScoreFunctionOP scorefxn_centroid_;
	std::ofstream output_;
	std::ofstream output2_;
	bool ssinput_;
	core::pose::Pose native_;

};

typedef utility::pointer::shared_ptr< Foldptn > FoldptnOP;


int
main( int argc, char * argv [] )
{
	try{
		ThisApplication::register_options();

		// init
		devel::init(argc, argv);

		// mover
		protocols::moves::MoverOP protocol;
		protocol = protocols::moves::MoverOP( new Foldptn() );

		// run
		protocols::jd2::JobDistributor::get_instance()->go( protocol );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
