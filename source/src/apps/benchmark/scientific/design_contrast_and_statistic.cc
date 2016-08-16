// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/pilot/yiliu/sqc_test.cc
///
/// @brief
/// @author Yi Liu


// C++ headers
#include <fstream>
#include <iostream>
#include <cstdlib>

// Unit Headers
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>


// Project Headers
#include <core/io/sequence_comparation/DesignContrast.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace basic;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using utility::vector1;
using utility::file::FileName;
using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.yiliu.DC" );

///////////////////////////////////////////////////////////////////////////////
// YAML helper function
std::ostream & writeYamlValue(std::ostream & S, std::string name, core::Real value)
{
	S << "'" << name << "' : " << value << ", ";
	return S;
}

///////////////////////////////////////////////////////////////////////////////
// YAML helper function
std::ostream & writeYamlValue(std::ostream & S, std::string name, bool value)
{
	std::string sv;
	if ( value ) sv = "True";
	else sv = "False";

	S << "'" << name << "' : " << sv << ", ";
	return S;
}


///////////////////////////////////////////////////////////////////////////////


int isHydroNonPolar(std::string const & resname ){
	if ( resname == "VAL" || resname == "ILE" || resname == "LEU" ||
			resname == "MET" || resname == "PHE" || resname == "GLY" ||
			resname == "ALA" || resname == "PRO" ||
			resname == "TRP" || resname == "TYR" ) {
		return 1;
	} else return 0;
}


int isPolarUncharge(std::string const & resname ){
	if ( resname == "SER" || resname == "THR" || resname == "GLN" || resname == "ASN" ) {
		return 1;
	} else return 0;
}

int isNegative(std::string const & resname ){
	if ( resname == "ASP" || resname == "GLU" ) {
		return 1;
	} else return 0;
}

int isPositive(std::string const & resname ){
	if ( resname == "ARG" || resname == "LYS" || resname == "HIS" ) {
		return 1;
	} else return 0;
}


int isBoundary(Size const & neighbor ) {
	if ( (neighbor>13)&&(neighbor<=18) ) {
		return 1;
	} else return 0;
}

int isBuried(Size const & neighbor ) {
	if ( (neighbor>18) ) {
		return 1;
	} else return 0;
}

// Statistics
void statistics( std::string filename ) {
	//open the redesign file
	TR << "debug file name: " << filename << std::endl;
	std::ifstream redesign_file(filename.c_str());
	std::string line;
	std::string nres, dres;
	unsigned int isTer;
	std::string pdbCode="";
	Real BoundaryNum(0), SurfaceNum(0), BuriedNum(0);
	Real hydroNonPolarTotal(0), polarUnchargeTotal(0), polarChargeTotal(0), basicTotal(0), cysTotal(0);
	Size neighbor(0),totalNo(0),lineNo(0);

	// Final calculation variables
	//Size actualTotal(0), RedesignhydroNonPolarTotal(0), RedesignpolarUnchargeTotal(0), RedesignpolarChargeTotal(0), RedesignbasicTotal(0);
	Real sum(0), basicSum(0), hydroNonPolarSum(0), polarUnchargeSum(0), polarChargeSum(0), cysSum(0);
	Real idTotal(0.0), idHydroNonPolar(0.0), idpolarUncharge(0.0), idNegative(0.0), idPositive(0.0);
	//std::map<std::string,Size> NativeRes, DesignRes, NativeBoundary, DesignBoundary, NativeBuried, DesignBuried, NativeSurface,DesignSurface, unchangedRes;
	//std::string aas[] = {"ALA", "VAL", "ILE", "LEU", "PHE",
	//           "MET", "GLY", "PRO", "TRP", "TYR",
	//           "THR", "GLN", "ASN", "SER", "GLU",
	//           "ASP", "ARG", "LYS", "HIS", "TRP"};
	//for (Size i = 0; i < 20; i++) {
	// NativeRes.insert(std::make_pair(aas[i],0));
	// DesignRes.insert(std::make_pair(aas[i],0));
	// NativeBoundary.insert(std::make_pair(aas[i],0));
	// DesignBoundary.insert(std::make_pair(aas[i],0));
	// NativeBuried.insert(std::make_pair(aas[i],0));
	// DesignBuried.insert(std::make_pair(aas[i],0));
	// NativeSurface.insert(std::make_pair(aas[i],0));
	// DesignSurface.insert(std::make_pair(aas[i],0));
	// unchangedRes.insert(std::make_pair(aas[i],0));
	//}

	while ( getline(redesign_file,line) ) {
		nres = line.substr(13,3);
		dres = line.substr(17,3);
		neighbor = std::atol(line.substr(24,2).c_str());
		if ( line.substr(0,4) != pdbCode ) {
			//dont count those N terminus MET
			if ( ( nres == "MET") || (dres == "MET") ) {
				isTer = 1;
			} else {
				isTer = 0;
			}
		} else {
			isTer = 0;
		}
		pdbCode = line.substr(0,4);
		if ( isTer == 0 ) {
			totalNo++;
			TR <<"hit" <<std::endl;
			//NativeRes[nres]++;
			//DesignRes[dres]++;
		}
		//NatBoundary[i] keeps track of the numbers in boundary for each residue in native state
		//DesBoundary[i] keeps track of the numbers in boundary for each residue in redesign state
		if ( isBoundary(neighbor) ) {
			//NativeBoundary[nres]++;
			//DesignBoundary[dres]++;
			BoundaryNum++;
		} else if ( isBuried(neighbor) ) {
			//NatBuried[i] keeps track of the numbers in buried for each residue in native state
			//DesBuried[i] keeps track of the numbers in buried for each residues in redesign state
			//NativeBuried[nres]++;
			//DesignBuried[dres]++;
			BuriedNum++;
		} else {
			//NatSurface[i] keeps track of the numbers in surface for each residue in native state
			//DesSurface[i] keeps track of the numbers in surface for each residue in redesign state
			//NativeSurface[nres]++;
			//DesignSurface[dres]++;
			SurfaceNum++;
		}
		//matrix[a,b] keeps all residues
		//--- matrix[a,b]++
		//group[i,j] is grouped by four categories
		if ( isHydroNonPolar(nres) ) {
			hydroNonPolarTotal++;
		} else if ( isPolarUncharge(nres) ) {
			polarUnchargeTotal++;
		} else if ( isNegative(nres) ) {
			polarChargeTotal++;
		} else if ( isPositive(nres) ) {
			basicTotal++;
		} else cysTotal++;
		if ( nres == dres ) {
			//unchangedRes[nres]++;
			sum++;
			if ( isHydroNonPolar(nres) ) {
				hydroNonPolarSum++;
			} else if ( isPolarUncharge(nres) ) {
				polarUnchargeSum++;
			} else if ( isNegative(nres) ) {
				polarChargeSum++;
			} else if ( isPositive(nres) ) {
				basicSum++;
			} else cysSum++;
		}
		++lineNo;
	}
	// Final calculations
	//actualTotal = totalNo-cysTotal;
	TR << sum << totalNo << std::endl;
	idTotal = (sum-cysSum)/(totalNo-cysTotal);
	idHydroNonPolar = hydroNonPolarSum/hydroNonPolarTotal;
	idpolarUncharge = polarUnchargeSum/polarUnchargeTotal;
	idNegative = polarChargeSum/polarChargeTotal;
	idPositive = basicSum/basicTotal;

	std::string results_fname( ".results.log" );
	std::ofstream staResult( results_fname.c_str() );
	if ( !staResult ) {
		TR.Error << "Can not open file " << results_fname;
	} else {

		char idT[100], idH[100], idPU[100],idN[100], idP[100];
		sprintf(idT,"%.1f%%",idTotal*100);
		sprintf(idH,"%.1f%%",idHydroNonPolar*100);
		sprintf(idPU,"%.1f%%",idpolarUncharge*100);
		sprintf(idN,"%.1f%%",idNegative*100);
		sprintf(idP,"%.1f%%",idPositive*100);
		staResult
			<< "the identity for all positions is: " << idT << "    [should be more than 0.32]" << " \n"
			<< "hydrophobic non-polar amino acids (VAL,ILE,LEU,MET,PHE,GLY,ALA,PRO,TRP,TYR): " << idH << " \n"
			<< "polar uncharged amino acid except CYS (SER,THR,ASN,GLN): " << idPU << " \n"
			<< "negative charged amino acid(ASP,GLU): " << idN << " \n"
			<< "positive charged amino acid(ARG,LYS,HIS): " << idP << " \n";
	}

	std::string yaml_fname( ".results.yaml" );
	std::ofstream yaml( yaml_fname.c_str() );
	if ( !yaml ) {
		TR.Error << "Can not open file " << yaml_fname;
	} else {
		yaml << "{ ";

		writeYamlValue(yaml, "IdentityForAllPositions", idTotal);
		writeYamlValue(yaml, "HydrophobicNon_PolarAminoAcids", idHydroNonPolar);
		writeYamlValue(yaml, "PolarUnchargedAmino", idpolarUncharge);
		writeYamlValue(yaml, "NegativeChargedAminoAcid", idNegative);
		writeYamlValue(yaml, "PositiveChargedAminoAcid", idPositive);

		writeYamlValue(yaml, "_isTestPassed", idTotal > 0.32 );

		yaml << "}\n";
	}

	TR << "the identity for all positions is: " << idTotal << "    [should be more than 0.32]" << std::endl;
	TR << "hydrophobic non-polar amino acids (VAL,ILE,LEU,MET,PHE,GLY,ALA,PRO,TRP,TYR): " << idHydroNonPolar << std::endl;
	TR << "polar uncharged amino acid except CYS (SER,THR,ASN,GLN): " << idpolarUncharge << std::endl;
	TR << "negative charged amino acid(ASP,GLU): " << idNegative << std::endl;
	TR << "positive charged amino acid(ARG,LYS,HIS): " << idPositive << std::endl;
}

int main( int argc, char * argv [] )
{
	try {

		using namespace core;
		using namespace core::io;
		devel::init(argc, argv);

		TR << "in the main" << std::endl;
		std::string out_path, redesign_name;
		vector1 <std::string> in_pdb_names, pdb_codes;
		sequence_comparation::DesignContrast dc;
		vector1 <pose::Pose> native_poses, decoy_poses;
		std::string sqc_file;

		dc.setNames();
		dc.setPdbCodes();
		in_pdb_names = dc.getPdbNames();
		pdb_codes = dc.getPdbCodes();
		// Initialize and use the mover
		using namespace core::pack::task;
		using namespace core::pack::task::operation;
		TaskFactoryOP main_task_factory( new TaskFactory );
		main_task_factory->push_back( TaskOperationCOP( new operation::InitializeFromCommandline ) );
		if ( option[ packing::resfile ].user() ) {
			main_task_factory->push_back( TaskOperationCOP( new operation::ReadResfile ) );
		}
		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();
		protocols::simple_moves::PackRotamersMoverOP pack_mover( new protocols::simple_moves::PackRotamersMover );
		pack_mover->task_factory( main_task_factory );
		pack_mover->score_function( score_fxn );
		TR << "before the loop" << std::endl;
		for ( Size i=1; i <= pdb_codes.size(); ++i ) {
			pose::Pose single_in_pose, single_out_pose;
			core::import_pose::pose_from_file(single_in_pose, in_pdb_names[i], core::import_pose::PDB_file);
			single_out_pose = single_in_pose;
			// Fixbb run
			pack_mover->apply(single_out_pose);
			dc.setNeighbors(single_in_pose);
			native_poses.push_back(single_in_pose);
			decoy_poses.push_back(single_out_pose);
		}

		if ( option [ out::file::design_contrast ].active() ) {
			//  std::string sqc_file = option [ out::file::design_contrast ].default_value();
			sqc_file = option [ out::file::design_contrast ].value();
			std::ofstream sqc;
			sqc.open(sqc_file.c_str());
			for ( Size j=1; j <=  decoy_poses.size(); ++j ) {
				dc.output_sqc_file(native_poses[j], decoy_poses[j], pdb_codes[j], sqc);
			}
			sqc.close();
			statistics(sqc_file);
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

