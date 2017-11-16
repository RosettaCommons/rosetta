// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
#include <core/types.hh>
#include <devel/init.hh>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iomanip>

#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>


#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/rms_util.hh>

#include <utility/io/izstream.hh>

#include <utility/excn/Exceptions.hh>

#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>

OPT_1GRP_KEY( File, out, pdb )
OPT_1GRP_KEY( Real, tnm, scaling )
OPT_1GRP_KEY( Integer, tnm, n_modes )
OPT_1GRP_KEY( File, tnm, input_modes )

static basic::Tracer TR( "protocols::TNMSampler" );

using namespace std;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;

void register_options() {

	OPT(in::file::s);
	OPT(in::file::native);
	OPT(in::file::residue_type_set);
	NEW_OPT( out::pdb, "provides a file name for the output PDB file","out_pose.pdb" );
	NEW_OPT( tnm::input_modes, "provides input normalo modes in dihedral space","");
	NEW_OPT( tnm::scaling, "a factor used to multiply a normal mode vector",0.1);
	NEW_OPT( tnm::n_modes, "how many modes should be used",10);
}

void tokenize_string(const std::string& sentence, utility::vector1<std::string>& words) {

	std::istringstream iss(sentence);
	words.clear();
	std::copy(std::istream_iterator<std::string>(iss),
		std::istream_iterator<std::string>(),
		std::back_inserter<utility::vector1<std::string> >(words));
}

class Mode {
public:
	Mode(core::Size id,core::Size n_res, core::Real e_value) : _id(id), _n_res(n_res), evalue(e_value),mode(2*n_res) { std::cerr << id << " "<<n_res<<" "<<evalue<<"\n"; } // No omega so far!
	core::Real e_value() { return evalue;}
	core::Real n_res() { return _n_res;}
	void write_modes(std::ostream& out);
	void insert_vector(utility::vector1<core::Real>& src);
	void set( core::Size dof_id, core::Real  value) { mode[dof_id] = value; }
	core::Size size() const {return mode.size();}
	inline core::Real phi(const core::Size i_pos) { if ( i_pos>1 ) return mode[i_pos*2-2];else return 0.0;}
	inline core::Real psi(const core::Size i_pos) { if ( i_pos<=_n_res ) return mode[i_pos*2-1];else return 0.0;}

private:
	core::Size _id;
	core::Size _n_res;
	core::Real evalue;
	utility::vector1<core::Real> mode;
};

class TNM {
public :
	TNM(const std::string& fname) {load_normal_modes(fname);}
	void write_modes(std::ostream& out);
	void apply(core::Size, core:: Real, core::pose::Pose& pose);
private:
	utility::vector1<Mode> modes;
	void load_normal_modes(const std::string& fname);
};

void TNM::write_modes(std::ostream& out) {

	// Print evalues as the first line, make it a comment
	out<<"# evalues:";
	for ( core::Size i_mode=1; i_mode<=modes.size(); ++i_mode ) {
		out<<" "<< setw(10) << setprecision(7)<< modes[i_mode].e_value();
	}
	out<<"\n";

	// Print Psi for the very first residue
	out<<"         1";
	for ( core::Size i_mode=1; i_mode<=modes.size(); ++i_mode ) {
		out<<" "<< setw(10) << setprecision(7)<< modes[i_mode].psi(1);
	}
	out<<" Psi1\n";

	// Print Phi
	for ( core::Size i_pos=2; i_pos<modes[1].n_res(); ++i_pos ) {
		out << setw(10)  <<  i_pos;
		for ( core::Size i_mode=1; i_mode<=modes.size(); ++i_mode ) {
			out<<" "<< setw(10) << setprecision(7)<< modes[i_mode].phi(i_pos);
		}
		out <<" phi"<<i_pos<<endl;

		out<<setw(10)  <<  i_pos;
		for ( core::Size i_mode=1; i_mode<=modes.size(); ++i_mode ) {
			out<<" "<< setw(10) << setprecision(7)<< modes[i_mode].psi(i_pos);
		}
		out<<" psi"<<i_pos<<endl;
	}
}

void TNM::apply(core::Size which_mode, core::Real lambda, core::pose::Pose& pose) {

	for ( core::Size i=2; i<pose.size(); ++i ) {
		pose.set_phi( i, pose.phi(i) + lambda * modes[which_mode].phi(i) );
		pose.set_psi( i, pose.psi(i) + lambda * modes[which_mode].psi(i) );
	}
}

void TNM::load_normal_modes(const std::string& fname) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	utility::io::izstream data(fname);

	std::string line;
	utility::vector1< core::Real > evalues;
	utility::vector1< std::string > tokens;
	utility::vector1< utility::vector1<std::string> > lines;

	std::cerr <<"Reading TNM from "<<fname<<"\n";

	while ( !data.fail() ) {
		char c = data.peek();
		if ( c == '#' || c == '\n' ) {
			getline(data, line); //comment, junk or e-values; the latter are quite important
			tokens.clear();
			tokenize_string(line,tokens);
			if ( tokens[2] == "e_values:" )  {
				for ( core::Size i=3; i<=tokens.size(); ++i ) {
					evalues.push_back( atof( tokens[i].c_str() ) );
				}
			}
			continue;
		}
		getline(data, line);
		utility::vector1< std::string > tokenss;
		tokenize_string(line,tokenss);
		lines.push_back(tokenss);
	}
	std::cerr <<"at "<<lines.size()<<" "<<lines[1].size()<<"\n";
	core::Size n_modes = lines[1].size()-2;
	core::Size n_res = lines.size()/2+1;
	std::cerr <<"n_res, n_modes: "<<n_res<<" "<<n_modes<<"\n";
	for ( core::Size i_modes=1; i_modes<=n_modes; ++i_modes ) {
		Mode m(i_modes,n_res,evalues[i_modes]);
		for ( core::Size i=1; i <=lines.size(); ++i ) {
			m.set(i, atof( lines[i][i_modes+1].c_str() ) );
		}
		modes.push_back(m);
	}
}

int main(int argc, char * argv[]) {

	try {
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using std::string;

		register_options();
		devel::init(argc, argv);

		core::pose::Pose extended_pose;
		std::string out_file_name = (option[out::pdb]());

		core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
			option[ in::file::residue_type_set ]()
		);
		string sequence;

		if ( !option[in::file::s].user() ) {
			std::cerr << "Please provide input structure with -in::file::s option" << std::endl;
			return 0;
		}
		core::pose::PoseOP tmp_pose(new core::pose::Pose);
		std::string fn = option[in::file::s](1);
		core::import_pose::pose_from_file(*tmp_pose, fn, core::import_pose::PDB_file);
		std::cerr <<"Input structure with "<<tmp_pose->size()<<" residues\n";

		core::pose::PoseOP reference_pose(new core::pose::Pose);
		std::string fnn = option[in::file::native];
		core::import_pose::pose_from_file(*reference_pose, fnn, core::import_pose::PDB_file);


		if ( !option[tnm::input_modes].user() ) {
			std::cerr << "Please provide input normal modes with -tnm::input_modes option" << std::endl;
			return 0;
		}

		TNM tnm(option[tnm::input_modes]().name());
		//        tnm.write_modes(std::cout);

		core::scoring::ScoreFunction scorefxn = *( core::scoring::get_score_function(false) );

		//        core::Size i_mode = 1;
		cout << scorefxn(*tmp_pose) << endl;

		int n_steps = 100;
		core::Real step = option[tnm::scaling]();
		core::Size n_modes = option[tnm::n_modes]();
		for ( core::Size i_mode=1; i_mode<=n_modes; ++i_mode ) {
			core::import_pose::pose_from_file(*tmp_pose, fn, core::import_pose::PDB_file);
			cout << "Loaded pose from: "<<fn<<", en = "<<" "<<scorefxn(*tmp_pose) << endl;

			for ( int i=-n_steps; i<=n_steps; ++i ) {
				core::Real lambda = i*step;
				tnm.apply(i_mode,lambda,*tmp_pose);
				core::Real rms = core::scoring::CA_rmsd(*tmp_pose, *reference_pose);
				cout << setw(2) << i_mode<<" "<<F(6,3,lambda)<<" "<<F(8,3,scorefxn(*tmp_pose)) << " " << F(8,3,rms)<<endl;
				tnm.apply(i_mode,-lambda,*tmp_pose);
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
