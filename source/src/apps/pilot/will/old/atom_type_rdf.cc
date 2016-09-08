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

#include <core/types.hh>
#include <devel/init.hh>

#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>

#include <utility/io/ozstream.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/jobdist/not_universal_main.hh>

std::string tag_from_pose( core::pose::Pose & pose ) {
	std::string tag( "empty_tag" );
	if ( pose.data().has( basic::JOBDIST_OUTPUT_TAG ) ) {
		tag =
			static_cast< basic::CacheableString const & >
			( pose.data().get( basic::JOBDIST_OUTPUT_TAG ) ).str();
	}
	return tag;
}

class	RDF_Mover : public protocols::moves::Mover {

public:

	RDF_Mover( core::Size natypes = 26, utility::vector1<bool> sep_ss=utility::vector1<bool>(26,false) )
	 : max_atype_(natypes) , sep_ss_(sep_ss)
	{
		using core::Size;
		using utility::vector1;
		using std::string;
		using std::map;
		hist_.resize(natypes);
		for( Size i = 1; i <= max_atype_; ++i ) {
			vector1<string> SSi;
			SSi.push_back("");
			if( sep_ss_[i] ) {	SSi.push_back("E"); SSi.push_back("H"); SSi.push_back("L");	}
			for( vector1<string>::iterator ssi = SSi.begin(); ssi != SSi.end(); ++ssi ) {
				hist_[i][*ssi] = vector1< map< string, vector1<Size> > >();
				hist_[i][*ssi].resize(natypes);
				for( Size j = 1; j <= max_atype_; ++j ) {
					vector1<string> SSj;
					SSj.push_back("");
					if( sep_ss_[j] ) {	SSj.push_back("E"); SSj.push_back("H"); SSj.push_back("L");	}
					for( vector1<string>::iterator ssj = SSj.begin(); ssj != SSj.end(); ++ssj ) {
						hist_[i][*ssi][j][*ssj].resize(200,0);
					}
				}
			}
		}
	}

	utility::vector1<core::Size> const &
	hist(
		core::Size atype1, std::string ss1,
		core::Size atype2, std::string ss2
	) const
	{
		return hist_[atype1].find(ss1)->second[atype2].find(ss2)->second;
	}

	utility::vector1<core::Size> &
	hist(
		core::Size atype1, std::string ss1,
		core::Size atype2, std::string ss2
	) {
		return hist_[atype1][ss1][atype2][ss2];
	}

	void
	apply(
		core::pose::Pose & pose
	) {
		using namespace std;
		using namespace core;
		using core::Size;
		using namespace scoring::etable::count_pair;

		std::cout << "counting pairs for " << tag_from_pose(pose) << std::endl;

		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose( pose );

		if( 0==type_names_.size() ) {
			std::cout << "initializing atom types from " << tag_from_pose(pose) << " (all must match!) " << std::endl;
			core::chemical::AtomTypeSet const & atset( pose.residue(1).atom_type_set() );
			for( Size i = 1; i <= atset.n_atomtypes(); ++i ) {
				type_names_.push_back(atset[i].name());
			}
		} else {
			core::chemical::AtomTypeSet const & atset( pose.residue(1).atom_type_set() );
			for( Size i = 1; i <= atset.n_atomtypes(); ++i ) {
				if( type_names_[i] != atset[i].name() ) {
					std::cerr << "RDF_Hist: atom type sets don't match!!!" << std::endl;
				}
			}
		}

		for( Size ir = 1; ir <= pose.size(); ++ir ) {
			for( Size ia = 1; ia <= pose.residue(ir).natoms(); ++ia ) {
				numeric::xyzVector<Real> const xyz1( pose.residue(ir).xyz(ia) );
				Size atype1 = pose.residue(ir).atom(ia).type();
				if( atype1 > max_atype_ ) continue;
				string ss1 = "";
				if( sep_ss_[atype1] ) ss1 = pose.secstruct(ir);
				assert(""==ss1|"E"==ss1||"H"==ss1||"L"==ss1);
				for( Size jr = ir; jr <= pose.size(); ++jr ) {
					CountPairFunctionOP cp = CountPairFactory::create_count_pair_function( pose.residue(ir), pose.residue(jr), CP_CROSSOVER_4 );
					for( Size ja = 1; ja <= pose.residue(jr).natoms(); ++ja ) {
						if( ir==jr && ia >= ja ) continue; // no dups & upper tri
						if( ir==jr ) continue;
						numeric::xyzVector<Real> const xyz2( pose.residue(jr).xyz(ja) );
						Size atype2 = pose.residue(jr).atom(ja).type();
						if( atype2 > max_atype_ ) continue;
						Real WT = 1.0; cp->count(ia,ja,WT);	if( WT != 1.0 ) continue;
						string ss2 = "";
						if( sep_ss_[atype2] ) ss2 = pose.secstruct(jr);
						assert(""==ss2|"E"==ss2||"H"==ss2||"L"==ss2);
						core::Real d2 = xyz1.distance_squared(xyz2);
						if( d2 > 100.0 ) continue;
						core::Real d = sqrt(d2);
						Size bin = (Size)ceil(d*20.0);
						++hist(atype1,ss1,atype2,ss2)[bin];
						++hist(atype2,ss2,atype1,ss1)[bin];
						if( ss1 != "" ) ++hist(atype1,"" ,atype2,ss2)[bin];
						if( ss1 != "" ) ++hist(atype2,ss2,atype1,"" )[bin];
						if( ss2 != "" ) ++hist(atype1,ss1,atype2,"" )[bin];
						if( ss2 != "" ) ++hist(atype2,"" ,atype1,ss1)[bin];
						if( ss1 != "" && ss2 != "" ) ++hist(atype1,"",atype2,"")[bin];
						if( ss1 != "" && ss2 != "" ) ++hist(atype2,"",atype1,"")[bin];
					}
				}
			}
		}
	}

	core::Size const max_atype() const {
		return max_atype_;
	}

	std::string const atom_type_name( core::Size i ) const {
		return type_names_[i];
	}

	void
	write_R_data_file( std::ostream & out ) {
		using core::Size;
		using utility::vector1;
		using std::string;
		if( max_atype_ > type_names_.size() ) {
			out << "NO DATA IN HISTOGRAM!!!!!" << std::endl;
			return;
		}
		for( Size i = 1; i <= max_atype_; ++i ) {
			vector1<string> SSi;
			SSi.push_back("");
			if( sep_ss_[i] ) {	SSi.push_back("E"); SSi.push_back("H"); SSi.push_back("L");	}
			for( vector1<string>::iterator ssi = SSi.begin(); ssi != SSi.end(); ++ssi ) {
				string lab1 = atom_type_name(i) + *ssi;
				for( Size j = i; j <= max_atype_; ++j ) {
					vector1<string> SSj;
					SSj.push_back("");
					if( sep_ss_[j] ) {	SSj.push_back("E"); SSj.push_back("H"); SSj.push_back("L");	}
					for( vector1<string>::iterator ssj = SSj.begin(); ssj != SSj.end(); ++ssj ) {
						utility::vector1<Size> const & counts( hist(i,*ssi,j,*ssj) );
						for( Size k = 1; k <= 200; ++k ) {
							string lab2 = atom_type_name(j) + *ssj;
							out << lab1 << " " << lab2;
							out << " " << ((core::Real)k)/20.0-0.025 << " " << counts[k] << std::endl;
						}
					}
				}
			}
		}
	}

private:

	utility::vector1< std::map< std::string, utility::vector1< std::map< std::string, utility::vector1<core::Size> > > > > hist_;
	core::Size max_atype_;
	utility::vector1<std::string> type_names_;
	utility::vector1<bool> sep_ss_;

};


int
main (int argc, char *argv[])
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility;

	devel::init( argc, argv );

	utility::vector1<bool> sep_ss(26,false);

	if( option[ in::rdf::sep_bb_ss ]() ) {
		for( core::Size i = 18; i <= 21; i++ ) sep_ss[i] = true;
		sep_ss[26] = true;
	}
	RDF_Mover hist(26,sep_ss);

	protocols::jobdist::not_universal_main( hist );

	std::string outname = "OUT_FILE_O_UNSPECIFIED.rdf";
	if( option[ out::file::o ].user() ) outname = option[ out::file::o ]();
	std::cout << "write RDF to " << outname << std::endl;
	utility::io::ozstream out(outname.c_str());
	out << "data generated by minirosetta atom_type_rdf.cc" << std::endl;
	hist.write_R_data_file(out);
	out.close();

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
