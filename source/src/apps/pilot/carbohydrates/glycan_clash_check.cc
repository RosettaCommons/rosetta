// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/carbohydrates/glycan_clash_check.cc
/// @brief A small app to calculate clashes between carbohydrate branches and other carbohydrates/chains.  Outputs info into result scorefile.  Does not echo structures.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// devel headers
#include <devel/init.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/selection.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.hh>

// protocol headers
// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/moves/Mover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <math.h>

static THREAD_LOCAL basic::Tracer TR("glycan_clash_check");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( carbohydrates::clash_check::glycan_branches );
	option.add_relevant( carbohydrates::clash_check::check_chains );
	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}

// Class Definitions //////////////////////////////////////////////////////////

/// @brief  Class to print out info specific to glycan poses.  Will be expanded as needed.
/// Hope to eventually add it as both a mover and a features reporter.
///  Needs to be completely rewritten.
class GlycanClashCheckMover : public protocols::moves::Mover {

	typedef std::map< std::string, utility::vector1< core::Size > >::const_iterator branch_it_type;


public:  // Standard methods
	/// @brief  Default constructor.
	GlycanClashCheckMover() : protocols::moves::Mover()
	{
		init();
	}

	/// @brief  Copy constructor.
	GlycanClashCheckMover( GlycanClashCheckMover const & src ) :
		Mover( src ),
		branches_(src.branches_),
		branch_point_resnums_(src.branch_point_resnums_),
		chains_(src.chains_),
		chain_nums_(src.chain_nums_),
		branch_residues_(src.branch_residues_),
		glycan_glycan_res_atomic_clashes_(src.glycan_glycan_res_atomic_clashes_),
		glycan_chain_res_atomic_clashes_(src.glycan_chain_res_atomic_clashes_),
		glycan_glycan_res_soft_atomic_clashes_(src.glycan_glycan_res_soft_atomic_clashes_),
		glycan_chain_res_soft_atomic_clashes_(src.glycan_chain_res_soft_atomic_clashes_),
		soft_clash_percent_(src.soft_clash_percent_),
		ignore_hydrogens_(src.ignore_hydrogens_),
		ignore_full_res_output_(src.ignore_full_res_output_),
		output_per_glycan_data_(src.output_per_glycan_data_),
		glycan_resnums_(src.glycan_resnums_),
		chain_resnums_(src.chain_resnums_),
		score_data_(src.score_data_),
		res_data_(src.res_data_),
		pymol_print_(src.pymol_print_),
		lj_n_(src.lj_n_)
	{

	}

	// Destructor
	virtual ~GlycanClashCheckMover() {}


public:  // Standard Rosetta methods



	/// @brief  Generate string representation of DockGlycansProtocol for debugging purposes.
	virtual
	void
	show( std::ostream & output=std::cout ) const
	{
		protocols::moves::Mover::show( output );  // name, type, tag
	}


	// Mover methods
	/// @brief  Return the name of the Mover.
	virtual
	std::string
	get_name() const
	{
		return type();
	}

	virtual
	protocols::moves::MoverOP
	clone() const
	{
		return protocols::moves::MoverOP( new GlycanClashCheckMover( *this ) );
	}

	virtual
	protocols::moves::MoverOP
	fresh_instance() const
	{
		return protocols::moves::MoverOP( new GlycanClashCheckMover );
	}

	virtual bool
	reinitialize_for_each_job() const {
		return true;
	}


	//////////////////////////////////// Actual Code ////////////////////////////////////////////////////



	/// @brief  Apply the corresponding protocol to <pose>.
	virtual
	void
	apply( core::pose::Pose & pose )
	{
		branch_residues_.clear();
		branch_point_resnums_.clear();
		chain_nums_.clear();
		res_data_.clear();
		score_data_.clear();

		using namespace core::pose::carbohydrates;

		//Setup branches
		for ( core::Size i = 1; i <= branches_.size(); ++i ) {
			core::Size resnum = core::pose::parse_resnum( branches_[ i ], pose );
			branch_point_resnums_.push_back(resnum);

			std::pair< utility::vector1< core::Size >, utility::vector1< core::Size > > res_and_tips;
			res_and_tips = get_carbohydrate_residues_and_tips_of_branch( pose, resnum );
			branch_residues_[ branches_[ i ] ] = res_and_tips.first; //These are the residue numbers of the branch
		}
		TR << "Branch point resnums: "<< utility::to_string( branch_point_resnums_ ) << std::endl;
		//Setup Chains
		for ( core::Size i = 1; i <= chains_.size(); ++i ) {
			utility::vector1< core::Size > chains = core::pose::get_chain_ids_from_chain( chains_[ i ], pose) ;
			chain_nums_.insert(chain_nums_.end(), chains.begin(), chains.end() );
		}

		//Setup Residues
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue( i ).is_carbohydrate() ) {
				glycan_resnums_.push_back( i );
			} else if ( chain_nums_.has_value( pose.residue( i ).chain() ) ) {
				chain_resnums_.push_back( i );
			} else {
				continue;
			}
		}

		//Measure distances, count.  Can be made much faster, I know.  I need this data yesterday.
		for ( branch_it_type it = branch_residues_.begin(); it != branch_residues_.end(); ++it ) {
			std::string branch_start = it->first;

			//Initialize data
			glycan_glycan_res_atomic_clashes_[     branch_start ].clear();
			glycan_chain_res_atomic_clashes_[      branch_start ].clear();
			glycan_glycan_res_soft_atomic_clashes_[ branch_start ].clear();
			glycan_chain_res_soft_atomic_clashes_[  branch_start ].clear();

			glycan_glycan_res_atomic_clashes_[     branch_start ].resize( pose.size(), 0 );
			glycan_chain_res_atomic_clashes_[      branch_start ].resize( pose.size(), 0 );
			glycan_glycan_res_soft_atomic_clashes_[ branch_start ].resize( pose.size(), 0 );
			glycan_chain_res_soft_atomic_clashes_[  branch_start ].resize( pose.size(), 0 );

			for ( core::Size index = 1; index <= it->second.size(); ++index ) {

				core::Size glycan_resnum = it->second[ index ];
				//TR << "Branch: " << branch_start << " Glycan Res: " << glycan_resnum << std::endl;
				calculate_atom_clashes("glycan", pose, branch_start, glycan_resnum, glycan_resnums_, glycan_glycan_res_atomic_clashes_, glycan_glycan_res_soft_atomic_clashes_);


				calculate_atom_clashes("chain", pose, branch_start, glycan_resnum,  chain_resnums_, glycan_chain_res_atomic_clashes_, glycan_chain_res_soft_atomic_clashes_);

			}
		}

		tabulate_individual_data( "glycan", "norm", pose, glycan_glycan_res_atomic_clashes_);
		tabulate_individual_data( "glycan", "soft", pose, glycan_glycan_res_soft_atomic_clashes_);

		tabulate_individual_data( "chain", "norm", pose, glycan_chain_res_atomic_clashes_);
		tabulate_individual_data( "chain", "soft", pose, glycan_chain_res_soft_atomic_clashes_);

		//Sum glycan and chain data.
		//std::cout << "Summing data " << std::endl;
		std::map< std::string, utility::vector1< core::Size > > summed_norm;
		std::map< std::string, utility::vector1< core::Size > > summed_soft;

		for ( core::Size index = 1; index <= branches_.size(); ++index ) {
			std::string branch = branches_[ index ];

			summed_norm[ branch ].resize( pose.size(), 0);
			summed_soft[ branch ].resize( pose.size(), 0);
			std::cout << "Summing " << branch << std::endl;
			for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {

				summed_norm[branch][resnum] = glycan_glycan_res_atomic_clashes_[branch][resnum] +
					glycan_chain_res_atomic_clashes_[branch][resnum];

				//TR << "Norm: "<< resnum << " "<< glycan_glycan_res_atomic_clashes_[branch][resnum] ;
				//TR << "Norm: "<< resnum << " "<< glycan_chain_res_atomic_clashes_[branch][resnum];

				summed_soft[branch][resnum] =  glycan_glycan_res_soft_atomic_clashes_[branch][resnum] +
					glycan_chain_res_soft_atomic_clashes_[branch][resnum];


				//TR << "Bad: "<< resnum << " "<< summed_soft[branch][resnum] << std::endl;
			}
		}
		//TR << "Summed " << std::endl;
		tabulate_individual_data( "combined", "norm", pose, summed_norm);
		tabulate_individual_data( "combined", "soft", pose, summed_soft);

		report_data( pose );

	}



private:  // Private methods



	// Initialize data members from arguments.
	void
	init()
	{
		using namespace basic::options;

		type( "GlycanClashCheckMover" );

		branches_ =            option[ OptionKeys::carbohydrates::clash_check::glycan_branches ]();
		chains_ =              option[ OptionKeys::carbohydrates::clash_check::check_chains]();
		soft_clash_percent_ =   option[ OptionKeys::carbohydrates::clash_check::soft_clash ]();
		ignore_hydrogens_ =   option[ OptionKeys::carbohydrates::clash_check::ignore_hydrogens]();
		ignore_full_res_output_ = option[ OptionKeys::carbohydrates::clash_check::ignore_full_res_output]();
		output_per_glycan_data_ = option[ OptionKeys::carbohydrates::clash_check::output_per_glycan_data]();

		lj_n_ = pow( 2, 1.0/6 );

	}

	//This should really be stored in AtomType.
	core::Real lj_radius_to_zero_e_radius( core::Real num) const {
		//TR << num << std::endl;
		return num/lj_n_;
	}

	///Return a pair of booleans denoting a clash and a soft clash.  gx is the glycan atom num
	std::pair< bool, bool >
	is_clashing( core::pose::Pose const & pose, core::Size const glycan_resnum, core::Size const other_resnum, core::Size const gx){

		bool clash = false;
		bool soft_clash = false;

		for ( core::Size ox = 1; ox <= pose.residue( other_resnum ).natoms(); ++ ox ) {
			if ( pose.residue( other_resnum ).atom_type( ox ).is_virtual() ) continue;
			if ( ignore_hydrogens_ && pose.residue( other_resnum ).atom_type( ox ).element() == "H" ) continue;



			numeric::xyzVector<core::Real> gx_xyz= pose.residue( glycan_resnum ).xyz(gx);
			numeric::xyzVector<core::Real> ox_xyz= pose.residue( other_resnum ).xyz(ox);
			core::Real gx_ox_dis = gx_xyz.distance(ox_xyz);

			core::Real glycan_vdw_radii = lj_radius_to_zero_e_radius( pose.residue( glycan_resnum ).atom_type( gx ).lj_radius() );

			core::Real other_vdw_radii = lj_radius_to_zero_e_radius( pose.residue( other_resnum ).atom_type( ox ).lj_radius() );

			core::Real normal_cutoff = glycan_vdw_radii + other_vdw_radii;
			core::Real soft_cutoff = normal_cutoff*(1 - soft_clash_percent_ );
			//std::cout << gx_ox_dis << " " << normal_cutoff << " " << soft_cutoff << std::endl;
			if ( gx_ox_dis < soft_cutoff ) {
				soft_clash = true;
				clash = true; //Bad Clash implies clash!
				break;
			}
			if ( gx_ox_dis < normal_cutoff ) {
				clash = true;
				//TR << std::endl;
				//TR << pose.residue( glycan_resnum ).atom_name( gx )  <<"  -  "<< pose.residue( other_resnum ).atom_name( ox ) << std::endl;
				//TR << glycan_vdw_radii <<"  -  "<<other_vdw_radii <<std::endl;
				//TR << "Clash" << glycan_resnum << " to " << other_resnum<<" "<< gx_ox_dis << " " << normal_cutoff << " " << soft_cutoff << std::endl;

			}
		}
		return std::make_pair( clash, soft_clash );

	}
	/// @brief Calculate regular and soft clashes.  Again, this could be made better. Expecially in regard to pre-computing the soft cutoffs, lj_radius at zero, etc.
	///
	/// @details
	///     N Clashes: Number of atoms making at least one clash to other resnum
	///     N Bad Clashes: Number of atoms making at least one soft clash to other resnum.
	///
	void
	calculate_atom_clashes( std::string const & clash_type, const core::pose::Pose & pose, std::string branch, core::Size const & glycan_resnum, utility::vector1< core::Size > other_resnums, std::map< std::string, utility::vector1< core::Size > > & atom_clashes, std::map< std::string, utility::vector1< core::Size > > & soft_atom_clashes) {

		TR << "calculating atom clashes: " << clash_type << " for branch " << branch << std::endl;
		for ( core::Size gx = 1; gx <= pose.residue( glycan_resnum ).natoms(); ++gx ) {
			if ( pose.residue( glycan_resnum ).atom_type( gx ).is_virtual() ) continue;
			if ( ignore_hydrogens_ && pose.residue( glycan_resnum ).atom_type( gx ).element() == "H" ) continue;

			core::Size other_resnum = 0;
			std::pair < bool, bool > clashing;

			bool normal_clash = false;
			bool soft_clash = false;
			for ( core::Size index = 1; index <= other_resnums.size(); ++ index ) {

				other_resnum = other_resnums[ index ];

				if ( other_resnum == glycan_resnum || branch_point_resnums_.has_value(other_resnum) ) continue;

				//Don't calculate glycan - self clashing for now.
				if ( clash_type == "glycan" && branch_residues_[ branch ].has_value( other_resnum ) ) {
					continue;
				}
				//std::cout << "Glycan: "<< glycan_resnum << " Other: " << other_resnum << " " << gx << std::endl;
				clashing  = is_clashing(pose, glycan_resnum, other_resnum, gx );

				if ( clashing.first ) {
					normal_clash = true;
				}
				//If we have found a clash and a soft clash, we are done here.
				if ( clashing.second ) {
					soft_clash = true;
					break;
				}

			}

			//Now we report the data.
			if ( normal_clash ) {
				atom_clashes[ branch ][ glycan_resnum ]+=1;
			}
			if ( soft_clash ) {
				soft_atom_clashes[ branch ][ glycan_resnum ]+=1;
			}

		}

	} // End Function

	//Calculate and report all the data.
	void
	tabulate_individual_data(std::string const & clash_type, std::string const & data_type, core::pose::Pose const & pose,
		std::map< std::string, utility::vector1< core::Size > > const & data
	){

		std::string main_prefix = clash_type+"_"+data_type;
		std::string prefix;

		//Residue information
		core::Size overall_atom_clashes = 0;
		core::Size overall_res_clashes_1 = 0;
		core::Size overall_res_clashes_3 = 0;
		core::Size overall_res_clashes_5 = 0;

		for ( branch_it_type it = data.begin(); it != data.end(); ++it ) {
			prefix = main_prefix + "_Branch-" + it->first;

			//TR << "Reporting " << clash_type << " data_type: "<< data_type << " branch_name: " << it->first << std::endl;

			//Total - just imagine how easy this is in python.  Then be sad.
			core::Size atom_clashes = total_atom_clashes( it->second);
			overall_atom_clashes += atom_clashes;

			core::Size res_clashes_1 = total_residue_clashes( it->second, 1);
			overall_res_clashes_1 += res_clashes_1;

			core::Size res_clashes_3 = total_residue_clashes( it->second, 3);
			overall_res_clashes_3 += res_clashes_3;

			core::Size res_clashes_5 = total_residue_clashes( it->second, 5);
			overall_res_clashes_5 += res_clashes_5;

			if ( output_per_glycan_data_ ) {

				score_data_[prefix + "_" + "total_atom_clashes"] = atom_clashes;
				score_data_[prefix + "_" + "total_residue_clashes_1"] = res_clashes_1;
				score_data_[prefix + "_" + "total_residue_clashes_5"] = res_clashes_5;
			}

			if ( ! ignore_full_res_output_ ) {
				for ( core::Size resnum = 1; resnum <= it->second.size(); ++resnum ) {
					if ( ! branch_residues_[it->first].has_value( resnum ) ) continue;
					std::string chain = utility::to_string(pose.pdb_info()->chain( resnum ));
					std::string pdb_num = utility::to_string(pose.pdb_info()->number( resnum ));
					std::string pdb_name =  pdb_num + chain;
					std::string icode = utility::to_string( pose.pdb_info()->icode( resnum ) );
					if ( icode != " " ) {
						pdb_name = pdb_name+"-"+icode;
					}
					//TR << main_prefix+"_"+pdb_name+"_total_atom_clashes " << it->second[ resnum ] << std::endl;
					res_data_[main_prefix+"_"+pdb_name+"_total_atom_clashes"] = it->second[ resnum ];

					std::string pymol_entry = "( resid " +pdb_num+" and chain "+chain+")";

					if ( it->second[resnum] > 0 && ! pymol_print_[main_prefix].has_value( pymol_entry ) ) {
						pymol_print_[main_prefix].push_back( pymol_entry );
					}
				}
			}



		}
		score_data_[main_prefix+"_"+"total_atom_clashes"] = overall_atom_clashes;
		score_data_[main_prefix+"_"+"total_residue_clashes_1"] = overall_res_clashes_1;
		//score_data_[main_prefix+"_"+"total_residue_clashes_3"] = overall_res_clashes_3;
		score_data_[main_prefix+"_"+"total_residue_clashes_5"] = overall_res_clashes_5;
		//TR << "reported " << main_prefix << std::endl;

	}

	core::Size
	total_residue_clashes(utility::vector1< core::Size > const clashes, core::Size at_least_x){
		core::Size total_residue_clashes = 0;
		for ( core::Size resnum = 1; resnum <= clashes.size(); ++resnum ) {
			if ( clashes[resnum] >= at_least_x ) {
				total_residue_clashes+=1;
			}
		}
		return total_residue_clashes;
	}

	core::Size
	total_atom_clashes(utility::vector1< core::Size > const clashes ){
		core::Size total_atom_clashes = 0;
		for ( core::Size resnum = 1; resnum <= clashes.size(); ++resnum ) {
			total_atom_clashes += clashes[resnum];
		}
		return total_atom_clashes;
	}

	void
	report_data(core::pose::Pose & pose){
		typedef std::map< std::string , core::Real >::const_iterator it_type;
		typedef std::map< std::string, utility::vector1< std::string > >::const_iterator pymol_it_type;

		TR << "Printing only those scores greater than 1" << std::endl;

		TR << "Res Data" << std::endl;
		if ( res_data_.size() > 0 ) {
			for ( it_type it = res_data_.begin(); it != res_data_.end(); ++it ) {
				core::pose::setPoseExtraScore(pose, it->first, it->second);
				if ( it->second > 0 ) {
					TR << it->first << " " << it->second << std::endl;
				}
			}
		}

		TR << std::endl;
		TR << " Selections with more than 1 clash: " << std::endl;
		if ( pymol_print_.size() > 0 ) {
			for ( pymol_it_type pymol_it = pymol_print_.begin(); pymol_it != pymol_print_.end(); ++pymol_it ) {

				if ( pymol_it->second.size() == 0 ) continue;
				TR << std::endl;
				TR << "Selections "<< pymol_it->first << std::endl;
				if ( pymol_it->second.size() == 1 ) {
					TR << pymol_it->second[ 1 ] << std::endl;
				} else {
					std::string pymol_string = pymol_it->second[ 1 ];
					for ( core::Size i = 2; i <= pymol_it->second.size(); ++ i ) {
						pymol_string = pymol_string + " or " + pymol_it->second[ i ];
					}
					TR << pymol_string << std::endl;
				}

			}

		}

		TR << std::endl;
		for ( it_type it2 = score_data_.begin(); it2 != score_data_.end(); ++it2 ) {
			core::pose::setPoseExtraScore(pose, it2->first, it2->second);
			if ( it2->second > 0 ) {
				TR << it2->first << " " << it2->second << std::endl;
			}
		}
	}


private:  // Private data

	/// @brief Starting glycan branches.  Strings set from command line.
	utility::vector1< std::string > branches_;

	utility::vector1< core::Size > branch_point_resnums_;

	/// @brief Chains we are getting info on.
	utility::vector1< std::string > chains_;

	/// @brief Rosetta chain numberings.
	utility::vector1< core::Size > chain_nums_;

	/// @brief Branching residue numbers for each branch start.
	std::map< std::string, utility::vector1< core::Size > >
		branch_residues_;

	/// @brief Glycan to Glycan total number of clashes defined by atom-atom fa_rep radii
	std::map< std::string, utility::vector1< core::Size > >
		glycan_glycan_res_atomic_clashes_;

	/// @brief Glycan to chain atomic clashes.
	std::map< std::string, utility::vector1< core::Size > >
		glycan_chain_res_atomic_clashes_;


	/// @brief Glycan to Glycan total number of clashes defined by atom-atom fa_rep radii
	std::map< std::string, utility::vector1< core::Size > >
		glycan_glycan_res_soft_atomic_clashes_;

	/// @brief Glycan to chain atomic clashes.
	std::map< std::string, utility::vector1< core::Size > >
		glycan_chain_res_soft_atomic_clashes_;


	/// @brief When we calculate atom-atom distances using VDW, clash is softly if distance < (atomI_vdw + atomJ_vdw)*(1 - soft_clash)
	core::Real soft_clash_percent_;

	bool ignore_hydrogens_;
	bool ignore_full_res_output_;
	bool output_per_glycan_data_;

	utility::vector1< core::Size > glycan_resnums_; // All carbohydrate residues.
	utility::vector1< core::Size > chain_resnums_; // All resnums that match the chains we are searching again.

	std::map< std::string, core::Real > score_data_; //Final Score Data
	std::map< std::string, core::Real > res_data_; //Final Score Data for each residue.
	std::map< std::string, utility::vector1< std::string > > pymol_print_;

	core::Real lj_n_; //pow( 2, 1/6 );

};

typedef utility::pointer::shared_ptr< GlycanClashCheckMover > GlycanClashCheckMoverOP;




























///////////////////////////////////////////////////////// Main //////////////////////////////////////////////////////


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );
		register_options();


		if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}

		// Make sure the default JobOutputter is SilentJobOutputter to ensure that when this
		// is called with default arguments is prints a proper scorefile and not the hacky thing that
		// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		// Copied from score_jd2.
		protocols::jd2::SilentFileJobOutputterOP jobout( new protocols::jd2::SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);
		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( protocols::jd2::JobDistributorFactory::create_job_outputter( jobout ));

		GlycanClashCheckMoverOP mover_protocol( new GlycanClashCheckMover() );

		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
