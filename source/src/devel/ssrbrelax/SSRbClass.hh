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
/// @author Vatsan Raman

#ifndef INCLUDED_devel_ssrbrelax_SSRbClass_hh
#define INCLUDED_devel_ssrbrelax_SSRbClass_hh

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.fwd.hh>

// C++ Headers
#include <iosfwd>
#include <map>


// Utility Headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

///////////////////////////////////////////////////////////////////////////////
namespace devel {
namespace ssrbrelax {
	class Rb {

		int index_;
		int seg_begin_;
		int seg_end_;
		int anchor_pos_;
		int loop_stub_n_;
		int loop_stub_c_;

	public:
	Rb():
		index_(0),
		seg_begin_(0),
			seg_end_(0),
		anchor_pos_(0),
		loop_stub_n_(0),
		loop_stub_c_(0)
	{}

		Rb(
			 int const index, int const seg_begin, int const seg_end,
			 int const anchor_pos, int loop_stub_n, int loop_stub_c
			 ):
			index_( index ),
		seg_begin_( seg_begin ),
			seg_end_( seg_end ),
		anchor_pos_( anchor_pos ),
		loop_stub_n_( loop_stub_n ),
		loop_stub_c_( loop_stub_c )
			{}

		inline int index() const { return index_; }
		inline int seg_begin() const { return seg_begin_; }
		inline int seg_end() const { return seg_end_; }
		inline int anchor_pos() const { return anchor_pos_; }
		inline int loop_stub_n() const { return loop_stub_n_; }
		inline int loop_stub_c() const { return loop_stub_c_; }

		friend std::ostream & operator<<(std::ostream & is, const Rb & rb);

	};

////////////////////////////////////////////////////////////////
//sraman
//A list of rigid body segments.
//Methods to add, delete, pick segments
//set up one segment fold tree
/////////////////////////////////////////////////////////////////

	class RbSegments {
	public:
		typedef std::vector< Rb > Rb_list;
		typedef Rb_list::iterator iterator;
		typedef Rb_list::const_iterator const_iterator;

	private:
		Rb_list rb_list;

	public:
		//constructor
		RbSegments(){};

		//copy constructor
		RbSegments( const RbSegments & src ):rb_list( src.rb_list ){}
		//operator
		RbSegments & operator =( RbSegments const & src ) {
			rb_list = src.rb_list;
			return *this;
		}
		//		std::map < std::pair < int, int >, int > segment_anchor_map; //stores segments and their corresponding anchor positions

		typedef std::map< std::pair< int, int >, float > RbFunc;
		typedef RbFunc::iterator RbFunc_it;
		typedef RbFunc::const_iterator RbFunc_const_it;
		RbFunc rbfunc;

		std::map< std::pair< int, int >, int > FlexJumpRes;
	 	typedef numeric::xyzVector< Real > Vector;
		///		typedef std::vector< numeric::xyzVector < float > > Perturbation_Axes;

		friend std::ostream & operator<<( std::ostream & os, const RbSegments & rbsegments );

		inline int num_rb() const { return rb_list.size(); }

		//Following functions assume that there is only one segment in the rb_list. This is needed because rb_list is declared as private

		inline int index() const;
		inline int seg_start() const;
		inline int seg_stop() const;
		inline int anchor_res() const;
		inline int n_term_loop() const;
		inline int c_term_loop() const;

		//		inline int get_anchor_position( int seg_begin, int seg_end ) { return segment_anchor_map[ std::pair < int, int >( seg_begin, seg_end ) ]; }

	void
	one_segment_fold_tree(
		core::kinematics::FoldTree & f,
		int const total_residue
	);

	void
	read_segments_from_file();

	void
	add_segment(
	int const index,
	int const seg_begin,
	int const seg_end,
	int const anchor_pos,
	int const loop_stub_n,
	int const loop_stub_c
	);

	void
	add_segment(
	const RbSegments::iterator & it
	);

	bool
	verify_segment(
	int const seg_begin,
	int const seg_end,
	int const anchor_pos,
	int const loop_stub_n,
	int const loop_stub_c
	);

	void
	delete_segment(
	 int const seg_begin,
	 int const seg_end
	 );

	void
	delete_segment(
	 const RbSegments::iterator & it
	 );

	iterator one_random_segment();

	void
	read_param_file();

		float
		get_gaussian_parameters(
			int const & index,
			int const & dof
		);

		bool
		distribution_exists(
												int const & index,
												int const & dof
												);

		int
		get_flex_jump_pos(
											int const & seg_begin,
											int const & seg_end
											);

		Vector const
		X_axis(
					 core::pose::Pose & pose
					 );

		Vector const
		Z_axis(
					 core::pose::Pose & pose
					 );

		Vector const
		Y_axis(
					 core::pose::Pose & pose
					 );

		Vector const
		alt_X_axis(
			core::pose::Pose & pose
		);

		Vector const
		alt_Y_axis(
			core::pose::Pose & pose
		);

		Vector const
		alt_Z_axis(
			core::pose::Pose & pose
		);

	};

	///////////////////////////////////////////////////////////////////////////
	//inline function definitions

	inline int RbSegments::index() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.index();
	}


	inline int RbSegments::seg_start() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.seg_begin();
	}


	inline int RbSegments::seg_stop() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.seg_end();
	}

	inline int RbSegments::anchor_res() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.anchor_pos();
	}

	inline int RbSegments::n_term_loop() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.loop_stub_n();
	}

	inline int RbSegments::c_term_loop() const
	{
		Rb rb( rb_list.at( 0 ) );
		return rb.loop_stub_c();
	}

} //namespace ssrbrelax
} //namespace devel


#endif


