////////////////////////////////////////////////////////////////////////////////
//     pathways_planners.h
//     This class is for local planners of motion, i.e. methods to
//     densely sample a conformational range between two conformations
//     such as a linear planner that linearly iterates over the DOFs
//     values in between two conformations
//
//	   Author: Barak Raveh & Angela Enosh
//     Created: 17/12/2007
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_planners_hh
#define INCLUDED_devel_path_rover_pathrover_planners_hh


// ObjexxFCL Headers
#include <ObjexxFCL/ObjexxFCL.hh>
#include <stdio.h>
#include <algorithm>
#include <string>
#include <vector>
#include "pose.h"


namespace pathways {

  // ********** Functors for addition and subtraction: ********

  // binary functor
  class Periodic_plus
    {
    public:
      Periodic_plus(double range = 360.0){
    	  _range = range;
      }

      double operator()(double s1, double s2)
	{
	  return periodic_range(s1 + s2, _range);
	}

    private:
      double _range;
    };

  // binary functor
  class Periodic_minus
    {
    public:
      Periodic_minus(double range = 360.0){
	_range = range;
      }

      double operator()(double s1, double s2)
	{
	  return periodic_range(s1 - s2, _range);
	}

    private:
      double _range;
    };


  // *********** Linear_planner_iterator class: **********

  // an iterator between "from" and "to" vectors of DOFs, according to specified step parameters
  // for use by Linear planners of the RRT algorithm
  class Linear_planner_iterator
  {
  public:
    Linear_planner_iterator(std::vector<double> v_from, std::vector<double> v_to, double max_step_size = 1)
    {
      _v_from = v_from;
      _v_to = v_to;
      _step_num = 0;
      _v_cur = v_from; // step #0 = from
      calc_step_vector(max_step_size);
    }


    // go to the first element of the iteration
    void go_begin()
    {  _step_num = 0; _v_cur = _v_from;   }

    // go beyond the last element of the iteration
    // i.e. invoking operator-- right after go_end() will
    //       move the iterator to the last element exactly
    void go_end()
    {
      _step_num = _total_steps_num;
      _v_cur = _v_to;
      (*this)++; // go beyond last valid place
    }

    bool is_begin()
    { return _step_num == 0; }

    bool is_end()
    { return _step_num > _total_steps_num; }

    int step_num()
      {	return _step_num; }

    void operator++() // prefix i.e. ++iterator
    {
      // std::cout << "Linear_planner_iterator::prefix++" << std::endl;
      (*this)++; // prefix and postfix are same...
    }

    void operator++(int)  // postfix i.e. iterator++
    {
      //      std::cout << "Linear_planner_iterator::postfix++" << std::endl;
      assert(!is_end());
      _step_num++;
      // add up step_delta to _v_cur
      std::transform(
    		  _v_cur.begin(), _v_cur.end(),
    		  _v_step_delta.begin(),
    		  _v_cur.begin(),
    		  Periodic_plus(360.0));
    }

    void operator--() // prefix i.e. ++iterator
    {
      (*this)--; // prefix and postfix are same...
    }

    void operator--(int)  // postfix i.e. iterator++
    {
      assert(!is_begin());
      _step_num--;
      // add up step_delta to _v_cur
      std::transform(
    		  _v_cur.begin(), _v_cur.end(),
    		  _v_step_delta.begin(),
    		  _v_cur.begin(),
    		  Periodic_minus(360.0));
    }

    std::vector<double>& operator*() // TODO: const iterator?
      {
    	return _v_cur;
      }

  private:

    // returns the angle differnce of to-from in the range [-180,180]
    double diff_angle(double from, double to) const;

    // calculates a vector of the delta-values for each step in the iteration,
    // s.t. the maximal step size of any given DOF is max_step_size
    void calc_step_vector(double max_step_size);

  private:
    int _step_num, _total_steps_num;
    std::vector<double> _v_from, _v_to; // vectores of "from" and "to" DOF assignments
    std::vector<double> _v_step_delta; // a vector for the differences in each step
    std::vector<double> _v_cur;

  };


} // namespace pathways

#endif

