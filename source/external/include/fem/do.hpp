#ifndef FEM_DO_HPP
#define FEM_DO_HPP

#define FEM_DO(i, f, l) for(i=f; i<=l; i++)

namespace fem {

  class do_safe
  {
    int* i_;
    int l_;

    public:

      do_safe(int* i, int f, int l) : i_(i), l_(l) { (*i) = f; }

      operator
      bool() { return (*i_) <= l_; }

      void
      incr() { (*i_)++; }
  };

} // namespace fem

#define FEM_DO_SAFE(i, f, l) \
  for(fem::do_safe fem_do(&(i), f, l);fem_do;fem_do.incr())

namespace fem {

  class dostep
  {
    int* i_;
    int l_;
    int s_;

    public:

      dostep(int* i, int f, int l, int s) : i_(i), l_(l), s_(s) { (*i) = f; }

      operator bool()
      {
        if (s_ < 0) return (*i_) >= l_;
        return (*i_) <= l_;
      }

      void
      incr() { (*i_) += s_; }
  };

} // namespace fem

#define FEM_DOSTEP(i, f, l, s) \
  for(fem::dostep fem_do(&(i), f, l, s);fem_do;fem_do.incr())

#endif // GUARD
