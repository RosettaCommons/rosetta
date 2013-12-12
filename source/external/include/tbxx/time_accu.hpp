#ifndef TBXX_TIME_ACCU_HPP
#define TBXX_TIME_ACCU_HPP

#include <sys/time.h>

namespace tbxx {

  struct time_accu
  {
    timeval accu;
    timeval mark;

    time_accu()
    {
      timerclear(&accu);
      timerclear(&mark);
    }

    void
    set_mark()
    {
      gettimeofday(&mark, 0);
    }

    time_accu&
    accumulate()
    {
      timeval curr;
      gettimeofday(&curr, 0);
      timeval diff;
      timersub(&curr, &mark, &diff);
      timerclear(&mark);
      timeval sum;
      timeradd(&accu, &diff, &sum);
      accu.tv_sec = sum.tv_sec;
      accu.tv_usec = sum.tv_usec;
      return *this;
    }

    double
    as_double() const
    {
      return static_cast<double>(accu.tv_sec)
           + static_cast<double>(accu.tv_usec) * 1.e-6;
    }
  };

} // namespace tbxx

#endif // GUARD
