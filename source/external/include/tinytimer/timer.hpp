// -*- C++ -*-
// timer.hpp
// A library to calculate elapsed times and easily write chrono::duration
// into an iostream
//
//
// MIT License
//
// Copyright (c) 2017 Roman Saldygashev
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#pragma once

#include <chrono>
#include <iostream>
#include <iomanip>
#include <ios>			// std::streamsize
#include <ratio>
#include <string>
#include <sstream>
#include <utility>


namespace tinytimer
{
    enum FF		// Format Flag
    { DAY, HOUR, MINUTE, SECOND };

    class Manip
    {
	FF _flag = DAY;
	std::streamsize _precision = 3;	// digits after point

    public:
	explicit
	Manip(const FF & flag,
	      const std::streamsize & precision) : _flag(flag), _precision(precision) {};

	static int
	getIndex()
	{
	    static int index = std::ios_base::xalloc();
	    return index;
	}

	FF
	getFF() const noexcept
	{ return _flag; }

	std::streamsize
	getPrecision() const noexcept
	{ return _precision; }

	friend
	std::ostream &
	operator<<(std::ostream & os, Manip && m)
	{
	    os.pword(getIndex()) = static_cast<void*>(&m);
	    return os;
	}
    };

    Manip
    set_duration_format(const FF flag = DAY,
			const std::streamsize precision = 3)
    { return Manip(flag, precision); }


    template <typename Clock = std::chrono::steady_clock>
    class ElapsedTimer
    {
	using ClockDuration = std::chrono::duration<typename Clock::rep, typename Clock::period>;
	using TimePoint = typename Clock::time_point;

	TimePoint startpoint = Clock::now();
    public:
	ElapsedTimer() {}

	ElapsedTimer(const ElapsedTimer & rhs) : startpoint(rhs.startpoint) {}

	explicit
	ElapsedTimer(const TimePoint & rhs) : startpoint(rhs) {}

	bool
	operator==(const ElapsedTimer & rhs) const noexcept
	{ return startpoint == rhs.startpoint; }
	bool
	operator!=(const ElapsedTimer & rhs) const noexcept
	{ return startpoint != rhs.startpoint; }

	// is younger
	bool
	operator<(const ElapsedTimer & rhs) const noexcept
	{ return startpoint > rhs.startpoint; }
	bool
	operator<=(const ElapsedTimer & rhs) const noexcept
	{ return startpoint >= rhs.startpoint; }

	// is older
	bool
	operator>(const ElapsedTimer & rhs) const noexcept
	{ return startpoint < rhs.startpoint; }
	bool
	operator>=(const ElapsedTimer & rhs) const noexcept
	{ return startpoint <= rhs.startpoint; }

	ClockDuration
	operator-(const ElapsedTimer & rhs) const noexcept
	{ return startpoint - rhs.startpoint; }

	ClockDuration
	elapsed() const noexcept
	{ return Clock::now() - startpoint; }

	ClockDuration
	restart() noexcept
	{
	    const ClockDuration ret = elapsed();
	    startpoint = Clock::now();
	    return ret;
	}

	template <typename Rep, typename Period>
	bool
	has_expired(const std::chrono::duration<Rep, Period> & timeout) const noexcept
	{ return elapsed() > timeout; }
    };

    template <typename F, typename Clock = std::chrono::steady_clock, typename ... Args>
    std::chrono::duration<typename Clock::rep, typename Clock::period>
    measure(F && callable, Args && ...args)
    {
	ElapsedTimer<Clock> timer;

	callable(std::forward<Args>(args)...);

	return timer.elapsed();
    }
}


namespace std
{
    namespace chrono
    {
	typedef duration<int64_t, ratio<86400>> days;
    }


    template <typename Rep, intmax_t Num, intmax_t Denom>
    ostream &
    operator<<(ostream & os,
	       chrono::duration<Rep, ratio<Num, Denom>> duration)
    {
	tinytimer::FF flag = tinytimer::SECOND;
	streamsize newprecision = 3;

	auto word = os.pword(tinytimer::Manip::getIndex());

	if (word) {
	    tinytimer::Manip *info = static_cast<tinytimer::Manip*>(word);
	    os.iword(tinytimer::Manip::getIndex()) = 0;
	    flag = info->getFF();
	    newprecision = info->getPrecision();
	}

	ios state(nullptr);
	state.copyfmt(os);

	os << setfill('0');

	switch (flag) {
	case tinytimer::DAY:
	    {
		auto days = chrono::duration_cast<chrono::days>(duration);
		os << days.count() << ':';
		duration -= days;
	    }
	case tinytimer::HOUR:
	    {
		auto hours = chrono::duration_cast<chrono::hours>(duration);
		os << setw(2) << hours.count() << ':';
		duration -= hours;
	    }
	case tinytimer::MINUTE:
	    {
		auto minutes = chrono::duration_cast<chrono::minutes>(duration);
		os << setw(2) << minutes.count() << ':';
		duration -= minutes;
	    }
	case tinytimer::SECOND:
	    break;
	}

	auto dseconds = chrono::duration_cast<chrono::seconds>(duration);
	duration -= dseconds;
	double fseconds = static_cast<double>(dseconds.count()) + static_cast<double>(duration.count()) / Denom;
	os << setprecision(newprecision) << setw(3+newprecision) << fixed << fseconds;
	os.copyfmt(state);

	return os;
    }


    template <typename Rep, typename Period>
    string
    to_string(const chrono::duration<Rep, Period> & duration,
	      const tinytimer::FF flag = tinytimer::DAY,
	      const streamsize precision = 3)
    {
	stringstream ss;
	ss << tinytimer::set_duration_format(flag, precision) << duration;
	return ss.str();
    }
}
