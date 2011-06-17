/*

Eigen3 version by Alexey Isavnin (2011)
Original version by Daniel Stahlke (2009)

*/

#ifndef GNUPLOT_IOSTREAM_H
#define GNUPLOT_IOSTREAM_H

#define BOOST_THREAD_USE_LIB

#include "curve.h"

// C system includes
#include <stdio.h>
#ifdef GNUPLOT_ENABLE_PTY
#include <termios.h>
#include <unistd.h>
#include <pty.h>
#endif // GNUPLOT_ENABLE_PTY

// C++ system includes
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

// library includes: double quotes make cpplint not complain
#include "boost/iostreams/device/file_descriptor.hpp"
#include "boost/iostreams/stream.hpp"
#include "boost/noncopyable.hpp"
#ifdef GNUPLOT_ENABLE_BLITZ
#include "blitz/array.h"
#endif

// Patch for Windows by Damien Loison
#ifdef WIN32
#define PCLOSE _pclose
#define POPEN  _popen
#define FILENO _fileno
#else
#define PCLOSE pclose
#define POPEN  popen
#define FILENO fileno
#endif

#ifdef GNUPLOT_ENABLE_PTY
// this is a private class
class GnuplotPty {
public:
	explicit GnuplotPty(bool debug_messages);
	~GnuplotPty();

	std::string pty_fn;
	FILE *pty_fh;
	int master_fd, slave_fd;
};
#else // GNUPLOT_ENABLE_PTY
class GnuplotPty { };
#endif // GNUPLOT_ENABLE_PTY

///////////////////////////////////////////////////////////

class Gnuplot : public boost::iostreams::stream<
	boost::iostreams::file_descriptor_sink>, private boost::noncopyable
{
public:
	explicit Gnuplot(std::string cmd = "gnuplot");
	~Gnuplot();

#ifdef GNUPLOT_ENABLE_PTY
	void getMouse(
		double &mx, double &my, int &mb,
		std::string msg="Click Mouse!"
	);
#endif // GNUPLOT_ENABLE_PTY

	// used for one STL container
	template <class T>
	Gnuplot &send(T p, T last) {
		while(p != last) {
			sendEntry(*p);
			*this << "\n";
			++p;
		}
		*this << "e" << std::endl;
		return *this;
	}

	// used for two STL containers
	template <class T, class U>
	Gnuplot &send(T x, T x_last, U y, U y_last) {
		while(x != x_last && y != y_last) {
			sendEntry(*x, *y);
			*this << "\n";
			++x;
			++y;
		}
		// assert inputs same size
		assert(x==x_last && y==y_last);
		*this << "e" << std::endl;
		return *this;
	}

	// this handles STL containers as well as blitz::Array<T, 1> and
	// blitz::Array<blitz::TinyVector<T, N>, 1>
	template <class Iter>
	Gnuplot &send(Iter arr) {
		send(arr.begin(), arr.end());
		return *this;
	}

  // used to plot the Curve
	Gnuplot &send(const Curve& curve) {
	  for (int i = 0; i < curve.size(); i++) {
	    *this << curve.cols().x(i) << " " << curve.cols().y(i) << "\n";
	  }
	  *this << "e" << std::endl;
	  return *this;
	}

#ifdef GNUPLOT_ENABLE_BLITZ
	// Note: T could be either a scalar or a blitz::TinyVector.
	template <class T>
	Gnuplot &send(const blitz::Array<T, 2> &a) {
		for(int i=a.lbound(0); i<=a.ubound(0); i++) {
			for(int j=a.lbound(1); j<=a.ubound(1); j++) {
				sendEntry(a(i, j));
				*this << "\n";
			}
			*this << "\n";
		}
		*this << "e" << std::endl;
		return *this;
	}

private:
	template <class T, int N>
	void sendEntry(blitz::TinyVector<T, N> v) {
		for(int i=0; i<N; i++) {
			sendEntry(v[i]);
		}
	}
#endif // GNUPLOT_ENABLE_BLITZ

private:
	template <class T>
	void sendEntry(T v) {
		*this << v << " ";
	}

	template <class T, class U>
	void sendEntry(std::pair<T, U> v) {
		sendEntry(v.first, v.second);
	}

	template <class T, class U>
	void sendEntry(T t, U u) {
		sendEntry(t);
		sendEntry(u);
	}

#ifdef GNUPLOT_ENABLE_PTY
	void allocPty() {
		if(!gp_pty) {
			gp_pty = new GnuplotPty(debug_messages);
		}
	}
#endif // GNUPLOT_ENABLE_PTY

private:
	FILE *pout;
	// this is included even in the absense of GNUPLOT_ENABLE_PTY, to
	// protect binary compatibility
	GnuplotPty *gp_pty;

public:
	bool debug_messages;
};

#endif // GNUPLOT_IOSTREAM_H

