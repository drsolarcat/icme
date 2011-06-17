
// project headers
#include "gnuplot.h"

#ifdef GNUPLOT_ENABLE_PTY

GnuplotPty::GnuplotPty(bool debug_messages) :
	pty_fh(NULL),
	master_fd(-1),
	slave_fd(-1)
{
// adapted from http://www.gnuplot.info/files/gpReadMouseTest.c
	if(0 > openpty(&master_fd, &slave_fd, NULL, NULL, NULL)) {
		perror("openpty");
		throw std::runtime_error("openpty failed");
	}
	char pty_fn_buf[1024];
	if(ttyname_r(slave_fd, pty_fn_buf, 1024)) {
		perror("ttyname_r");
		throw std::runtime_error("ttyname failed");
	}
	pty_fn = std::string(pty_fn_buf);
	if(debug_messages) {
		std::cerr << "fn=" << pty_fn << std::endl;
	}

	// disable echo
	struct termios tios;
	if(tcgetattr(slave_fd, &tios) < 0) {
		perror("tcgetattr");
		throw std::runtime_error("tcgetattr failed");
	}
	tios.c_lflag &= ~(ECHO | ECHONL);
	if(tcsetattr(slave_fd, TCSAFLUSH, &tios) < 0) {
		perror("tcsetattr");
		throw std::runtime_error("tcsetattr failed");
	}

	pty_fh = fdopen(master_fd, "r");
	if(!pty_fh) {
		throw std::runtime_error("fdopen failed");
	}
}

GnuplotPty::~GnuplotPty() {
	if(pty_fh) fclose(pty_fh);
	if(master_fd > 0) ::close(master_fd);
	if(slave_fd  > 0) ::close(slave_fd);
}
#endif // GNUPLOT_ENABLE_PTY

///////////////////////////////////////////////////////////

Gnuplot::Gnuplot(std::string cmd) :
	boost::iostreams::stream<boost::iostreams::file_descriptor_sink>(
		FILENO(pout = POPEN(cmd.c_str(), "w")),
		boost::iostreams::never_close_handle),
	gp_pty(NULL),
	debug_messages(false)
{
	setf(std::ios::scientific, std::ios::floatfield);
	precision(18);
}

Gnuplot::~Gnuplot() {
	if(debug_messages) {
		std::cerr << "closing gnuplot" << std::endl;
	}

	// FIXME - boost's close method calls close() on the file descriptor, but
	// we need to use pclose instead.  For now, just skip calling boost's close
	// and use flush just in case.
	flush();
	//close();

	if(PCLOSE(pout)) {
		std::cerr << "pclose returned error" << std::endl;
	}

	if(gp_pty) delete(gp_pty);
}

#ifdef GNUPLOT_ENABLE_PTY
void Gnuplot::getMouse(
	double &mx, double &my, int &mb,
	std::string msg
) {
	allocPty();

	*this << "set mouse; set print \"" << gp_pty->pty_fn << "\"" << std::endl;
	*this << "pause mouse \"" << msg << "\\n\"" << std::endl;
	*this << "print MOUSE_X, MOUSE_Y, MOUSE_BUTTON" << std::endl;
	if(debug_messages) {
		std::cerr << "begin scanf" << std::endl;
	}
	if(3 != fscanf(gp_pty->pty_fh, "%lf %lf %d", &mx, &my, &mb)) {
		throw std::runtime_error("could not parse reply");
	}
	if(debug_messages) {
		std::cerr << "end scanf" << std::endl;
	}
}
#endif // GNUPLOT_ENABLE_PTY

