#include "stopwatch.h"

#include <iostream>
#include <sstream>

Stopwatch::Stopwatch( const std::experimental::filesystem::path& protocolFile ) : protocol_( protocolFile ) {
	if ( !protocol_.is_open() ) {
		std::cerr << "Stopwatch protocol cannot be opened. Continue without procolling" << std::endl;
	}
};

Stopwatch::~Stopwatch() {
	if ( protocol_.is_open() ) {
		protocol_.close();
	}
}

Stopwatch::TimeStamp Stopwatch::start() {
	start_ = std::chrono::high_resolution_clock::now();
	return start_;
}

Stopwatch::TimeStamp Stopwatch::stop() {
	stop_ = std::chrono::high_resolution_clock::now();
	return stop_;
}