#pragma once
#include <string>
#include <fstream>
#include <experimental\filesystem>
#include <chrono>
#include <type_traits>

class Stopwatch
{
	using TimeStamp = decltype( std::chrono::high_resolution_clock::now() );
	TimeStamp start_;
	TimeStamp stop_;
	std::ofstream protocol_;

public:

	using Milliseconds = std::chrono::milliseconds;
	using Seconds = std::chrono::seconds;
	using Nanoseconds = std::chrono::nanoseconds;
	using Microseconds = std::chrono::microseconds;
	using Minutes = std::chrono::minutes;

	Stopwatch() {};
	Stopwatch( const std::experimental::filesystem::path& protocolFile );

	~Stopwatch();
	
	TimeStamp start();
	TimeStamp stop();
	template<class T> TimeStamp stop(const std::string& tag);
	template<class T> auto runtime();
	template<class T> std::string runtimeStr();
};

template<class T>
Stopwatch::TimeStamp Stopwatch::stop( const std::string& tag ) {
	stop_ = std::chrono::high_resolution_clock::now();
	if(protocol_.is_open() )
		protocol_ << tag << ": " << runtimeStr<T>() << std::endl;
	return stop_;
}

template<class T>
auto Stopwatch::runtime() {
	return std::chrono::duration_cast<T>( stop_ - start_ ).count();
}

template<class T>
std::string Stopwatch::runtimeStr() {
	
	std::string s = std::to_string( runtime<T>() );

	if ( std::is_same<T, Milliseconds>::value ) {
		s += " ms";
	}
	else if ( std::is_same<T, Seconds>::value ) {
		s += " s";
	}
	else if ( std::is_same<T, Nanoseconds>::value ) {
		s += " ns";
	}
	else if ( std::is_same<T, Microseconds>::value ) {
		s += " mus";
	}
	else if ( std::is_same<T, Minutes>::value ) {
		s += " min";
	}

	return s;
}