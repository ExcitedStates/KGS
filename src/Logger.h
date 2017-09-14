/*

Excited States software: KGS
Contributors: See CONTRIBUTORS.txt
Contact: kgs-contact@simtk.org

Copyright (C) 2009-2017 Stanford University

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

This entire text, including the above copyright notice and this permission notice
shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
IN THE SOFTWARE.

*/
#ifndef LOGGER_H
#define LOGGER_H

#include <map>
#include <string>
#include <iostream>


//From http://stackoverflow.com/questions/760301/implementing-a-no-op-stdostream/3261439#3261439
//Only used to define onullstream which is a drain
//template <class cT, class traits = std::char_traits<cT> >
//class basic_nullbuf: public std::basic_streambuf<cT, traits> {
//    typename traits::int_type overflow(typename traits::int_type c)
//    {
//        return traits::not_eof(c); // indicate success
//    }
//};
//
//template <class cT, class traits = std::char_traits<cT> >
//class basic_onullstream: public std::basic_ostream<cT, traits> {
//public:
//    basic_onullstream():
//        std::basic_ios<cT, traits>(&m_sbuf),
//        std::basic_ostream<cT, traits>(&m_sbuf)
//    {
//        init(&m_sbuf);
//    }
//
//private:
//    basic_nullbuf<cT, traits> m_sbuf;
//};
//
//typedef basic_onullstream<char> onullstream;

////From http://asmodehn.wordpress.com/2010/06/20/busy-c-coding-and-testing/
////Only used to define onullstream which is a drain
//template <class cT, class traits = std::char_traits<cT> >
//class basic_nullstreambuf: public std::basic_streambuf<cT, traits>
//{
//public:
//    basic_nullstreambuf();
//    ~basic_nullstreambuf();
//    typename traits::int_type overflow(typename traits::int_type c)
//    {
//        return traits::not_eof(c); // indicate success
//    }
//};
//
//template <class _CharT, class _Traits>
//basic_nullstreambuf<_CharT, _Traits>::basic_nullstreambuf() : std::basic_streambuf<_CharT, _Traits>(){}
//
//template <class _CharT, class _Traits>
//basic_nullstreambuf<_CharT, _Traits>::~basic_nullstreambuf(){}
//
//template <class cT, class traits = std::char_traits<cT> >
//class basic_onullstream: public std::basic_ostream<cT, traits> {
//public:
//    basic_onullstream():
//        std::basic_ios<cT, traits>(),
//        std::basic_ostream<cT, traits>(0),
//        m_sbuf()
//    {
//        this->init(&m_sbuf);
//    }
//private:
//    basic_nullstreambuf<cT, traits> m_sbuf;
//};
//
//typedef basic_onullstream<char> onullstream;

//From: http://stackoverflow.com/questions/11826554/standard-no-op-output-stream
class NullBuffer : public std::streambuf
{
 public:
  int overflow(int c) { return c; }
};

class onullstream : public std::ostream {
 public:
  onullstream() :
      m_sb(),
      std::ostream(&m_sb)
  {}
 private:
  NullBuffer m_sb;
};


/* A logger class for maintaining log output. Each call to log returns an ostream associated with the
   supplied logger-name. Depending on the logger-name the message can be suppressed or propagated.

   A set of global methods ease the use of this class. Example of use:

   enableLogger("test");
   enableLogger("debug", cerr);
   log("default")<<"Message 1"<<endl; 	//Will NOT be printed
   log("test")<<"Message 2"<<endl; 	//Will be printed to cout
   log("debug")<<"Message 3"<<endl; 	//Will be printed to cerr
   disableLogger("debug");
   log("debug")<<"Message 4"<<endl; 	//Will NOT be printed

   log()<<"Message 5"<<endl;			//Will NOT be printed
   enableLogger("default");
   log()<<"Message 6"<<endl;			//Will be printed to cout

   If no logger-name is supplied "default" is assumed. The "default" logger is DISABLED by default.
*/
class Logger{
public:

    /** Get the stream associated with the logger-name. If no stream is associated with the name a
      drain is returned. */
    std::ostream& log(const std::string& loggerName);

    /** Set the stream associated with the logger-name. */
    void enableLogger(const std::string& loggerName, std::ostream& stream);

    /** Disable a logger-name by disassociating it with any stream. */
    void disableLogger(const std::string& loggerName);

    /** Indicate if the loggerName is enabled */
    bool loggerEnabled(const std::string& loggerName);

    static Logger* getInstance();
private:
    Logger();
    std::map<std::string, std::ostream*> activeLoggers;
    onullstream drain;
    static Logger* instance;
};


/** Get the "default" logger stream. */
std::ostream& log();

/** Get the stream associated with the logger-name. If no stream is associated with the name a
      drain is returned. */
std::ostream& log(const std::string& loggerName);

/** Set the stream associated with the logger-name to cout. */
void enableLogger(const std::string& loggerName);

/** Set the stream associated with the logger-name. */
void enableLogger(const std::string& loggerName, std::ostream& stream);

/** Disable a logger-name by disassociating it with any stream. */
void disableLogger(const std::string& loggerName);

/** Indicate if the loggerName is enabled */
bool loggerEnabled(const std::string& loggerName);

/** Turn the logger-name restrictions on/off. */
void setAllLogsEnabled(bool enable);

#endif // LOGGER_H
