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


#include "CTKTimer.h"
#include <stdlib.h>

#ifdef WIN32
#include <windows.h>
// #pragma comment(lib,"")
#define GETCURRENTTIME(x) x=timeGetTime()
#else
#define GETCURRENTTIME(x) gettimeofday(&x,nullptr)
#endif //WIN32

// Sadly, timersub isn't defined in Solaris. :(
// So we use this instead. (added by Ryan)

#if defined (__SVR4) && defined (__sun)
#include "timersub.h"
#endif

CTKTimer::CTKTimer()
{
  Reset();
}

void CTKTimer::Reset()
{
  GETCURRENTTIME(start);
  current=start;
}

long long CTKTimer::ElapsedTicks()
{
  GETCURRENTTIME(current);
  return LastElapsedTicks();
}

long long CTKTimer::LastElapsedTicks() const
{
#ifdef WIN32
  return current-start;
#else
  timeval delta;
  timersub(&current,&start,&delta);
  long long ticks = delta.tv_sec*1000 + delta.tv_usec/1000;
  return ticks;
#endif //WIN32
}
    
double CTKTimer::ElapsedTime()
{
  GETCURRENTTIME(current);
  return LastElapsedTime();
}

double CTKTimer::LastElapsedTime() const
{
#ifdef WIN32
  return double(current-start)/1000.0;
#else
  timeval delta;
  timersub(&current,&start,&delta);
  double secs=double(delta.tv_sec);
  secs += double(delta.tv_usec)/1000000.0;
  return secs;
#endif //WIN32
}


double CTKTimer::getTimeNow()
{
  timeval tod;

  gettimeofday(&tod, nullptr);
  double time_seconds = (double) tod.tv_sec + ((double) tod.tv_usec / 1000000.0);
  return time_seconds;
}


/*
clock_t CTKTimer::ElapsedTicks()
{
  current = clock();
  return (current-start);
}

double CTKTimer::ElapsedTime()
{
  current = clock();
  return double(current-start)/CLOCKS_PER_SEC;
}

clock_t CTKTimer::LastElapsedTicks() const
{
  return current-start;
}

double CTKTimer::LastElapsedTime() const
{
  return double(current-start)/CLOCKS_PER_SEC;
}
*/
