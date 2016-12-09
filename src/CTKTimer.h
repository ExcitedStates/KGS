/*
    LoopTK: Protein Loop Kinematic Toolkit
    Copyright (C) 2007 Stanford University

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#ifndef MY_TIMER_H
#define MY_TIMER_H

#ifdef WIN32
#include <windows.h>
typedef DWORD TimerCounterType;
#else
#include <sys/time.h>
typedef timeval TimerCounterType;
#endif //WIN32

class CTKTimer
{
 public:
  CTKTimer();
  void Reset();
  double getTimeNow(); 

  // Returns elapsed time in milliseconds,seconds respectively
  long long ElapsedTicks();
  double ElapsedTime();

  // Doesn't refresh the current time
  long long LastElapsedTicks() const;
  double LastElapsedTime() const;

 private:
  TimerCounterType start;
  TimerCounterType current;
};

#endif
