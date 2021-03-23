/*
    Simple Sky-Dome Projector Library
    Copyright (C) 2021  B. E. Pieters, 
    IEK-5 Photovoltaik, Forschunszentrum Juelich

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef _UTIL_H
#define _UTIL_H
typedef enum {QUIET, VERBOSE, VVERBOSE} VERB;
extern VERB ssdp_verbosity;

void Print(VERB v, const char *format_str, ...);
void Fatal( const char *format_str, ...);
void Warning( const char *format_str, ...);
int ProgressBar(int pcn, int pco, int len, int tics);
#endif
