/*
    KGSX: Biomolecular Kino-geometric Sampling and Fitting of Experimental Data
    Yao et al, Proteins. 2012 Jan;80(1):25-43
    e-mail: latombe@cs.stanford.edu, vdbedem@slac.stanford.edu, julie.bernauer@inria.fr

        Copyright (C) 2011-2013 Stanford University

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
#include <iostream>
#include <sstream>
#include <cmath>

#include "Util.h"

using namespace std;
 
string Util::i2s(int x, int length) {
	string s = i2s(x);
	while (s.length() < length) {
		s = "0"+s;
	}
	return s;
}

string Util::d2s (double x) {
	std::ostringstream o;
	o << x;
	return o.str();
} 

string Util::f2s (float x) {
	std::ostringstream o;
	o << x;
	return o.str();
} 

string Util::i2s (int x) {
        std::ostringstream o;
        o << x;
        return o.str();
}

bool Util::stob (string s) {
        return ( s == "true" || s == "1" );
}

string Util::trim(string s,char c) { 
	int start= 0, end = s.size() -1; 
	int i = 0; 
	while (s[i++]==c) ++start; 
	i = s.size() -1; 
	while (s[i--]==c) --end; 
	return s.substr(start, end-start+1);
}

string Util::getPath (string s) {
	string path = "./";
	int index = s.find_last_of("/");
	if ( index!=string::npos )
		path = s.substr(0,index+1);
	return path;
}

string Util::getBaseName (string s) {
	string base_name = "";
	int name_index = s.find_first_of(".");
	int path_index = s.find_last_of("/");
	if ( name_index!=string::npos && path_index!=string::npos && name_index>path_index )
		base_name = s.substr(path_index+1,name_index-path_index-1);
	return base_name;
}

string Util::formatNumber(string number, int digits_num) {
        while ( number.length()<digits_num ) {
                number = "0" + number;
        }
        return number;
}

float Util::round(float number, int precision) {
	const int prec = std::pow(10.,precision);
	return std::floor( number*prec + 0.5 )/prec;
}

double Util::round(double number, int precision) {
	const int prec = std::pow(10.,precision);
	return std::floor( number*prec + 0.5 )/prec;
}

string Util::cutDecimal(float number, int decimals) {
	return f2s(round(number,decimals));
}

vector<string>& Util::split( const string& s, char delim, vector<string>& words ) {
  words.clear();
  stringstream ss(s);
  string item;
  while( getline(ss, item, delim) )
    words.push_back(item);
  return words;
}

vector<int>& Util::split( const string& s, char delim, vector<int>& numbers ) {
//See outcommented section of Selection
  return numbers;
}

vector<string> Util::split( const string& s, char delim ) {
  vector<string> words;
  split(s, delim, words);
  return words;
}

vector<string> Util::split( const string& s, const string& delim ) {
  vector<string> words;
  auto i = 0, pos = 0;
  do{
    pos = s.find(delim, i);
    words.push_back(s.substr(i, pos));
    i = pos+delim.size();
  }while(pos!=string::npos);

  return words;
}

bool Util::contains( const string& s, const string& substring ) {
  return s.find(substring)!=string::npos;
}

bool Util::startsWith( const string& s, const string& substring ) {
  return s.compare(0, substring.length(), substring)==0;
}

