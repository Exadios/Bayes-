#ifndef _BAYES_FILTER_EXCEPTION
#define _BAYES_FILTER_EXCEPTION

/*
 * Bayes++ the Bayesian Filtering Library
 * Copyright (c) 2002 Michael Stevens
 * Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * See Bayes++.htm for copyright license details
 *
 * $Header$
 * $NoKeywords: $
 */

/*
 * Exception types: Exception heirarchy for Bayesian filtering 
 */
 
// Common headers required for declerations
#include <exception>

/* Filter namespace */
namespace Bayesian_filter
{


class Filter_exception : virtual public std::exception
/*
 *	Base class for all exception produced by filter heirachy
 */
{
public:
	const char *what() const throw()
	{	return error_description;
	}
protected:
	Filter_exception (const char* description)
	{	error_description = description;
	};
private:
	const char* error_description;
};

class Logic_exception : virtual public Filter_exception
/*
 * Logic Exception
 */
{
public:
	Logic_exception (const char* description) :
		Filter_exception (description)
	{};
};

class Numeric_exception : virtual public Filter_exception
/*
 * Numeric Exception
 */
{
public:
	Numeric_exception (const char* description) :
		Filter_exception (description)
	{};
};


}//namespace
#endif
