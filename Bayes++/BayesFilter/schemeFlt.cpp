/*
 * Bayesian Filtering Library
 * (c) Michael Stevens, Australian Centre for Field Robotics 2000
 * $Header$
 * $NoKeywords: $
 */

/*
 * Filter Scheme implementation
 *  Sadly the template specialisation definitions and where it resides varies greatly with compilers
 */
#include "matSup.h"
#include "UDFlt.h"
#include "unsFlt.h"

#include "schemeFlt.h"

/* Filter namespace */
namespace Bayesian_filter
{
	namespace FM = Bayesian_filter_matrix;


template <>
Filter_scheme<UD_filter>::Filter_scheme(FM::Subscript x_size, FM::Subscript q_size, FM::Subscript z_size) :
	UD_filter (x_size, q_size, z_size)
{}

}//namespace

