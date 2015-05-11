/*
 * Open-source multi-Physics Phenomena Analyzer (OP2A) ver. 0.1
 *
 * 		Copyright (c) 2015 MINKWAN KIM
 *
 * 	Initial Developed Date: May 11, 2015
 *      			Author: Minkwan Kim
 *
 * NonCopyable.hpp
 * 			-  
 *  
 */
#ifndef NONCOPYABLE_HPP_
#define NONCOPYABLE_HPP_

namespace OP2A{
namespace Common{

template < typename TYPE >
class NonCopyable {
public:
  /// Default inline constructor
  NonCopyable () {}

  /// Default inline destructor
  virtual ~NonCopyable () {}

private:
  /// private (non defined) copy constructor to prevent
  /// copying of the object
  NonCopyable (const NonCopyable & Source);

  /// private (non defined) assignment operator to prevent
  /// copy assignment of the object
  const NonCopyable & operator = (const NonCopyable & Source);
}; // end class NonCopyable


}
}



#endif /* NONCOPYABLE_HPP_ */
