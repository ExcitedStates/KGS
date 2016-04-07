//
// Created by Rasmus Fonseca on 14/03/16.
//

#ifndef KGS_KGSEXCEPTIONS_H
#define KGS_KGSEXCEPTIONS_H

#include <exception>
#include <stdexcept>
#include <string>

class KGSRuntimeException: public std::runtime_error {
 public:
  KGSRuntimeException(const std::string& classSource, const std::string& msg);
};


#endif //KGS_KGSEXCEPTIONS_H
