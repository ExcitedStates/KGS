//
// Created by Rasmus Fonseca on 14/03/16.
//

#include "KGSExceptions.h"

KGSRuntimeException::KGSRuntimeException(const std::string& classSource, const std::string& msg):
    runtime_error(classSource+" - "+msg)
{}
