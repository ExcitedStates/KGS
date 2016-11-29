
#ifndef DIHEDRAL_H_
#define DIHEDRAL_H_

#include <string>

#include "Metric.h"
#include "core/Configuration.h"

namespace metrics {

class Dihedral : public Metric {
 public:
  Dihedral(Selection &selection);

  double distance(Configuration *, Configuration *);

};

}
#endif
