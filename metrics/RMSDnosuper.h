#ifndef RMSD_H_
#define RMSD_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <string>
#include <list>

#include "Metric.h"
#include "core/Configuration.h"
#include "core/Molecule.h"

namespace metrics{


class RMSDnosuper: public Metric{
 public:
  RMSDnosuper(Selection& selection);

  double distance(Configuration*, Configuration*);

};


}
#endif
