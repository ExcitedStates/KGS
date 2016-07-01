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


class RMSD: public Metric{
 public:
  RMSD();
  RMSD(const std::vector<Atom*>* atomsRMSD);

  double distance(Configuration*, Configuration*);

 private:
  const std::vector<Atom*>* m_atomsRMSD;

};


}
#endif
