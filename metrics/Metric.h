
#ifndef METRIC_H_
#define METRIC_H_

#include "core/Configuration.h"
#include "Selection.h"

namespace metrics{

class Metric{
 public:
  Metric(Selection&);

  virtual double distance(Configuration*, Configuration*) = 0;

  Selection& getSelection();

 protected:
  Selection& m_selection;
};

}
#endif

