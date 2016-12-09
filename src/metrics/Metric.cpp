#include "Metric.h"

namespace metrics {

Metric::Metric( Selection &selection ) :
    m_selection(selection)
{ }

Selection& Metric::getSelection(){
  return m_selection;
}

}
