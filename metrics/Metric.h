
#ifndef METRIC_H_
#define METRIC_H_

#include "core/Configuration.h"

namespace metrics{

	class Metric{
		public:
			virtual double distance(Configuration*, Configuration*) = 0;
	};

}
#endif

