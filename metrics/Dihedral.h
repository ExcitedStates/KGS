
#ifndef DIHEDRAL_H_
#define DIHEDRAL_H_

#include <string>

#include "metrics/Metric.h"
#include "core/Configuration.h"

namespace metrics{

	class Dihedral: public Metric{
		public:
			Dihedral(): Dihedral("ALL") {}
			Dihedral(const std::string& atom_selection);

			double distance(Configuration*, Configuration*);

		private:
			const std::string atom_selection;

	};
}
#endif
