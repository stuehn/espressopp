/*
  Copyright (C) 2014
      Pierre de Buyl
  Copyright (C) 2012,2013
      Max Planck Institute for Polymer Research
  Copyright (C) 2008,2009,2010,2011
      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI

  This file is part of ESPResSo++.

  ESPResSo++ is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "python.hpp"
#include "LennardJones_sphwall30.hpp"
#include "Real3D.hpp"

namespace espressopp {
  namespace interaction {

    typedef class SingleParticleInteractionTemplate <LennardJones_sphwall30>
    SingleParticleLennardJones_sphwall30;

    //////////////////////////////////////////////////
    // REGISTRATION WITH PYTHON
    //////////////////////////////////////////////////
    void
    LennardJones_sphwall30::registerPython() {
      using namespace espressopp::python;

      class_< LennardJones_sphwall30, bases< SingleParticlePotential > >
        ("interaction_LennardJones_sphwall30", init<>())
        .def("setParams", &LennardJones_sphwall30::setParams)
        .def("getParams", &LennardJones_sphwall30::getParams)
        ;

      class_< SingleParticleLennardJones_sphwall30, bases< Interaction > >
        ("interaction_SingleParticleLennardJones_sphwall30", init< shared_ptr<System>, shared_ptr<LennardJones_sphwall30> >())
        .def("setPotential", &SingleParticleLennardJones_sphwall30::setPotential)
        .def("getPotential", &SingleParticleLennardJones_sphwall30::getPotential)
       ;
    }
  }
}
