#  Copyright (C) 2014
#      Pierre de Buyl
#  Copyright (C) 2012,2013
#      Max Planck Institute for Polymer Research
#  Copyright (C) 2008,2009,2010,2011
#      Max-Planck-Institute for Polymer Research & Fraunhofer SCAI
#
#  This file is part of ESPResSo++.
#
#  ESPResSo++ is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  ESPResSo++ is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

r"""
*****************************************
espressopp.interaction.nLennardJones104Wall30
*****************************************

This class defines a nnLennard-Jones 10-4 SingleParticlePotential in the direction x.

.. math:: V(r) = \epsilon \left( 0.2*\left(\frac{\sigma}{r}\right)^10 - 0.5* \left(\frac{\sigma}{r}\right)^4 \right)

where :math:`r` is the distance from the lower or upper wall in the x
direction. :math:`V(r)=0` after a distance `sigmaCutoff`.

The parameters have to be defined for every species present in the system with
`setParams` and can be retrieved with `getParams`.

# Slab.BC does not work properly, so size of simulation box is set to Lx0+3 but the wall is located at Lx0

Example:

    >>> nLJ104 = espressopp.interaction.nnLennardJones104Wall30()
    >>> nLJ104.setParams(0, eps_wall, 1., wall_cutoff, r0=0.0,Lx0)   
    >>> SPnLJ104 = espressopp.interaction.SingleParticlennLennardJones104Wall30(system, nLJ104)
    >>> system.addInteraction(SPnLJ104)


.. function:: espressopp.interaction.nLennardJones104Wall30()


.. function:: espressopp.interaction.nLennardJones104Wall30.getParams(type_var)

                :param type_var:
                :type type_var:
                :rtype:

.. function:: espressopp.interaction.nnLennardJones104Wall30.setParams(type_var, epsilon, sigma, sigmaCutoff, r0,Lx0)

                :param type_var:
                :param epsilon:
                :param sigma:
                :param sigmaCutoff:
                :param r0
                :param Lx0:
                :type type_var:
                :type epsilon:
                :type sigma:
                :type sigmaCutoff:
                :type r0:
                :type Lx0:


.. function:: espressopp.interaction.SingleParticlenLennardJones104Wall30(system, potential)

                :param system:
                :param potential:
                :type system:
                :type potential:

.. function:: espressopp.interaction.SingleParticlenLennardJones104Wall30.setPotential(potential)

                :param potential:
                :type potential:
"""
from espressopp import pmi
from espressopp.esutil import *

from espressopp.interaction.SingleParticlePotential import *
from espressopp.interaction.Interaction import *
from _espressopp import interaction_nLennardJones104Wall30, interaction_SingleParticlenLennardJones104Wall30


class nLennardJones104Wall30Local(SingleParticlePotentialLocal, interaction_nLennardJones104Wall30):

    def __init__(self):

        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_nLennardJones104Wall30)
    def getParams(self, type_var):


        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            return self.cxxclass.getParams(self, type_var)

    def setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0, Lx0):


        self.cxxclass.setParams(self, type_var, epsilon, sigma, sigmaCutoff, r0, Lx0)

class SingleParticlenLennardJones104Wall30Local(InteractionLocal, interaction_SingleParticlenLennardJones104Wall30):

    def __init__(self, system, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            cxxinit(self, interaction_SingleParticlenLennardJones104Wall30, system, potential)

    def setPotential(self, potential):
        if not (pmi._PMIComm and pmi._PMIComm.isActive()) or pmi._MPIcomm.rank in pmi._PMIComm.getMPIcpugroup():
            self.cxxclass.setPotential(self, potential)

if pmi.isController:
    class nLennardJones104Wall30(SingleParticlePotential):
        'The nLennardJones104Wall30 potential.'
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.nLennardJones104Wall30Local',
            pmicall = ['setParams', 'getParams']
            )

    class SingleParticlenLennardJones104Wall30(Interaction, metaclass=pmi.Proxy):
        pmiproxydefs = dict(
            cls = 'espressopp.interaction.SingleParticlenLennardJones104Wall30Local',
            pmicall = ['setPotential']
            )
