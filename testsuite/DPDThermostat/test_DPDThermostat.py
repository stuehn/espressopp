import sys
import time
import espressopp
import numpy as np
import mpi4py.MPI as MPI

import unittest

class TestDPDThermostat(unittest.TestCase):
    def setUp(self):
        # set up system
        system = espressopp.System()
        rng = espressopp.esutil.RNG()
        rng.seed(1)
        system.rng = rng
        box = (10, 10, 10)
        system.bc = espressopp.bc.OrthorhombicBC(system.rng, box)
        system.skin = 0.3
        system.comm = MPI.COMM_WORLD
        self.system = system

    def test_normal(self):
        # set up normal domain decomposition
        nodeGrid = espressopp.tools.decomp.nodeGrid(espressopp.MPI.COMM_WORLD.size)
        cellGrid = espressopp.tools.decomp.cellGrid((10, 10, 10), nodeGrid, 1.5, self.system.skin)
        self.system.storage = espressopp.storage.DomainDecomposition(self.system, nodeGrid, cellGrid)

        # add some particles (normal, coarse-grained particles only)
        x = [ espressopp.Real3D(5.5, 5.0, 5.0),
              espressopp.Real3D(6.0, 5.0, 5.0) ]

        v = [ espressopp.Real3D(0.5,0.25,0.25),
              espressopp.Real3D(-0.25,0.5,-0.25) ]

        particle_list = [
            (1, 1, x[0], v[0], 1.0),
            (2, 1, x[1], v[1], 1.0)
        ]
        self.system.storage.addParticles(particle_list, 'id', 'type', 'pos', 'v', 'mass')
        self.system.storage.decompose()

        # generate a verlet list
        vl = espressopp.VerletList(self.system, cutoff=1.5)

        # integrator
        integrator = espressopp.integrator.VelocityVerlet(self.system)
        integrator.dt = 0.01

        # DPD Thermostat
        dpd = espressopp.integrator.DPDThermostat(self.system,vl)
        dpd.gamma = 2.0
        dpd.tgamma = 5.0
        dpd.temperature = 2.0
        integrator.addExtension(dpd)
        
        #This is the first step, which is always a heat-up step, hence the factor of sqrt(3)
        noise_pref = np.sqrt( 3.0 * 24.0 * dpd.temperature / integrator.dt )

        # coordinates of particles before integration
        dist = espressopp.Real3D( x[0] - x[1] )
        dist_norm = np.sqrt( dist*dist )
        dist_unitv = dist / dist_norm
        veldiff = espressopp.Real3D( v[0] - v[1] )

        omega = 1.0 - dist_norm / 1.5
        
        f_expected = [ espressopp.Real3D(0.0),
                       espressopp.Real3D(0.0) ]
        
       
        #longitudinal contribution
        r0 = self.system.rng()-0.5
        
        f_damp = ( dist_unitv*veldiff ) * dpd.gamma * omega * omega 
        f_noise =  noise_pref * np.sqrt(dpd.gamma) * omega * r0 
        
        f_expected[0] += (f_noise - f_damp) * dist_unitv
        f_expected[1] -= (f_noise - f_damp) * dist_unitv
       
        #transversal contribution
        tf_damp = espressopp.Real3D(0.0)

        tf_damp[0] =   (1.0 - dist_unitv[0]*dist_unitv[0])*veldiff[0] \
                     - dist_unitv[0]*dist_unitv[1]*veldiff[1]         \
                     - dist_unitv[0]*dist_unitv[2]*veldiff[2] 
        tf_damp[1] =   (1.0 - dist_unitv[1]*dist_unitv[1])*veldiff[1] \
                     - dist_unitv[1]*dist_unitv[0]*veldiff[0]         \
                     - dist_unitv[1]*dist_unitv[2]*veldiff[2] 
        tf_damp[2] =   (1.0 - dist_unitv[2]*dist_unitv[2])*veldiff[2] \
                     - dist_unitv[2]*dist_unitv[0]*veldiff[0]         \
                     - dist_unitv[2]*dist_unitv[1]*veldiff[1] 

        
        tf_noise = espressopp.Real3D(0.0)
        
        r1 = self.system.rng()-.5
        r2 = self.system.rng()-.5
        r3 = self.system.rng()-.5

        #the order of evaluation is left to right in python but probably right to left in the c++ implementation
        #therefore we switch here to get equal results
        randvec = espressopp.Real3D(r3,r2,r1)

        tf_noise[0] =   (1.0 - dist_unitv[0]*dist_unitv[0])*randvec[0] \
                      - dist_unitv[0]*dist_unitv[1]*randvec[1]         \
                      - dist_unitv[0]*dist_unitv[2]*randvec[2] 
        tf_noise[1] =   (1.0 - dist_unitv[1]*dist_unitv[1])*randvec[1] \
                      - dist_unitv[1]*dist_unitv[0]*randvec[0]         \
                      - dist_unitv[1]*dist_unitv[2]*randvec[2] 
        tf_noise[2] =   (1.0 - dist_unitv[2]*dist_unitv[2])*randvec[2] \
                      - dist_unitv[2]*dist_unitv[0]*randvec[0]         \
                      - dist_unitv[2]*dist_unitv[1]*randvec[1] 

        f_expected[0] += (  noise_pref * np.sqrt(dpd.tgamma) * omega * tf_noise\
                          - dpd.tgamma * omega * omega * tf_damp )
        f_expected[1] -= (  noise_pref * np.sqrt(dpd.tgamma) * omega * tf_noise\
                          - dpd.tgamma * omega * omega * tf_damp )

        # calculate forces
        self.system.rng.seed(1)
        integrator.run(0)

        f_result = [ self.system.storage.getParticle(1).f,
                     self.system.storage.getParticle(2).f ]


        #the checks follow

        #Do we obey Sir Isaac?
        self.assertEqual(f_result[0],-1.0*f_result[1])

        #Check if the forces are (almost) equal
        self.assertAlmostEqual(f_expected[0][0],f_result[0][0],places=5)
        self.assertAlmostEqual(f_expected[0][1],f_result[0][1],places=5)
        self.assertAlmostEqual(f_expected[0][2],f_result[0][2],places=5)
        
        self.assertAlmostEqual(f_expected[1][0],f_result[1][0],places=5)
        self.assertAlmostEqual(f_expected[1][1],f_result[1][1],places=5)
        self.assertAlmostEqual(f_expected[1][2],f_result[1][2],places=5)
        

if __name__ == '__main__':
    unittest.main()
