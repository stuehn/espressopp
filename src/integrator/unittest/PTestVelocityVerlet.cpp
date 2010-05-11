#define PARALLEL_TEST_MODULE VelocityVerlet
#include "ut.hpp"

#include "types.hpp"
#include "mpi.hpp"
#include "logging.hpp"
#include "esutil/RNG.hpp"
#include "storage/DomainDecomposition.hpp"
#include "interaction/LennardJones.hpp"
#include "integrator/VelocityVerlet.hpp"
#include "iterator/CellListIterator.hpp"
#include "bc/OrthorhombicBC.hpp"
#include "System.hpp"
#include "VerletList.hpp"
#include "Real3D.hpp"

using namespace espresso;
using namespace storage;
using namespace interaction;
using namespace integrator;
using namespace iterator;
using namespace bc;
using namespace esutil;

struct LoggingFixture {  
  LoggingFixture() { 
    LOG4ESPP_CONFIGURE();
    log4espp::Logger::getRoot().setLevel(log4espp::Logger::WARN);
    // log4espp::Logger::getInstance("VerletList").setLevel(log4espp::Logger::WARN);
    // log4espp::Logger::getInstance("Storage").setLevel(log4espp::Logger::DEBUG);
    log4espp::Logger::getInstance("MDIntegrator").setLevel(log4espp::Logger::INFO);
  }
};

static real cutoff1 = 1.4;
static real cutoff2 = 1.0;
static int N = 3;
static real size = N;

BOOST_GLOBAL_FIXTURE(LoggingFixture);

struct Fixture {
  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  Fixture() {

    int N = 3;  // number of particles in each dimension

    real density = 1.0;

    real SIZE = N * density;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = " << density);

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;

    ConstReal3DRef boxLRef(boxL);

    int nodeGrid[3] = { 1, 1, mpiWorld->size() };
    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {
      int ncells = 1;
       while (SIZE / (ncells * nodeGrid[i]) >= cutoff1) {
         ncells ++;
       }
       ncells--;
       cellGrid[i] = ncells;
    }

    BOOST_MESSAGE("ncells in each dim / proc: " << cellGrid[0] << " x " <<
                                 cellGrid[1] << " x " << cellGrid[2]);

    system = make_shared< System >();
    system->rng = make_shared< RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxLRef);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid);
    esutil::RNG rng;

    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          int m = (i + 2*j + 3*k) % 11;
          real r = 0.45 + m * 0.01;
          real x = (i + r) / N * SIZE;
          real y = (j + r) / N * SIZE; 
          real z = (k + r) / N * SIZE;
          real pos[3] = { x, y, z };
      
          if (mpiWorld->rank() == 0) {
          
            printf("add particle at %f %f %f\n", x, y, z);
            domdec->addParticle(id, pos);

          }

          id++;
        }
      }
    }

    system->storage = domdec;

    BOOST_MESSAGE("number of particles in storage =  " << 
                  domdec->getNRealParticles());
  }
};

// sets up a storage with particles in a lattice

struct DomainFixture {

  shared_ptr<DomainDecomposition> domdec;
  shared_ptr<System> system;

  real SIZE;

  DomainFixture(double density) {

    SIZE = pow(N * N * N / density, 1.0/3.0) ;

    BOOST_MESSAGE("box SIZE = " << SIZE << ", density = 1.0");

    Real3D boxL(SIZE, SIZE, SIZE);

    real skin   = 0.3;

    int nodeGrid[3] = { 1, 1, mpiWorld->size() };

    int cellGrid[3] = { 1, 1, 1 };

    for (int i = 0; i < 3; i++) {

       // take largest number of cells, but cell size must be >= cutoff
       // so cutoff used here should be the largest cutoff

       int ncells = 1;

       while (SIZE / (ncells * nodeGrid[i]) >= cutoff2) {
          ncells ++;
       }

       ncells--;

       cellGrid[i] = ncells;
    }

    BOOST_MESSAGE("ncells in each dim / proc: " << cellGrid[0] << " x " <<
                                 cellGrid[1] << " x " << cellGrid[2]);

    system = make_shared< System >();
    system->rng = make_shared< RNG >();
    system->bc = make_shared< OrthorhombicBC >(system->rng, boxL);
    system->skin = skin;
    domdec = make_shared< DomainDecomposition >(system,
                                                mpiWorld,
                                                nodeGrid,
                                                cellGrid);
    system->storage = domdec;

  }
};

// sets up a storage with particles in a lattice

struct LatticeFixture : DomainFixture {

  LatticeFixture() : DomainFixture(1.0) {

    int id = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        for (int k = 0; k < N; k++) {

          real r = 0.5;
          real x = (i + r) / N * SIZE;
          real y = (j + r) / N * SIZE;
          real z = (k + r) / N * SIZE;
          real pos[3] = { x, y, z };

          if (mpiWorld->rank() == 0) {
            BOOST_MESSAGE("add particle at pos " << x << " " << y << " " << z);
            domdec->addParticle(id, pos);
          }

          id++;
        }
      }
    }

    CellList realCells = domdec->getRealCells();

    // loop over all particles of the real cells and set velocity in x - direction

    for (CellListIterator cit(realCells); !cit.isDone(); ++cit) {
      cit->m.v[0] = 1.0;
    }

    BOOST_MESSAGE("number of lattice particles in storage = " <<
                   domdec->getNRealParticles());

    domdec->resortParticles();
  }
};

BOOST_FIXTURE_TEST_CASE(moveParticles, LatticeFixture)
{
  // define a 0 potential, but with a cutoff

  BOOST_MESSAGE("starting to build verlet lists");

  shared_ptr<MDIntegrator> integrator = 
     make_shared<VelocityVerlet>(system);

  integrator->setTimeStep(0.01);

  int niter = N * 100;

  BOOST_MESSAGE("run " << N << " iterations with timestep = 0.01");

  integrator->run(niter);

  // Due to velocity 1.0 particles should be around 

  CellList realCells = system->storage->getRealCells();

  for(CellListIterator cit(realCells); !cit.isDone(); ++cit) {

    int id = cit->p.id;

    int z = id % N;
    int y = (id - z) / N % N;
    int x = (id - N * y - z);  x /= N * N; x %= N;

    BOOST_CHECK_SMALL((x + 0.5) / N * size - cit->r.p[0], 1e-6);
    BOOST_CHECK_SMALL((y + 0.5) / N * size - cit->r.p[1], 1e-6);
    BOOST_CHECK_SMALL((z + 0.5) / N * size - cit->r.p[2], 1e-6);
  }

  // make a final reduction to synchronize the processors

  int val = 0;

  boost::mpi::all_reduce(*mpiWorld, val, val, std::plus<int>());

  printf("ready\n");
}