#define RAND_MAX
#define NUMBER_OF_PARTICLES 33000
#define BOX_SIZE 200

#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>

#include "Particle.hpp"
#include "LennardJonesInteraction.hpp"
#include "FullPairIterator.hpp"

template<class Interaction, class PairIterator>
void compute (std::vector<Particle> *pc, Interaction lj, PairIterator* it)

{
  double rsq;
  double en;

  PairIterator it1 = *it;

  en = 0.0;

  //loop over all pairs assuming non-periodic BCs

  for (it1.reset(); !it1.done(); it1.next()) {
     Particle* Pi = &(*pc)[it1.first()];
     Particle* Pj = &(*pc)[it1.second()];
     rsq   = pow(Pi->getx() - Pj->getx(), 2);
     rsq  += pow(Pi->gety() - Pj->gety(), 2);
     rsq  += pow(Pi->getz() - Pj->getz(), 2);
     en += lj.computeLennardJonesEnergy(rsq);
  }

  //write out the total LJ energy
  std::cout << "en = " << en << std::endl;

}

int main() {

  //variables for storing random numbers
  double rx;
  double ry;
  double rz;

  //create a vector to store the particles
  std::vector<Particle> pc(NUMBER_OF_PARTICLES);

  //assign random positions to the particles on r[0, BOX_SIZE]
  for(int i = 0; i < pc.size(); i++) {
    rx = BOX_SIZE * double(rand()) / RAND_MAX;
    ry = BOX_SIZE * double(rand()) / RAND_MAX;
    rz = BOX_SIZE * double(rand()) / RAND_MAX;
    pc[i] = Particle(rx, ry, rz);
  }

  //print a particle to standard output as a test
  std::cout << pc[1].toString() << std::endl;

  //create a LJ interaction and set its cutoff
  LennardJonesInteraction lj = LennardJonesInteraction();
  lj.setCutoff(2.5);

  FullPairIterator it = FullPairIterator(pc.size());

  compute<LennardJonesInteraction, FullPairIterator>(&pc, lj, &it);

  return 0;
}

