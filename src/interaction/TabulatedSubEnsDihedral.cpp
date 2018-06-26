/*
  Copyright (C) 2016
      Jakub Krajniak (jkrajniak at gmail.com)
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
#include "TabulatedSubEnsDihedral.hpp"
#include "InterpolationLinear.hpp"
#include "InterpolationAkima.hpp"
#include "InterpolationCubic.hpp"
#include "FixedTripleListInteractionTemplate.hpp"
#include "FixedQuadrupleListInteractionTemplate.hpp"
#include "FixedQuadrupleListTypesInteractionTemplate.hpp"

namespace espressopp {
    namespace interaction {

        void TabulatedSubEnsDihedral::setFilenames(int dim,
            int itype, boost::python::list _filenames) {
            boost::mpi::communicator world;
            filenames.resize(dim);
            colVarRef.setDimension(dim);
            numInteractions = dim;
            for (int i=0; i<dim; ++i) {
              filenames[i] = boost::python::extract<std::string>(_filenames[i]);
              colVarRef[i].setDimension(4);
              if (itype == 1) { // create a new InterpolationLinear
                  tables[i] = make_shared <InterpolationLinear> ();
                  tables[i]->read(world, filenames[i].c_str());
              }

              else if (itype == 2) { // create a new InterpolationAkima
                  tables[i] = make_shared <InterpolationAkima> ();
                  tables[i]->read(world, filenames[i].c_str());
              }

              else if (itype == 3) { // create a new InterpolationCubic
                  tables[i] = make_shared <InterpolationCubic> ();
                  tables[i]->read(world, filenames[i].c_str());
              }
            }
        }

        void TabulatedSubEnsDihedral::addInteraction(int itype,
            boost::python::str fname, const RealND& _cvref) {
            boost::mpi::communicator world;
            int i = numInteractions;
            numInteractions += 1;
            colVarRef.setDimension(numInteractions);
            // Dimension 6: angle, bond, dihed, sd_angle, sd_bond, sd_dihed
            colVarRef[i].setDimension(6);
            colVarRef[i] = _cvref;
            filenames.push_back(boost::python::extract<std::string>(fname));
            weights.push_back(0.);
            weightSum.push_back(0.);
            targetProb.push_back(0.);
            if (itype == 1) { // create a new InterpolationLinear
                  tables.push_back(make_shared <InterpolationLinear> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
              else if (itype == 2) { // create a new InterpolationAkima
                  tables.push_back(make_shared <InterpolationAkima> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
              else if (itype == 3) { // create a new InterpolationCubic
                  tables.push_back(make_shared <InterpolationCubic> ());
                  tables[i]->read(world, filenames[i].c_str());
              }
        }

        void TabulatedSubEnsDihedral::setColVarRef(
            const RealNDs& cvRefs){
            // Set the reference values of the collective variables
            // aka cluster centers
            for (int i=0; i<numInteractions; ++i)
                colVarRef[i] = cvRefs[i];
        }

        void TabulatedSubEnsDihedral::computeColVarWeights(
            const Real3D& dist21, const Real3D& dist32,
            const Real3D& dist43, const bc::BC& bc){
            // Compute the weights for each force field
            // given the reference and instantaneous values of ColVars
            setColVar(dist21, dist32, dist43, bc);
            // Compute weights up to next to last FF
            real maxWeight = 0.;
            int maxWeightI = 0;
            // Check first whether we're stuck in a surface
            bool stuck = false;
            for (int i=0; i<numInteractions; ++i) {
                if (weights[i] > maxWeight) {
                    maxWeight = weights[i];
                    maxWeightI = i;
                }
            }
            if (weightCounts > 0 &&
                maxWeightI < numInteractions-1 &&
                weightSum[maxWeightI]/weightCounts < 0.98*targetProb[maxWeightI])
                stuck = true;
            if (!stuck) {
                maxWeight = 0.;
                maxWeightI = numInteractions-1;
                for (int i=0; i<numInteractions-1; ++i) {
                    weights[i]    = 1.;
                    real norm_d_i = 0.;
                    real norm_l_i = 0.;
                    for (int j=0; j<colVar.getDimension(); ++j) {
                        int k = 0;
                        // Choose between dihed, bond, and angle
                        if (j == 0) k = 0;
                        else if (j>0 && j<1+colVarBondList->size()) k = 1;
                        else k = 2;
                        norm_d_i += pow((colVar[j] -  colVarRef[i][k]) / colVarSd[k], 2);
                        norm_l_i += pow(colVarRef[i][3+k], 2);
                    }
                    if (norm_d_i > norm_l_i)
                      weights[i] = exp(- (sqrt(norm_d_i) - sqrt(norm_l_i)) / alpha);
                    if (weights[i] > maxWeight) {
                      maxWeight = weights[i];
                      maxWeightI = i;
                    }
                }
                for (int i=0; i<numInteractions-1; ++i) {
                    if (i != maxWeightI)
                        weights[i] = 0.;
                    else {
                        if (weightCounts > 0 &&
                            weights[i] > 0.01 &&
                            weightSum[i]/weightCounts < 0.98*targetProb[i]) {
                            weights[i] = 1.;
                            maxWeight = 1.;
                        }
                    }
                }
                if (maxWeightI == numInteractions-1)
                    maxWeight = 1.;
                weights[numInteractions-1] = 1. - maxWeight;
            }

            // Update weightSum
            for (int i=0; i<numInteractions; ++i)
                weightSum[i] += weights[i];
            weightCounts += 1;
        }

        // Collective variables
        void TabulatedSubEnsDihedral::setColVar(
            const Real3D& dist21, const Real3D& dist32,
            const Real3D& dist43, const bc::BC& bc) {
            colVar.setDimension(1+colVarBondList->size()+colVarAngleList->size());
            // compute phi
            real dist21_sqr = dist21 * dist21;
            real dist32_sqr = dist32 * dist32;
            real dist43_sqr = dist43 * dist43;
            real dist21_magn = sqrt(dist21_sqr);
            real dist32_magn = sqrt(dist32_sqr);
            real dist43_magn = sqrt(dist43_sqr);

            // cos0
            real sb1 = 1.0 / dist21_sqr;
            real sb2 = 1.0 / dist32_sqr;
            real sb3 = 1.0 / dist43_sqr;
            real rb1 = sqrt(sb1);
            real rb3 = sqrt(sb3);
            real c0 = dist21 * dist43 * rb1 * rb3;


            // 1st and 2nd angle
            real ctmp = dist21 * dist32;
            real r12c1 = 1.0 / (dist21_magn * dist32_magn);
            real c1mag = ctmp * r12c1;

            ctmp = (-1.0 * dist32) * dist43;
            real r12c2 = 1.0 / (dist32_magn * dist43_magn);
            real c2mag = ctmp * r12c2;


            //cos and sin of 2 angles and final cos
            real sin2 = 1.0 - c1mag * c1mag;
            if (sin2 < 0) sin2 = 0.0;
            real sc1 = sqrt(sin2);
            sc1 = 1.0 / sc1;

            sin2 = 1.0 - c2mag * c2mag;
            if (sin2 < 0) sin2 = 0.0;
            real sc2 = sqrt(sin2);
            sc2 = 1.0 / sc2;

            real s1 = sc1 * sc1;
            real s2 = sc2 * sc2;
            real s12 = sc1 * sc2;
            real c = (c0 + c1mag * c2mag) * s12;

            Real3D cc = dist21.cross(dist32);
            real cmag = sqrt(cc * cc);
            real dx = cc * dist43 / cmag / dist43_magn;

            if (c > 1.0) c = 1.0;
            else if (c < -1.0) c = -1.0;

            // phi
            real phi = acos(c);
            if (dx < 0.0) phi *= -1.0;
            colVar[0] = phi;
            // Now all bonds in colVarBondList
            int i=1;
            for (FixedPairList::PairList::Iterator it(*colVarBondList); it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              Real3D dist12;
              bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
              colVar[i] = sqrt(dist12 * dist12);
              i+=1;
            }
            // Now all angles in colVarAngleList
            for (FixedTripleList::TripleList::Iterator it(*colVarAngleList); it.isValid(); ++it) {
              Particle &p1 = *it->first;
              Particle &p2 = *it->second;
              Particle &p3 = *it->third;
              Real3D dist12, dist32;
              bc.getMinimumImageVectorBox(dist12, p1.position(), p2.position());
              bc.getMinimumImageVectorBox(dist32, p3.position(), p2.position());
              real dist12_sqr = dist12 * dist12;
              real dist32_sqr = dist32 * dist32;
              real dist1232 = sqrt(dist12_sqr) * sqrt(dist32_sqr);
              real cos_theta = dist12 * dist32 / dist1232;
              colVar[i] = acos(cos_theta);
              i+=1;
            }
        }


        typedef class FixedQuadrupleListInteractionTemplate <TabulatedSubEnsDihedral>
            FixedQuadrupleListTabulatedSubEnsDihedral;

        typedef class FixedQuadrupleListTypesInteractionTemplate<TabulatedSubEnsDihedral>
            FixedQuadrupleListTypesTabulatedSubEnsDihedral;

        //////////////////////////////////////////////////
        // REGISTRATION WITH PYTHON
        //////////////////////////////////////////////////
        void TabulatedSubEnsDihedral::registerPython() {
            using namespace espressopp::python;

            class_ <TabulatedSubEnsDihedral, bases <DihedralPotential> >
                ("interaction_TabulatedSubEnsDihedral", init <>())
                .def("dimension_get", &TabulatedSubEnsDihedral::getDimension)
                .def("filenames_get", &TabulatedSubEnsDihedral::getFilenames)
                .def("filename_get", &TabulatedSubEnsDihedral::getFilename)
                .def("filename_set", &TabulatedSubEnsDihedral::setFilename)
                .def("targetProb_get", &TabulatedSubEnsDihedral::getTargetProb)
                .def("targetProb_set", &TabulatedSubEnsDihedral::setTargetProb)
                .def("colVarMu_get", &TabulatedSubEnsDihedral::getColVarMus)
                .def("colVarMu_set", &TabulatedSubEnsDihedral::setColVarMu)
                .def("colVarSd_get", &TabulatedSubEnsDihedral::getColVarSds)
                .def("colVarSd_set", &TabulatedSubEnsDihedral::setColVarSd)
                .def("weight_get", &TabulatedSubEnsDihedral::getWeights)
                .def("weight_set", &TabulatedSubEnsDihedral::setWeight)
                .def("alpha_get", &TabulatedSubEnsDihedral::getAlpha)
                .def("alpha_set", &TabulatedSubEnsDihedral::setAlpha)
                .def("addInteraction", &TabulatedSubEnsDihedral::addInteraction)
                .def("colVarRefs_get", &TabulatedSubEnsDihedral::getColVarRefs)
                .def("colVarRef_get", &TabulatedSubEnsDihedral::getColVarRef)
                .def_pickle(TabulatedSubEnsDihedral_pickle())
                ;

            class_ <FixedQuadrupleListTabulatedSubEnsDihedral, bases <Interaction> >
                ("interaction_FixedQuadrupleListTabulatedSubEnsDihedral",
                        init <shared_ptr<System>,
                              shared_ptr<FixedQuadrupleList>,
                              shared_ptr<TabulatedSubEnsDihedral> >())
                .def("setPotential", &FixedQuadrupleListTabulatedSubEnsDihedral::setPotential)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTabulatedSubEnsDihedral::getFixedQuadrupleList);

            class_< FixedQuadrupleListTypesTabulatedSubEnsDihedral, bases< Interaction > >
                ("interaction_FixedQuadrupleListTypesTabulatedSubEnsDihedral",
                 init< shared_ptr<System>, shared_ptr<FixedQuadrupleList> >())
                .def("setPotential", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::setPotential)
                .def("getPotential", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::getPotentialPtr)
                .def("setFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::setFixedQuadrupleList)
                .def("getFixedQuadrupleList", &FixedQuadrupleListTypesTabulatedSubEnsDihedral::getFixedQuadrupleList);

        }

    } // ns interaction
} // ns espressopp
