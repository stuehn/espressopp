/*
  Copyright (C) 2019-2021
      Max Planck Institute for Polymer Research & JGU Mainz
  Copyright (C) 2012,2013,2014,2015,2016,2017,2018
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

// ESPP_CLASS
#ifndef VEC_INTERACTION_VERLETLISTLENNARDJONES_HPP
#define VEC_INTERACTION_VERLETLISTLENNARDJONES_HPP

#include "types.hpp"
#include "interaction/Interaction.hpp"
#include "Real3D.hpp"
#include "Tensor.hpp"
#include "Particle.hpp"
#include "esutil/Array2D.hpp"
#include "bc/BC.hpp"
#include "storage/Storage.hpp"

#include "LennardJones.hpp"
#include "VerletListInteractionTemplate.hpp"
#include "vec/VerletList.hpp"
#include "vec/Vectorization.hpp"

namespace espressopp { namespace vec {
  namespace interaction {

    class VerletListLennardJones
      : public VerletListInteractionTemplate<LennardJones>
    {
    protected:
      typedef LennardJones _Potential;
      typedef _Potential Potential;
      typedef VerletListInteractionTemplate<LennardJones> base;

      struct LJCoefficients
      {
        LJCoefficients(real const& ff1, real const& ff2): ff1(ff1), ff2(ff2){}
        LJCoefficients(){}
        real ff1, ff2;
      };

    public:
      VerletListLennardJones(shared_ptr<VerletList> _verletList)
        : base(_verletList)
        , np_types(0)
        , p_types(0)
      {

      }

      void rebuildPotential()
      {
        np_types = potentialArray.size_n();
        p_types = potentialArray.size_m();
        ffs = AlignedVector<LJCoefficients>(np_types*p_types);
        cutoffSqr = AlignedVector<real>(np_types*p_types);
        AlignedVector<LJCoefficients>::iterator it1 = ffs.begin();
        AlignedVector<real>::iterator it3 = cutoffSqr.begin();
        for(auto& p: potentialArray){
          *(it1++) = LJCoefficients(p.getff1(),p.getff2());
          *(it3++) = p.getCutoffSqr();
        }
        needRebuildPotential = false;
      }
      virtual void addForces();

    protected:

      template <bool ONETYPE, bool VEC_MODE_AOS>
      void addForces_impl(
        ParticleArray & particleArray,
        VerletList::NeighborList const& neighborList);

      size_t np_types, p_types;
      AlignedVector<LJCoefficients> ffs;
      AlignedVector<real> cutoffSqr;
      bool needRebuildPotential = true;
    };

    //////////////////////////////////////////////////
    // INLINE IMPLEMENTATION
    //////////////////////////////////////////////////
    inline void
    VerletListLennardJones::
    addForces()
    {
      LOG4ESPP_DEBUG(_Potential::theLogger, "loop over verlet list pairs and add forces");

      // lookup table for LJ variables
      // ideal for low number of types vs number of particle pairs
      // trigger rebuild on setPotential and on size modification from getPotential
      int vlmaxtype = verletList->getMaxType();
      Potential max_pot = getPotential(vlmaxtype,vlmaxtype);
      if(needRebuildPotential) rebuildPotential();
      bool VEC_MODE_AOS = verletList->getVectorization()->modeAOS();

      auto& pa = verletList->getVectorization()->particles;
      const auto& nl = verletList->getNeighborList();

      if(np_types==1 && p_types==1){
        if(VEC_MODE_AOS){
          addForces_impl<true,true>(pa, nl);
        }
        else{
          addForces_impl<true,false>(pa, nl);
        }
      }
      else {
        if(VEC_MODE_AOS){
          addForces_impl<false,true>(pa, nl);
        }
        else{
          addForces_impl<false,false>(pa, nl);
        }
      }
    }

    template <bool ONETYPE, bool VEC_MODE_AOS>
    inline void
    VerletListLennardJones::
    addForces_impl(ParticleArray & particleArray, VerletList::NeighborList const& neighborList)
    {
      using namespace vec::storage;

      {
        real ff1_, ff2_, cutoffSqr_;
        if(ONETYPE){
          ff1_ = ffs[0].ff1;
          ff2_ = ffs[0].ff2;
          cutoffSqr_ = cutoffSqr[0];
        }

        const Real3DInt *pa_pos  = particleArray.position.data();
        Real4D *pa_force         = particleArray.force.data();

        const lint* __restrict pa_type = particleArray.type.data();
        const real* __restrict pa_p_x  = particleArray.p_x.data();
        const real* __restrict pa_p_y  = particleArray.p_y.data();
        const real* __restrict pa_p_z  = particleArray.p_z.data();
        real* __restrict pa_f_x        = particleArray.f_x.data();
        real* __restrict pa_f_y        = particleArray.f_y.data();
        real* __restrict pa_f_z        = particleArray.f_z.data();

        const auto* __restrict plist   = neighborList.plist.data();
        const auto* __restrict prange  = neighborList.prange.data();
        const auto* __restrict nplist  = neighborList.nplist.data();
        const int ip_max = neighborList.plist.size();

        for(int ip=0; ip<ip_max; ip++)
        {
          int p = plist[ip];
          int p_lookup;
          real p_x, p_y, p_z;

          if(VEC_MODE_AOS)
          {
            p_x = pa_pos[p].x;
            p_y = pa_pos[p].y;
            p_z = pa_pos[p].z;
            if(!ONETYPE){
              p_lookup  = pa_pos[p].t*np_types;
            }
          }
          else
          {
            if(!ONETYPE){
              p_lookup  = pa_type[p]*np_types;
            }
            p_x = pa_p_x[p];
            p_y = pa_p_y[p];
            p_z = pa_p_z[p];
          }

          real f_x = 0.0;
          real f_y = 0.0;
          real f_z = 0.0;

          const int in_min=prange[ip].first;
          const int in_max=prange[ip].second;

          #pragma vector always
          #pragma vector aligned
          #pragma ivdep
          for(int in=in_min; in<in_max; in++)
          {
            auto np_ii = nplist[in];
            {
              int np_lookup;
              real dist_x, dist_y, dist_z;
              if(VEC_MODE_AOS)
              {
                dist_x   = p_x - pa_pos[np_ii].x;
                dist_y   = p_y - pa_pos[np_ii].y;
                dist_z   = p_z - pa_pos[np_ii].z;
                if(!ONETYPE) np_lookup = pa_pos[np_ii].t+p_lookup;
              }
              else
              {
                dist_x   = p_x - pa_p_x[np_ii];
                dist_y   = p_y - pa_p_y[np_ii];
                dist_z   = p_z - pa_p_z[np_ii];
                if(!ONETYPE) np_lookup = pa_type[np_ii]+p_lookup;
              }

              real distSqr  = dist_x*dist_x + dist_y*dist_y + dist_z*dist_z;
              if(!ONETYPE){
                cutoffSqr_ = cutoffSqr[np_lookup];
              }

              #if defined(ESPP_VECTOR_MASK)
              if(distSqr <= cutoffSqr_)
              #endif
              {
                real frac2 = 1.0 / distSqr;
                real frac6 = frac2 * frac2 * frac2;
                real ffactor;

                if(ONETYPE)
                  ffactor = ff1_ * frac6 - ff2_;
                else
                  ffactor = ffs[np_lookup].ff1 * frac6 - ffs[np_lookup].ff2;

                #if !defined(ESPP_VECTOR_MASK)
                if(distSqr > cutoffSqr_) ffactor = 0.0;
                #endif

                ffactor = frac6 * ffactor * frac2;

                f_x += dist_x * ffactor;
                f_y += dist_y * ffactor;
                f_z += dist_z * ffactor;

                if(VEC_MODE_AOS)
                {
                  auto& np_force = pa_force[np_ii];
                  np_force.x -= dist_x * ffactor;
                  np_force.y -= dist_y * ffactor;
                  np_force.z -= dist_z * ffactor;
                }
                else
                {
                  pa_f_x[np_ii] -= dist_x * ffactor;
                  pa_f_y[np_ii] -= dist_y * ffactor;
                  pa_f_z[np_ii] -= dist_z * ffactor;
                }
              }
            }
          }
          if(VEC_MODE_AOS)
          {
            auto& p_force = pa_force[p];
            p_force.x += f_x;
            p_force.y += f_y;
            p_force.z += f_z;
          }
          else
          {
            pa_f_x[p] += f_x;
            pa_f_y[p] += f_y;
            pa_f_z[p] += f_z;
          }
        }
      }
    }
  }
}}

#endif//VEC_INTERACTION_VERLETLISTLENNARDJONES_HPP
