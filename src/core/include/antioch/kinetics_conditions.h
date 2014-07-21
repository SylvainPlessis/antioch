//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_KINETICS_CONDITIONS_H
#define ANTIOCH_KINETICS_CONDITIONS_H

// Antioch
#include "antioch/antioch_asserts.h"
#include "antioch/particle_flux.h"

// C++
#include <vector>
#include <map>

namespace Antioch{

  /*! This class contains the conditions of the chemistry

      Conditions are:
        - Temperature(s)
        - particle flux
        - any other not yet implemented conditions (catalysis for instance)

      Pressure is given by the molecular composition of the mixture,
      we are in a gaz phase, so the state equation will give pressure (typically
      ideal gas).

      The idea is to avoid any copy of anything as much as we can, 
      we store pointers but deal only with references. Nothing
      belongs to this object, if there is any cleaning to do,
      it must be done elsewhere.

      Might be interesting to see about the "double" case
      for temperature: faster to copy instead of copying
      the reference, and it's expected to be the classic
      case.
   */
  template <typename StateType, 
            typename VectorStateType = std::vector<StateType> >
  class KineticsConditions
  {
        public:

          KineticsConditions();
          KineticsConditions(const StateType & temperature);
          ~KineticsConditions();

          void add_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr);

          void change_temperature(const StateType & temperature);

          void change_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr);

          const StateType & T() const;

          const ParticleFlux<VectorStateType> & particle_flux(int nr) const;

        private:
        // pointer's not const, temperature is
          StateType const * _temperature; 
        // pointer's not const, particle flux is
          std::map<unsigned int,ParticleFlux<VectorStateType> const * > _map_pf; 

  };

  template <typename StateType, typename VectorStateType>
  inline
  KineticsConditions<StateType,VectorStateType>::KineticsConditions():
        _temperature(NULL)
  {
    return;
  }

  template <typename StateType, typename VectorStateType>
  inline
  KineticsConditions<StateType,VectorStateType>::KineticsConditions(const StateType & temperature):
        _temperature(&temperature)
  {
    return;
  }

  template <typename StateType, typename VectorStateType>
  inline
  KineticsConditions<StateType,VectorStateType>::~KineticsConditions()
  {
    return;
  }

  template <typename StateType, typename VectorStateType>
  inline
  void KineticsConditions<StateType,VectorStateType>::add_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr)
  {
     _map_pf.insert(std::make_pair(nr, &pf));
  }

  template <typename StateType, typename VectorStateType>
  inline
  void KineticsConditions<StateType,VectorStateType>::change_temperature(const StateType & temperature)
  {
      _temperature = &temperature;
  }

  template <typename StateType, typename VectorStateType>
  inline
  void KineticsConditions<StateType,VectorStateType>::change_particle_flux(const ParticleFlux<VectorStateType> & pf, unsigned int nr)
  {
     antioch_assert(_map_pf.count(nr));
     _map_pf[nr] = &pf;
  }

  template <typename StateType, typename VectorStateType>
  inline
  const StateType & KineticsConditions<StateType,VectorStateType>::T() const
  {
     antioch_assert(_temperature);
     return *_temperature;
  }

  template <typename StateType, typename VectorStateType>
  inline
  const ParticleFlux<VectorStateType> & KineticsConditions<StateType,VectorStateType>::particle_flux(int nr) const
  {
     antioch_assert(_map_pf.count(nr));
     return *(_map_pf.at(nr));
  }


} //end namespace Antioch

#endif
