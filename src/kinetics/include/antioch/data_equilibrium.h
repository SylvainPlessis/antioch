//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-

#ifndef ANTIOCH_DATA_EQUILIBRIUM_H
#define ANTIOCH_DATA_EQUILIBRIUM_H

// antioch
#include "antioch/chemical_mixture.h"

//C++
#include <vector>
#include <string>

namespace Antioch
{

//forward declaration
  template<typename CoeffType>
  class ReactionSet;

  template<typename CoeffType = double>
  class DataEquilibrium
  {
     public:
          DataEquilibrium(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set, const std::string &key = std::string());
          virtual ~DataEquilibrium();

          void fill_constrain(const std::vector<CoeffType> &molar_densities,
                                        std::vector<CoeffType> &F,std::vector<std::vector<CoeffType> > &jacob);

          void set_constrain(const std::string &key,const std::vector<CoeffType> &molar_densities);

          const CoeffType T() const {return _T;}
          const CoeffType P() const {return _P;}
          const CoeffType local_pressure(const std::vector<CoeffType> &molar_densities) const; 
          const ReactionSet<CoeffType> &reaction_set() const {return _reac_set;}

          const unsigned int n_constrain() const {return _n_constrain;}
          void update_constrain(unsigned int icstr, const CoeffType &delta_cstr);

          void set_mass(const CoeffType &m) {_fixed_mass = m;}
          void set_T(const CoeffType &t)    {_T = t;}
          void set_P(const CoeffType &p)    {_P = p;}
          int stoi(unsigned int nreac,unsigned int nspec) const {return reac_species_stoi[nreac][nspec];}
          CoeffType max_xi(unsigned int nreac) const {return xi_max[nreac];}
          CoeffType min_xi(unsigned int nreac) const {return xi_min[nreac];}


     protected:
     const CoeffType &_T;
     const CoeffType &_P;
     const ReactionSet<CoeffType> &_reac_set;
     std::vector<std::vector<int> > reac_species_stoi;
     std::vector<CoeffType> xi_max;
     std::vector<CoeffType> xi_min;

     private:
//I don't want it to be used, change it if you want
     DataEquilibrium(){return;}
     //!\brief fill the matrix reac / species
     void fill_stoi_matrix(const std::vector<CoeffType> &molar_densities);


     CoeffType _fixed_pressure;
     CoeffType _fixed_mass;
     unsigned int m_constrain;
     unsigned int a_constrain;
     unsigned int p_constrain;
     unsigned int _n_constrain;

  };


  template<typename CoeffType>
  inline
  DataEquilibrium<CoeffType>::DataEquilibrium(const CoeffType &T_mix, const CoeffType &P_mix, 
                    const ReactionSet<CoeffType> &reac_set, const std::string &key):
     _T(T_mix),_P(P_mix),_reac_set(reac_set),m_constrain(1),a_constrain(0),p_constrain(0),_n_constrain(1)
  {
     _fixed_pressure = this->_P/(Constants::R_universal<CoeffType>() * this->_T);
     return;
  }

  template<typename CoeffType>
  inline
  DataEquilibrium<CoeffType>::~DataEquilibrium()
  {
    return;
  }


  template<typename CoeffType>
  inline
  void DataEquilibrium<CoeffType>::set_constrain(const std::string &key,const std::vector<CoeffType> &molar_densities)
  {
     if(key.find("m") != std::string::npos)m_constrain = 1;
     if(key.find("p") != std::string::npos){p_constrain = 1;m_constrain = 0;}
     fill_stoi_matrix(molar_densities);
     return;
  }

  template<typename CoeffType>
  inline
  void DataEquilibrium<CoeffType>::fill_stoi_matrix(const std::vector<CoeffType> &molar_densities)
  {
     reac_species_stoi.resize(this->_reac_set.n_reactions());
     xi_max.resize(this->_reac_set.n_reactions());
     xi_min.resize(this->_reac_set.n_reactions());
     for(unsigned int rxn = 0; rxn < this->_reac_set.n_reactions(); rxn++)
     {
        reac_species_stoi[rxn].resize(this->_reac_set.n_species());
//n_a A + n_b B <=> n_c C + n_d D, min([A]/n_a,[B]/n_b) limits the forward direction (extent of reaction > 0)
        xi_max[rxn] = 1e20;
        for (unsigned int r=0; r < this->_reac_set.reaction(rxn).n_reactants(); r++)
        {
          reac_species_stoi[rxn][this->_reac_set.reaction(rxn).reactant_id(r)] = -this->_reac_set.reaction(rxn).reactant_stoichiometric_coefficient(r);
          CoeffType loc_xi = molar_densities[this->_reac_set.reaction(rxn).reactant_id(r)] / this->_reac_set.reaction(rxn).reactant_stoichiometric_coefficient(r);
          if(loc_xi < xi_max[rxn])xi_max[rxn] = loc_xi;
        }
//n_a A + n_b B <=> n_c C + n_d D, max(-[C]/n_c,-[D]/n_d) limits the backward direction (extent of reaction < 0)
        xi_min[rxn] = -1e20;
        for (unsigned int p=0; p < this->_reac_set.reaction(rxn).n_products(); p++)
        {
          reac_species_stoi[rxn][this->_reac_set.reaction(rxn).product_id(p)] = this->_reac_set.reaction(rxn).product_stoichiometric_coefficient(p);
          CoeffType loc_xi = molar_densities[this->_reac_set.reaction(rxn).product_id(p)] / -this->_reac_set.reaction(rxn).product_stoichiometric_coefficient(p);
          if(loc_xi > xi_min[rxn])xi_min[rxn] = loc_xi;
        }
     }
  }

  template<typename CoeffType>
  inline
  const CoeffType DataEquilibrium<CoeffType>::local_pressure(const std::vector<CoeffType> &molar_densities) const
  {
     antioch_assert_equal_to(molar_densities.size(),this->_reac_set.n_species());
     CoeffType sum_mol;
     Antioch::set_zero(sum_mol);
     for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
     {
        sum_mol += molar_densities[i];
     }
     return (sum_mol * this->_T * Constants::R_universal<CoeffType>());
  }

  template<typename CoeffType>
  inline
  void DataEquilibrium<CoeffType>::fill_constrain(const std::vector<CoeffType> &molar_densities,std::vector<CoeffType> &F, std::vector<std::vector<CoeffType> > &jacob)
  {
     if(m_constrain)
     {
       antioch_assert_equal_to(F.size(),jacob.size());
       for(unsigned int i = 0; i < F.size(); i++)
       {
         antioch_assert_equal_to(jacob[i].size(),F.size());
       }

       CoeffType mass_sum  = Antioch::zero_clone(this->_P);
       for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
       {
          mass_sum += molar_densities[i] * this->_reac_set.chemical_mixture().M(i);
          jacob.back()[i] = 1.;
       }
       F.back() = mass_sum - _fixed_mass;
     }
     if(p_constrain)
     {
       antioch_assert_equal_to(F.size(),jacob.size());
       for(unsigned int i = 0; i < F.size(); i++)
       {
         antioch_assert_equal_to(jacob[i].size(),F.size());
       }

       CoeffType molar_sum  = Antioch::zero_clone(this->_P);
       for(unsigned int i = 0; i < this->_reac_set.n_species(); i++)
       {
          molar_sum += molar_densities[i];
          jacob.back()[i] = 1./this->_reac_set.chemical_mixture().M(i);
       }
       F.back() = (molar_sum - _fixed_pressure);
     }
  }

} // end namespace Antioch

#endif // ANTIOCH_DATA_EQUILIBRIUM_H
