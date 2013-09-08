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
#ifndef ANTIOCH_EQUILIBRIUM_EVALUATOR_H
#define ANTIOCH_EQUILIBRIUM_EVALUATOR_H

#include "antioch_config.h"

#ifdef ANTIOCH_HAVE_EIGEN
#include "antioch/eigen_utils_decl.h"
#include "antioch/eigen_utils.h"

// antioch
#include "antioch/kinetics_evaluator.h"
#include "antioch/cea_thermo.h"
#include "antioch/data_equilibrium.h"

//C++
#include <vector>

namespace Antioch
{

  template<typename CoeffType = double,typename StateType=CoeffType>
  class EquilibriumEvaluator
  {
  public:
    EquilibriumEvaluator(DataEquilibrium<CoeffType> &data_eq,
                         KineticsEvaluator<CoeffType,StateType>&kin,
                         const StateType tol = std::numeric_limits<StateType>::epsilon() * 100);

    ~EquilibriumEvaluator();

    void first_guess_mass_fraction(const std::vector<CoeffType> &first) {eq_mass_fraction = first;}
    void first_guess_molar_fraction(const std::vector<CoeffType> &first);


    const std::vector<CoeffType> mass_fraction_equilibrium()   const {return eq_mass_fraction;}
    const std::vector<CoeffType> molar_densities_equilibrium() const {return eq_molar_densities;}
    const CoeffType Peq() const;

    void equilibrium();

    const bool success()              const {return !over_threshold;}
    const StateType residual()        const {return thres;}
    const StateType conv_threshold()  const {return threshold;}
    const unsigned int max_loop_tol() const {return max_loop;}

    void provide_first_guess(const std::vector<CoeffType> &first_mass_fraction)
    {eq_mass_fraction = first_mass_fraction;}


  private:

    void first_guess_iso();
    void init_from_mass();
    void calculate_function_and_jacobian(std::vector<StateType> &F,
                                         std::vector<std::vector<StateType> > &jacob);

    std::vector<CoeffType> eq_molar_densities;
    std::vector<CoeffType> eq_mass_fraction;
    std::vector<CoeffType> eq_mass;
    CoeffType mass_tot;
    CoeffType mass_tot_ini;
    KineticsEvaluator<CoeffType,StateType> &kinetics_eval;
    DataEquilibrium<CoeffType> & data_storage_and_constrain;

    void solve_decomp(const std::vector<std::vector<CoeffType> > &jacob,
                      const std::vector<CoeffType> &target,
                      Eigen::Matrix<CoeffType,Eigen::Dynamic,1> &x);

    Eigen::Matrix<CoeffType,Eigen::Dynamic,Eigen::Dynamic> A;
    Eigen::Matrix<CoeffType,Eigen::Dynamic,1> b;


    //useless but needed
    std::vector<CoeffType> dGibbs_RT_dT;
    std::vector<CoeffType> dmass_dT;
    std::vector<std::vector<StateType> > xistoi;
    std::vector<CoeffType> eq_xi;

    //local variables
    const StateType threshold;
    StateType thres;
    const unsigned int max_loop;
    bool over_threshold;

    //! 
    EquilibriumEvaluator();
  };

  template<typename CoeffType, typename StateType>
  inline
  EquilibriumEvaluator<CoeffType,StateType>::~EquilibriumEvaluator()
  {
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  EquilibriumEvaluator<CoeffType,StateType>::EquilibriumEvaluator( DataEquilibrium<CoeffType>& data_eq,
                                                                   KineticsEvaluator<CoeffType,StateType>& kin,
                                                                   const StateType tol )
    : mass_tot(0.),
      kinetics_eval(kin),
      data_storage_and_constrain(data_eq),
      A(data_eq.reaction_set().n_reactions(),data_eq.reaction_set().n_reactions()),
      b(data_eq.reaction_set().n_reactions()),
      threshold(tol),max_loop(100),
      over_threshold(true)
  {
    eq_mass_fraction.resize(data_storage_and_constrain.reaction_set().n_species(),0.);
    eq_mass.resize(data_storage_and_constrain.reaction_set().n_species(),0.);
    eq_molar_densities.resize(data_storage_and_constrain.reaction_set().n_species(),0.);
    eq_xi.resize(data_storage_and_constrain.reaction_set().n_reactions());
    xistoi.resize(data_storage_and_constrain.reaction_set().n_reactions());
    for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions(); rxn++)
    {
      xistoi[rxn].resize(data_storage_and_constrain.reaction_set().n_species());
    }
    first_guess_iso();
    return;
  }

  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::equilibrium()
  {

    init_from_mass();

    std::vector<StateType> F;
    std::vector<std::vector<StateType> > jacob;
    F.resize(data_storage_and_constrain.reaction_set().n_reactions());
    jacob.resize(data_storage_and_constrain.reaction_set().n_reactions());
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_reactions(); i++)
    {
      jacob[i].resize(data_storage_and_constrain.reaction_set().n_reactions());
    }

    over_threshold = true;
    unsigned int nloop(0);
    while(over_threshold)
    {
      nloop++;
      calculate_function_and_jacobian(F,jacob);
      Antioch::set_zero(thres);
      antioch_assert_equal_to(F.size(),data_storage_and_constrain.reaction_set().n_reactions());
      antioch_assert_equal_to(jacob.size(),data_storage_and_constrain.reaction_set().n_reactions());
      for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_reactions(); i++)
      {
        thres += (F[i] > 0.)?F[i]:-F[i];
        antioch_assert_equal_to((int)jacob[i].size(),b.innerSize());
        b(i) = - F[i];
        for(unsigned int j = 0; j < jacob[i].size(); j++)
        {
          A(i,j) = jacob[i][j];
        }
      }  
      over_threshold = thres > threshold;
std::cout << "thresh " << thres << " v/s " << threshold << std::endl;
(over_threshold)?std::cout << "pre: loop again" << std::endl:std::cout << "pre: out" << std::endl;
      if(!over_threshold)break;

std::cout << A << " and\n" << b << std::endl;
      Eigen::PartialPivLU<Eigen::Matrix<StateType,Eigen::Dynamic,Eigen::Dynamic> > mypartialPivLu(A);
      Eigen::Matrix<StateType,Eigen::Dynamic,1> x(data_storage_and_constrain.reaction_set().n_reactions());
      x = mypartialPivLu.solve(b);
std::cout << "gives\n" << x << std::endl;

      Antioch::set_zero(mass_tot);
      Antioch::set_zero(thres);
std::cout << "before update" << std::endl;
      for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_reactions(); i++)
      {
        std::cout << eq_xi[i] << std::endl;
      }
      for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_reactions(); i++)
      {
        thres += (x(i) > 0.)?x(i):-x(i);
        eq_xi[i] += x(i);
        if(eq_xi[i] > data_storage_and_constrain.max_xi(i))eq_xi[i] = data_storage_and_constrain.max_xi(i);
        if(eq_xi[i] < data_storage_and_constrain.min_xi(i))eq_xi[i] = data_storage_and_constrain.min_xi(i);
      }
std::cout << "after update" << std::endl;
      for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_reactions(); i++)
      {
        std::cout << eq_xi[i] << std::endl;
      }
//        for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
//        {
            //std::cout << data_storage_and_constrain.reaction_set().chemical_mixture().chemical_species()[i]->species() << " (" << minus_function(i) << "): " 
            //          << eq_mass[i] << " ** ";
/*            eq_mass[i] += x(i);
            if(eq_mass[i] < 0.)eq_mass[i] = -eq_mass[i];
            eq_molar_densities[i] = eq_mass[i]/data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
*/            //std::cout << eq_mass[i] <<  std::endl;
//            mass_tot += eq_mass[i];
//          }

std::cout << "mass tot & ini " << mass_tot << "  " << mass_tot_ini << std::endl;
std::cout << "thresh " << thres << " v/s " << threshold << std::endl;
std::cout << "nloop " << nloop << "/" << max_loop << std::endl;

      if(thres != thres)break; //!TODO quick nan fix in equilibrium_evaluator.h, need better treatment
      over_threshold = thres > threshold;
      (over_threshold)?std::cout << "loop again" << std::endl:std::cout << "out" << std::endl;
      if(nloop > max_loop)break;
    }
    Antioch::set_zero(mass_tot);
std::cout << "out of loop conv" << std::endl;
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
    {
      for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions(); rxn++)
      {
        eq_molar_densities[i] += eq_xi[rxn]*data_storage_and_constrain.stoi(rxn,i);
      }
      eq_mass[i] = eq_molar_densities[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      mass_tot += eq_mass[i];
      std::cout << eq_molar_densities[i] << std::endl;
    }
std::cout << "mass tot = " << mass_tot << std::endl;

    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
    {
      eq_mass_fraction[i] = eq_mass[i]/mass_tot;
    }

  }


  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::calculate_function_and_jacobian(std::vector<StateType>& F, std::vector<std::vector<StateType> > &jacob)
  {
     Antioch::set_zero(F);//first fill with pure kinetics
     Antioch::set_zero(jacob);//first fill with pure kinetics
    ///kinetics part only first
    //updating
/*    std::cout << std::endl;
    for(unsigned int i = 0; i < eq_mass_fraction.size(); i++)
      {
        std::cout << eq_molar_densities[i] << std::endl;
      }
    std::cout << std::endl;
*/


// energy matrix, d2G/dx_idx_j * Delta_x = dG/dx
    std::vector<StateType> tmp;
    tmp.resize(data_storage_and_constrain.reaction_set().n_species(),0.);
    for(unsigned int nsp = 0; nsp < data_storage_and_constrain.reaction_set().n_species(); nsp++)
    {
      for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions(); rxn++)
      {
         tmp[nsp] += data_storage_and_constrain.stoi(rxn,nsp) * eq_xi[rxn];
      }
      tmp[nsp] += eq_molar_densities[nsp];
std::cout << "denominateur " << tmp[nsp] << std::endl;
    }
    for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions()-1; rxn++)
    {
      F[rxn] = 0.;
      for(unsigned int nsp = 0; nsp < data_storage_and_constrain.reaction_set().n_species(); nsp++)
      {
          if(tmp[nsp] == 0.)continue;
          F[rxn] += data_storage_and_constrain.stoi(rxn,nsp)/tmp[nsp];
      }
      F[rxn] *= Constants::R_universal<CoeffType>() * data_storage_and_constrain.T();
      for(unsigned int rxm = 0; rxm < data_storage_and_constrain.reaction_set().n_reactions(); rxm++)
      {
        jacob[rxn][rxm] = 0.;
        for(unsigned int nsp = 0; nsp < data_storage_and_constrain.reaction_set().n_species(); nsp++)
        {
            if(tmp[nsp] == 0.)continue;
            jacob[rxn][rxm] += data_storage_and_constrain.stoi(rxn,nsp)*data_storage_and_constrain.stoi(rxm,nsp)/
                                (tmp[nsp] * tmp[nsp]);
        }
      }
    }
    for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions(); rxn++)
    {
        jacob.back()[rxn] = 1.;
    }
    F.back() = 0;
    
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        eq_mass_fraction[i] = eq_mass[i]/mass_tot;
      }

    return;

//kinetics matrix: d_omega_dot/d_rho * Delta_rho = omega_dot
    //setting the system
    std::vector<StateType> h_RT_minus_s_R(data_storage_and_constrain.reaction_set().n_species());

    Antioch::CEAThermodynamics<StateType> thermo(data_storage_and_constrain.reaction_set().chemical_mixture() );
    typedef typename Antioch::CEAThermodynamics<StateType>::template Cache<StateType> Cache;
    thermo.h_RT_minus_s_R(Cache(data_storage_and_constrain.T()),h_RT_minus_s_R);

    //omega_dot & jacobian calculations
    /*std::cout << "T: " << data_storage_and_constrain.T()<< std::endl;
      std::cout << "local P: " << cur_P << std::endl;
      std::cout << "rho: " << rho << std::endl;
      std::cout << "R_mix: " << R_mix << std::endl;
      std::cout << "n_species: " << data_storage_and_constrain.reaction_set().n_species() << std::endl;
      std::cout << "n_reaction: " << data_storage_and_constrain.reaction_set().n_reactions() << std::endl;
    */   kinetics_eval.compute_mass_sources_and_derivs(data_storage_and_constrain.T(),
                                                       eq_molar_densities,
                                                       h_RT_minus_s_R,
                                                       dGibbs_RT_dT,
                                                       F,
                                                       dmass_dT,
                                                       jacob);
/*        std::vector<StateType> scale(data_storage_and_constrain.reaction_set().n_species());
        kinetics_eval.compute_mass_sources(1000.,
                                           eq_molar_densities,
                                           h_RT_minus_s_R,
                                           scale);
*/
    /*  std::cout << "F size " << F.size() << std::endl << "\t";
        for(unsigned int i = 0; i < F.size(); i++)std::cout << F[i] << " ";
        std::cout << std::endl;
    */
    
 /*   for(unsigned int i = 0; i < jacob.size()-1; i++)
    {
      F[i] /= scale[i];
      for(unsigned int j= 0; j < jacob[i].size(); j++)
      {
        jacob[i][j] /= scale[i];
      }
    }
    
*/
    data_storage_and_constrain.fill_constrain(eq_molar_densities,F,jacob); //constrain added
  }

  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::first_guess_iso()
  {
    CoeffType iso = 1./(double)(data_storage_and_constrain.reaction_set().n_species());
    std::fill(eq_mass_fraction.begin(),eq_mass_fraction.end(),iso);
  }


  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::init_from_mass()
  {
    eq_molar_densities = Antioch::zero_clone(eq_mass_fraction);
    eq_mass = Antioch::zero_clone(eq_mass_fraction);
    dGibbs_RT_dT = Antioch::zero_clone(eq_mass_fraction);
    dmass_dT = Antioch::zero_clone(eq_mass_fraction);
    mass_tot_ini = data_storage_and_constrain.reaction_set().chemical_mixture().M(eq_mass_fraction) *
      data_storage_and_constrain.P() /
      (Constants::R_universal<CoeffType>() * data_storage_and_constrain.T());
    mass_tot = mass_tot_ini;
    data_storage_and_constrain.set_mass(mass_tot);

    for(unsigned int i = 0; i < eq_mass_fraction.size(); i++)
      {
        eq_mass[i] = mass_tot * eq_mass_fraction[i];
        eq_molar_densities[i] = eq_mass[i] / data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      }
    Antioch::set_zero(eq_xi);

    data_storage_and_constrain.set_constrain("",eq_molar_densities);

    return;

/// first approx prod/loss*Cs
    std::vector<std::vector<StateType> > lossM,prodM,netM;
    //setting the system
    std::vector<StateType> h_RT_minus_s_R(data_storage_and_constrain.reaction_set().n_species());

    Antioch::CEAThermodynamics<StateType> thermo(data_storage_and_constrain.reaction_set().chemical_mixture() );
    typedef typename Antioch::CEAThermodynamics<StateType>::template Cache<StateType> Cache;
    thermo.h_RT_minus_s_R(Cache(data_storage_and_constrain.T()),h_RT_minus_s_R);
    data_storage_and_constrain.reaction_set().print_chemical_scheme(std::cout,data_storage_and_constrain.T(),
                                                    eq_molar_densities,h_RT_minus_s_R,
                                                    lossM,prodM,netM);
    Antioch::set_zero(mass_tot);
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        StateType prod;
        StateType loss;
        Antioch::set_zero(prod);
        Antioch::set_zero(loss);
        for(unsigned int rxn = 0; rxn < data_storage_and_constrain.reaction_set().n_reactions(); rxn++)
        {
           prod += prodM[i][rxn];
           loss += lossM[i][rxn];
        }
        eq_molar_densities[i] *= prod/loss;
        eq_mass[i] = eq_molar_densities[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
        mass_tot += eq_mass[i];
      }
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
      {
        eq_mass_fraction[i] = eq_mass[i] / mass_tot;
      }

  }

  template<typename CoeffType, typename StateType>
  inline
  const CoeffType EquilibriumEvaluator<CoeffType, StateType>::Peq() const
  {
    return mass_tot / (data_storage_and_constrain.reaction_set().chemical_mixture().M(eq_mass_fraction)) *
      (Constants::R_universal<CoeffType>()*data_storage_and_constrain.T());
  }



  template<typename CoeffType, typename StateType>
  inline
  void EquilibriumEvaluator<CoeffType, StateType>::first_guess_molar_fraction(const std::vector<CoeffType> &first)
  {
    antioch_assert_equal_to(first.size(),data_storage_and_constrain.reaction_set().n_species() );
    eq_molar_densities = Antioch::zero_clone(eq_mass_fraction);
    eq_mass = Antioch::zero_clone(eq_mass_fraction);
    Antioch::set_zero(mass_tot);
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
    {
      eq_molar_densities[i] = first[i] * data_storage_and_constrain.P()/(Constants::R_universal<CoeffType>()*data_storage_and_constrain.T());
      eq_mass[i] = eq_molar_densities[i] * data_storage_and_constrain.reaction_set().chemical_mixture().M(i);
      mass_tot += eq_mass[i];
    }
    for(unsigned int i = 0; i < data_storage_and_constrain.reaction_set().n_species(); i++)
    {
      eq_mass_fraction[i] = eq_mass[i]/mass_tot;
    }
    return;
  }

} // end namespace Antioch

#endif // ANTIOCH_EQUILIBRIUM_EVALUATOR_H
#endif //ANTIOCH_HAVE_EIGEN
