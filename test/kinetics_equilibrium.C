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
//

// C++
#include <limits>
#include <string>
#include <vector>

// Antioch
#include "antioch/eigen_utils_decl.h"
#include "antioch/vector_utils.h"

#include "antioch/antioch_asserts.h"
#include "antioch/equilibrium_evaluator.h"
#include "antioch/data_equilibrium.h"
#include "antioch/chemical_species.h"
#include "antioch/chemical_mixture.h"
#include "antioch/reaction_set.h"
#include "antioch/read_reaction_set_data_xml.h"
#include "antioch/cea_thermo.h"
#include "antioch/kinetics_evaluator.h"

#ifdef ANTIOCH_HAVE_EIGEN

template <typename Scalar>
int tester(const std::string& input_name)
{
  using std::abs;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 5;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );
  species_str_list.push_back( "O2" );
  species_str_list.push_back( "NO" );
  species_str_list.push_back( "N" );
  species_str_list.push_back( "O" );

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );
  Antioch::ReactionSet<Scalar> reaction_set( chem_mixture );

  Antioch::read_reaction_set_data_xml<Scalar>( input_name, true, reaction_set );

  const Scalar T = 2500.0;
  const Scalar P = 1.0e5;

  Antioch::KineticsEvaluator<Scalar> kinetics( reaction_set, 0 );


  Antioch::DataEquilibrium<Scalar> equil(T, P,reaction_set);

  std::vector<Scalar> first;
  first.push_back(0.7);
  first.push_back(0.2);
  first.push_back(0.03);
  first.push_back(0.03);
  first.push_back(0.04);

  Antioch::EquilibriumEvaluator<Scalar> eq_solver(equil,kinetics);

  eq_solver.first_guess_molar_fraction(first);
  eq_solver.equilibrium();

  // solution testing => omega_dot = 0
  std::vector<Scalar> mass_fraction_eq = eq_solver.mass_fraction_equilibrium();
  std::vector<Scalar> molar_densities_eq = eq_solver.molar_densities_equilibrium();

  const Scalar R_mix = chem_mixture.R(mass_fraction_eq); // get R_tot in J.kg-1.K-1
  const Scalar rho = P/(R_mix*T); // kg.m-3

  std::vector<Scalar> molar_densities(n_species,0.);
  chem_mixture.molar_densities(rho,mass_fraction_eq,molar_densities);

  std::vector<Scalar> h_RT_minus_s_R(n_species);
  std::vector<Scalar> dh_RT_minus_s_R_dT(n_species);

  Antioch::CEAThermodynamics<Scalar> thermo( chem_mixture );
  typedef typename Antioch::CEAThermodynamics<Scalar>::template Cache<Scalar> Cache;
  thermo.h_RT_minus_s_R(Cache(T),h_RT_minus_s_R);
  thermo.dh_RT_minus_s_R_dT(Cache(T),dh_RT_minus_s_R_dT);

  std::vector<Scalar> omega_dot(n_species);
  
  kinetics.compute_mass_sources( T, molar_densities, h_RT_minus_s_R, omega_dot );

  int return_flag = 0;

  if(eq_solver.success())
  {

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 100;

//first test consistency
  Scalar cons;
  Antioch::set_zero(cons);
  for(unsigned int isp = 0; isp < n_species; isp++)
  {
      cons += (molar_densities[isp] - molar_densities_eq[isp] > 0.)?
                molar_densities[isp] - molar_densities_eq[isp]:
                molar_densities_eq[isp] - molar_densities[isp];
  }
   if(cons > tol)
   {
      return_flag = 1;
       std::cout << "tolerance is " << tol << " and diff in concentration is " << cons << std::endl;
      for(unsigned int isp = 0; isp < n_species; isp++)
      {
         std::cout << chem_mixture.chemical_species()[isp]->species() 
                   << " eq_solver = " << molar_densities_eq[isp] 
                   << " , recomputed = " << molar_densities[isp] << std::endl;
      }
   }


  Scalar sum_dot(0.L);
  for( unsigned int s = 0; s < n_species; s++)
    {
      sum_dot += (omega_dot[s] > 0.)?omega_dot[s]:-omega_dot[s];
    }

   if(sum_dot > tol)
   {
      return_flag = 1;
       std::cout << "tolerance is " << tol << " and it is " << sum_dot << std::endl;
   }

   }//else
   {
     return_flag = 1;
     std::cout << "equilibrium failed (" << eq_solver.max_loop_tol() << "):\n\tsum " 
               << eq_solver.residual() << ", threshold: " << eq_solver.conv_threshold() << std::endl;
   }

  if(return_flag == 1)
  {
    for( unsigned int s = 0; s < n_species; s++)
      {
        std::cout << std::scientific << std::setprecision(16)
  		  << "omega_dot(" << chem_mixture.chemical_species()[s]->species() << ") = "
		  << omega_dot[s] << std::endl;
      }
   }


  return return_flag;
}
#endif

int main(int argc, char* argv[])
{
#ifdef ANTIOCH_HAVE_EIGEN
  // Check command line count.
  if( argc < 2 )
    {
      // TODO: Need more consistent error handling.
      std::cerr << "Error: Must specify reaction set XML input file." << std::endl;
      antioch_error();
    }

  int fl = (tester<float>(std::string(argv[1])));// ||
      fl = tester<double>(std::string(argv[1])); /*||
          tester<long double>(std::string(argv[1])) || */
        //  );
  return fl;

#else
 return 77;
#endif
}
