//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// Antioch - A Gas Dynamics Thermochemistry Library
//
// Copyright (C) 2014 Paul T. Bauman, Benjamin S. Kirk, Sylvain Plessis,
//                    Roy H. Stonger
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

// C++
#include <iostream>
#include <cmath>

// Antioch
#include "antioch/stat_mech_thermo.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_curve_fit.h"
#include "antioch/pure_species_thermal_conductivity.h"
#include "antioch/chemical_mixture.h"
#include "antioch/metaprogramming.h"

template <typename Scalar>
int test_k( const Scalar k, const Scalar k_exact, const Scalar tol )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity" << std::endl
		<< "k       = " << k << std::endl
		<< "k_exact = " << k_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename Scalar>
int tester()
{
  std::vector<std::string> species_str_list;
  const unsigned int n_species = 1;
  species_str_list.reserve(n_species);
  species_str_list.push_back( "N2" );

  const Scalar LJ_depth_N2 = 97.53L;
  const Scalar Z_298 = 4.0L;
  std::vector<Scalar> coeffs_N2(14,0.);
  std::vector<Scalar> temps_N2(3,0.);

  coeffs_N2[0]  = 0.03298677E+02;
  coeffs_N2[1]  = 0.14082404E-02;;
  coeffs_N2[2]  = -0.03963222E-04;
  coeffs_N2[3]  =  0.05641515E-07;
  coeffs_N2[4]  = -0.02444854E-10;
  coeffs_N2[5]  = -0.10208999E+04;
  coeffs_N2[6]  =  0.03950372E+02;
  coeffs_N2[7]  =  0.02926640E+02;
  coeffs_N2[8]  =  0.14879768E-02;
  coeffs_N2[9]  = -0.05684760E-05;
  coeffs_N2[10] =  0.10097038E-09;
  coeffs_N2[11] = -0.06753351E-13;
  coeffs_N2[12] = -0.09227977E+04;
  coeffs_N2[13] =  0.05980528E+02;

  temps_N2[0] = 300.L;
  temps_N2[1] = 1000.L;
  temps_N2[2] = 5000.L;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::NASAThermoMixture<Scalar, Antioch::NASACurveFit<Scalar> > thermo_mix( chem_mixture );
  thermo_mix.add_curve_fit("N2",coeffs_N2,temps_N2);
  Antioch::NASAEvaluator<Scalar, Antioch::NASACurveFit<Scalar> > thermo( thermo_mix );

  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k( thermo, Z_298, LJ_depth_N2);

  const Scalar mu  = 3.14e-3;
  const Scalar dss = 5.23e-5;
  const Scalar rho = 1.4;
  const Scalar T = 1500.1L;

  // from bc
  const Scalar k_N2_exact = 123.296800605309501138839307257794302255577891781560672779553;

  int return_flag = 0;

  const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 2;

  int return_flag_temp = 0;
  return_flag_temp = test_k( k(0,mu,T,rho,dss), k_N2_exact, tol );
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
//         tester<long double>() ||
          tester<float>());
}
