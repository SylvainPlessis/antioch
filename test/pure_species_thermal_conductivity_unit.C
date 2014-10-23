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
int test_k( const Scalar k, const Scalar k_exact, const Scalar tol, const std::string & words )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << "Error: Mismatch in thermal conductivity " << words << std::endl
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
  const unsigned int n_species = 3;
  species_str_list.reserve(n_species);
//linear
  species_str_list.push_back( "N2" );
//atom
  species_str_list.push_back( "H" );
//other
  species_str_list.push_back( "H2O" );

  const Scalar LJ_depth_N2 = 97.53L;
  const Scalar Z_298_N2 = 4.0L;
  std::vector<Scalar> coeffs_N2(14,0.);
  std::vector<Scalar> temps_N2(3,0.);

  const Scalar LJ_depth_H = 145.0L;
  const Scalar Z_298_H = 0.0L;
  std::vector<Scalar> coeffs_H(14,0.);
  std::vector<Scalar> temps_H(3,0.);

  const Scalar LJ_depth_H2O = 572.400L;
  const Scalar Z_298_H2O = 4.0L;
  std::vector<Scalar> coeffs_H2O(14,0.);
  std::vector<Scalar> temps_H2O(3,0.);

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

  coeffs_H[0]  =  2.50000000E+00;
  coeffs_H[1]  =  7.05332819E-13;
  coeffs_H[2]  = -1.99591964E-15;
  coeffs_H[3]  =  2.30081632E-18;
  coeffs_H[4]  = -9.27732332E-22;
  coeffs_H[5]  =  2.54736599E+04;
  coeffs_H[6]  = -4.46682853E-01;
  coeffs_H[7]  =  2.50000001E+00;
  coeffs_H[8]  = -2.30842973E-11;
  coeffs_H[9]  =  1.61561948E-14;
  coeffs_H[10] = -4.73515235E-18;
  coeffs_H[11] =  4.98197357E-22;
  coeffs_H[12] =  2.54736599E+04;
  coeffs_H[13] = -4.46682914E-01;

  temps_H[0] = 200.L;
  temps_H[1] = 1000.L;
  temps_H[2] = 3500.L;

  coeffs_H2O[0]  =  4.19864056E+00;
  coeffs_H2O[1]  = -2.03643410E-03;
  coeffs_H2O[2]  =  6.52040211E-06;
  coeffs_H2O[3]  = -5.48797062E-09;
  coeffs_H2O[4]  =  1.77197817E-12;
  coeffs_H2O[5]  = -3.02937267E+04;
  coeffs_H2O[6]  = -8.49032208E-01;
  coeffs_H2O[7]  =  3.03399249E+00;
  coeffs_H2O[8]  =  2.17691804E-03;
  coeffs_H2O[9]  = -1.64072518E-07;
  coeffs_H2O[10] = -9.70419870E-11;
  coeffs_H2O[11] =  1.68200992E-14;
  coeffs_H2O[12] = -3.00042971E+04;
  coeffs_H2O[13] =  4.96677010E+00;

  temps_H2O[0] = 200.L;
  temps_H2O[1] = 1000.L;
  temps_H2O[2] = 3500.L;

  Antioch::ChemicalMixture<Scalar> chem_mixture( species_str_list );

  Antioch::NASAThermoMixture<Scalar, Antioch::NASACurveFit<Scalar> > thermo_mix( chem_mixture );
  thermo_mix.add_curve_fit("N2", coeffs_N2, temps_N2);
  thermo_mix.add_curve_fit("H",  coeffs_H,  temps_H);
  thermo_mix.add_curve_fit("H2O",coeffs_H2O,temps_H2O);
  Antioch::NASAEvaluator<Scalar, Antioch::NASACurveFit<Scalar> > thermo( thermo_mix );

  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_N2( thermo, Z_298_N2, LJ_depth_N2);
  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_H( thermo, Z_298_H, LJ_depth_H);
  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_H2O( thermo, Z_298_H2O, LJ_depth_H2O);

//made up values
  const Scalar mu  = 3.14e-3;
  const Scalar dss = 5.23e-5;
  const Scalar rho = 1.4;
  const Scalar T_low = 800.1L;
  const Scalar T_high = 1500.1L;

  // from bc
  const Scalar k_N2_exact_low   = 3.093167156837462471870738790881329983241829808672669033364;
  const Scalar k_H_exact_low    = 97.12578494791666666666666666666666666666666666666666666666;
  const Scalar k_H2O_exact_low  = 4.854337293026230683188193320063326358735232984125875255551;

  const Scalar k_N2_exact_high  = 3.194831520584533016571929766390779999673621653723027368129;
  const Scalar k_H_exact_high   = 97.12578494791666666666666666666666666666666666666666666666;
  const Scalar k_H2O_exact_high = 5.143662825303267444216513226471740441147579013327958436906;
  int return_flag = 0;

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 7e-17)?
                        7e-17:
                        std::numeric_limits<Scalar>::epsilon()*10;
  Scalar tolH = tol*50000.;
  if(tolH < 2e-11)tolH = 2e-11;

  int return_flag_temp = 0;
  return_flag_temp = test_k( k_N2(0,mu,T_low,rho,dss),   k_N2_exact_low,   tol,  "N2 low" )   || return_flag_temp;
  return_flag_temp = test_k( k_N2(0,mu,T_high,rho,dss),  k_N2_exact_high,  tol,  "N2 high" )  || return_flag_temp;
  return_flag_temp = test_k( k_H(1,mu,T_low,rho,dss),    k_H_exact_low,    tolH, "H low" )    || return_flag_temp; //NASA polynomials cv = cp - R
  return_flag_temp = test_k( k_H(1,mu,T_high,rho,dss),   k_H_exact_high,   tolH, "H high" )   || return_flag_temp; //NASA polynomials cv = cp - R
  return_flag_temp = test_k( k_H2O(2,mu,T_low,rho,dss),  k_H2O_exact_low,  tol,  "H2O low" )  || return_flag_temp;
  return_flag_temp = test_k( k_H2O(2,mu,T_high,rho,dss), k_H2O_exact_high, tol,  "H2O high" ) || return_flag_temp;
  if( return_flag_temp != 0 ) return_flag = 1;

  return return_flag;
}

int main()
{
  return (tester<double>() ||
          tester<long double>() ||
          tester<float>());
}
