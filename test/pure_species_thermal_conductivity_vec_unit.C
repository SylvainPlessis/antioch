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
//--------------------------------------------------------------------------

// valarray has to be declared before Antioch or gcc can't find the
// right versions of exp() and pow() to use??

#include "antioch_config.h"

#include <valarray>

#ifdef ANTIOCH_HAVE_EIGEN
#include "Eigen/Dense"
#endif

#ifdef ANTIOCH_HAVE_METAPHYSICL
#include "metaphysicl/numberarray.h"
#endif

#ifdef ANTIOCH_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

#include "antioch/eigen_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/stat_mech_thermo.h"
#include "antioch/nasa_evaluator.h"
#include "antioch/nasa_mixture.h"
#include "antioch/nasa_curve_fit.h"
#include "antioch/chemical_mixture.h"
#include "antioch/pure_species_thermal_conductivity.h"

#include "antioch/eigen_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/vexcl_utils.h"

#ifdef ANTIOCH_HAVE_GRVY
#include "grvy.h"

GRVY::GRVY_Timer_Class gt;
#endif

#include <cmath>
#include <limits>

template <typename Scalar, typename Element>
int test_k( const Element & k, const Element & k_exact, const Scalar & tol, const std::string & words )
{
  using std::abs;

  int return_flag = 0;

  const Scalar rel_error = abs( (k - k_exact)/k_exact);

  if( rel_error  > tol )
    {
      std::cerr << std::setprecision(20) << std::scientific
                << "Error: Mismatch in thermal conductivity " << words << std::endl
		<< "k       = " << k << std::endl
		<< "k_exact = " << k_exact << std::endl
		<< "rel_error = " << rel_error << std::endl
		<< "tol = " << tol << std::endl;
      return_flag = 1;
    }

  return return_flag;
}

template <typename PairScalars>
int vectester(const PairScalars& example, const std::string& testname)
{

  typedef typename Antioch::value_type<PairScalars>::type Scalar;

  std::vector<std::string> species_str_list;
  const unsigned int n_species = 1;
  species_str_list.reserve(n_species);
//linear
  species_str_list.push_back( "N2" );
//atom
  species_str_list.push_back( "H" );
//other
  species_str_list.push_back( "H2O" );

  const Scalar LJ_depth_N2 = 97.53L; // in K
  const Scalar Z_298_N2 = 4.0L;         // dimensionless
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

 // Antioch::StatMechThermodynamics<Scalar> thermo( chem_mixture );
  Antioch::NASAThermoMixture<Scalar, Antioch::NASACurveFit<Scalar> > thermo_mix( chem_mixture );
  thermo_mix.add_curve_fit("N2", coeffs_N2, temps_N2);
  thermo_mix.add_curve_fit("H",  coeffs_H,  temps_H);
  thermo_mix.add_curve_fit("H2O",coeffs_H2O,temps_H2O);
  Antioch::NASAEvaluator<Scalar, Antioch::NASACurveFit<Scalar> > thermo( thermo_mix );

  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_N2( thermo, Z_298_N2, LJ_depth_N2);
  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_H( thermo, Z_298_H, LJ_depth_H);
  Antioch::PureSpeciesThermalConductivity<Antioch::NASAEvaluator<Scalar,Antioch::NASACurveFit<Scalar> >, Scalar > k_H2O( thermo, Z_298_H2O, LJ_depth_H2O);

  // Construct from example to avoid resizing issues
  PairScalars mu  = example;
  PairScalars dss = example;
  PairScalars rho = example;
  PairScalars T = example;
  PairScalars k_exact_N2 = example;
  PairScalars k_exact_H = example;
  PairScalars k_exact_H2O = example;
  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {
      T[2*tuple]     = 800.1L;
      T[2*tuple+1]   = 1500.1L;
      for(unsigned int i=0; i < 2; i++)
      {
        mu[2*tuple+i]    = 3.14e-3L;
        dss[2*tuple+i]   = 5.23e-5L;
        rho[2*tuple+i]   = 1.4L;
      }

      k_exact_N2[2*tuple]    = 3.093167156837462471870738790881329983241829808672669033364;
      k_exact_H[2*tuple]     = 97.12578494791666666666666666666666666666666666666666666666;
      k_exact_H2O[2*tuple]   = 4.854337293026230683188193320063326358735232984125875255551;
      k_exact_N2[2*tuple+1]  = 3.194831520584533016571929766390779999673621653723027368129;
      k_exact_H[2*tuple+1]   = 97.12578494791666666666666666666666666666666666666666666666;
      k_exact_H2O[2*tuple+1] = 5.143662825303267444216513226471740441147579013327958436906;
    }
  

  int return_flag = 0;

#ifdef ANTIOCH_HAVE_GRVY
  gt.BeginTimer(testname);
#endif

  const PairScalars k_ps_N2 = k_N2(0,mu,T,rho,dss);
  const PairScalars k_ps_H = k_H(1,mu,T,rho,dss);
  const PairScalars k_ps_H2O = k_H2O(2,mu,T,rho,dss);

#ifdef ANTIOCH_HAVE_GRVY
  gt.EndTimer(testname);
#endif

  const Scalar tol = (std::numeric_limits<Scalar>::epsilon()*10 < 7e-17)?
                        7e-17:
                        std::numeric_limits<Scalar>::epsilon()*10;
  Scalar tolH = tol*50000.;
  if(tolH < 2e-11)tolH = 2e-11;

  for (unsigned int tuple=0; tuple != ANTIOCH_N_TUPLES; ++tuple)
    {

      return_flag = test_k( k_ps_N2[2*tuple],    k_exact_N2[2*tuple],    tol,  "N2 low" )   || return_flag;
      return_flag = test_k( k_ps_H[2*tuple],     k_exact_H[2*tuple],     tolH, "H low" )    || return_flag;
      return_flag = test_k( k_ps_H2O[2*tuple],   k_exact_H2O[2*tuple],   tol,  "H2O low" )  || return_flag;
      return_flag = test_k( k_ps_N2[2*tuple+1],  k_exact_N2[2*tuple+1],  tol,  "N2 high" )  || return_flag;
      return_flag = test_k( k_ps_H[2*tuple+1],   k_exact_H[2*tuple+1],   tolH, "H high" )   || return_flag;
      return_flag = test_k( k_ps_H2O[2*tuple+1], k_exact_H2O[2*tuple+1], tol,  "H2O high" ) || return_flag;
    }

  return return_flag;
}


int main()
{
  int returnval = 0;

  returnval = returnval ||
    vectester (std::valarray<float>(2*ANTIOCH_N_TUPLES), "valarray<float>");
  returnval = returnval ||
    vectester (std::valarray<double>(2*ANTIOCH_N_TUPLES), "valarray<double>");
  returnval = returnval ||
    vectester (std::valarray<long double>(2*ANTIOCH_N_TUPLES), "valarray<ld>");
#ifdef ANTIOCH_HAVE_EIGEN
  returnval = returnval ||
    vectester (Eigen::Array<float, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXf");
  returnval = returnval ||
    vectester (Eigen::Array<double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXd");
  returnval = returnval ||
    vectester (Eigen::Array<long double, 2*ANTIOCH_N_TUPLES, 1>(), "Eigen::ArrayXld");
#endif
#ifdef ANTIOCH_HAVE_METAPHYSICL
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, float> (0), "NumberArray<float>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, double> (0), "NumberArray<double>");
  returnval = returnval ||
    vectester (MetaPhysicL::NumberArray<2*ANTIOCH_N_TUPLES, long double> (0), "NumberArray<ld>");
#endif
#ifdef ANTIOCH_HAVE_VEXCL
  vex::Context ctx_f (vex::Filter::All);
  if (!ctx_f.empty())
    returnval = returnval ||
      vectester (vex::vector<float> (ctx_f, 2*ANTIOCH_N_TUPLES), "vex::vector<float>");

  vex::Context ctx_d (vex::Filter::DoublePrecision);
  if (!ctx_d.empty())
    returnval = returnval ||
      vectester (vex::vector<double> (ctx_d, 2*ANTIOCH_N_TUPLES), "vex::vector<double>");
#endif

#ifdef ANTIOCH_HAVE_GRVY
  gt.Finalize();
  gt.Summarize();
#endif

  return returnval;
}
