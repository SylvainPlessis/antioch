//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

/*! \file generic_program.C
 * \brief Example program fully commented
 *
 * This example shows how to set up kinetics
 * and transport using Antioch.
 * The models we have to define here
 * concerns the transport:
 *  - viscosity:
 *     * [Sutherland]{@ref SutherlandViscosity},
 *     * [Blottner]{@ref BlottnerViscosity},
 *     * [kinetics theory]{@ref KineticsTheoryKinetics},
 *  - diffusion:
 *     * [constant Lewis]{@ref ConstantLewisDiffusivity},
 *     * [kinetics theory]{@ref MolecularBinaryDiffusion},
 *  - thermal conduction:
 *     * [Eucken]{@ref EuckenThermalConductivity},
 *     * [kinetics theory]{@ref KineticsTheoryThermalConductivity}.
 * The kinetics is all set up at run
 * time through input files.
 * 
 * For this example, we will choose
 * kinetics theory for everyone.
 *
 * First, we need to include the necessary
 * header files.
*/

/*! vectorization,
 * first the declarations, at the
 * end the implementation
 */
#include "antioch/eigen_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"


/*! kinetics,
 *  those files contains all that
 *  is necessary to compute the kinetics.
 */ 
#include "antioch/kinetics_evaluator.h"
#include "antioch/reaction_set.h"
#include "antioch/reaction.h"
#include "antioch/kinetics_parsing.h"
/*! transport,
 * 
 * 
 * 
 * 
 * 
 * 
 */

/*! viscosity,
 * the models, here we choose
 * kinetics theory
*/
//#include "antioch/sutherland_viscosity.h"
//#include "antioch/blottner_viscosity.h"
#include "antioch/kinetics_theory_viscosity.h"

/*! diffusion,
 * the models, here we choose
 * kinetics theory
 */
//#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/molecular_binary_diffusion.h"

/*! thermal conduction
 * the models, here we choose
 * kinetics theory
 */
//#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"

/*! Setting the transport
 */

#include "antioch/wilke_transport_evaluator.h"
#include "antioch/wilke_transport_mixture.h"
#include "antioch/physical_set.h"
#include "antioch/transport_mixture.h"

/*! The parsing,
 *  Antioch provides three kind of parsers:
 *    - ASCIIParser
 *    - XMLParser
 *    - ChemKinParser
 * 
 * 
 * 
 */




/*! metaprogrammation,
 * Antioch uses heavily generic
 * programming.  Antioch keeps a
 * declaration/implementation scheme.
 * Thus now the implementation for
 * both general generic programming and
 * the chosen vectorization, if any.
 */
#include "antioch/metaprogramming.h"
#include "antioch/eigen_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/vexcl_utils.h"
