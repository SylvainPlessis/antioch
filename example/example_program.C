//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//-----------------------------------------------------------------------el-

/*! \file example_program.C
 * \brief Example program fully commented
 *
 * This example shows how to set up kinetics
 * and transport using Antioch.
 * The models we have to define here
 * concerns the transport:
 *  - viscosity:
 *     * [Sutherland](@ref Antioch::SutherlandViscosity),
 *     * [Blottner](@ref Antioch::BlottnerViscosity),
 *     * [kinetics theory](@ref Antioch::KineticsTheoryViscosity),
 *  - diffusion:
 *     * [constant Lewis](@ref Antioch::ConstantLewisDiffusivity),
 *     * [kinetics theory](@ref Antioch::MolecularBinaryDiffusion),
 *  - thermal conduction:
 *     * [Eucken](@ref Antioch::EuckenThermalConductivity),
 *     * [kinetics theory](@ref Antioch::KineticsTheoryThermalConductivity).
 *
 * The kinetics is all set up at run
 * time through input files.
 * 
 * For this example, we will choose
 * kinetics theory for everyone.
 *
 * # Header files
 *
 * ## Vectorization
 * 
 * A first choice is the vectorization (or not)
 * of the program to be run. Antioch provides
 * several possibilities:
 *  - standard vector    (\c std::vector<>)
 *  - eigen vector       (e.g. \c Eigen::Array<scalar type, size , 1>)
 *  - valarray vector    (\c std::valarray<>)
 *  - metaphysicl vector (e.g. \c MetaPhysicL::NumberArray<size, scalar type>)
 *  - vexcl              (\c vex::vector<>)
 *
 * As Antioch uses heavily metaprogrammation,
 * she requires to know the necessary generic
 * methods to be used with the chosen vectorization.
 * The header files are designed ``classicaly'': a
 * declaration file and an implementation file.
 * Here comes the declaration files, named on the
 * design:\n
 * \c VECTORIZATION_CHOICE_utils_decl.h
 *
 * The file ``metaprogramming_decl.h'' is not necessary
 * here as it is included in every vectorization declaration,
 * but is shown for completeness. This file is required though
 * if you choose no vectorization but yet require (which will
 * certainly be the case) the generic capabilities with non-vectorized
 * types.
 *
 * ~~~~~~~~
#include "antioch/metaprogramming_decl.h" //general declarations
#include "antioch/vector_utils_decl.h"
#include "antioch/eigen_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"
 * ~~~~~~~~
 *
 * ## Kinetics
 * 
 * those files contains all that
 * is necessary to compute the kinetics.
 * ~~~~~~~~
#include "antioch/kinetics_evaluator.h"
#include "antioch/reaction_set.h"
#include "antioch/reaction.h"
#include "antioch/kinetics_parsing.h"
 * ~~~~~~~~
 *
 * ## Transport
 *
 * We need here to include the viscosity,
 * diffusion and thermal diffusion
 * header files.
 *
 * ### Viscosity
 *
 * the models, here we choose
 * kinetics theory
 * ~~~~~~~~~~~~~~~
//#include "antioch/sutherland_viscosity.h"
//#include "antioch/blottner_viscosity.h"
#include "antioch/kinetics_theory_viscosity.h"
 * ~~~~~~~~~~~~~~~
 *
 * ### Diffusion
 *
 * the models, here we choose
 * kinetics theory
 * ~~~~~~~~~~~~~~~~~~~
//#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/molecular_binary_diffusion.h"
 * ~~~~~~~~~~~~~~~~~~~
 *
 * ### Thermal conduction
 *
 * the models, here we choose
 * kinetics theory
 *
 * ~~~~~~~~~~~~~~~~~~~~
//#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"
 * ~~~~~~~~~~~~~~~~~~~~
 *
 * ## Setting the transport
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~
#include "antioch/wilke_transport_evaluator.h"
#include "antioch/wilke_transport_mixture.h"
#include "antioch/physical_set.h"
#include "antioch/transport_mixture.h"
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ## The parsing
 *
 * Antioch provides three kind of parsers:
 *  - [ASCIIParser](@ref Antioch::ASCIIParser)
 *  - [XMLParser](@ref Antioch::XMLParser)
 *  - [ChemKinParser](@ref Antioch::ChemKinParser)
 * 
 * We need to include only the necessary format
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "antioch/ascii_parser.h"
#include "antioch/xml_parser.h"
//#include "antioch/chemkin_parser.h"
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * ## Metaprogrammation
 *
 * Antioch uses heavily generic
 * programming.  Antioch keeps a
 * declaration/implementation scheme.
 * Thus now the implementation for
 * both general generic programming and
 * the chosen vectorization, if any.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "antioch/metaprogramming.h"
#include "antioch/eigen_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/vexcl_utils.h"
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * # The actual code
 *
 * ~~~~~~~~~~~~~~~~~
int main()
{

   typedef double CoeffScalar
   typedef double StateNumeric
 * ~~~~~~~~~~~~~~~~~
 * The composition of the chemical mixture
 * is either in a file, here noted by `std::string` _chemical_species_file_ or 
 * directly in `std::vector<std::string>` _species_.
 * The species' data in `std::string` _species_data_file_,
 * the species' vibrational data in `std::string` _species_vib_data_file_,
 * the species' electronic data in `std::string` _species_elec_data_file_.
 *
 * In this case, we consider the definition of the mixture in a file,
 * formatted in the ASCII format, as defined by the [ASCIIParser](@ref Antioch::ASCIIParser) class.
 *
 *
 *
 *
 *
 *
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~
   Antioch::ASCIIParser<CoeffScalar> chemical_parser(chemical_species_file,true);
   Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(&chemical_parser,
                                                          species_data_file,
                                                          species_vib_data_file,
                                                          species_elec_data_file);
 * ~~~~~~~~~~~~~~~~~~~~
 * The other possibilities, commented out here:
 * ~~~~~~~~~~~~~~~~~~~~~~~~
  // Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(chemical_species_file, true,
  //                                                        species_data_file,
  //                                                        species_vib_data_file,
  //                                                        species_elec_data_file);
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 * and
 * ~~~~~~~~~~~~~~~~~~~~~~~~
  // Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(species, true,
  //                                                        species_data_file,
  //                                                        species_vib_data_file,
  //                                                        species_elec_data_file);
 * ~~~~~~~~~~~~~~~~~~~~~~~~
 * Now we define the reaction set:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
   Antioch::ReactionSet<CoeffScalar> reaction_set(chemical_mixture);
 * ~~~~~~~~~~~~~~~~~~~~~~~~~
 * the reaction set is defined in 'std::string' _reaction_set_file_, this
 * file is formatted in XML 'Antioch::XML' or as defined by the
 * ChemKin program 'Antioch::CHEMKIN'
 * ~~~~~~~~~~~~~~~~~~
   Antioch::read_reaction_set_data(reaction_set_file, true, reaction_set,Antioch::XML);

   StateNumeric StateNumeric_exemple;

   Antioch::KineticsEvaluator<CoeffScalar,StateNumeric> kinetics_evaluator(reaction_set, StateNumeric_exemple);


   return 0;
}
 * ~~~~~~~~~~~~~~~~~~
 */

#include "antioch/metaprogramming_decl.h" //general declarations
#include "antioch/eigen_utils_decl.h"
#include "antioch/valarray_utils_decl.h"
#include "antioch/metaphysicl_utils_decl.h"
#include "antioch/vexcl_utils_decl.h"

#include "antioch/kinetics_evaluator.h"
#include "antioch/reaction_set.h"
#include "antioch/reaction.h"
#include "antioch/kinetics_parsing.h"

//#include "antioch/sutherland_viscosity.h"
//#include "antioch/blottner_viscosity.h"
#include "antioch/kinetics_theory_viscosity.h"
//#include "antioch/constant_lewis_diffusivity.h"
#include "antioch/molecular_binary_diffusion.h"
//#include "antioch/eucken_thermal_conductivity.h"
#include "antioch/kinetics_theory_thermal_conductivity.h"

#include "antioch/wilke_transport_evaluator.h"
#include "antioch/wilke_transport_mixture.h"
#include "antioch/physical_set.h"
#include "antioch/transport_mixture.h"

#include "antioch/ascii_parser.h"
#include "antioch/xml_parser.h"
//#include "antioch/chemkin_parser.h"

#include "antioch/metaprogramming.h"  //define generic programming for non-vectorized type
#include "antioch/eigen_utils.h"
#include "antioch/valarray_utils.h"
#include "antioch/metaphysicl_utils.h"
#include "antioch/vexcl_utils.h"

int main()
{

   typedef double CoeffScalar
   typedef double StateNumeric

   Antioch::ASCIIParser<CoeffScalar> chemical_parser(chemical_species_file,true);
   Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(&chemical_parser,
                                                          species_data_file,
                                                          species_vib_data_file,
                                                          species_elec_data_file);
 /*
   Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(chemical_species_file, true,
                                                          species_data_file,
                                                          species_vib_data_file,
                                                          species_elec_data_file);

   Antioch::ChemicalMixture<CoeffScalar> chemical_mixture(species, true,
                                                          species_data_file,
                                                          species_vib_data_file,
                                                          species_elec_data_file);

*/

   Antioch::ReactionSet<CoeffScalar> reaction_set(chemical_mixture);

   Antioch::read_reaction_set_data(reaction_set_file, true, reaction_set,Antioch::XML);

   StateNumeric StateNumeric_exemple;

   Antioch::KineticsEvaluator<CoeffScalar,StateNumeric> kinetics_evaluator(reaction_set, StateNumeric_exemple);


   return 0;
}
