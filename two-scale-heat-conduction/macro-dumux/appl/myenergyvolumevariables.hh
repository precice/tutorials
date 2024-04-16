// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NIModel
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 */

#ifndef DUMUX_MY_ENERGY_VOLUME_VARIABLES_HH
#define DUMUX_MY_ENERGY_VOLUME_VARIABLES_HH

#include <dune/common/std/type_traits.hh>
#include <type_traits>

#include <dumux/material/solidsystems/1csolid.hh>
#include <dumux/porousmediumflow/volumevariables.hh>

namespace Dumux {

// forward declaration
template <class IsothermalTraits, class Impl, bool enableEnergyBalance>
class MyEnergyVolumeVariablesImplementation;

/*!
 * \ingroup NIModel
 * \brief Base class for the model specific class which provides
 *        access to all volume averaged quantities.
 *
 * The volume variables base class is specialized for isothermal and
 * non-isothermal models.
 */
template <class IsothermalTraits, class Impl>
using MyEnergyVolumeVariables = MyEnergyVolumeVariablesImplementation<
    IsothermalTraits, Impl,
    IsothermalTraits::ModelTraits::enableEnergyBalance()>;

/*!
 * \ingroup NIModel
 * \brief The isothermal base class
 */
template <class IsothermalTraits, class Impl>
class MyEnergyVolumeVariablesImplementation<IsothermalTraits, Impl, false> {
  using Scalar                  = typename IsothermalTraits::PrimaryVariables::value_type;
  static constexpr int dimWorld = 2; // hardcoded for now
  using DimWorldMatrix          = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
  using FluidState  = typename IsothermalTraits::FluidState;
  using SolidState  = typename IsothermalTraits::SolidState;
  using FluidSystem = typename IsothermalTraits::FluidSystem;

  //! The temperature is obtained from the problem as a constant for isothermal
  //! models
  template <class ElemSol, class Problem, class Element, class Scv>
  void updateTemperature(
      const ElemSol &elemSol, const Problem &problem, const Element &element,
      const Scv &scv, FluidState &fluidState,
      SolidState &solidState)
  { // retrieve temperature from solution vector,
    // all phases have the same temperature
    Scalar T = problem.spatialParams().temperature(element, scv, elemSol);
    for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
      fluidState.setTemperature(phaseIdx, T);
    }
    solidState.setTemperature(T);
  }

  template <class ElemSol, class Problem, class Element, class Scv>
  void updateSolidEnergyParams(const ElemSol &elemSol, const Problem &problem,
                               const Element &element, const Scv &scv,
                               SolidState &solidState) {}

  //! The phase enthalpy is zero for isothermal models
  //! This is needed for completing the fluid state
  template <class FluidState, class ParameterCache>
  static Scalar enthalpy(const FluidState     &fluidState,
                         const ParameterCache &paramCache, const int phaseIdx)
  {
    return 0;
  }

  //! The effective thermal conductivity is zero for isothermal models
  template <class ElemSol, class Problem, class Element, class Scv>
  void updateEffectiveThermalConductivity(const ElemSol &elemSol,
                                          const Problem &problem,
                                          const Element &element,
                                          const Scv     &scv,
                                          SolidState    &solidState) {}
};

//! The non-isothermal implicit volume variables base class
template <class Traits, class Impl>
class MyEnergyVolumeVariablesImplementation<Traits, Impl, true> {
  using Scalar                  = typename Traits::PrimaryVariables::value_type;
  static constexpr int dimWorld = 2; // hardcoded for now
  using DimWorldMatrix          = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
  using Idx                     = typename Traits::ModelTraits::Indices;
  using ParentType              = PorousMediumFlowVolumeVariables<Traits>;

  static constexpr int temperatureIdx = Idx::temperatureIdx;
  static constexpr int numEnergyEq    = Traits::ModelTraits::numEnergyEq();

  static constexpr bool fullThermalEquilibrium  = (numEnergyEq == 1);
  static constexpr bool fluidThermalEquilibrium = (numEnergyEq == 2);

public:
  // export the fluidstate
  using FluidState = typename Traits::FluidState;
  //! export the underlying fluid system
  using FluidSystem = typename Traits::FluidSystem;
  //! Export the indices
  using Indices = Idx;
  // export the solidstate
  using SolidState = typename Traits::SolidState;
  //! export the underlying solid system
  using SolidSystem = typename Traits::SolidSystem;

  //! The temperature is obtained from the problem as a constant for isothermal
  //! models
  template <class ElemSol, class Problem, class Element, class Scv>
  void updateTemperature(const ElemSol &elemSol, const Problem &problem,
                         const Element &element, const Scv &scv,
                         FluidState &fluidState, SolidState &solidState)
  {
    if constexpr (fullThermalEquilibrium) {
      // retrieve temperature from solution vector, all phases have the same
      // temperature
      const Scalar T = elemSol[scv.localDofIndex()][temperatureIdx];
      for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
        fluidState.setTemperature(phaseIdx, T);
      }
      solidState.setTemperature(T);
    }

    else {
      // this means we have 1 temp for fluid phase, one for solid
      if constexpr (fluidThermalEquilibrium) {
        const Scalar T = elemSol[scv.localDofIndex()][temperatureIdx];
        for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
          fluidState.setTemperature(phaseIdx, T);
        }
      }
      // this is for numEnergyEqFluid > 1
      else {
        for (int phaseIdx = 0; phaseIdx < FluidSystem::numPhases; ++phaseIdx) {
          // retrieve temperatures from solution vector, phases might have
          // different temperature
          const Scalar T =
              elemSol[scv.localDofIndex()][temperatureIdx + phaseIdx];
          fluidState.setTemperature(phaseIdx, T);
        }
      }
      const Scalar solidTemperature =
          elemSol[scv.localDofIndex()][temperatureIdx + numEnergyEq - 1];
      solidState.setTemperature(solidTemperature);
    }
  }

  template <class ElemSol, class Problem, class Element, class Scv>
  void updateSolidEnergyParams(const ElemSol &elemSol, const Problem &problem,
                               const Element &element, const Scv &scv,
                               SolidState &solidState)
  {
    Scalar cs = solidHeatCapacity_(elemSol, problem, element, scv, solidState);
    solidState.setHeatCapacity(cs);

    Scalar rhos = solidDensity_(elemSol, problem, element, scv, solidState);
    solidState.setDensity(rhos);
  }

  template <class ElemSol, class Problem, class Element, class Scv>
  void updateEffectiveThermalConductivity(const ElemSol &elemSol,
                                          const Problem &problem,
                                          const Element &element,
                                          const Scv     &scv,
                                          SolidState    &solidState)
  {
    lambdaEff_ =
        solidThermalConductivity_(elemSol, problem, element, scv, solidState);
  }

  /*!
   * \brief Returns the total internal energy of a phase in the
   *        sub-control volume.
   *
   * \param phaseIdx The phase index
   */
  Scalar internalEnergy(const int phaseIdx) const
  {
    return asImp_().fluidState().internalEnergy(phaseIdx);
  }

  /*!
   * \brief Returns the total enthalpy of a phase in the sub-control
   *        volume.
   *
   * \param phaseIdx The phase index
   */
  Scalar enthalpy(const int phaseIdx) const
  {
    return asImp_().fluidState().enthalpy(phaseIdx);
  }

  /*!
   * \brief Returns the temperature in fluid / solid phase(s)
   *        the sub-control volume.
   */
  Scalar temperatureSolid() const
  {
    return asImp_().solidState().temperature();
  }

  /*!
   * \brief Returns the temperature of a fluid phase assuming thermal
   * nonequilibrium the sub-control volume. \param phaseIdx The local index of
   * the phases
   */
  Scalar temperatureFluid(const int phaseIdx) const
  {
    return asImp_().fluidState().temperature(phaseIdx);
  }

  /*!
   * \brief Returns the total heat capacity \f$\mathrm{[J/(kg K)]}\f$ of the
   * rock matrix in the sub-control volume.
   */
  Scalar solidHeatCapacity() const
  {
    return asImp_().solidState().heatCapacity();
  }

  /*!
   * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of the rock matrix
   * in the sub-control volume.
   */
  Scalar solidDensity() const
  {
    return asImp_().solidState().density();
  }

  /*!
   * \brief Returns the effective thermal conductivity \f$\mathrm{[W/(m*K)]}\f$
   * in the sub-control volume. Specific to equilibirum models (case
   * fullThermalEquilibrium).
   */
  template <bool enable                   = fullThermalEquilibrium,
            std::enable_if_t<enable, int> = 0>
  DimWorldMatrix effectiveThermalConductivity() const
  {
    return lambdaEff_;
  }

  //! The phase enthalpy is zero for isothermal models
  //! This is needed for completing the fluid state
  template <class ParameterCache>
  static Scalar enthalpy(const FluidState     &fluidState,
                         const ParameterCache &paramCache, const int phaseIdx)
  {
    return FluidSystem::enthalpy(fluidState, paramCache, phaseIdx);
  }

protected:
  const Impl &asImp_() const
  {
    return *static_cast<const Impl *>(this);
  }
  Impl &asImp_()
  {
    return *static_cast<Impl *>(this);
  }

private:
  /*!
   * It has to be decided if the full solid system / solid state interface is
   * used (general option, but more complicated), or the simple nonisothermal
   * spatial params interface (simpler but less general). In the simple
   * nonisothermal spatial params interface the functions solidHeatCapacity,
   * solidDensity, and solidThermalConductivity in the spatial params overwrite
   * the parameters given in the solid system. This only makes sense in
   * combination with the simplest solid system InertSolidPhase, and can be used
   * to quickly change parameters in certain domain regions. For setups with
   * more general solids with several components these functions should not
   * exist. Instead, the solid system determines the values for
   * solidHeatCapacity, solidDensity, and solidThermalConductivity depending on
   * the given composition.
   */

  /*!
   * \name Access functions for the solidsystem / solidstate interface
   */
  // \{

  /*!
   * \brief Gets the solid heat capacity in an scv.
   *
   * \param elemSol the element solution vector
   * \param problem the problem to solve
   * \param element the element (codim-0-entity) the scv belongs to
   * \param scv the sub control volume
   * \param solidState the solid state
   * \note this gets selected if the user uses the solidsystem / solidstate
   * interface
   */
  template <class ElemSol, class Problem, class Element, class Scv,
            std::enable_if_t<!Detail::hasSolidHeatCapacity<
                                 typename Problem::SpatialParams, Element, Scv,
                                 ElemSol, SolidState>(),
                             int> = 0>
  Scalar solidHeatCapacity_(const ElemSol &elemSol, const Problem &problem,
                            const Element &element, const Scv &scv,
                            const SolidState &solidState)
  {
    return SolidSystem::heatCapacity(solidState);
  }

  /*!
   * \brief Gets the solid density in an scv.
   *
   * \param elemSol the element solution vector
   * \param problem the problem to solve
   * \param element the element (codim-0-entity) the scv belongs to
   * \param scv the sub control volume
   * \param solidState the solid state
   * \note this gets selected if the user uses the solidsystem / solidstate
   * interface
   */
  template <class ElemSol, class Problem, class Element, class Scv,
            std::enable_if_t<
                !Detail::hasSolidDensity<typename Problem::SpatialParams,
                                         Element, Scv, ElemSol, SolidState>(),
                int> = 0>
  Scalar solidDensity_(const ElemSol &elemSol, const Problem &problem,
                       const Element &element, const Scv &scv,
                       const SolidState &solidState)
  {
    return SolidSystem::density(solidState);
  }

  // \}

  /*!
   * \name Access functions for the simple nonisothermal spatial params
   * interface in combination with an InertSolidPhase as solid system
   */
  // \{

  /*!
   * \brief Gets the solid heat capacity in an scv.
   *
   * \param elemSol the element solution vector
   * \param problem the problem to solve
   * \param element the element (codim-0-entity) the scv belongs to
   * \param scv the sub control volume
   * \param solidState the solid state
   * \note this gets selected if the user uses the simple spatial params
   * interface in combination with an InertSolidPhase as solid system
   */
  template <class ElemSol, class Problem, class Element, class Scv,
            std::enable_if_t<Detail::hasSolidHeatCapacity<
                                 typename Problem::SpatialParams, Element, Scv,
                                 ElemSol, SolidState>(),
                             int> = 0>
  Scalar solidHeatCapacity_(const ElemSol &elemSol, const Problem &problem,
                            const Element &element, const Scv &scv,
                            const SolidState &solidState)
  {
    static_assert(Detail::isInertSolidPhase<SolidSystem>::value,
                  "solidHeatCapacity can only be overwritten in the spatial "
                  "params when the solid system is a simple InertSolidPhase\n"
                  "If you select a proper solid system, the solid heat "
                  "capacity will be computed as stated in the solid system!");
    std::cout << "sHC called in myenergyvolumevariables" << std::endl; // DEBUG
    return problem.spatialParams().solidHeatCapacity(element, scv, elemSol,
                                                     solidState);
  }

  /*!
   * \brief Gets the solid density in an scv.
   *
   * \param elemSol the element solution vector
   * \param problem the problem to solve
   * \param element the element (codim-0-entity) the scv belongs to
   * \param scv the sub control volume
   * \param solidState the solid state
   * \note this gets selected if the user uses the simple spatial params
   * interface in combination with an InertSolidPhase as solid system
   */
  template <class ElemSol, class Problem, class Element, class Scv,
            std::enable_if_t<
                Detail::hasSolidDensity<typename Problem::SpatialParams,
                                        Element, Scv, ElemSol, SolidState>(),
                int> = 0>
  Scalar solidDensity_(const ElemSol &elemSol, const Problem &problem,
                       const Element &element, const Scv &scv,
                       const SolidState &solidState)
  {
    static_assert(Detail::isInertSolidPhase<SolidSystem>::value,
                  "solidDensity can only be overwritten in the spatial params "
                  "when the solid system is a simple InertSolidPhase\n"
                  "If you select a proper solid system, the solid density will "
                  "be computed as stated in the solid system!");
    return problem.spatialParams().solidDensity(element, scv, elemSol,
                                                solidState);
  }

  /*!
   * \brief Gets the solid's thermal conductivity in an scv.
   *
   * \param elemSol the element solution vector
   * \param problem the problem to solve
   * \param element the element (codim-0-entity) the scv belongs to
   * \param scv the sub control volume
   * \param solidState the solid state
   * \note this gets selected if the user uses the solidsystem / solidstate
   * interface
   */
  template <class ElemSol, class Problem, class Element, class Scv,
            std::enable_if_t<!Detail::hasSolidThermalConductivity<
                                 typename Problem::SpatialParams, Element, Scv,
                                 ElemSol, SolidState>(),
                             int> = 0>
  DimWorldMatrix
  solidThermalConductivity_(const ElemSol &elemSol, const Problem &problem,
                            const Element &element, const Scv &scv,
                            const SolidState &solidState)
  {
    static_assert(
        Detail::isInertSolidPhase<SolidSystem>::value,
        "solidThermalConductivity can only be overwritten in the spatial "
        "params when the solid system is a simple InertSolidPhase\n"
        "If you select a proper solid system, the solid thermal conductivity "
        "will be computed as stated in the solid system!");
    return problem.spatialParams().solidThermalConductivity(element, scv);
  }

  DimWorldMatrix lambdaEff_;
};

} // end namespace Dumux

#endif
