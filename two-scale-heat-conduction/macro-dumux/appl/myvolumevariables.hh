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
 * \ingroup My OnePModel
 * \brief Quantities required by the one-phase fully implicit model defined on a
 * vertex. Copy of OnePModel, inherits from MyEnergyVolumeVariables
 */

#ifndef DUMUX_MY_1P_VOLUME_VARIABLES_HH
#define DUMUX_MY_1P_VOLUME_VARIABLES_HH

#include <dumux/material/fluidstates/immiscible.hh>
#include <dumux/material/solidstates/updatesolidvolumefractions.hh>
#include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
#include <dumux/porousmediumflow/volumevariables.hh>
#include "myenergyvolumevariables.hh"

namespace Dumux {

/*!
 * \ingroup OnePModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the one-phase model.
 *
 * \tparam Traits Class encapsulating types to be used by the vol vars
 */
template <class Traits>
class MyOnePVolumeVariables
    : public PorousMediumFlowVolumeVariables<Traits>,
      public MyEnergyVolumeVariables<Traits, MyOnePVolumeVariables<Traits>> {
  using ThisType   = MyOnePVolumeVariables<Traits>;
  using ParentType = PorousMediumFlowVolumeVariables<Traits>;
  using EnergyVolVars =
      MyEnergyVolumeVariables<Traits, MyOnePVolumeVariables<Traits>>;

  using Scalar                       = typename Traits::PrimaryVariables::value_type;
  using PermeabilityType             = typename Traits::PermeabilityType;
  static constexpr int numFluidComps = ParentType::numFluidComponents();
  static constexpr int dimWorld      = 2; // hardcoded for now
  using DimWorldMatrix               = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
  //! Export the underlying fluid system
  using FluidSystem = typename Traits::FluidSystem;
  //! Export the fluid state type
  using FluidState = typename Traits::FluidState;
  //! Export the indices
  using Indices = typename Traits::ModelTraits::Indices;
  //! Export type of solid state
  using SolidState = typename Traits::SolidState;
  //! Export type of solid system
  using SolidSystem = typename Traits::SolidSystem;

  /*!
   * \brief Updates all quantities for a given control volume.
   *
   * \param elemSol A vector containing all primary variables connected to the
   * element \param problem The object specifying the problem which ought to be
   * simulated \param element An element which contains part of the control
   * volume \param scv The sub-control volume
   */
  template <class ElemSol, class Problem, class Element, class Scv>
  void update(const ElemSol &elemSol, const Problem &problem,
              const Element &element, const Scv &scv)
  {
    ParentType::update(elemSol, problem, element, scv);

    // porosity
    completeFluidState(elemSol, problem, element, scv, fluidState_,
                       solidState_);

    // porosity and permeability
    updateSolidVolumeFractions(elemSol, problem, element, scv, solidState_,
                               numFluidComps);
    EnergyVolVars::updateSolidEnergyParams(elemSol, problem, element, scv,
                                           solidState_);
    permeability_ = problem.spatialParams().permeability(element, scv, elemSol);
    EnergyVolVars::updateEffectiveThermalConductivity(elemSol, problem, element,
                                                      scv, solidState_);
  }

  /*!
   * \brief Sets complete fluid state
   *
   * \param elemSol A vector containing all primary variables connected to the
   * element \param problem The object specifying the problem which ought to be
   * simulated \param element An element which contains part of the control
   * volume \param scv The sub-control volume \param fluidState A container with
   * the current (physical) state of the fluid \param solidState A container
   * with the current (physical) state of the solid
   */
  template <class ElemSol, class Problem, class Element, class Scv>
  void completeFluidState(const ElemSol &elemSol, const Problem &problem,
                          const Element &element, const Scv &scv,
                          FluidState &fluidState, SolidState &solidState)
  {
    EnergyVolVars::updateTemperature(elemSol, problem, element, scv, fluidState,
                                     solidState);
    fluidState.setSaturation(/*phaseIdx=*/0, 1.);

    const auto &priVars = elemSol[scv.localDofIndex()];
    fluidState.setPressure(/*phaseIdx=*/0, priVars[Indices::pressureIdx]);

    // saturation in a single phase is always 1 and thus redundant
    // to set. But since we use the fluid state shared by the
    // immiscible multi-phase models, so we have to set it here...
    fluidState.setSaturation(/*phaseIdx=*/0, 1.0);

    typename FluidSystem::ParameterCache paramCache;
    paramCache.updatePhase(fluidState, /*phaseIdx=*/0);

    Scalar value = FluidSystem::density(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setDensity(/*phaseIdx=*/0, value);

    value = FluidSystem::viscosity(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setViscosity(/*phaseIdx=*/0, value);

    // compute and set the enthalpy
    value = EnergyVolVars::enthalpy(fluidState, paramCache, /*phaseIdx=*/0);
    fluidState.setEnthalpy(/*phaseIdx=*/0, value);
  }

  /*!
   * \brief Returns the temperature \f$\mathrm{[K]}\f$ inside the sub-control
   * volume.
   *
   * Note that we assume thermodynamic equilibrium, i.e. the
   * temperatures of the rock matrix and of all fluid phases are
   * identical.
   */
  Scalar temperature() const
  {
    return fluidState_.temperature();
  }

  /*!
   * \brief Returns the phase state for the control volume.
   */
  const SolidState &solidState() const
  {
    return solidState_;
  }

  /*!
   * \brief Returns the effective pressure \f$\mathrm{[Pa]}\f$ of a given phase
   * within the control volume.
   */
  Scalar pressure(int phaseIdx = 0) const
  {
    return fluidState_.pressure(phaseIdx);
  }

  /*!
   * \brief Returns the saturation.
   */
  Scalar saturation(int phaseIdx = 0) const
  {
    return 1.0;
  }

  /*!
   * \brief Returns the mass density \f$\mathrm{[kg/m^3]}\f$ of a given phase
   * within the control volume.
   */
  Scalar density(int phaseIdx = 0) const
  {
    return fluidState_.density(phaseIdx);
  }

  /*!
   * \brief Returns the dynamic viscosity \f$\mathrm{[Pa s]}\f$ of the fluid
   * within the control volume.
   */
  Scalar viscosity(int phaseIdx = 0) const
  {
    return fluidState_.viscosity(phaseIdx);
  }

  /*!
   * \brief Returns the mobility \f$\mathrm{[1/(Pa s)]}\f$.
   *
   * This function enables the use of ImplicitDarcyFluxVariables
   * with the 1p fully implicit model, ALTHOUGH the term mobility is
   * usually not employed in the one phase context.
   *
   * \param phaseIdx The phase index
   */
  Scalar mobility(int phaseIdx = 0) const
  {
    return 1.0 / fluidState_.viscosity(phaseIdx);
  }

  /*!
   * \brief Returns the average porosity \f$\mathrm{[-]}\f$ within the control
   * volume.
   */
  Scalar porosity() const
  {
    return solidState_.porosity();
  }

  /*!
   * \brief Returns the permeability within the control volume in \f$[m^2]\f$.
   */
  const PermeabilityType &permeability() const
  {
    return permeability_;
  }

  /*!
   * \brief Returns the fluid state of the control volume.
   */
  const FluidState &fluidState() const
  {
    return fluidState_;
  }

  DimWorldMatrix effectiveThermalConductivity() const
  {
    return EnergyVolVars::effectiveThermalConductivity();
  }

protected:
  FluidState       fluidState_;
  SolidState       solidState_;
  PermeabilityType permeability_;
};

} // namespace Dumux

#endif
