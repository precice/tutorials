// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief The free flow sub problem
 */
#ifndef DUMUX_STOKES_SUBPROPERTIES_HH
#define DUMUX_STOKES_SUBPROPERTIES_HH

#ifndef DIMWORLD
#define DIMWORLD 2
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/freeflow/navierstokes/mass/1p/model.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>
#include <dumux/multidomain/freeflow/couplingmanager.hh>
#include <dumux/multidomain/traits.hh>

#include "problem.hh"

namespace Dumux {

namespace Properties {
// Create new type tags
namespace TTag {
struct FreeFlowModel {
};
struct FreeFlowSubMomentum {
  using InheritsFrom = std::
      tuple<FreeFlowModel, NavierStokesMomentum, FaceCenteredStaggeredModel>;
};
struct FreeFlowSubMass {
  using InheritsFrom =
      std::tuple<FreeFlowModel, NavierStokesMassOneP, CCTpfaModel>;
};
} // end namespace TTag

// the fluid system
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::FreeFlowModel> {
  using Scalar = GetPropType<TypeTag, Properties::Scalar>;
  using type =
      FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};

// Set the grid type
template <class TypeTag>
struct Grid<TypeTag, TTag::FreeFlowModel> {
  using type = Dune::YaspGrid<DIMWORLD,
                              Dune::EquidistantOffsetCoordinates<
                                  GetPropType<TypeTag, Properties::Scalar>,
                                  DIMWORLD>>;
};

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowSubMomentum> {
  using type =
      Dumux::StokesSubProblem<TypeTag,
                              Dumux::NavierStokesMomentumProblem<TypeTag>>;
};
template <class TypeTag>
struct Problem<TypeTag, TTag::FreeFlowSubMass> {
  using type =
      Dumux::StokesSubProblem<TypeTag,
                              Dumux::NavierStokesMassProblem<TypeTag>>;
};

template <class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::FreeFlowModel> {
  static constexpr bool value = true;
};
template <class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::FreeFlowModel> {
  static constexpr bool value = true;
};
template <class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::FreeFlowModel> {
  static constexpr bool value = true;
};

// Define the DuMux coupling manager to couple the momentum and mass subproblems
// of the freeflow participant
template <class TypeTag>
struct CouplingManager<TypeTag, TTag::FreeFlowModel> {
  using Traits =
      MultiDomainTraits<TTag::FreeFlowSubMomentum, TTag::FreeFlowSubMass>;
  using type = FreeFlowCouplingManager<Traits>;
};
} // end namespace Properties

} // namespace Dumux

#endif // DUMUX_STOKES_SUBPROPERTIES_HH
