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

#ifndef DUMUX_1PNI_CONDUCTION_PROBLEM_PROPERTIES_HH
#define DUMUX_1PNI_CONDUCTION_PROBLEM_PROPERTIES_HH

#include <cmath>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/material/fluidmatrixinteractions/1p/thermalconductivityaverage.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dune/grid/yaspgrid.hh>

#include "mysimpleliquid.hh"
#include "myvolumevariables.hh"
#include "problem.hh"
#include "spatialparams.hh"

namespace Dumux::Properties {
// Create new type tags
namespace TTag {
struct OnePNIConduction {
  using InheritsFrom = std::tuple<OnePNI>;
};
struct OnePNIConductionCCTpfa {
  using InheritsFrom = std::tuple<OnePNIConduction, CCTpfaModel>;
};
} // end namespace TTag

// Set the grid type
template <class TypeTag>
struct Grid<TypeTag, TTag::OnePNIConduction> {
  using type = Dune::YaspGrid<2>;
}; // structured parallel 2D grid
// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::OnePNIConduction> {
  using type = OnePNIConductionProblem<TypeTag>;
};

// Set the fluid system
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePNIConduction> {
  using type = FluidSystems::OnePLiquid<
      GetPropType<TypeTag, Properties::Scalar>,
      Components::MySimpleLiquid<GetPropType<TypeTag, Properties::Scalar>>>;
};

//! Set the volume variables property
template <class TypeTag>
struct VolumeVariables<TypeTag, TTag::OnePNIConduction> {
private:
  using PV  = GetPropType<TypeTag, Properties::PrimaryVariables>;
  using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
  using FST = GetPropType<TypeTag, Properties::FluidState>;
  using SSY = GetPropType<TypeTag, Properties::SolidSystem>;
  using SST = GetPropType<TypeTag, Properties::SolidState>;
  using PT  = typename GetPropType<TypeTag,
                                  Properties::SpatialParams>::PermeabilityType;
  using MT  = GetPropType<TypeTag, Properties::ModelTraits>;

  using Traits = OnePVolumeVariablesTraits<PV, FSY, FST, SSY, SST, PT, MT>;

public:
  using type = MyOnePVolumeVariables<Traits>;
};

// Set the spatial parameters
template <class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePNIConduction> {
  using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
  using Scalar       = GetPropType<TypeTag, Properties::Scalar>;
  using type         = OnePNISpatialParams<GridGeometry, Scalar>;
};

} // namespace Dumux::Properties
// end namespace Dumux
#endif
