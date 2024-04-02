// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
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

#ifndef CELL_PROBLEM_MODEL_HH
#define CELL_PROBLEM_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/common/properties/model.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <string>

#include "../spatialparams_cellproblem.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"

namespace Dumux {

struct CellProblemModelTraits {
  using Indices = CellProblemIndices<>;
  static constexpr int numEq()
  {
    return 2;
  }
  static constexpr int numComponents()
  {
    return 1;
  }
};

template <class PV, class MT>
struct CellProblemVolumeVariablesTraits {
  using PrimaryVariables = PV;
  using ModelTraits      = MT;
};

namespace Properties {
namespace TTag {
struct CellModel {
  using InheritsFrom = std::tuple<ModelProperties>;
};
} // namespace TTag

template <class TypeTag>
struct SpatialParams<TypeTag, TTag::CellModel> {
private:
  using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
  using Scalar       = GetPropType<TypeTag, Properties::Scalar>;

public:
  using type = CellProblemSpatialParams<GridGeometry, Scalar>;
};

template <class TypeTag>
struct LocalResidual<TypeTag, TTag::CellModel> {
  using type = CellProblemLocalResidual<TypeTag>;
};

template <class TypeTag>
struct ModelTraits<TypeTag, TTag::CellModel> {
  using type = CellProblemModelTraits;
};

template <class TypeTag>
struct VolumeVariables<TypeTag, TTag::CellModel> {
private:
  using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
  using MT = GetPropType<TypeTag, Properties::ModelTraits>;

  using Traits = CellProblemVolumeVariablesTraits<PV, MT>;

public:
  using type = CellProblemVolumeVariables<Traits>;
};

template <class TypeTag>
struct FluxVariablesCache<TypeTag, TTag::CellModel> {
  using type = FluxVariablesCaching::EmptyCache<
      GetPropType<TypeTag, Properties::Scalar>>;
};

template <class TypeTag>
struct FluxVariablesCacheFiller<TypeTag, TTag::CellModel> {
  using type = FluxVariablesCaching::EmptyCacheFiller;
};

} // end namespace Properties
} // end namespace Dumux

#endif