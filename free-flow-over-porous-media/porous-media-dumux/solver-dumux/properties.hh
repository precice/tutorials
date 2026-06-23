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
 *
 * \brief The porous medium flow sub problem
 */
#ifndef DUMUX_DARCY_SUBPROPERTIES_HH
#define DUMUX_DARCY_SUBPROPERTIES_HH

#ifndef DIMWORLD
#define DIMWORLD 2
#endif

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cctpfa.hh>

#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DarcyOneP {
  using InheritsFrom = std::tuple<OneP, CCTpfaModel>;
};
} // end namespace TTag

// Set the problem property
template <class TypeTag>
struct Problem<TypeTag, TTag::DarcyOneP> {
  using type = Dumux::DarcySubProblem<TypeTag>;
};

// the fluid system
template <class TypeTag>
struct FluidSystem<TypeTag, TTag::DarcyOneP> {
  using Scalar = GetPropType<TypeTag, Properties::Scalar>;
  using type =
      FluidSystems::OnePLiquid<Scalar, Dumux::Components::SimpleH2O<Scalar>>;
};

// Set the grid type
template <class TypeTag>
struct Grid<TypeTag, TTag::DarcyOneP> {
  using type = Dune::YaspGrid<DIMWORLD>;
};

template <class TypeTag>
struct SpatialParams<TypeTag, TTag::DarcyOneP> {
  using type = OnePSpatialParams<GetPropType<TypeTag, GridGeometry>,
                                 GetPropType<TypeTag, Scalar>>;
};

} // end namespace Dumux::Properties

#endif // DUMUX_DARCY_SUBPROPERTIES_HH
