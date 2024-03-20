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

// adapted from dumux/examples/1ptracer/properties.hh

#ifndef DUMUX_CELL_PROBLEM_PROPERTIES_HH
#define DUMUX_CELL_PROBLEM_PROPERTIES_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dune/grid/yaspgrid.hh>

#include "cell_problem/model.hh"
#include "problem_cellproblem.hh"

namespace Dumux::Properties {

namespace TTag {
struct CellProblem {
  using InheritsFrom = std::tuple<CellModel, CCTpfaModel>;
};
} // namespace TTag

template <class TypeTag>
struct Grid<TypeTag, TTag::CellProblem> {
  using type = Dune::SPGrid<double, 2>;
};

template <class TypeTag>
struct Problem<TypeTag, TTag::CellProblem> {
  using type = CellProblemProblem<TypeTag>;
};

template <class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::CellProblem> {
  static constexpr bool value = true;
};

template <class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::CellProblem> {
  static constexpr bool value = true;
};

} // end namespace Dumux::Properties

#endif
