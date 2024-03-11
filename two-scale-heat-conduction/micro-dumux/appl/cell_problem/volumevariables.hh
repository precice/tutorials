// -*- mode: C++; tab-width: 3; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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

#ifndef CELL_PROBLEM_VOLUME_VARIABLES_HH
#define CELL_PROBLEM_VOLUME_VARIABLES_HH

#include <dumux/phasefield/volumevariables.hh>

namespace Dumux {

template <class Traits>
class CellProblemVolumeVariables : public PhasefieldVolumeVariables<Traits> {
  using Scalar = typename Traits::PrimaryVariables::value_type;

public:
  template <class Problem, class Scv>
  Scalar phi0delta(const Problem &problem, const Scv &scv) const
  {
    return problem.spatialParams().phi0delta(scv);
  }
};

} // end namespace Dumux

#endif
