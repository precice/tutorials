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

#ifndef DUMUX_CELL_PROBLEM_SPATIAL_PARAMS_HH
#define DUMUX_CELL_PROBLEM_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template <class GridGeometry, class Scalar>
class CellProblemSpatialParams
    : public FVPorousMediumFlowSpatialParamsOneP<
          GridGeometry, Scalar,
          CellProblemSpatialParams<GridGeometry, Scalar>> {
  using GridView          = typename GridGeometry::GridView;
  using FVElementGeometry = typename GridGeometry::LocalView;
  using SubControlVolume  = typename FVElementGeometry::SubControlVolume;
  using Element           = typename GridView::template Codim<0>::Entity;
  using ParentType        = FVPorousMediumFlowSpatialParamsOneP<
      GridGeometry, Scalar, CellProblemSpatialParams<GridGeometry, Scalar>>;
  static constexpr int dimWorld = GridView::dimensionworld;
  using GlobalPosition          = typename SubControlVolume::GlobalPosition;

  using AllenCahnTypeTag = Properties::TTag::AllenCahn;
  using ACSolutionVector =
      GetPropType<AllenCahnTypeTag, Properties::SolutionVector>;
  using Vector = Dune::FieldVector<Scalar, 1>;
  using PhasefieldIndices =
      typename GetPropType<AllenCahnTypeTag, Properties::ModelTraits>::Indices;

  enum { phiIdx = PhasefieldIndices::phiIdx };

public:
  using PermeabilityType = Scalar;

  CellProblemSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
      : ParentType(gridGeometry)
  {
    ks_      = getParam<Scalar>("Problem.ks");
    kg_      = getParam<Scalar>("Problem.kg");
    numDofs_ = gridGeometry->numDofs();
    phi_.resize(gridGeometry->numDofs());
    updatePhiCalled_ = false;
  }

  void updatePhi(ACSolutionVector &sol)
  {
    for (int i = 0; i < numDofs_; i++) {
      phi_[i] = sol[i];
    }
    updatePhiCalled_ = true;
  }

  Scalar phasefield(const SubControlVolume &scv) const
  {
    return phi_[scv.elementIndex()];
  }

  /*!
   * \brief Returns \f$ \Phi*k_s + (1-\Phi)*k_g\f$ for by scv.
   */
  Scalar phi0delta(const SubControlVolume &scv) const
  {
    return phasefield(scv) * ks_ + (1 - phasefield(scv)) * kg_;
  }

  /*!
   * \brief Returns \f$ \Phi*k_s + (1-\Phi)*k_g\f$ for by elementIndex.
   */
  Scalar phi0deltaIdx(int idx)
  {
    return ks_ * phi_[idx] + kg_ * (1 - phi_[idx]);
  }

private:
  std::vector<Scalar> phi_;
  Scalar              ks_;
  Scalar              kg_;
  Scalar              numDofs_;
  bool                updatePhiCalled_;
};
} // end namespace Dumux

#endif
