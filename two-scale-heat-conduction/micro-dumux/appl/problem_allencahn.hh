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

#ifndef DUMUX_ALLENCAHN_PROBLEM_HH
#define DUMUX_ALLENCAHN_PROBLEM_HH

#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/common/integrate.hh>
#include <dumux/common/timeloop.hh>
#include <dune/istl/matrix.hh>
#include <fstream>
#include <iostream>

namespace Dumux {

template <class TypeTag>
class AllenCahnProblem : public FVProblemWithSpatialParams<TypeTag> {
  using ParentType = FVProblemWithSpatialParams<TypeTag>;
  using GridView =
      typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
  using FVGridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
  using Scalar         = GetPropType<TypeTag, Properties::Scalar>;
  using Indices =
      typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
  using BoundaryTypes = Dumux::BoundaryTypes<
      GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
  using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
  using NumEqVector      = Dumux::NumEqVector<PrimaryVariables>;
  using ElementVolumeVariables =
      typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
  using FVElementGeometry =
      typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
  using SubControlVolume = typename FVElementGeometry::SubControlVolume;
  using SolutionVector   = GetPropType<TypeTag, Properties::SolutionVector>;

  static constexpr int phiIdx = Indices::phiIdx;

  using Element        = typename GridView::template Codim<0>::Entity;
  using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
  AllenCahnProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
      : ParentType(fvGridGeometry)
  {
    omega_   = getParam<Scalar>("Problem.omega");
    alpha_   = 1.0;
    xi_      = getParam<Scalar>("Problem.xi");
    kt_      = getParam<Scalar>("Problem.kt");
    eqconc_  = getParam<Scalar>("Problem.eqconc");
    centerX_ = (getParam<GlobalPosition>("Grid.UpperRight")[0] -
                getParam<GlobalPosition>("Grid.LowerLeft")[0]) /
               2;
    centerY_ = (getParam<GlobalPosition>("Grid.UpperRight")[1] -
                getParam<GlobalPosition>("Grid.LowerLeft")[1]) /
               2;
    radius_  = getParam<Scalar>("Problem.Radius");
    factor_  = getParam<Scalar>("Problem.PhasefieldICScaling");
    maxPoro_ = getParam<Scalar>("Problem.MaxPorosity");
  }

  BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
  {
    BoundaryTypes bcTypes;

    bcTypes.setAllNeumann();

    return bcTypes;
  }

  /*!
   * \brief Calculates the derivative P' of the double-well potential \f$ P = 8
   * \phi^2 (1-\phi)^2\f$.
   */
  Scalar pPrime(Scalar phi) const
  {
    return 16.0 * phi * (1.0 - phi) * (1.0 - 2.0 * phi);
  }

  /*!
   * \brief The source term is calculated as \f$ -\omega *P'(Phi) - 4*xi_ * Phi
   * * (1-Phi)*F(T)\f$.
   */
  NumEqVector source(const Element                &element,
                     const FVElementGeometry      &fvGeometry,
                     const ElementVolumeVariables &elemVolVars,
                     const SubControlVolume       &scv) const
  {
    NumEqVector source;

    const auto &priVars = elemVolVars[scv].priVars();

    source[phiIdx] = -omega_ * pPrime(priVars[phiIdx]);
    source +=
        4 * xi_ * priVars[phiIdx] * (1.0 - priVars[phiIdx]) * reactionRate();
    return source;
  }

  /*!
   * \brief Set the macro temperature/concentration
   */
  void updateConcentration(Scalar conc)
  {
    conc_ = conc;
  }

  /*!
   * \brief Returns the interfaceVelocity to use in the source term (-F(T))
   */
  Scalar reactionRate() const
  {
    return -kt_ *
           ((concentration() / eqconc_) * (concentration() / eqconc_) - 1);
  }

  /*!
   * \brief Returns the macro temperature/concentration
   */
  Scalar concentration() const
  {
    return conc_;
  }

  /*!
   * \brief Returns the initial analytic value of the phasefield
   * \note \f$ 1/(1+\exp(-4\lambda \sqrt{(y-y_0)² + (x-x_0)²}-R_0 )) \f$,
   *       where \f$|(x_0,y_0|\f$ is cell center
   *       and \f$ R_0\f$ is the initial radius of the grain.
   */
  PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
  {
    PrimaryVariables values;

    Scalar s =
        std::sqrt((globalPos[0] - centerX_) * (globalPos[0] - centerX_) +
                  (globalPos[1] - centerY_) * (globalPos[1] - centerY_)) -
        radius_;
    values[phiIdx] = 1.0 / (1.0 + std::exp(-factor_ * s / xi_));
    return values;
  }

  Scalar getAlpha() const
  {
    return alpha_;
  }

  Scalar getOmega() const
  {
    return omega_;
  }

  /*!
   * \brief Calculates the upscaled porosity by integrating phi
   */
  Scalar calculatePorosity(SolutionVector &sol) const
  {
    std::size_t order = 2;
    Scalar      poro  = integrateGridFunction(this->gridGeometry(), sol, order);
    if (poro <= maxPoro_)
      return poro;
    else
      return maxPoro_;
  }

private:
  Scalar              xi_;
  Scalar              omega_;
  Scalar              alpha_;
  Scalar              kt_;
  Scalar              eqconc_;
  std::vector<Scalar> poro_;
  Scalar              conc_;
  Scalar              maxPoro_;

  Scalar centerX_;
  Scalar centerY_;
  Scalar radius_;
  Scalar factor_;
};

} // end namespace Dumux

#endif
