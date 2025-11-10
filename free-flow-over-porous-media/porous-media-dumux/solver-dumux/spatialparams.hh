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
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the 1p cc model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {
/*!
 * \ingroup OnePModel
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p cc model
 */
template <class FVGridGeometry, class Scalar>
class OnePSpatialParams : public FVPorousMediumFlowSpatialParamsOneP<
                              FVGridGeometry,
                              Scalar,
                              OnePSpatialParams<FVGridGeometry, Scalar>> {
  using GridView   = typename FVGridGeometry::GridView;
  using ParentType = FVPorousMediumFlowSpatialParamsOneP<
      FVGridGeometry,
      Scalar,
      OnePSpatialParams<FVGridGeometry, Scalar>>;

  using Element        = typename GridView::template Codim<0>::Entity;
  using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
  // export permeability type
  using PermeabilityType = Scalar;

  OnePSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
      : ParentType(fvGridGeometry)
  {
    permeability_ = getParam<Scalar>("SpatialParams.Permeability");
    porosity_     = getParam<Scalar>("SpatialParams.Porosity");
    alphaBJ_      = getParam<Scalar>("SpatialParams.AlphaBeaversJoseph");
  }

  /*!
   * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
   *
   * \param globalPos The global position
   * \return the intrinsic permeability
   */
  PermeabilityType permeabilityAtPos(const GlobalPosition &globalPos) const
  {
    return permeability_;
  }

  /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position
   */
  Scalar porosityAtPos(const GlobalPosition &globalPos) const
  {
    return porosity_;
  }

  /*! \brief Define the Beavers-Joseph coefficient in [-].
   *
   * \param globalPos The global position
   */
  Scalar beaversJosephCoeffAtPos(const GlobalPosition &globalPos) const
  {
    return alphaBJ_;
  }

  /*!
   * \brief Return the temperature within the domain in [K].
   *
   * This problem assumes a temperature of 10 degrees Celsius.
   */
  Scalar temperatureAtPos(const GlobalPosition &globalPos) const
  {
    return 273.15 + 10; // 10°C
  }

private:
  Scalar permeability_;
  Scalar porosity_;
  Scalar alphaBJ_;
};

} // namespace Dumux

#endif
