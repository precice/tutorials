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

#ifndef DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH
#define DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH

#include <dune/grid/common/partitionset.hh>

#include <dumux/parallel/vectorcommdatahandle.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/porousmediumflow/properties.hh>

#include <dumux-precice/couplingadapter.hh>

namespace Dumux {

/*!
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
template <class GridGeometry, class Scalar>
class OnePNISpatialParams
    : public FVPorousMediumFlowSpatialParamsOneP<
          GridGeometry, Scalar, OnePNISpatialParams<GridGeometry, Scalar>> {
  using GridView = typename GridGeometry::GridView;
  using ThisType = OnePNISpatialParams<GridGeometry, Scalar>;
  using ParentType =
      FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
  static const int dimWorld = GridView::dimensionworld;
  using DimWorldMatrix      = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
  using Element             = typename GridView::template Codim<0>::Entity;
  using GlobalPosition      = typename Element::Geometry::GlobalCoordinate;
  using FVElementGeometry   = typename GridGeometry::LocalView;
  using SubControlVolume    = typename FVElementGeometry::SubControlVolume;

public:
  // export permeability type
  using PermeabilityType = Scalar;

  OnePNISpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
      : ParentType(gridGeometry), couplingParticipant_(Dumux::Precice::CouplingAdapter::getInstance()), couplingData_(gridGeometry->numDofs()), couplingDataHandle_(this->gridGeometry().elementMapper(), couplingData_)
  {
  }

  /*!
   * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
   *
   * \param globalPos The global position
   */
  PermeabilityType permeabilityAtPos(const GlobalPosition &globalPos) const
  {
    return getParam<Scalar>("Problem.Permeability");
  }

  /*!
   * \brief Defines the porosity \f$\mathrm{[-]}\f$.
   *
   * \param globalPos The global position
   */
  template <class ElementSolution>
  Scalar porosity(const Element &element, const SubControlVolume &scv,
                  const ElementSolution &elemSol) const
  {
    if (getParam<bool>("Precice.RunWithCoupling") == true)
      return couplingData_[scv.elementIndex()][0];
    else
      return getParam<Scalar>("Problem.DefaultPorosity");
  }

  /*!
   * \brief Defines the conductivity tensor \f$ K \f$.
   */

  DimWorldMatrix solidThermalConductivity(const Element          &element,
                                          const SubControlVolume &scv) const
  {
    DimWorldMatrix K;

    if (getParam<bool>("Precice.RunWithCoupling") == true) {
      K[0][0] = couplingData_[scv.elementIndex()][1];
      K[0][1] = couplingData_[scv.elementIndex()][2];
      K[1][0] = couplingData_[scv.elementIndex()][3];
      K[1][1] = couplingData_[scv.elementIndex()][4];
    } else {
      K[0][0] = getParam<Scalar>("Component.SolidThermalConductivity");
      K[0][1] = 0.0;
      K[1][0] = 0.0;
      K[1][1] = getParam<Scalar>("Component.SolidThermalConductivity");
    }
    return K;
  }

  void updateCouplingData()
  {
    for (const auto &element : elements(this->gridGeometry().gridView(), Dune::Partitions::interior)) {
      auto fvGeometry = localView(this->gridGeometry());
      fvGeometry.bindElement(element);
      for (const auto &scv : scvs(fvGeometry)) {
        const auto elementIdx = scv.elementIndex();
        couplingData_[elementIdx][0] =
            couplingParticipant_.getScalarQuantityOnFace("Macro-Mesh", "Porosity", elementIdx);
        couplingData_[elementIdx][1] =
            couplingParticipant_.getScalarQuantityOnFace("Macro-Mesh", "K00", elementIdx);
        couplingData_[elementIdx][2] =
            couplingParticipant_.getScalarQuantityOnFace("Macro-Mesh", "K01", elementIdx);
        couplingData_[elementIdx][3] =
            couplingParticipant_.getScalarQuantityOnFace("Macro-Mesh", "K10", elementIdx);
        couplingData_[elementIdx][4] =
            couplingParticipant_.getScalarQuantityOnFace("Macro-Mesh", "K11", elementIdx);
      }
    }
    // Trigger exchange of coupling data between neighboring ranks, if the domain is partitioned
    if (this->gridGeometry().gridView().comm().size() > 1) {
      this->gridGeometry().gridView().communicate(couplingDataHandle_,
                                                  Dune::InteriorBorder_All_Interface, Dune::ForwardCommunication);
    }
  }

private:
  Dumux::Precice::CouplingAdapter                &couplingParticipant_;
  Dune::BlockVector<Dune::FieldVector<double, 5>> couplingData_;
  Dumux::VectorCommDataHandleEqual<
      typename GridGeometry::ElementMapper,
      Dune::BlockVector<Dune::FieldVector<double, 5>>,
      /* Entity codimension = */ 0>
      couplingDataHandle_;
};

} // end namespace Dumux

#endif
