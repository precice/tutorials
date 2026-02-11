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
#ifndef DUMUX_STOKES_SUBPROBLEM_HH
#define DUMUX_STOKES_SUBPROBLEM_HH

#include <dune/common/fvector.hh>

#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/mass/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/problem.hh>

#include <dumux/freeflow/navierstokes/mass/1p/advectiveflux.hh>
#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

#include <dumux-precice/couplingadapter.hh>

namespace Dumux {

/*!
 * \brief The free flow sub problem
 */
template <class TypeTag, class BaseProblem>
class StokesSubProblem : public BaseProblem {
  using ParentType = BaseProblem;

  using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
  using GridView     = typename GridGeometry::GridView;
  using Scalar       = GetPropType<TypeTag, Properties::Scalar>;

  using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
  using Indices     = typename ModelTraits::Indices;

  using BoundaryTypes   = typename ParentType::BoundaryTypes;
  using InitialValues   = typename ParentType::InitialValues;
  using Sources         = typename ParentType::Sources;
  using DirichletValues = typename ParentType::DirichletValues;
  using BoundaryFluxes  = typename ParentType::BoundaryFluxes;

  using FVElementGeometry = typename GridGeometry::LocalView;
  using SubControlVolumeFace =
      typename FVElementGeometry::SubControlVolumeFace;
  using Element = typename GridView::template Codim<0>::Entity;

  using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

  using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

  using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

  using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

  static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
  using VelocityVector           = Dune::FieldVector<Scalar, dimWorld>;

public:
  StokesSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<CouplingManager>    couplingManager)
      : ParentType(gridGeometry, couplingManager, "FreeFlow"),
        eps_(1e-6),
        couplingParticipant_(Dumux::Precice::CouplingAdapter::getInstance())
  {
    deltaP_ = getParamFromGroup<Scalar>(this->paramGroup(),
                                        "Problem.PressureDifference");
  }

  /*!
   * \name Problem parameters
   */
  // \{

  /*!
   * \brief Return the sources within the domain.
   *
   * \param globalPos The global position
   */
  Sources sourceAtPos(const GlobalPosition &globalPos) const
  {
    return Sources(0.0);
  }
  // \}

  /*!
   * \name Boundary conditions
   */
  // \{

  /*!
   * \brief Specifies which kind of boundary condition should be
   *        used for which equation on a given boundary segment.
   *
   * \param element The finite element
   * \param scvf The sub control volume face
   */
  BoundaryTypes boundaryTypes(const Element              &element,
                              const SubControlVolumeFace &scvf) const
  {
    BoundaryTypes values;

    const auto &globalPos = scvf.center();

    if constexpr (ParentType::isMomentumProblem()) {
      if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
        values.setAllNeumann();
      }
      // slip boundary with coupling interface
      else if (onLowerBoundary_(globalPos)) {
        values.setAllNeumann();
        // TODO: Check the handling of the corners
        if (!onLeftBoundary_(scvf.ipGlobal()) &&
            !onRightBoundary_(scvf.ipGlobal()))
          values.setDirichlet(Indices::velocityYIdx);
      } else {
        values.setAllDirichlet();
      }
    } else // mass subproblem
    {
      if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
        values.setNeumann(Indices::conti0EqIdx);
      } else {
        values.setNeumann(Indices::conti0EqIdx);
      }
    }

    return values;
  }

  /*!
   * \brief Evaluate the boundary conditions for a Dirichlet control volume.
   *
   * \param globalPos The global position
   */
  DirichletValues dirichlet(const Element              &element,
                            const SubControlVolumeFace &scvf) const
  {
    DirichletValues values(0.0);
    values = initialAtPos(scvf.ipGlobal());

    if constexpr (ParentType::isMomentumProblem()) {
      const auto faceId = scvf.index();
      if (couplingParticipant_.isCoupledEntity(faceId)) {
        values[Indices::velocityYIdx] =
            couplingParticipant_.getScalarQuantityOnFace(
                "Free-Flow-Mesh", "Velocity", faceId);
      }
    } else {
      auto pressure                = onLeftBoundary_(scvf.ipGlobal()) ? deltaP_ : 0.0;
      values[Indices::pressureIdx] = pressure;
    }

    return values;
  }

  /*!
   * \brief Evaluate the boundary conditions for a Neumann control volume.
   *
   * \param element The element for which the Neumann boundary condition is set
   * \param fvGeomentry The fvGeometry
   * \param elemVolVars The element volume variables
   * \param elemFaceVars The element face variables
   * \param scvf The boundary sub control volume face
   */
  template <class ElementVolumeVariables, class ElementFluxVariablesCache>
  BoundaryFluxes neumann(const Element                   &element,
                         const FVElementGeometry         &fvGeometry,
                         const ElementVolumeVariables    &elemVolVars,
                         const ElementFluxVariablesCache &elemFluxVarsCache,
                         const SubControlVolumeFace      &scvf) const
  {
    BoundaryFluxes values(0.0);

    const auto &globalPos = scvf.ipGlobal();
    if constexpr (ParentType::isMomentumProblem()) {
#if DUMUX_VERSION_MAJOR >= 3 & DUMUX_VERSION_MINOR >= 9
      using SlipVelocityPolicy = NavierStokesSlipVelocity<
          typename GridGeometry::DiscretizationMethod,
          NavierStokes::SlipConditions::BJ>;
      using FluxHelper = NavierStokesMomentumBoundaryFlux<
          typename GridGeometry::DiscretizationMethod,
          SlipVelocityPolicy>;
#else
      using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;
#endif
      if (onSlipBoundary(fvGeometry, scvf)) {
        values += FluxHelper::slipVelocityMomentumFlux(
            *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache);
      } else if (onLeftBoundary_(globalPos) ||
                 onRightBoundary_(globalPos)) {
        auto pressure = onLeftBoundary_(globalPos) ? deltaP_ : 0.0;
        values        = FluxHelper::fixedPressureMomentumFlux(
                   *this, fvGeometry, scvf, elemVolVars, elemFluxVarsCache,
                   pressure, /* zeroNormalVelocityGradient = */ true);
      }
    } else {
      using FluxHelper = NavierStokesScalarBoundaryFluxHelper<
          AdvectiveFlux<ModelTraits>>;
      if (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)) {
        values = FluxHelper::scalarOutflowFlux(
            *this, element, fvGeometry, scvf, elemVolVars);
      } else if (onSlipBoundary(fvGeometry, scvf)) {
        const Scalar density =
            1000; // TODO how to handle compressible fluids?
        // TODO: Use flux helper with outside data?
        // TODO: remove hard-coded values index y=1 and dirsign = -1.0
        values[Indices::conti0EqIdx] =
            density * this->faceVelocity(element, fvGeometry, scvf)[1] *
            -1.0; // scvf.directionSign();
      }
    }
    return values;
  }

  bool onSlipBoundary(const FVElementGeometry    &fvGeometry,
                      const SubControlVolumeFace &scvf) const
  {
    return onLowerBoundary_(scvf.ipGlobal());
  }

  // \}

  /*!
   * \name Volume terms
   */
  // \{

  /*!
   * \brief Evaluate the initial value for a control volume.
   *
   * \param globalPos The global position
   */
  InitialValues initialAtPos(const GlobalPosition &globalPos) const
  {
    InitialValues values(0.0);

    if constexpr (ParentType::isMomentumProblem()) {
      // values[Indices::velocityYIdx] = -1e-6 * globalPos[0] * (this->gridGeometry().bBoxMax()[0] - globalPos[0]);
    } else {
      if (onLeftBoundary_(globalPos))
        values[Indices::pressureIdx] = deltaP_;
      if (onRightBoundary_(globalPos))
        values[Indices::pressureIdx] = 0.0;
    }

    return values;
  }

  /*!
   * \brief Returns the intrinsic permeability of required as input parameter for the Beavers-Joseph-Saffman boundary condition
   */
  Scalar permeability(const Element              &element,
                      const SubControlVolumeFace &scvf) const
  {
    return 1e-10; // TODO transfer information or just use constant value
  }

  /*!
   * \brief Returns the alpha value required as input parameter for the Beavers-Joseph-Saffman boundary condition
   */
  Scalar alphaBJ(const FVElementGeometry    &fvGeometry,
                 const SubControlVolumeFace &scvf) const
  {
    return 1.0; // TODO transfer information or just use constant value
  }

  /*!
   * \brief Returns the beta value required as input parameter for the Beavers-Joseph-Saffman boundary condition
   */
  Scalar betaBJ(const FVElementGeometry    &fvGeometry,
                const SubControlVolumeFace &scvf,
                const GlobalPosition       &tangentialVector) const
  {
    return 1e+5; // TODO transfer information or just use constant value
  }

  /*!
   * \brief Returns the velocity in the porous medium (which is 0 by default according to Saffman).
   */
  VelocityVector porousMediumVelocity(const FVElementGeometry    &fvGeometry,
                                      const SubControlVolumeFace &scvf) const
  {
    VelocityVector velocity(0.0);
    if constexpr (ParentType::isMomentumProblem()) {
      const auto faceId = scvf.index();
      if (couplingParticipant_.isCoupledEntity(faceId)) {
        velocity[Indices::velocityYIdx] =
            couplingParticipant_.getScalarQuantityOnFace(
                "Free-Flow-Mesh", "Velocity", faceId);
      }
    }
    return velocity;
  }

  /*!
   * \brief calculate the analytical velocity in x direction based on Beavers & Joseph (1967)
   */
  void calculateAnalyticalVelocityX() const
  {
    analyticalVelocityX_.resize(this->gridGeometry().gridView().size(0));

    using std::sqrt;
    const Scalar        dPdX = -deltaP_ / (this->gridGeometry().bBoxMax()[0] -
                                    this->gridGeometry().bBoxMin()[0]);
    static const Scalar mu   = FluidSystem::viscosity(273.15 + 10, 1e5);
    static const Scalar alpha =
        getParam<Scalar>("Darcy.SpatialParams.AlphaBeaversJoseph");
    static const Scalar K =
        getParam<Scalar>("Darcy.SpatialParams.Permeability");
    static const Scalar sqrtK = sqrt(K);
    const Scalar        sigma = (this->gridGeometry().bBoxMax()[1] -
                          this->gridGeometry().bBoxMin()[1]) /
                         sqrtK;

    const Scalar uB =
        -K / (2.0 * mu) *
        ((sigma * sigma + 2.0 * alpha * sigma) / (1.0 + alpha * sigma)) *
        dPdX;

    for (const auto &element : elements(this->gridGeometry().gridView())) {
      const auto eIdx =
          this->gridGeometry().gridView().indexSet().index(element);
      const Scalar y = element.geometry().center()[1] -
                       this->gridGeometry().bBoxMin()[1];

      const Scalar u =
          uB * (1.0 + alpha / sqrtK * y) +
          1.0 / (2.0 * mu) * (y * y + 2 * alpha * y * sqrtK) * dPdX;
      analyticalVelocityX_[eIdx] = u;
    }
  }

  /*!
   * \brief Get the analytical velocity in x direction
   */
  const std::vector<Scalar> &getAnalyticalVelocityX() const
  {
    if (analyticalVelocityX_.empty())
      calculateAnalyticalVelocityX();
    return analyticalVelocityX_;
  }

  // \}

private:
  bool onLeftBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_;
  }

  bool onRightBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_;
  }

  bool onLowerBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_;
  }

  bool onUpperBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_;
  }

  Scalar eps_;
  Scalar deltaP_;

  Dumux::Precice::CouplingAdapter &couplingParticipant_;

  mutable std::vector<Scalar> analyticalVelocityX_;
};
} // namespace Dumux

#endif // DUMUX_STOKES_SUBPROBLEM_HH
