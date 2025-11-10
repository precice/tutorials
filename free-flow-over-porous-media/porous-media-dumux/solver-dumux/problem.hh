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
#ifndef DUMUX_DARCY_SUBPROBLEM_HH
#define DUMUX_DARCY_SUBPROBLEM_HH

#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

#include <dumux-precice/couplingadapter.hh>

namespace Dumux {

/*!
 * \brief The porous medium flow sub problem
 */
template <class TypeTag>
class DarcySubProblem : public PorousMediumFlowProblem<TypeTag> {
  using ParentType = PorousMediumFlowProblem<TypeTag>;
  using GridView =
      typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
  using Scalar           = GetPropType<TypeTag, Properties::Scalar>;
  using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

  using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

  using BoundaryTypes = Dumux::BoundaryTypes<
      GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
  using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
  using FVElementGeometry =
      typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
  using SubControlVolume = typename FVElementGeometry::SubControlVolume;
  using SubControlVolumeFace =
      typename FVElementGeometry::SubControlVolumeFace;
  using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

  using Indices =
      typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

  using Element        = typename GridView::template Codim<0>::Entity;
  using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
  DarcySubProblem(std::shared_ptr<const GridGeometry> fvGridGeometry)
      : ParentType(fvGridGeometry, "Darcy"),
        eps_(1e-7),
        couplingParticipant_(Dumux::Precice::CouplingAdapter::getInstance())
  {
  }

  /*!
   * \brief Specifies which kind of boundary condition should be
   *        used for which equation on a given boundary control volume.
   *
   * \param element The element
   * \param scvf The boundary sub control volume face
   */
  BoundaryTypes boundaryTypes(const Element              &element,
                              const SubControlVolumeFace &scvf) const
  {
    BoundaryTypes values;

    // set Neumann BCs to all boundaries first
    values.setAllNeumann();

    const auto faceId = scvf.index();
    if (couplingParticipant_.isCoupledEntity(faceId))
      values.setAllDirichlet();
    return values;
  }

  /*!
   * \brief Evaluate the boundary conditions for a Dirichlet control volume.
   *
   * \param element The element for which the Dirichlet boundary condition is set
   * \param scvf The boundary subcontrolvolumeface
   *
   * For this method, the \a values parameter stores primary variables.
   */
  PrimaryVariables dirichlet(const Element              &element,
                             const SubControlVolumeFace &scvf) const
  {
    // set p = 0 at the bottom
    PrimaryVariables values(0.0);
    values = initial(element);

    const auto faceId = scvf.index();
    if (couplingParticipant_.isCoupledEntity(faceId))
      values = couplingParticipant_.getScalarQuantityOnFace(
          "Porous-Media-Mesh", "Pressure", faceId);

    return values;
  }

  /*!
   * \brief Evaluate the boundary conditions for a Neumann control volume.
   *
   * \param element The element for which the Neumann boundary condition is set
   * \param fvGeomentry The fvGeometry
   * \param elemVolVars The element volume variables
   * \param scvf The boundary sub control volume face
   *
   * For this method, the \a values variable stores primary variables.
   */
  template <class ElementVolumeVariables, class ElementFluxVarsCache>
  NumEqVector neumann(const Element                &element,
                      const FVElementGeometry      &fvGeometry,
                      const ElementVolumeVariables &elemVolVars,
                      const ElementFluxVarsCache   &elemFluxVarsCache,
                      const SubControlVolumeFace   &scvf) const
  {
    // no-flow everywhere ...
    NumEqVector values(0.0);
    return values;
  }

  // \}

  /*!
   * \name Volume terms
   */
  // \{
  /*!
   * \brief Evaluate the source term for all phases within a given
   *        sub-control-volume.
   *
   * \param element The element for which the source term is set
   * \param fvGeomentry The fvGeometry
   * \param elemVolVars The element volume variables
   * \param scv The subcontrolvolume
   */
  template <class ElementVolumeVariables>
  NumEqVector source(const Element                &element,
                     const FVElementGeometry      &fvGeometry,
                     const ElementVolumeVariables &elemVolVars,
                     const SubControlVolume       &scv) const
  {
    return NumEqVector(0.0);
  }

  // \}

  /*!
   * \brief Evaluate the initial value for a control volume.
   *
   * \param element The element
   *
   * For this method, the \a priVars parameter stores primary
   * variables.
   */
  PrimaryVariables initial(const Element &element) const
  {
    static const Scalar p =
        getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialP");
    return PrimaryVariables(p);
  }

  // \}

private:
  bool onLeftBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[0] < this->fvGridGeometry().bBoxMin()[0] + eps_;
  }

  bool onRightBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
  }

  bool onLowerBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + eps_;
  }

  bool onUpperBoundary_(const GlobalPosition &globalPos) const
  {
    return globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - eps_;
  }

  Scalar eps_;

  Dumux::Precice::CouplingAdapter &couplingParticipant_;
};
} // namespace Dumux

#endif // DUMUX_DARCY_SUBPROBLEM_HH
