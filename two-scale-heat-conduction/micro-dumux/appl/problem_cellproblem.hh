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

#ifndef DUMUX_CELL_PROBLEM_HH
#define DUMUX_CELL_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class CellProblemProblem : public FVProblemWithSpatialParams<TypeTag> {
  using ParentType = FVProblemWithSpatialParams<TypeTag>;
  using GridView =
      typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
  using Element        = typename GridView::template Codim<0>::Entity;
  using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

  using Scalar           = GetPropType<TypeTag, Properties::Scalar>;
  using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
  using FVElementGeometry =
      typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
  using SubControlVolume     = typename FVElementGeometry::SubControlVolume;
  using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
  using GridGeometry         = GetPropType<TypeTag, Properties::GridGeometry>;
  using BoundaryTypes        = Dumux::BoundaryTypes<PrimaryVariables::size()>;
  using SolutionVector       = GetPropType<TypeTag, Properties::SolutionVector>;
  using DimWorldVector       = Dune::FieldVector<Scalar, GridView::dimensionworld>;
  using Extrusion            = Extrusion_t<GridGeometry>;
  using ModelTraits          = GetPropType<TypeTag, Properties::ModelTraits>;
  using Indices =
      typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
  static constexpr auto numEq = ModelTraits::numEq();

public:
  CellProblemProblem(std::shared_ptr<const GridGeometry> gridGeometry)
      : ParentType(gridGeometry)
  {
    // the partial derivatives of Psi
    d0psi1_.resize(gridGeometry->numDofs());
    d1psi1_.resize(gridGeometry->numDofs());
    d0psi2_.resize(gridGeometry->numDofs());
    d1psi2_.resize(gridGeometry->numDofs());

    // some memory for internal convenience
    d_.resize(gridGeometry->numDofs());
    dPsi_.resize(gridGeometry->numDofs());
    delta_ij_.resize(gridGeometry->numDofs());
  }

  BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
  {
    BoundaryTypes bcTypes;

    bcTypes.setAllNeumann();

    return bcTypes;
  }

  //! Enable internal Dirichlet constraints
  static constexpr bool enableInternalDirichletConstraints()
  {
    return true;
  }

  /*!
   * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
   *        If true is returned for a dof, the equation for this dof is replaced
   *        by the constraint that its primary variable values must match the
   *        user-defined values obtained from the function internalDirichlet(),
   *        which must be defined in the problem.
   *
   * \param element The finite element
   * \param scv The sub-control volume
   */
  std::bitset<numEq>
  hasInternalDirichletConstraint(const Element          &element,
                                 const SubControlVolume &scv) const
  {
    // the pure Neumann problem is only defined up to a constant
    // we create a well-posed problem by fixing the pressure at one dof in the
    // middle of the domain
    std::bitset<numEq> values;
    if (scv.dofIndex() ==
        static_cast<std::size_t>(this->gridGeometry().numDofs() - 1)) {
      values.set(Indices::psi1Idx);
      values.set(Indices::psi2Idx);
    }
    return values;
  }

  /*!
   * \brief Define the values of internal Dirichlet constraints for a degree of
   * freedom. \param element The finite element \param scv The sub-control
   * volume
   */
  PrimaryVariables internalDirichlet(const Element          &element,
                                     const SubControlVolume &scv) const
  {
    return PrimaryVariables(1.0);
  }

  /*!
   * \brief Calculates the upscaled conductivity tensor components by
   * integrating over the effective conductivity field \param psiIdx The
   * index/dimension of the cell problem (0 or 1) \param derivIdx The
   * index/dimension of the partial derivative (0 or 1)
   */
  Scalar calculateConductivityTensorComponent(int psiIdx, int derivIdx)
  {
    std::size_t order = 4;
    return integrateGridFunction(this->gridGeometry(),
                                 effectiveConductivityField(psiIdx, derivIdx),
                                 order);
  }

  /*!
   * \brief Calculates the effective conductivity field
   * \note  \f$ (\Phi*k_s + (1-\Phi)*kg_)*(\delta_{ij} + \partial y_i \Psi^{j})
   * \f$. \param psiIdx The index/dimension of the cell problem (0 or 1) = j
   * \param derivIdx The index/dimension of the partial derivative (0 or 1) = i
   */
  std::vector<Scalar> &effectiveConductivityField(int psiIdx, int derivIdx)
  {

    if (psiIdx == derivIdx)
      std::fill(delta_ij_.begin(), delta_ij_.end(), 1.0);
    else
      std::fill(delta_ij_.begin(), delta_ij_.end(), 0.0);

    dPsi_ = partialDerivativePsi(psiIdx, derivIdx);

    for (int elemIdx = 0; elemIdx < dPsi_.size(); ++elemIdx) {
      d_[elemIdx] = this->spatialParams().phi0deltaIdx(elemIdx) *
                    (delta_ij_[elemIdx] + dPsi_[elemIdx]);
    }
    return d_;
  }

  /*!
   * \brief Returns the previously calculated partial derivative
   */
  std::vector<Scalar> &partialDerivativePsi(int psiIdx, int derivIdx)
  {
    assert((psiIdx == 0) || (psiIdx == 1));
    assert((derivIdx == 0) || (derivIdx == 1));
    if ((psiIdx == 0) && (derivIdx == 0))
      return d0psi1_;
    else if ((psiIdx == 0) && (derivIdx == 1))
      return d1psi1_;
    else if ((psiIdx == 1) && (derivIdx == 0))
      return d0psi2_;
    else
      return d1psi2_;
  }

  /*!
   * \brief Computes the Psi Derivatives.
   */
  template <class Problem, class Assembler, class GridVariables,
            class SolutionVector>
  void computePsiDerivatives(const Problem &problem, const Assembler &assembler,
                             const GridVariables  &gridVars,
                             const SolutionVector &psi)
  {
    const auto &gridGeometry = this->gridGeometry();
    auto        fvGeometry   = localView(gridGeometry);

    const auto gridVolVars = assembler.gridVariables().curGridVolVars();
    auto       elemVolVars = localView(gridVolVars);

    for (const auto &element : elements(gridGeometry.gridView())) {

      fvGeometry.bindElement(element);
      elemVolVars.bindElement(element, fvGeometry, psi);

      for (int k = 0; k < Indices::numIdx; k++) {

        DimWorldVector cellDeriv(0.0);
        Scalar         scvVolume(0.0);

        for (const auto &scvf : scvfs(fvGeometry)) {
          if (!scvf.boundary()) {
            // Get the inside and outside volume variables
            const auto &insideScv  = fvGeometry.scv(scvf.insideScvIdx());
            const auto &outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

            const auto &insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto  valInside     = insideVolVars.priVar(k);

            const Scalar ti =
                computeTpfaTransmissibility(fvGeometry, scvf, insideScv, 1.0,
                                            insideVolVars.extrusionFactor());

            // faces might lie on the periodic boundary, requiring the matching
            // scvf of the scv on the other side of the periodic boundary.
            auto        outsideFvGeometry = localView(gridGeometry);
            const auto &periodicElement =
                gridGeometry.element(outsideScv.elementIndex());
            outsideFvGeometry.bind(periodicElement);
            auto outsideElemVolVars = localView(gridVolVars);
            outsideElemVolVars.bindElement(periodicElement, outsideFvGeometry,
                                           psi);

            Scalar tij        = 0.0;
            Scalar valOutside = 0.0;
            for (const auto &outsideScvf : scvfs(outsideFvGeometry)) {
              if (outsideScvf.unitOuterNormal() * scvf.unitOuterNormal() <
                  -1 + 1e-6) {
                const auto &outsideVolVars =
                    outsideElemVolVars[outsideScvf.insideScvIdx()];

                valOutside      = outsideVolVars.priVar(k);
                const Scalar tj = computeTpfaTransmissibility(
                    fvGeometry, outsideScvf, outsideScv, 1.0,
                    outsideVolVars.extrusionFactor());
                tij = scvf.area() * (tj) / (ti + tj);
                break;
              }
            }

            cellDeriv += scvf.area() * (valOutside - valInside) * tij *
                         scvf.unitOuterNormal();
            scvVolume = insideScv.volume();
          }
        }

        if (scvVolume > 0.0)
          cellDeriv /= scvVolume;
        const int eIdxGlobal = gridGeometry.elementMapper().index(element);

        if (k == 0) {
          d0psi1_[eIdxGlobal] = cellDeriv[0];
          d1psi1_[eIdxGlobal] = cellDeriv[1];
        } else if (k == 1) {
          d0psi2_[eIdxGlobal] = cellDeriv[0];
          d1psi2_[eIdxGlobal] = cellDeriv[1];
        }
      } // end iteration over psiIdx
    }   // end iteration over elements
  }

private:
  std::vector<Scalar> d0psi1_;
  std::vector<Scalar> d1psi1_;
  std::vector<Scalar> d0psi2_;
  std::vector<Scalar> d1psi2_;

  std::vector<Scalar> d_;
  std::vector<Scalar> delta_ij_;
  std::vector<Scalar> dPsi_;
};
} // end namespace Dumux

#endif
