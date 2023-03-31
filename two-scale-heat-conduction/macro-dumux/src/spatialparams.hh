// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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

#ifndef DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH
#define DUMUX_TEST_1PNI_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

#include <dumux-precice/couplingadapter.hh>

namespace Dumux {

/*!
 * \brief Definition of the spatial parameters for the 1pni problems.
 */
template<class GridGeometry, class Scalar>
class OnePNISpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         OnePNISpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using ThisType = OnePNISpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    static const int dimWorld = GridView::dimensionworld;
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePNISpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry),
      couplingInterface_(Dumux::Precice::CouplingAdapter::getInstance())
     {}

    /*!
     * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param globalPos The global position
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {return getParam<Scalar>("Problem.Permeability"); } 

    /*!
     * \brief Defines the porosity \f$\mathrm{[-]}\f$.
     *
     * \param globalPos The global position
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element, 
                    const SubControlVolume& scv, 
                    const ElementSolution& elemSol) const
    {   
        if (getParam<bool>("Precice.RunWithCoupling") == true)
            return couplingInterface_.getScalarQuantityOnFace(porosityId_,scv.elementIndex());
        else
            return getParam<Scalar>("Problem.DefaultPorosity"); 
    } 

    /*!
     * \brief Defines the conductivity tensor \f$ K \f$.
     */

    DimWorldMatrix solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv) const
    {   DimWorldMatrix K;
        if (getParam<bool>("Precice.RunWithCoupling") == true)
        {   
            K[0][0] = couplingInterface_.getScalarQuantityOnFace(k00Id_,scv.elementIndex());
            K[0][1] = couplingInterface_.getScalarQuantityOnFace(k01Id_,scv.elementIndex());
            K[1][0] = couplingInterface_.getScalarQuantityOnFace(k10Id_,scv.elementIndex());
            K[1][1] = couplingInterface_.getScalarQuantityOnFace(k11Id_,scv.elementIndex());
        } 
        else
        { 
            K[0][0] = getParam<Scalar>("Component.SolidThermalConductivity");
            K[0][1] = 0.0;
            K[1][0] = 0.0;
            K[1][1] = getParam<Scalar>("Component.SolidThermalConductivity");
        }
        return K; 
    }

    void updatePreciceDataIds() 
    {   
        porosityId_ = couplingInterface_.getIdFromName("porosity");
        k00Id_ = couplingInterface_.getIdFromName("k_00");
        k01Id_ = couplingInterface_.getIdFromName("k_01");
        k10Id_ = couplingInterface_.getIdFromName("k_10");
        k11Id_ = couplingInterface_.getIdFromName("k_11");
    }

private:
    Dumux::Precice::CouplingAdapter &couplingInterface_;
    int porosityId_; 
    int k00Id_; 
    int k01Id_; 
    int k10Id_; 
    int k11Id_; 
};

} // end namespace Dumux

#endif
