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
      couplingParticipant_(Dumux::Precice::CouplingAdapter::getInstance())
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
        const precice::string_view meshNameView("macro-mesh", 9);
        const precice::string_view readDataPorosity("porosity", 8);
        if (getParam<bool>("Precice.RunWithCoupling") == true)
            return couplingParticipant_.getScalarQuantityOnFace(meshNameView,readDataPorosity,scv.elementIndex());
        else
            return getParam<Scalar>("Problem.DefaultPorosity"); 
    } 

    /*!
     * \brief Defines the conductivity tensor \f$ K \f$.
     */

    DimWorldMatrix solidThermalConductivity(const Element &element,
                                    const SubControlVolume& scv) const
    {   DimWorldMatrix K;
        const precice::string_view meshNameView("macro-mesh", 9);
        const precice::string_view readDatak00("k_00", 4);
        const precice::string_view readDatak01("k_01", 4);   
        const precice::string_view readDatak10("k_10", 4);   
        const precice::string_view readDatak11("k_11", 4);   
        if (getParam<bool>("Precice.RunWithCoupling") == true)
        {   
            K[0][0] = couplingParticipant_.getScalarQuantityOnFace(meshNameView,readDatak00,scv.elementIndex());
            K[0][1] = couplingParticipant_.getScalarQuantityOnFace(meshNameView,readDatak01,scv.elementIndex());
            K[1][0] = couplingParticipant_.getScalarQuantityOnFace(meshNameView,readDatak10,scv.elementIndex());
            K[1][1] = couplingParticipant_.getScalarQuantityOnFace(meshNameView,readDatak11,scv.elementIndex());
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

private:
    Dumux::Precice::CouplingAdapter &couplingParticipant_;
};

} // end namespace Dumux

#endif
