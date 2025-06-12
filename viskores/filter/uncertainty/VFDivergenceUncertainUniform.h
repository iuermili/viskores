//============================================================================
//  The contents of this file are covered by the Viskores license. See
//  LICENSE.txt for details.
//
//  By contributing to this file, all contributors agree to the Developer
//  Certificate of Origin Version 1.1 (DCO 1.1) as stated in DCO.txt.
//============================================================================

//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//============================================================================

#ifndef viskores_filter_uncertainty_VFDivergenceUncertainUniform_h
#define viskores_filter_uncertainty_VFDivergenceUncertainUniform_h

#include <viskores/filter/Filter.h>
#include <viskores/filter/uncertainty/viskores_filter_uncertainty_export.h>

namespace viskores {
namespace filter {
namespace uncertainty {

class VISKORES_FILTER_UNCERTAINTY_EXPORT VFDivergenceUncertainUniform : public viskores::filter::Filter
{
public:
    VISKORES_CONT VFDivergenceUncertainUniform() = default;
    
}

} // namespace uncertainty
} // namespace filter
} // namespace viskores

#endif