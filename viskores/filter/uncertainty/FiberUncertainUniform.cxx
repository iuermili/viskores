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

//  This code is based on the algorithm presented in the following papers:
//  Hari, G., Joshi, N., Wang, Z., Gong, Q., Pugmire, D., Moreland, K.,
//  Johnson, C. R., Klasky, S., Podhorszki, N., & Athawale, T.
//  (2024). FunM^2C: A Filter for Uncertainty Visualization of MultiVariate Data
//  on Multi-Core Devices. Oak Ridge National Laboratory (ORNL),
//  Oak Ridge, TN (United States).

#include <viskores/cont/DataSet.h>
#include <viskores/cont/Timer.h>
#include <viskores/filter/uncertainty/FiberUncertainUniform.h>
#include <viskores/filter/uncertainty/worklet/FiberUncertainUniform.h>


namespace viskores
{
namespace filter
{
namespace uncertainty
{
VISKORES_CONT viskores::cont::DataSet FiberUncertainUniform::DoExecute(const viskores::cont::DataSet& input)
{
  std::string FieldName;

  
  viskores::cont::Field EnsembleMinX = this->GetFieldFromDataSet(0, input);
  viskores::cont::Field EnsembleMaxX = this->GetFieldFromDataSet(1, input);
  viskores::cont::Field EnsembleMinY = this->GetFieldFromDataSet(2, input);
  viskores::cont::Field EnsembleMaxY = this->GetFieldFromDataSet(3, input);

  // Output Field
  viskores::cont::UnknownArrayHandle OutputProbability;


  //  For Invoker
  auto resolveType = [&](auto ConcreteEnsembleMinX) {
    //  Obtaining Type
    using ArrayType = std::decay_t<decltype(ConcreteEnsembleMinX)>;
    using ValueType = typename ArrayType::ValueType;

    ArrayType ConcreteEnsembleMaxX;
    ArrayType ConcreteEnsembleMinY;
    ArrayType ConcreteEnsembleMaxY;

    viskores::cont::ArrayCopyShallowIfPossible(EnsembleMaxX.GetData(), ConcreteEnsembleMaxX);
    viskores::cont::ArrayCopyShallowIfPossible(EnsembleMinY.GetData(), ConcreteEnsembleMinY);
    viskores::cont::ArrayCopyShallowIfPossible(EnsembleMaxY.GetData(), ConcreteEnsembleMaxY);

    viskores::cont::ArrayHandle<ValueType> Probability;
    // Invoker

    if (this->Approach == ApproachEnum::MonteCarlo)
    {
      FieldName = "MonteCarlo";
      VISKORES_LOG_S(viskores::cont::LogLevel::Info, "Adopt Monte Carlo with numsamples " << this->NumSamples);
      this->Invoke(viskores::worklet::detail::MultiVariateMonteCarlo{ this->minAxis,
                                                                  this->maxAxis,
                                                                  this->NumSamples },
                   ConcreteEnsembleMinX,
                   ConcreteEnsembleMaxX,
                   ConcreteEnsembleMinY,
                   ConcreteEnsembleMaxY,
                   Probability);
    }
    else if (this->Approach == ApproachEnum::ClosedForm)
    {
      FieldName = "ClosedForm";
      std::cout << "Adopt ClosedForm" << std::endl;
      this->Invoke(viskores::worklet::detail::MultiVariateClosedForm{ this->minAxis, this->maxAxis },
                   ConcreteEnsembleMinX,
                   ConcreteEnsembleMaxX,
                   ConcreteEnsembleMinY,
                   ConcreteEnsembleMaxY,
                   Probability);
    }
    else if (this->Approach == ApproachEnum::Mean)
    {
      FieldName = "Mean";
      std::cout << "Adopt Mean" << std::endl;
      this->Invoke(viskores::worklet::detail::MultiVariateMean{ this->minAxis, this->maxAxis },
                   ConcreteEnsembleMinX,
                   ConcreteEnsembleMaxX,
                   ConcreteEnsembleMinY,
                   ConcreteEnsembleMaxY,
                   Probability);
    }
    else if (this->Approach == ApproachEnum::Truth)
    {
      FieldName = "Truth";
      std::cout << "Adopt Truth" << std::endl;
      this->Invoke(viskores::worklet::detail::MultiVariateTruth{ this->minAxis, this->maxAxis },
                   ConcreteEnsembleMinX,
                   ConcreteEnsembleMaxX,
                   ConcreteEnsembleMinY,
                   ConcreteEnsembleMaxY,
                   Probability);
    }
    else
    {
      throw std::runtime_error(std::string("Unsupported approach."));
    }


    OutputProbability = Probability;
  };
  this->CastAndCallScalarField(EnsembleMinX, resolveType);


  viskores::cont::DataSet result = this->CreateResult(input);
  result.AddPointField(FieldName, OutputProbability);


  return result;
}
}
}
}
