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
#include <viskores/Math.h>
#include <viskores/Types.h>
#include <viskores/cont/ArrayHandle.h>
#include <viskores/cont/ArrayHandleView.h>
#include <viskores/cont/DataSet.h>
#include<viskores/cont/Field.h>
#include <viskores/cont/Initialize.h>
#include <viskores/cont/Invoker.h>
#include <viskores/filter/Filter.h>
#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>
#include <viskores/worklet/WorkletMapField.h>
#include <viskores/worklet/WorkletPointNeighborhood.h>
#include <viskores/worklet/WorkletMapTopology.h>
#include <viskores/CellShape.h>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

// Forward Declaration for Analytical and Sampling Approaches
namespace viskores
{
namespace filter
{
namespace uncertainty
{

class VFDivergenceAnalytical;

namespace sampling
{
class VFDivergenceSampling;
} // namespace sampling

} // namespace uncertainty
} // namespace filter
} // namespace viskores

// Filter Declaration for Analytical Approach
namespace viskores
{
namespace filter
{
namespace uncertainty
{

class VFDivergenceAnalytical : public viskores::filter::Filter
{
public:
    VFDivergenceAnalytical();

    void SetIsovalue(viskores::Float64 isovalue)
    { 
        this->Isovalue = isovalue;
    }

    private:
        viskores::cont::DataSet DoExecute(const viskores::cont::DataSet& input) override;
        viskores::Float64 Isovalue;
};
        
// Analytical Approach
namespace
{

// Worklet 1: computes divergence mean and variance for a 2D vector field.
struct ComputeDivergence : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood meanX,
                                    FieldInNeighborhood varX,
                                    FieldInNeighborhood meanY,
                                    FieldInNeighborhood varY,
                                    FieldOut divMean,
                                    FieldOut divVariance);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5, _6, _7);
    using InputDomain = _1;

    template <typename BoundaryType, typename NeighborhoodType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                    const NeighborhoodType& uMean,
                                    const NeighborhoodType& uVar,
                                    const NeighborhoodType& vMean,
                                    const NeighborhoodType& vVar,
                                    OutType& divergenceMean,
                                    OutType& divergenceVariance) const
    {
        OutType divUMean, varSquaredU;
        OutType divVMean, varSquaredV;

        if (boundary.MinNeighborIndices(1)[1] == 0)
        {
            divUMean = uMean.Get(0, 1, 0) - uMean.Get(0, 0, 0);
            varSquaredU = uVar.Get(0, 1, 0) + uVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[1] == 0) {
            divUMean = uMean.Get(0, 0, 0) - uMean.Get(0, -1, 0);
            varSquaredU = uVar.Get(0, 0, 0) + uVar.Get(0, -1, 0);
        }
        else {
            divUMean = (uMean.Get(0, 1, 0) - uMean.Get(0, -1, 0)) / 2.0;
            varSquaredU = uVar.Get(0, 1, 0) + uVar.Get(0, -1, 0);
        }

        if (boundary.MinNeighborIndices(1)[0] == 0)
        {
            divVMean = vMean.Get(1, 0, 0) - vMean.Get(0, 0, 0);
            varSquaredV = vVar.Get(1, 0, 0) + vVar.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[0] == 0)
        {
            divVMean = vMean.Get(0, 0, 0) - vMean.Get(-1, 0, 0);
            varSquaredV = vVar.Get(0, 0, 0) + vVar.Get(-1, 0, 0);
        }
        else
        {
            divVMean = (vMean.Get(1, 0, 0) - vMean.Get(-1, 0, 0)) / 2.0;
            varSquaredV = vVar.Get(1, 0, 0) + vVar.Get(-1, 0, 0);
        }

        divergenceMean = divUMean + divVMean;
        OutType totalVariance = varSquaredU + varSquaredV;
        divergenceVariance = viskores::Sqrt(viskores::Max(totalVariance, static_cast<OutType>(0.0)));
    }
};

// Worklet 2 Helper Function: computes Gaussian CDF.
template <typename T>
VISKORES_EXEC T GaussianCDF(T x, T mu, T sigma)
{
    return static_cast<T>(0.5) * (static_cast<T>(1.0) + erf((x - mu) / (sigma * viskores::Sqrt(2.0))));
}

// Worklet 2: computes crossing probability for a 2D vector field based on divergence mean and variance.
struct CrossingProbability : public viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInPoint divMean,
                                    FieldInPoint divVariance,
                                    FieldOutCell crossingProb);
    using ExecutionSignature = _4(_2, _3);
    using InputDomain = _1;

    viskores::Float64 Isovalue;

    template <typename MeanVecType, typename VarianceVecType>
    VISKORES_EXEC viskores::Float64 operator()(const MeanVecType& divMeanAtCorners,
                                    const VarianceVecType& divVarianceAtCorners) const
    {
        viskores::Float64 probNeg = 1.0;
        viskores::Float64 probPos = 1.0;

        for (viskores::IdComponent i = 0; i < divMeanAtCorners.GetNumberOfComponents(); ++i)
        {
            auto mean = static_cast<viskores::Float64>(divMeanAtCorners[i]);
            auto variance = static_cast<viskores::Float64>(divVarianceAtCorners[i]);
            
            auto cdfValue = viskores::Sqrt(viskores::Max(variance, static_cast<viskores::Float64>(0.0)));

            auto probBelow = GaussianCDF(this->Isovalue, mean, cdfValue);
            auto probAbove = 1.0 - probBelow;

            probNeg *= probBelow;
            probPos *= probAbove;
        }

        return 1.0 - probPos - probNeg;
    }
};

} // namespace

// Filter Implementation
VISKORES_CONT
VFDivergenceAnalytical::VFDivergenceAnalytical()
{
    this->Isovalue = 0.0;
}

VISKORES_CONT viskores::cont::DataSet
VFDivergenceAnalytical::DoExecute(const viskores::cont::DataSet& input)
{
    const auto& meanXField = input.GetField("meanX");
    const auto& varXField = input.GetField("varX");
    const auto& meanYField = input.GetField("meanY");
    const auto& varYField = input.GetField("varY");

    viskores::cont::ArrayHandle<viskores::FloatDefault> divMeanHandle;
    viskores::cont::ArrayHandle<viskores::FloatDefault> divVarianceHandle;
    viskores::cont::ArrayHandle<viskores::FloatDefault> crossingProbabilityHandle;

    auto resolveType = [&](const auto& concreteMeanX)
    {
        using ArrayType = std::decay_t<decltype(concreteMeanX)>;

        ArrayType concreteVarX, concreteMeanY, concreteVarY;
        varXField.GetData().AsArrayHandle(concreteVarX);
        meanYField.GetData().AsArrayHandle(concreteMeanY);
        varYField.GetData().AsArrayHandle(concreteVarY);

        this->Invoke(ComputeDivergence{},
                     input.GetCellSet(),
                     concreteMeanX,
                     concreteVarX,
                     concreteMeanY,
                     concreteVarY,
                     divMeanHandle,
                     divVarianceHandle);

        CrossingProbability crossingProbWorklet;
        crossingProbWorklet.Isovalue = this->Isovalue;
        this->Invoke(crossingProbWorklet,
                     input.GetCellSet(),
                     divMeanHandle,
                     divVarianceHandle,
                     crossingProbabilityHandle);
    };

    this->CastAndCallScalarField(meanXField, resolveType);

    viskores::cont::DataSet result;
    result.AddCoordinateSystem(input.GetCoordinateSystem());
    result.SetCellSet(input.GetCellSet());

    result.AddPointField("divergenceMean_analytical", divMeanHandle);
    result.AddPointField("divergenceVariance_analytical", divVarianceHandle);
    result.AddCellField("crossingProbability_analytical", crossingProbabilityHandle);

    return result;
}

} // namespace uncertainty
} // namespace filter
} // namespace viskores

// Filter Declaration for Sampling Approach
namespace viskores
{
namespace filter
{
namespace uncertainty
{
namespace sampling
{

class VFDivergenceSampling : public viskores::filter::Filter
{
public:
    VFDivergenceSampling();
    void SetIsovalue(viskores::Float64 isovalue) { this->Isovalue = isovalue; }
    void SetNumSamples(viskores::Id numSamples) { this->NumSamples = numSamples; }
    
    private:
        viskores::cont::DataSet DoExecute(const viskores::cont::DataSet& input) override;
        viskores::Float64 Isovalue;
        viskores::Id NumSamples;
};

// Sampling Approach
namespace
{

// Worklet 2: computes divergence mean and variance of a 2D vector field using sampling.
struct ComputeDivergenceSampling: public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood meanX,
                                    FieldInNeighborhood varX,
                                    FieldInNeighborhood meanY,
                                    FieldInNeighborhood varY,
                                    FieldOut divMean,
                                    FieldOut divVariance);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5, _6, _7);
    using InputDomain = _1;

    template <typename BoundaryType, typename NeighborhoodType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                        const NeighborhoodType& uMean,
                                        const NeighborhoodType& uVar,
                                        const NeighborhoodType& vMean,
                                        const NeighborhoodType& vVar,
                                        OutType& divergenceMean,
                                        outType& divergenceVariance) const
    {
        OutType divUMean, varSquaredU;
        OutTyp divVMean, varSquaredV;

        if (boundary.MinNeighborIndices(1)[1] == 0) {
            divUMean = uMean.Get(0, 1, 0) - uMean.Get(0, 0, 0);
            varSquaredU = uVar.Get(0, 1, 0) + uVar.Get(0, 0, 0);
        } else if (boundary.MaxNeighborIndices(1)[1] == 0) {
            divUMean = uMean.Get(0, 0, 0) - uMean.Get(0, -1, 0);
            varSquaredU = uVar.Get(0, 0, 0) + uVar.Get(0, -1, 0);
        } else {
            divUMean = (uMean.Get(0, 1, 0) - uMean.Get(0, -1, 0)) / 2.0;
            varSquaredU = (uVar.Get(0, 1, 0) + uVar.Get(0, -1, 0))
        }

        if (boundary.MinNeighborIndices(1)[0] == 0) {
            divVMean = vMean.Get(1, 0, 0) - vMean.Get(0, 0, 0);
            varSquaredV = vVar.Get(1, 0, 0) + vVar.Get(0, 0, 0);
        } else if (boundary.MaxNeighborIndices(1)[0] == 0) {
            divVMean = vMean.Get(0, 0, 0) - vMean.Get(-1, 0, 0);
            varSquaredV = vVar.Get(0, 0, 0) + vVar.Get(-1, 0, 0);
        } else {
            divVMean = (vMean.Get(1, 0, 0) - vMean.Get(-1, 0, 0)) / 2.0;
            varSquaredV = (vVar.Get(1, 0, 0) + vVar.Get(-1, 0, 0)) / 4.0; // Variance of (A-B)/2 is (Var(A)+Var(B))/4
        }
        
        divergenceMean = divUMean + divVMean;
        // The final variance is the sum of the individual variances
        divergenceVariance = varSquaredU + varSquaredV;
    }
};


// Worklet 3: computes crossing probability for a 2D vector field based on divergence mean and variance using sampling.
struct CrossingProbabilitySampling : public viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInPoint divMean,
                                    FieldInPoint divVar,
                                    FieldInPoint, seeds
                                    FieldOutCell crossingProb);
    using ExecutionSignature = _5(CellShape, _2, _3, _4);
    using InputDomain = _1;

    viskores::Float64 Isovalue;
    viskores::Id NumSamples;

    template <typename MeanVecType, typename VarVecType, typename SeedVecType>
    VISKORES_EXEC viskores::Float64 operator()(viskores::CellShapeTagQuad,
                                                const MeanVecType& divMeanAtCorners,
                                                const VarVecType& divVarAtCorners,
                                                const SeedVecType& seedsAtCorners) const
    {
        viskores::Id crossings = 0;

        for (viskores::Id sampleId = 0; sampleId < this->NumSamples; ++sampleId)
        {
            boost::random::mt19937 rng0(12345)(seedsAtCorners[0] + sampleId);
            boost::random::normal_distribution<viskores::Float64> dist0(divMeanAtCorners[0], varMeanAtCorners[0]);
            viskores::FloatDefault v0 = dist0(rng0);

            boost::random::mt19937 rng1(12345)(seedsAtCorners[1] + sampleId);
            boost::random::normal_distribution<viskores::Float64> dist1(divMeanAtCorners[1], varMeanAtCorners[1]);
            viskores::FloatDefault v1 = dist1(rng1);

            boost::random::mt19937 rng2(12345)(seedsAtCorners[2] + sampleId);
            boost::random::normal_distribution<viskores::Float64> dist2(divMeanAtCorners[2], varMeanAtCorners[2]);
            viskores::FloatDefault v2 = dist2(rng2);
            
            boost::random::mt19937 rng3(12345)(seedsAtCorners[3] + sampleId);
            boost::random::normal_distribution<viskores::Float64> dist3(divMeanAtCorners[3], varMeanAtCorners[3]);
            viskores::FloatDefault v3 = dist3(rng3);

            viskores::Float64 minVal = viskores::Min(viskores::Min(v0, v1), viskores::Min(v2, v3));
            viskores::Float64 maxVal = viskores::Max(viskores::Max(v0, v1), viskores::Max(v2, v3));

            if ((this->Isovalue > minVal) && (this->Isovalue < maxVal))
            {
                crossings++;
            }
        }

        return static_cast<viskores::Float64>(crossings) / this->NumSamples;
    }
};

} // namespace

// Filter Implementation
VISKORES_CONT
VFDivergenceSampling::VFDivergenceSampling()
{
    this->Isovalue = 0.0;
    this->NumSamples = 1000;
}

viskores::cont::DataSet VFDivergenceSampling::DoExecute(const viskores::cont::DataSet& input)
{
    const auto& meanXField = input.GetField("meanX").GetData();
    const auto& varXField = input.GetField("varX").GetData();
    const auto& meanYField = input.GetField("meanY").GetData();
    const auto& varYField = input.GetField("varY").GetData();

    viskores::Id numPoints = input.GetCoordinateSystem().GetNumberOfPoints();
    
    viskores::cont::ArrayHandle<viskores::FloatDefault> divMeanHandle;
    viskores::cont::ArrayHandle<viskores::FloatDefault> divVarianceHandle;
    viskores::cont::ArrayHandle<viskores::Float64> crossingProbabilityHandle; 

    auto resolveType = [&](const auto& concreteMeanX)
    {
        using ArrayType = std::decay_t<decltype(concreteMeanX)>;

        ArrayType concreteVarX, concreteMeanY, concreteVarY;
        varXField.GetData().AsArrayHandle(concreteVarX);
        meanYField.GetData().AsArrayHandle(concreteMeanY);
        varYField.GetData().AsArrayHandle(concreteVarY);

        this->Invoke(ComputeDivergenceSampling{},
                        input.GetCellSet(),
                        concreteMeanX,
                        concreteVarX,
                        concreteMeanY,
                        concreteVarY,
                        uSamples,
                        vSamples);

        viskores::cont::ArrayHandle<viskores::UInt32> seeds;
        seeds.Allocate(numPoints);
        auto seedsPortal = seeds.WritePortal();
        boost::random::mt19937 rng(12345);
        for (viskores::Id i = 0; i < numPoints; ++i)
        {
            seedsPortal.Set(i, rng());
        }

        CrossingProbabilitySampling crossingProbWorklet;
        crossingProbWorklet.Isovalue = this->Isovalue;
        crossingProbWorklet.NumSamples = this->NumSamples;
        this->Invoke(crossingProbWorklet,
                        input.GetCellSet(),
                        divMeanHandle,
                        divVarianceHandle,
                        seeds,
                        crossingProbabilityHandle);
    };

    this->CastAndCallScalarField(meanXField, resolveType);

    viskores::cont::DataSet result;
    result.AddCellField("crossingProbability_sampling", crossingProbabilityHandle);

    return result;
}

} // namespace sampling
} // namespace uncertainty
} // namespace filter
} // namespace viskores

// Main Function: reads VTK dataset, creates instance of filter, sets filter parameters, executes filter, and writes output to file.
int main(int argc, char* argv[])
{
    viskores::cont::Initialize(argc, argv);

    // reads VTK dataset
    viskores::io::VTKDataSetReader reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    // creates instance of analytical filter
    viskores::filter::uncertainty::VFDivergenceAnalytical analyticalFilter;

    // sets filter parameters
    analyticalFilter.SetIsovalue(-5.05);

    // executes filter
    viskores::cont::DataSet analyticalResult = analyticalFilter.Execute(ds);

    // creates instance of sampling filter
    viskores::filter::uncertainty::sampling::VFDivergenceSampling samplingFilter;

    // sets filter parameters
    samplingFilter.SetIsovalue(-5.05);
    samplingFilter.SetNumSamples(1000);

    // executes filter
    viskores::cont::DataSet samplingResult = samplingFilter.Execute(ds);

    // writes output to file
    viskores::io::VTKDataSetWriter analyticalWriter("out_vf_divergence_analytical.vtk");
    analyticalWriter.WriteDataSet(analyticalResult);

    viskores::io::VTKDataSetWriter samplingWriter("out_vf_divergence_sampling.vtk");
    samplingWriter.WriteDataSet(samplingResult);

    return 0;
}
