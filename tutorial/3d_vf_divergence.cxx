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

#include <iostream>
#include <chrono>
#include <random>
#include <cmath>

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

// Analytical Approach
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
        
namespace
{

// Worklet 1: computes divergence mean and variance for a 2D vector field.
struct ComputeDivergenceMeanAndVar : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood meanX,
                                    FieldInNeighborhood varX,
                                    FieldInNeighborhood meanY,
                                    FieldInNeighborhood varY,
                                    FieldInNeighborhood meanZ,
                                    FieldInNeighborhood varZ,
                                    FieldOut divMean,
                                    FieldOut divVar);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5, _6, _7, _8, _9);
    using InputDomain = _1;

    template <typename BoundaryType, typename NeighborhoodType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                    const NeighborhoodType& meanX,
                                    const NeighborhoodType& varX,
                                    const NeighborhoodType& meanY,
                                    const NeighborhoodType& varY,
                                    const NeighborhoodType& meanZ,
                                    const NeighborhoodType& varZ,
                                    OutType& divMean,
                                    OutType& divVar) const
    {
        OutType divmeanX, varsquaredX;
        OutType divmeanY, varsquaredY;
        OutType divmeanZ, varsquaredZ;

        if (boundary.MinNeighborIndices(1)[0] == 0)
        {
            divmeanX = meanX.Get(0, 1, 0) - meanX.Get(0, 0, 0);
            varsquaredX = varX.Get(0, 1, 0) + varX.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[0] == 0) {
            divmeanX = meanX.Get(0, 0, 0) - meanX.Get(0, -1, 0);
            varsquaredX = varX.Get(0, 0, 0) + varX.Get(0, -1, 0);
        }
        else {
            divmeanX = (meanX.Get(0, 1, 0) - meanX.Get(0, -1, 0)) / 2.0;
            varsquaredX = (varX.Get(0, 1, 0) + varX.Get(0, -1, 0)) / 4.0;
        }

        if (boundary.MinNeighborIndices(1)[1] == 0)
        {
            divmeanY = meanY.Get(1, 0, 0) - meanY.Get(0, 0, 0);
            varsquaredY = varY.Get(1, 0, 0) + varY.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[1] == 0)
        {
            divmeanY = meanY.Get(0, 0, 0) - meanY.Get(-1, 0, 0);
            varsquaredY = varY.Get(0, 0, 0) + varY.Get(-1, 0, 0);
        }
        else
        {
            divmeanY = (meanY.Get(1, 0, 0) - meanY.Get(-1, 0, 0)) / 2.0;
            varsquaredY = (varY.Get(1, 0, 0) + varY.Get(-1, 0, 0)) / 4.0;
        }

        if (boundary.MinNeighborIndices(1)[2] == 0)
        {
            divmeanZ = meanZ.Get(0, 0, 1) - meanZ.Get(0, 0, 0);
            varsquaredZ = varZ.Get(0, 0, 1) + varZ.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[2] == 0)
        {
            divmeanZ = meanZ.Get(0, 0, 0) - meanZ.Get(0, 0, -1);
            varsquaredZ = varZ.Get(0, 0, 0) + varZ.Get(0, 0, -1);
        }
        else
        {
            divmeanZ = (meanZ.Get(0, 0, 1) - meanZ.Get(0, 0, -1)) / 2.0;
            varsquaredZ = (varZ.Get(0, 0, 1) + varZ.Get(0, 0, -1)) / 4.0;
        }

        divMean = divmeanX + divmeanY + divmeanZ;
        divVar = varsquaredX + varsquaredY + varsquaredZ;
    }
};

// Worklet 2 Helper Function: computes Gaussian CDF.
template <typename T>
VISKORES_EXEC T GaussianCDF(T x, T mu, T sigma)
{
    if (sigma <= static_cast<T>(0.0)){
        if (x < mu) {
            return static_cast<T>(0.0);
        }
        else {
            return static_cast<T>(1.0);
        }
    }
    
    return static_cast<T>(0.5) * (static_cast<T>(1.0) + erf((x - mu) / (sigma * viskores::Sqrt(2.0))));
}

// Worklet 2: computes crossing probability for a 2D vector field based on divergence mean and variance.
struct CrossingProbabilityAnalytical : public viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInPoint divMean,
                                    FieldInPoint divVar,
                                    FieldOutCell crossingProb);
    using ExecutionSignature = _4(_2, _3);
    using InputDomain = _1;

    viskores::Float64 Isovalue;

    template <typename MeanVecType, typename VarVecType>
    VISKORES_EXEC viskores::Float64 operator()(const MeanVecType& divMean,
                                                const VarVecType& divVar) const
    {
        viskores::Float64 probNeg = 1.0;
        viskores::Float64 probPos = 1.0;

        for (viskores::IdComponent i = 0; i < divMean.GetNumberOfComponents(); ++i)
        {
            auto mean = static_cast<viskores::Float64>(divMean[i]);
            auto variance = static_cast<viskores::Float64>(divVar[i]);

            auto probBelow = GaussianCDF(this->Isovalue, mean, viskores::Sqrt(viskores::Max(variance, static_cast<viskores::Float64>(0.0))));
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
    const auto& meanX_Field = input.GetField("meanX");
    const auto& varX_Field = input.GetField("varX");
    const auto& meanY_Field = input.GetField("meanY");
    const auto& varY_Field = input.GetField("varY");
    const auto& meanZ_Field = input.GetField("meanZ");
    const auto& varZ_Field = input.GetField("varZ");

    viskores::cont::ArrayHandle<viskores::Float64> divMean_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> divVar_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> crossingProb_Handle;

    auto resolveType = [&](const auto& concrete_meanX)
    {
        using ArrayType = std::decay_t<decltype(concrete_meanX)>;

        ArrayType concrete_varX, concrete_meanY, concrete_varY, concrete_meanZ, concrete_varZ;
        varX_Field.GetData().AsArrayHandle(concrete_varX);
        meanY_Field.GetData().AsArrayHandle(concrete_meanY);
        varY_Field.GetData().AsArrayHandle(concrete_varY);
        meanZ_Field.GetData().AsArrayHandle(concrete_meanZ);
        varZ_Field.GetData().AsArrayHandle(concrete_varZ);

        this->Invoke(ComputeDivergenceMeanAndVar{},
                     input.GetCellSet(),
                     concrete_meanX,
                     concrete_varX,
                     concrete_meanY,
                     concrete_varY,
                     concrete_meanZ,
                     concrete_varZ,
                     divMean_Handle,
                     divVar_Handle);

        CrossingProbabilityAnalytical crossingProb_Worklet;
        crossingProb_Worklet.Isovalue = this->Isovalue;
        this->Invoke(crossingProb_Worklet,
                     input.GetCellSet(),
                     divMean_Handle,
                     divVar_Handle,
                     crossingProb_Handle);
    };

    this->CastAndCallScalarField(meanX_Field, resolveType);

    viskores::cont::DataSet result;
    result.AddCoordinateSystem(input.GetCoordinateSystem());
    result.SetCellSet(input.GetCellSet());

    result.AddPointField("divergenceMean_analytical", divMean_Handle);
    result.AddPointField("divergenceVariance_analytical", divVar_Handle);
    result.AddCellField("crossingProbability_analytical", crossingProb_Handle);

    return result;
}

} // namespace uncertainty
} // namespace filter
} // namespace viskores

// Sampling Approach
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
    void SetIsovalue(viskores::Float64 isovalue)
    { 
        this->Isovalue = isovalue; 
    }
    
    private:
        viskores::cont::DataSet DoExecute(const viskores::cont::DataSet& input) override;
        viskores::Float64 Isovalue;
};

namespace
{

constexpr viskores::Id NUM_SAMPLES = 1000;
using SampleVec = viskores::Vec<viskores::Float64, NUM_SAMPLES>;

// Worklet 1: sampling.
struct Sampling: public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood meanX,
                                    FieldInNeighborhood varX,
                                    FieldInNeighborhood meanY,
                                    FieldInNeighborhood varY,
                                    FieldInNeighborhood meanZ,
                                    FieldInNeighborhood varZ,
                                    FieldOut uSamples,
                                    FieldOut vSamples,
                                    FieldOut wSamples);
    using ExecutionSignature = void(_2, _3, _4, _5, _6, _7, _8, _9, _10);
    using InputDomain = _1;

    viskores::UInt32 Seed;

    template <typename MeanXType, typename VarXType, typename MeanYType, typename VarYType>
    VISKORES_EXEC void operator()(const MeanXType& meanX,
                                    const VarXType& varX,
                                    const MeanYType& meanY,
                                    const VarYType& varY,
                                    const VarXType& varZ,
                                    const MeanYType& meanZ,
                                    SampleVec& uSamples,
                                    SampleVec& vSamples,
                                    SampleVec& wSamples) const 
    {
        std::size_t hash = std::hash<double>()(meanX.Get(0,0,0)) ^ std::hash<double>()(meanY.Get(0,0,0)) ^ std::hash<double>(meanZ.Get(0,0,0,));
        std::mt19937 rng(Seed ^ hash);
        std::normal_distribution<double> distU(meanX.Get(0,0,0), std::sqrt(varX.Get(0,0,0)));
        std::normal_distribution<double> distV(meanY.Get(0,0,0), std::sqrt(varY.Get(0,0,0)));
        std::normal_distribution<double> distW(meanZ.Get(0,0,0), std::sqrt(varZ.Get(0,0,0)));
        for (viskores::Id i = 0; i < NUM_SAMPLES; ++i)
        {
            uSamples[i] = distU(rng);
            vSamples[i] = distV(rng);
            wSamples[i] = distW(rng);
        }
    }
};

// Worklet 2: computes divergence for a 2d vector field via sampling.
struct ComputeDivergence : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInNeighborhood uSamples,
                                    FieldInNeighborhood vSamples,
                                    FieldInNeighborhood wSamples,
                                    FieldOut divergenceSamples);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5);
    using InputDomain = _1;

    viskores::Id NumSamples;

    template <typename BoundaryType, typename USamplesType, typename VSamplesType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                    const USamplesType& uSamples,
                                    const VSamplesType& vSamples,
                                    const USamplesType& wSamples,
                                    OutType& divSamples) const
    {
        for (viskores::Id k = 0; k < NUM_SAMPLES; ++k)
        {
            double dudx, dvdy, dwdz;

            if (boundary.MinNeighborIndices(1)[1] == 0) {
                dudx = uSamples.Get(0, 1, 0)[k] - uSamples.Get(0, 0, 0)[k];
            } 
            else if (boundary.MaxNeighborIndices(1)[1] == 0) {
                dudx = uSamples.Get(0, 0, 0)[k] - uSamples.Get(0, -1, 0)[k];
            } 
            else {
                dudx = (uSamples.Get(0, 1, 0)[k] - uSamples.Get(0, -1, 0)[k]) / 2.0;
            }

            if (boundary.MinNeighborIndices(1)[0] == 0) {
                dvdy = vSamples.Get(1, 0, 0)[k] - vSamples.Get(0, 0, 0)[k];
            } 
            else if (boundary.MaxNeighborIndices(1)[0] == 0) {
                dvdy = vSamples.Get(0, 0, 0)[k] - vSamples.Get(-1, 0, 0)[k];
            } 
            else {
                dvdy = (vSamples.Get(1, 0, 0)[k] - vSamples.Get(-1, 0, 0)[k]) / 2.0;
            }

            if (boundary.MinNeighborIndices(1)[2] == 0) {
                dwdz = wSamples.Get(0, 0, 1)[k] - wSamples.Get(0, 0, 0)[k];
            } 
            else if (boundary.MaxNeighborIndices(1)[2] == 0) {
                dwdz = wSamples.Get(0, 0, 0)[k] - wSamples.Get(0, 0, -1)[k];
            } 
            else {
                dwdz = (wSamples.Get(0, 0, 1)[k] - wSamples.Get(0, 0, -1)[k]) / 2.0;
            }

            divSamples[k] = dudx + dvdy + dwdz;
        }
    }
};

// Worklet 3: computes crossing probability for a 2D vector field based on divergence using sampling.
struct CrossingProbabilitySampling : public viskores::worklet::WorkletVisitCellsWithPoints
{
    using ControlSignature = void(CellSetIn domain,
                                    FieldInPoint divSamples,
                                    FieldOutCell crossingProb);
    using ExecutionSignature = _3(CellShape, _2);
    using InputDomain = _1;

    viskores::Float64 Isovalue;
    viskores::Id NumSamples;

    template <typename CellShapeTag, typename DivSamplesType>
    VISKORES_EXEC viskores::Float64 operator()(CellShapeTag,
                                                const DivSamplesType& divSamples) const
    {
        viskores::Id crossings = 0;

        for (viskores::Id k = 0; k < NUM_SAMPLES; ++k)
        {
            const int numCorners = divSamples.GetNumberOfComponents();
            double minVal = divSamples[0][k];
            double maxVal = minVal;

            for (int c = 1; c < numCorners; ++c)
            {
                double val = divSamples[c][k];
                if (val < minVal) {
                    minVal = val;
                }
                if (val > maxVal) {
                    maxVal = val;
                }
            }

            if (this->Isovalue > minVal && this->Isovalue < maxVal){
                crossings++;
            }
        }

        return static_cast<double>(crossings) / NUM_SAMPLES;
    }
};

} // namespace

// Filter Implementation
VISKORES_CONT
VFDivergenceSampling::VFDivergenceSampling()
{
    this->Isovalue = 0.0;
}

viskores::cont::DataSet VFDivergenceSampling::DoExecute(const viskores::cont::DataSet& input)
{
    const auto& meanX_Field = input.GetField("meanX");
    const auto& varX_Field = input.GetField("varX");
    const auto& meanY_Field = input.GetField("meanY");
    const auto& varY_Field = input.GetField("varY");
    const auto& meanZ_Field = input.GetField("meanZ");
    const auto& varZ_Field = input.GetField("varZ");
    
    viskores::cont::ArrayHandle<SampleVec> uSamples_Handle, vSamples_Handle, wSamples_Handle;
    viskores::cont::ArrayHandle<SampleVec> divSamples_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> crossingProb_Handle;

    auto resolveType = [&](const auto& concrete_meanX)
    {
        using ArrayType = std::decay_t<decltype(concrete_meanX)>;

        ArrayType concrete_varX, concrete_meanY, concrete_varY, concrete_meanZ, concrete_varZ;
        varX_Field.GetData().AsArrayHandle(concrete_varX);
        meanY_Field.GetData().AsArrayHandle(concrete_meanY);
        varY_Field.GetData().AsArrayHandle(concrete_varY);
        meanZ_Field.GetData().AsArrayHandle(concrete_meanZ);
        varZ_Field.GetData().AsArrayHandle(concrete_varZ);

        Sampling sampling_Worklet;
        sampling_Worklet.Seed = std::random_device{}();
        this->Invoke(sampling_Worklet,
                        input.GetCellSet(),
                        concrete_meanX,
                        concrete_varX,
                        concrete_meanY,
                        concrete_varY,
                        concrete_meanZ,
                        concrete_varZ,
                        uSamples_Handle,
                        vSamples_Handle,
                        wSamples_Handle);
            
        this->Invoke(ComputeDivergence{},
                        input.GetCellSet(),
                        uSamples_Handle,
                        vSamples_Handle,
                        wSamples_Handle,
                        divSamples_Handle);

        CrossingProbabilitySampling crossingProb_Worklet;
        crossingProb_Worklet.Isovalue = this->Isovalue;
        this->Invoke(crossingProb_Worklet,
                        input.GetCellSet(),
                        divSamples_Handle,
                        crossingProb_Handle);
    };

    this->CastAndCallScalarField(meanX_Field, resolveType);

    viskores::cont::DataSet result;
    result.AddCoordinateSystem(input.GetCoordinateSystem());
    result.SetCellSet(input.GetCellSet());

    result.AddCellField("crossingProbability_sampling", crossingProb_Handle);

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

    // Uncertain Vector Field Dataset
    viskores::io::VTKDataSetReader vf_reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet vf_ds = vf_reader.ReadDataSet();

    auto vf_startAnalytical = std::chrono::high_resolution_clock::now();

    viskores::filter::uncertainty::VFDivergenceAnalytical vf_analyticalFilter;
    vf_analyticalFilter.SetIsovalue(-5.05);
    viskores::cont::DataSet vf_analyticalResult = vf_analyticalFilter.Execute(vf_ds);

    auto vf_endAnalytical = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> vf_analyticalDuration = vf_endAnalytical - vf_startAnalytical;
    std::cout << "(uncertainVectorField.vtk) Analytical Approach Computation Time: " << vf_analyticalDuration.count() << " seconds" << std::endl;

    viskores::io::VTKDataSetWriter vf_analyticalWriter("out_3d_uncertainVectorField_analytical.vtk");
    vf_analyticalWriter.WriteDataSet(vf_analyticalResult);

    auto vf_startSampling = std::chrono::high_resolution_clock::now();

    viskores::filter::uncertainty::sampling::VFDivergenceSampling vf_samplingFilter;
    vf_samplingFilter.SetIsovalue(-5.05);
    viskores::cont::DataSet vf_samplingResult = vf_samplingFilter.Execute(vf_ds);

    auto vf_endSampling = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> vf_samplingDuration = vf_endSampling - vf_startSampling;
    std::cout << "(uncertainVectorField.vtk) Sampling Approach Computation Time: " << vf_samplingDuration.count() << " seconds" << std::endl;

    viskores::io::VTKDataSetWriter vf_samplingWriter("out_3d_uncertainVectorField_sampling.vtk");
    vf_samplingWriter.WriteDataSet(vf_samplingResult);

    /*// Uncertain Red Sea Dataset
    viskores::io::VTKDataSetReader redSea_reader("data/uncertainVectorFieldRedSea.vtk");
    viskores::cont::DataSet redSea_ds = redSea_reader.ReadDataSet();

    auto redSea_startAnalytical = std::chrono::high_resolution_clock::now();

    viskores::filter::uncertainty::VFDivergenceAnalytical redSea_analyticalFilter;
    redSea_analyticalFilter.SetIsovalue(0);
    viskores::cont::DataSet redSea_analyticalResult = redSea_analyticalFilter.Execute(redSea_ds);

    auto redSea_endAnalytical = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> redSea_analyticalDuration = redSea_endAnalytical - redSea_startAnalytical;
    std::cout << "(uncertainVectorFieldRedSea.vtk) Analytical Approach Computation Time: " << redSea_analyticalDuration.count() << " seconds" << std::endl;

    viskores::io::VTKDataSetWriter redSea_analyticalWriter("out_uncertainVectorFieldRedSea_analytical.vtk");
    redSea_analyticalWriter.WriteDataSet(redSea_analyticalResult);

    auto redSea_startSampling = std::chrono::high_resolution_clock::now();

    viskores::filter::uncertainty::sampling::VFDivergenceSampling redSea_samplingFilter;
    redSea_samplingFilter.SetIsovalue(0.0);
    viskores::cont::DataSet redSea_samplingResult = redSea_samplingFilter.Execute(redSea_ds);

    auto redSea_endSampling = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> redSea_samplingDuration = redSea_endSampling - redSea_startSampling;
    std::cout << "(uncertainVectorFieldRedSea.vtk) Sampling Approach Computation Time: " << redSea_samplingDuration.count() << " seconds" << std::endl;

    viskores::io::VTKDataSetWriter redSea_samplingWriter("out_uncertainVectorFieldRedSea_sampling.vtk");
    redSea_samplingWriter.WriteDataSet(redSea_samplingResult);*/

    return 0;
}
