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
#include <viskores/cont/Field.h>
#include <viskores/cont/Initialize.h>
#include <viskores/cont/Invoker.h>
#include <viskores/cont/RuntimeDeviceTracker.h>
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
#include <string>
#include <vector>
#include <numeric>

#include <omp.h>

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
                                  FieldOut divMean,
                                  FieldOut divVar);
    using ExecutionSignature = void(Boundary, _2, _3, _4, _5, _6, _7);
    using InputDomain = _1;

    template <typename BoundaryType, typename NeighborhoodType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                  const NeighborhoodType& meanX,
                                  const NeighborhoodType& varX,
                                  const NeighborhoodType& meanY,
                                  const NeighborhoodType& varY,
                                  OutType& divMean,
                                  OutType& divVar) const
    {
        OutType divmeanX, varsquaredX;
        OutType divmeanY, varsquaredY;

        if (boundary.MinNeighborIndices(1)[0] == 0)
        {
            divmeanX = meanX.Get(0, 1, 0) - meanX.Get(0, 0, 0);
            varsquaredX = varX.Get(0, 1, 0) + varX.Get(0, 0, 0);
        }
        else if (boundary.MaxNeighborIndices(1)[0] == 0)
        {
            divmeanX = meanX.Get(0, 0, 0) - meanX.Get(0, -1, 0);
            varsquaredX = varX.Get(0, 0, 0) + varX.Get(0, -1, 0);
        }
        else
        {
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

        divMean = divmeanX + divmeanY;
        divVar = varsquaredX + varsquaredY;
    }
};

// Worklet 2 Helper Function: computes Gaussian CDF.
template <typename T>
VISKORES_EXEC T GaussianCDF(T x, T mu, T sigma)
{
    if (sigma <= static_cast<T>(0.0))
    {
        if (x < mu)
        {
            return static_cast<T>(0.0);
        }
        else
        {
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

    viskores::cont::ArrayHandle<viskores::Float64> divMean_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> divVar_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> crossingProb_Handle;

    auto resolveType = [&](const auto& concrete_meanX) {
        using ArrayType = std::decay_t<decltype(concrete_meanX)>;

        ArrayType concrete_varX, concrete_meanY, concrete_varY;
        varX_Field.GetData().AsArrayHandle(concrete_varX);
        meanY_Field.GetData().AsArrayHandle(concrete_meanY);
        varY_Field.GetData().AsArrayHandle(concrete_varY);

        this->Invoke(ComputeDivergenceMeanAndVar{},
                     input.GetCellSet(),
                     concrete_meanX,
                     concrete_varX,
                     concrete_meanY,
                     concrete_varY,
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
    void SetNumSamples(viskores::Id numSamples)
    {
        this->NumSamples = numSamples;
    }

private:
    viskores::cont::DataSet DoExecute(const viskores::cont::DataSet& input) override;
    viskores::Float64 Isovalue;
    viskores::Id NumSamples;
};

namespace
{
constexpr viskores::Id MAX_SAMPLES = 2000;
using SampleVec = viskores::Vec<viskores::Float64, MAX_SAMPLES>;

// Worklet 1: sampling.
struct Sampling : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                  FieldInNeighborhood meanX,
                                  FieldInNeighborhood varX,
                                  FieldInNeighborhood meanY,
                                  FieldInNeighborhood varY,
                                  FieldOut uSamples,
                                  FieldOut vSamples);
    using ExecutionSignature = void(_2, _3, _4, _5, _6, _7);
    using InputDomain = _1;

    viskores::UInt32 Seed;
    viskores::Id NumSamples; // Actual number of samples to generate

    template <typename MeanXType, typename VarXType, typename MeanYType, typename VarYType>
    VISKORES_EXEC void operator()(const MeanXType& meanX,
                                  const VarXType& varX,
                                  const MeanYType& meanY,
                                  const VarYType& varY,
                                  SampleVec& uSamples,
                                  SampleVec& vSamples) const
    {
        std::size_t hash = std::hash<double>()(meanX.Get(0, 0, 0)) ^ std::hash<double>()(meanY.Get(0, 0, 0));
        std::mt19937 rng(Seed ^ hash);
        std::normal_distribution<double> distU(meanX.Get(0, 0, 0), std::sqrt(varX.Get(0, 0, 0)));
        std::normal_distribution<double> distV(meanY.Get(0, 0, 0), std::sqrt(varY.Get(0, 0, 0)));
        
        // Generate the requested number of samples
        for (viskores::Id i = 0; i < this->NumSamples; ++i)
        {
            uSamples[i] = distU(rng);
            vSamples[i] = distV(rng);
        }
        // Zero out the rest of the vector to avoid garbage values
        for (viskores::Id i = this->NumSamples; i < MAX_SAMPLES; ++i)
        {
            uSamples[i] = 0.0;
            vSamples[i] = 0.0;
        }
    }
};

// Worklet 2: computes divergence for a 2d vector field via sampling.
struct ComputeDivergence : public viskores::worklet::WorkletPointNeighborhood
{
    using ControlSignature = void(CellSetIn domain,
                                  FieldInNeighborhood uSamples,
                                  FieldInNeighborhood vSamples,
                                  FieldOut divergenceSamples);
    using ExecutionSignature = void(Boundary, _2, _3, _4);
    using InputDomain = _1;

    viskores::Id NumSamples; // Actual number of samples to compute

    template <typename BoundaryType, typename USamplesType, typename VSamplesType, typename OutType>
    VISKORES_EXEC void operator()(const BoundaryType& boundary,
                                  const USamplesType& uSamples,
                                  const VSamplesType& vSamples,
                                  OutType& divSamples) const
    {
        for (viskores::Id k = 0; k < this->NumSamples; ++k)
        {
            double dudx, dvdy;

            if (boundary.MinNeighborIndices(1)[1] == 0)
            {
                dudx = uSamples.Get(0, 1, 0)[k] - uSamples.Get(0, 0, 0)[k];
            }
            else if (boundary.MaxNeighborIndices(1)[1] == 0)
            {
                dudx = uSamples.Get(0, 0, 0)[k] - uSamples.Get(0, -1, 0)[k];
            }
            else
            {
                dudx = (uSamples.Get(0, 1, 0)[k] - uSamples.Get(0, -1, 0)[k]) / 2.0;
            }

            if (boundary.MinNeighborIndices(1)[0] == 0)
            {
                dvdy = vSamples.Get(1, 0, 0)[k] - vSamples.Get(0, 0, 0)[k];
            }
            else if (boundary.MaxNeighborIndices(1)[0] == 0)
            {
                dvdy = vSamples.Get(0, 0, 0)[k] - vSamples.Get(-1, 0, 0)[k];
            }
            else
            {
                dvdy = (vSamples.Get(1, 0, 0)[k] - vSamples.Get(-1, 0, 0)[k]) / 2.0;
            }

            divSamples[k] = dudx + dvdy;
        }
        // Zero out the rest of the vector
        for (viskores::Id k = this->NumSamples; k < MAX_SAMPLES; ++k)
        {
            divSamples[k] = 0.0;
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
    viskores::Id NumSamples; // Actual number of samples to check

    template <typename CellShapeTag, typename DivSamplesType>
    VISKORES_EXEC viskores::Float64 operator()(CellShapeTag,
                                              const DivSamplesType& divSamples) const
    {
        if (this->NumSamples == 0)
        {
            return 0.0;
        }

        viskores::Id crossings = 0;

        for (viskores::Id k = 0; k < this->NumSamples; ++k)
        {
            const int numCorners = divSamples.GetNumberOfComponents();
            double minVal = divSamples[0][k];
            double maxVal = minVal;

            for (int c = 1; c < numCorners; ++c)
            {
                double val = divSamples[c][k];
                if (val < minVal)
                {
                    minVal = val;
                }
                if (val > maxVal)
                {
                    maxVal = val;
                }
            }

            if (this->Isovalue > minVal && this->Isovalue < maxVal)
            {
                crossings++;
            }
        }

        return static_cast<double>(crossings) / this->NumSamples;
    }
};

} // namespace

// Filter Implementation
VISKORES_CONT
VFDivergenceSampling::VFDivergenceSampling()
{
    this->Isovalue = 0.0;
    this->NumSamples = 100; // Default value
}

viskores::cont::DataSet VFDivergenceSampling::DoExecute(const viskores::cont::DataSet& input)
{
    const auto& meanX_Field = input.GetField("meanX");
    const auto& varX_Field = input.GetField("varX");
    const auto& meanY_Field = input.GetField("meanY");
    const auto& varY_Field = input.GetField("varY");

    viskores::cont::ArrayHandle<SampleVec> uSamples_Handle, vSamples_Handle;
    viskores::cont::ArrayHandle<SampleVec> divSamples_Handle;
    viskores::cont::ArrayHandle<viskores::Float64> crossingProb_Handle;

    auto resolveType = [&](const auto& concrete_meanX) {
        using ArrayType = std::decay_t<decltype(concrete_meanX)>;

        ArrayType concrete_varX, concrete_meanY, concrete_varY;
        varX_Field.GetData().AsArrayHandle(concrete_varX);
        meanY_Field.GetData().AsArrayHandle(concrete_meanY);
        varY_Field.GetData().AsArrayHandle(concrete_varY);

        Sampling sampling_Worklet;
        sampling_Worklet.Seed = std::random_device{}();
        sampling_Worklet.NumSamples = this->NumSamples;
        this->Invoke(sampling_Worklet,
                     input.GetCellSet(),
                     concrete_meanX,
                     concrete_varX,
                     concrete_meanY,
                     concrete_varY,
                     uSamples_Handle,
                     vSamples_Handle);

        ComputeDivergence computeDivergence_Worklet;
        computeDivergence_Worklet.NumSamples = this->NumSamples;
        this->Invoke(computeDivergence_Worklet,
                     input.GetCellSet(),
                     uSamples_Handle,
                     vSamples_Handle,
                     divSamples_Handle);

        CrossingProbabilitySampling crossingProb_Worklet;
        crossingProb_Worklet.Isovalue = this->Isovalue;
        crossingProb_Worklet.NumSamples = this->NumSamples;
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

int main(int argc, char* argv[])
{
    viskores::cont::Initialize(argc, argv);

    viskores::io::VTKDataSetReader reader("data/uncertainRedSea2D.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    const viskores::FloatDefault isovalue = 0.003;

    std::cout << "--- Running Analytical Approach ---\n";
    auto startAnalytical = std::chrono::high_resolution_clock::now();
    viskores::filter::uncertainty::VFDivergenceAnalytical analyticalFilter;
    analyticalFilter.SetIsovalue(isovalue);
    viskores::cont::DataSet analyticalResult = analyticalFilter.Execute(ds);
    auto endAnalytical = std::chrono::high_resolution_clock::now();
    double analyticalTime = std::chrono::duration<double>(endAnalytical - startAnalytical).count();

    std::cout << "Analytical Computation Time: " << analyticalTime << " Seconds\n";
    viskores::io::VTKDataSetWriter analyticalWriter("out_uncertainRedSea2D_analytical.vtk");
    analyticalWriter.WriteDataSet(analyticalResult);
    std::cout << "Wrote analytical results to out_uncertainRedSea2D_analytical.vtk\n";
    std::cout << "-----------------------------------\n\n";

    std::cout << "--- Benchmarking Sampling Approach ---\n";
    for (int samples = 100; samples <= 2000; samples += 100)
    {
        std::cout << "Running with " << samples << " samples..." << std::flush;

        auto startSampling = std::chrono::high_resolution_clock::now();

        viskores::filter::uncertainty::sampling::VFDivergenceSampling samplingFilter;
        samplingFilter.SetNumSamples(samples);
        samplingFilter.SetIsovalue(isovalue);
        viskores::cont::DataSet samplingResult = samplingFilter.Execute(ds);

        auto endSampling = std::chrono::high_resolution_clock::now();
        double samplingTime = std::chrono::duration<double>(endSampling - startSampling).count();

        std::cout << " Time: " << samplingTime << " seconds." << std::endl;

        // Create a unique filename for each run
        std::string output_filename = "result_data/red_sea/out_uncertainRedSea2D_" + std::to_string(samples) + ".vtk";
        
        try
        {
            viskores::io::VTKDataSetWriter samplingWriter(output_filename);
            samplingWriter.WriteDataSet(samplingResult);
            std::cout << "  -> Wrote results to " << output_filename << std::endl;
        }
        catch (const std::exception& e)
        {
            std::cerr << "Error writing file: " << output_filename << ". " << e.what() << std::endl;
            std::cerr << "Please ensure the directory 'result_data/red_sea/' exists." << std::endl;
        }
    }
    std::cout << "--------------------------------------\n";

    return 0;
}
