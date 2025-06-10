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
#include <viskores/cont/Initialize.h>
#include <viskores/io/VTKDataSetReader.h>
#include <viskores/io/VTKDataSetWriter.h>

#include <viskores/filter/vector_analysis/VectorMagnitude.h>

// Manual computation of vector magnitude
/* double ComputeVectorMagnitude(double x, double y, double z)
{
    return std::sqrt(x * x + y * y + z * z);
} */

int main(int argc, char** argv)
{
    viskores::cont::Initialize(argc, argv);

    // Read input file
    viskores::io::VTKDataSetReader reader("data/uncertainVectorField.vtk");
    viskores::cont::DataSet ds = reader.ReadDataSet();

    // Apply Vector Magnitude filter
    viskores::filter::vector_analysis::VectorMagnitude filter;
    filter.SetActiveField("uncertainty");
    filter.SetOutputFieldName("magnitude");
    viskores::cont::DataSet output = filter.Execute(ds);

    // Write output to file
    viskores::io::VTKDataSetWriter writer("out_vf_magnitude.vtk");
    writer.WriteDataSet(output);

    return 0;
}
