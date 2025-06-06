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

int main(int argc, char** argv)
{
    viskores::cont::Initialize(argc, argv);

    // Declare visualization file
    const char* input = "data/uncertainVectorField.vtk";

    // Read and write visualization file
    viskores::io::VTKDataSetReader reader(input);
    viskores::cont::DataSet ds = reader.ReadDataSet();
    viskores::io::VTKDataSetWriter writer("out_vf_divergence.vtk");
    writer.WriteDataSet(ds);

    return 0;
}
