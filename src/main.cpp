#include <io.hpp>
#include <example.hpp>
#include <measures.hpp>
#include <colors.hpp>
#include <smoothing.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <types.hpp>

int main(int argc, char *argv[]){
    
    // Print a greeting message to indicate the program has started
    std::cout << "Hello !" << std::endl;

    // Check if the program has received the required number of arguments
    if(argc < 2){
        // If not, throw an exception indicating that at least one argument (path to a mesh) is required
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    // Extract the mesh file path from the command line argument
    const std::string meshPath = std::string{argv[1]};
    
    // Declare a Polyhedron object to store the mesh data
    geomAlgoLib::Polyhedron myMesh;

    // Read the mesh from the .off file specified by the user
    geomAlgoLib::readOFF(meshPath, myMesh);

    // Compute the genus of the mesh and print the result
    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    // Declare maps to store data angles
    geomAlgoLib::Facet_double_map anglesNormalX, anglesNormalY, anglesNormalZ;    

    //calculate local properties
    geomAlgoLib::Facet_double_map areas = geomAlgoLib::computeArea(myMesh);                   // Area
    geomAlgoLib::Facet_double_map angles = geomAlgoLib::computeSmallestAngle(myMesh);         // Smallest angle 
    geomAlgoLib::computeAngleNormalAxes(myMesh, anglesNormalX, anglesNormalY, anglesNormalZ); // Angle between normal and X, Y, Z axes
    geomAlgoLib::Facet_double_map altitudes = geomAlgoLib::computeAltitudes(myMesh);          // Altitude
    geomAlgoLib::Facet_double_map roughness = geomAlgoLib::computeRoughness(myMesh);          // Roughness

    // Create a colored mesh based on the areas of the facets and save it to a file
    geomAlgoLib::computeOFFFile(myMesh, areas, "../data/coloredMesh.off");

    // Initialize iterators to loop through the maps of area, angles, altitudes, and roughness
    auto itAreas = areas.begin();
    auto itAngles = angles.begin();
    auto itAltitudes = altitudes.begin();
    auto itAnglesNormalX = anglesNormalX.begin(), itAnglesNormalY = anglesNormalY.begin(), itAnglesNormalZ = anglesNormalZ.begin();
    auto itRough = roughness.begin();

    // Loop through all faces of the mesh and print the corresponding values for each metric
    for(int i = 0; i < areas.size(); i++)
    {
        std::cout << "---Face " << i+1 << "---" << std::endl;
        std::cout << "Area = " << itAreas->second << std::endl;
        std::cout << "Smallest angle = " << itAngles->second << " (degrees)" << std::endl;
        std::cout << "Angles between normal and axes : X = " << itAnglesNormalX->second <<" Y = " << itAnglesNormalY->second << " Z = " << itAnglesNormalZ->second << " degrees" << std::endl;
        std::cout << "Altitude : " << itAltitudes->second << std::endl;
        std::cout << "Roughness : " << itRough->second << std::endl;
        std::cout << "------------" << std::endl;

        // Move to the next entry in each map
        itAreas++; 
        itAngles++;
        itAltitudes++;
        itAnglesNormalX++; 
        itAnglesNormalY++; 
        itAnglesNormalZ++;
        itRough++;
    }

    // Perform segmentation based on different criteria (area, slope, altitude, roughness)
    geomAlgoLib::Facet_string_map segAreas = geomAlgoLib::computeSegmentation(areas, 2.3, "Large face", "Small face");
    geomAlgoLib::Facet_string_map segSlope = geomAlgoLib::computeSegmentation(anglesNormalZ, 165, "Steep slope", "Gentle slope");
    geomAlgoLib::Facet_string_map segAltitudes = geomAlgoLib::computeSegmentation(altitudes, 20, "High altitude", "Low altitude");
    geomAlgoLib::Facet_string_map segRoughness = geomAlgoLib::computeSegmentation(roughness, 20, "High roughness", "Low roughness");

    // Generate colored meshes for each segmentation and save them to respective files
    geomAlgoLib::computeColoredMesh(myMesh, segAreas, "../data/segmentationAreas.off", purple, yellow, "Large face", "Small face");
    geomAlgoLib::computeColoredMesh(myMesh, segSlope, "../data/segmentationSlope.off", black, white, "Steep slope", "Gentle slope");
    geomAlgoLib::computeColoredMesh(myMesh, segAltitudes, "../data/segmentationAltitudes.off", red, green, "High altitude", "Low altitude");
    geomAlgoLib::computeColoredMesh(myMesh, segRoughness, "../data/segmentationRoughness.off", black, white, "High roughness", "Low roughness");

    // Combine all segmentations into a single multimap for multi-class segmentation
    geomAlgoLib::Facet_string_multimap allSegmentations;
    allSegmentations.insert(segSlope.begin(), segSlope.end());
    allSegmentations.insert(segAltitudes.begin(), segAltitudes.end());
    allSegmentations.insert(segRoughness.begin(), segRoughness.end());

    // Generate a colored mesh for the multi-class segmentation and save it to a file
    geomAlgoLib::computeColoredMeshMultiClass(myMesh, allSegmentations, "../data/allSegmentations.off", colors, classNames);
    
    // Write the final mesh (after processing) to an output file
    geomAlgoLib::writeOFF(myMesh,"output.off");

    // Smoothing
    geomAlgoLib::Polyhedron smoothedMeshLaplacian = geomAlgoLib::laplacianSmoothing(myMesh, 10);
    geomAlgoLib::writeOFF(smoothedMeshLaplacian, "../data/smoothedMeshLaplacian.off");

    geomAlgoLib::Polyhedron smoothedMeshGaussian = geomAlgoLib::gaussianSmoothing(myMesh, 10, 1);
    geomAlgoLib::writeOFF(smoothedMeshGaussian, "../data/smoothedMeshGaussian.off");

    // Print a message indicating that the program has finished executing
    std::cout << "The end..." << std::endl;

    return 0;
}
