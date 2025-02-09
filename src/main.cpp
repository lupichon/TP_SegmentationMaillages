#include <io.hpp>
#include <example.hpp>
#include <measures.hpp>
#include <colors.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <types.hpp>

int main(int argc, char *argv[]){

    
    std::cout << "Hello !" << std::endl;

    if(argc < 2){
        throw std::invalid_argument("This program expects at least 1 argument (path to a mesh).");
    }

    const std::string meshPath = std::string{argv[1]};
    
    geomAlgoLib::Polyhedron myMesh;

    geomAlgoLib::readOFF(meshPath, myMesh);

    auto genus = geomAlgoLib::computeGenus(myMesh);
    std::cout << "The Genus of [" << meshPath << "] is = " << std::to_string(genus) << std::endl;

    geomAlgoLib::Facet_double_map anglesNormalX, anglesNormalY, anglesNormalZ;

    geomAlgoLib::Facet_double_map areas = geomAlgoLib::computeArea(myMesh);
    geomAlgoLib::Facet_double_map angles = geomAlgoLib::computeSmallestAngle(myMesh);
    geomAlgoLib::computeAngleNormalAxes(myMesh, anglesNormalX, anglesNormalY, anglesNormalZ);
    geomAlgoLib::Facet_double_map altitudes = geomAlgoLib::computeAltitudes(myMesh); 
    geomAlgoLib::Facet_double_map roughness = geomAlgoLib::computeRoughness(myMesh);
    geomAlgoLib::computeOFFFile(myMesh, areas, "../data/coloredMesh.off");

    auto itAreas = areas.begin();
    auto itAngles = angles.begin();
    auto itAltitudes = altitudes.begin();
    auto itAnglesNormalX = anglesNormalX.begin(), itAnglesNormalY = anglesNormalY.begin(), itAnglesNormalZ = anglesNormalZ.begin();
    auto itRough = roughness.begin();

    for(int i = 0; i < areas.size(); i++)
    {
        std::cout << "---Face " << i+1 << "---" << std::endl;
        std::cout << "Area = " << itAreas->second << std::endl;
        std::cout << "Smallest angle = " << itAngles->second << " (degrees)" << std::endl;
        std::cout << "Angles between normal and axes : X = " << itAnglesNormalX->second <<" Y = " << itAnglesNormalY->second << " Z = " << itAnglesNormalZ->second << " degrees" << std::endl;
        std::cout << "Altitude : " << itAltitudes->second << std::endl;
        std::cout << "Roughness : " << itRough->second << std::endl;
        std::cout << "------------" << std::endl;

        itAreas ++; 
        itAngles ++;
        itAltitudes ++;
        itAnglesNormalX ++; 
        itAnglesNormalY ++; 
        itAnglesNormalZ ++;
        itRough ++;
    }

    geomAlgoLib::Facet_string_map segAreas = geomAlgoLib::computeSegmentation(areas, 2.3, "Grande face", "Petite face");
    geomAlgoLib::Facet_string_map segSlope = geomAlgoLib::computeSegmentation(anglesNormalZ, 165, "Grande pente", "Petite pente");
    geomAlgoLib::Facet_string_map segAltitudes = geomAlgoLib::computeSegmentation(altitudes, 20, "Haute altitude", "Basse altitude");
    geomAlgoLib::Facet_string_map segRoughness = geomAlgoLib::computeSegmentation(roughness, 20, "Grande rugosite", "Petite rugosite");

    geomAlgoLib::computeColoredMesh(myMesh, segAreas, "../data/segmentationAreas.off", purple, yellow, "Grande face", "Petite face");
    geomAlgoLib::computeColoredMesh(myMesh, segSlope, "../data/segmentationSlope.off", black, white, "Grande pente", "Petite pente");
    geomAlgoLib::computeColoredMesh(myMesh, segRoughness, "../data/segmentationRoughness.off", black, white, "Grande rugosite", "Petite rugosite");
    geomAlgoLib::computeColoredMesh(myMesh, segAltitudes, "../data/segmentationAltitudes.off", red, green, "Haute altitude", "Basse altitude");

    geomAlgoLib::Facet_string_multimap allSegmentations;
    allSegmentations.insert(segSlope.begin(), segSlope.end());
    allSegmentations.insert(segAltitudes.begin(), segAltitudes.end());
    allSegmentations.insert(segRoughness.begin(), segRoughness.end());

    geomAlgoLib::computeColoredMeshMultiClass(myMesh, allSegmentations, "../data/allSegmentations.off", colors, classNames);

    geomAlgoLib::writeOFF(myMesh,"output.off");

    std::cout << "The end..." << std::endl;

    return 0;
}