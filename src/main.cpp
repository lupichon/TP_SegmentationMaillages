#include <io.hpp>
#include <example.hpp>
#include <measures.hpp>

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

    geomAlgoLib::Facet_double_map areas = geomAlgoLib::computeArea(myMesh);
    geomAlgoLib::Facet_double_map angles = geomAlgoLib::computeSmallestAngle(myMesh);

    geomAlgoLib::Facet_double_map anglesNormalX, anglesNormalY, anglesNormalZ;
    geomAlgoLib::computeAngleNormalAxes(myMesh, anglesNormalX, anglesNormalY, anglesNormalZ);

    auto itAreas = areas.begin();
    auto itAngles = angles.begin();
    auto itAnglesNormalX = anglesNormalX.begin(), itAnglesNormalY = anglesNormalY.begin(), itAnglesNormalZ = anglesNormalZ.begin();

    for(int i = 0; i < areas.size(); i++)
    {
        std::cout << "---Face " << i+1 << "---" << std::endl;
        std::cout << "Area = " << itAreas->second << std::endl;
        std::cout << "Smallest angle = " << itAngles->second << " (degrees)" << std::endl;
        std::cout << "Angles between normal and axes : X = " << itAnglesNormalX->second <<" Y = " << itAnglesNormalY->second << " Z = " << itAnglesNormalZ->second << " degrees" << std::endl;
        std::cout << "------------" << std::endl;

        itAreas ++; 
        itAngles ++;
        itAnglesNormalX ++; 
        itAnglesNormalY ++; 
        itAnglesNormalZ ++;
    }

    geomAlgoLib::Facet_string_map segAreas = geomAlgoLib::computeSegmentation(areas, 2.3, "Grande face", "Petite face");
    geomAlgoLib::Facet_string_map segSlope = geomAlgoLib::computeSegmentation(anglesNormalZ, 165, "Grande pente", "Petite pente");

    geomAlgoLib::computeColoredMesh(myMesh, segAreas, "../data/segmentationAreas.off", CGAL::Color(0,255,0), CGAL::Color(255,0,0), "Grande face", "Petite face");
    geomAlgoLib::computeColoredMesh(myMesh, segSlope, "../data/segmentationSlope.off", CGAL::Color(0,0,0), CGAL::Color(255,255,255), "Grande pente", "Petite pente");


    geomAlgoLib::writeOFF(myMesh,"output.off");
    geomAlgoLib::computeSmallestAngle(myMesh);

    geomAlgoLib::computeOFFFile(myMesh, areas, "../data/coloredMesh.off");

    std::cout << "The end..." << std::endl;
    return 0;
}