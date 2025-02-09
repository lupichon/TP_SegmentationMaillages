#ifndef HPP_COLOR_HPP
#define HPP_COLOR_HPP

#include <types.hpp>

const geomAlgoLib::Color white = {1.0f, 1.0f, 1.0f};
const geomAlgoLib::Color black = {0.0f, 0.0f, 0.0f};
const geomAlgoLib::Color red = {0.8f, 0.2f, 0.2f};
const geomAlgoLib::Color green = {0.2f, 0.8f, 0.2f};
const geomAlgoLib::Color yellow = {0.8f, 0.8f, 0.2f};
const geomAlgoLib::Color purple = {0.4f, 0.2f, 0.8f};
const geomAlgoLib::Color orange = {0.8f, 0.4f, 0.0f};
const geomAlgoLib::Color blue = {0.0f, 0.0f, 1.0f};

const std::vector<geomAlgoLib::Color> colors = {
    red,  // Grande pente (Rouge)
    green,  // Petite pente (Vert)
    yellow,  // Grande rugosité (Jaune)
    white,  // Petite rugosité (Bleu clair)
    blue,  // Haute altitude (Jaune clair)
    purple,  // Basse altitude (Violet)
};

const std::vector<std::string> classNames{
    "Grande pente", "Petite pente", 
    "Grande rugosite", "Petite rugosite", 
    "Haute altitude", "Basse altitude"
};

#endif