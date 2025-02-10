#ifndef HPP_COLOR_HPP
#define HPP_COLOR_HPP

#include <types.hpp>

/**
 * @file color.hpp
 * @brief This header file defines color constants and related data for terrain feature segmentation.
 * 
 * The file includes color definitions using the RGB model, where each color is represented as a `geomAlgoLib::Color` structure 
 * containing three float values for red, green, and blue components. These colors are used for visualizing and classifying different 
 * terrain features based on characteristics like slope, roughness, and altitude in a mesh. It also defines lists of colors and their 
 * corresponding class names used for terrain classification.
 *
 * **Color Definitions:**
 * - `white`, `black`, `red`, `green`, `yellow`, `purple`, `orange`, `blue`
 * 
 * **Segmentation Classes:**
 * - Terrain features are segmented into classes such as "Steep slope", "Gentle slope", "High roughness", 
 *   "Low roughness", "High altitude", and "Low altitude". These classes are assigned specific colors for visualization.
 * 
 * **Usage:**
 * - These color definitions and class names are intended for use in terrain feature classification and segmentation tasks, 
 *   such as coloring the faces of a mesh based on certain properties (e.g., slope, altitude, roughness).
 */

const geomAlgoLib::Color white = {1.0f, 1.0f, 1.0f}; // white: RGB values for white color
const geomAlgoLib::Color black = {0.0f, 0.0f, 0.0f}; // black: RGB values for black color
const geomAlgoLib::Color red = {0.8f, 0.2f, 0.2f};   // red: RGB values for red color
const geomAlgoLib::Color green = {0.2f, 0.8f, 0.2f}; // green: RGB values for green color
const geomAlgoLib::Color yellow = {0.8f, 0.8f, 0.2f}; // yellow: RGB values for yellow color
const geomAlgoLib::Color purple = {0.4f, 0.2f, 0.8f}; // purple: RGB values for purple color
const geomAlgoLib::Color orange = {0.8f, 0.4f, 0.0f}; // orange: RGB values for orange color
const geomAlgoLib::Color blue = {0.0f, 0.0f, 1.0f};  // blue: RGB values for blue color

/**
 * @brief List of colors used for the global segmentation
 * The colors represent different classes based on terrain features like slope, roughness, and altitude.
 * Each color corresponds to a specific class of features.
 * - Steep slope (Red)
 * - Gentle slope (Green)
 * - High roughness (Yellow)
 * - Low roughness (White)
 * - High altitude (Blue)
 * - Low altitude (Purple)
 */
const std::vector<geomAlgoLib::Color> colors = {
    red,    // Steep slope (Red)
    green,  // Gentle slope (Green)
    yellow, // High roughness (Yellow)
    white,  // Low roughness (White)
    blue,   // High altitude (Blue)
    purple, // Low altitude (Purple)
};

/**
 * @brief List of class names for terrain features used in segmentation
 * Each class corresponds to a specific terrain characteristic.
 * - Steep slope
 * - Gentle slope
 * - High roughness
 * - Low roughness
 * - High altitude
 * - Low altitude
 */
const std::vector<std::string> classNames{
    "Steep slope", "Gentle slope", 
    "High roughness", "Low roughness", 
    "High altitude", "Low altitude"
};

#endif
