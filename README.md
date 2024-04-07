# Rasterizer

## Introduction
This project is a 3D rasterizer designed for rendering scenes described in XML files. It supports complex scene configurations, including multiple objects, camera views, color management, and geometric transformations (rotation and scaling).

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Features](#features)
- [Dependencies](#dependencies)
- [Configuration](#configuration)
- [Examples](#examples)
- [Contributors](#contributors)
- [License](#license)

## Installation
To compile the rasterizer, ensure you have `g++` installed and use the following command in the project root directory:
```
make all
```
This compiles the source files into an executable named `rasterizer`.

## Usage
To run the rasterizer, use the following command:
```
./rasterizer <input_file_name>
```
The `<input_file_name>` should be the path to an XML file describing the scene you wish to render.

## Features
- **XML Scene Configuration:** Define scenes with multiple objects, cameras, and lights in an XML format.
- **Camera Support:** Configurable camera properties including position, gaze direction, and projection parameters.
- **Color Management:** Define and manage object colors.
- **Geometric Transformations:** Apply rotation and scaling transformations to objects within the scene.

## Dependencies
- g++ (for compilation)
- TinyXML-2 library (for XML parsing)

## Configuration
The project expects scene descriptions in XML format. These files should specify the scene's cameras, objects, and their properties.

## Examples
An example XML configuration might look like this:
```xml
<Scene>
    <Cameras>
        <!-- Camera definitions go here -->
    </Cameras>
    <Rotations>
        <!-- Rotation definitions go here -->
    </Rotations>
    .
    .
    .
</Scene>
```
Refer to the `input_outputs.zip` for more XML scene files and expected results.

## Contributors
- Kaan Karaçanta
- Taha Emir Gökçegöz
