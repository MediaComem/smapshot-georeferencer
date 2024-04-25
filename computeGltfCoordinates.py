#!/usr/bin/env python
# coding: utf-8
import json
import sys

import georeferencerUtils as gu

"""Compute image coordinates to generate a gltf

Parameters:
azimuthDeg (float): azimuth in degrees (compatible with cesium)
tiltDeg (float): tilt in degrees (compatible with cesium)
rollDeg (float): roll in degrees (compatible with cesium)
width (float): image width
height (float): image height
focal (float): focal in pixels

Returns:
json: image coordinates

"""
# Georeferencer args
azimuthDeg = float(sys.argv[1])
tiltDeg = float(sys.argv[2])
rollDeg = float(sys.argv[3])
width = float(sys.argv[4])
height = float(sys.argv[5])
focal = float(sys.argv[6])

# Convert cesium angles to LM angles
azimuth, tilt, roll = gu.cesiumToLmAngles(azimuthDeg, tiltDeg, rollDeg)

# Create pose vector
p = [0.0, 0.0, 0.0, azimuth, tilt, roll, focal, 0, 0]

# Generate the images coordinates to be inserted in the collada file
imageCoordinatesForGltf = gu.computeImageCoordinatesForGltf(p, width, height)

# Store results
result = {}
result["imageCoordinates"] = imageCoordinatesForGltf

# Print results (node capture the printed parameters)
print(json.dumps(result))

sys.stdout.flush()
