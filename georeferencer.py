#!/usr/bin/env python
# coding: utf-8
import sys
import json
import georeferencerUtils as gu
import numpy as np

# Georeferencer args
lng = float(sys.argv[1])
lat = float(sys.argv[2])
alt = float(sys.argv[3])
azimuthDeg = float(sys.argv[4])
tiltDeg = float(sys.argv[5])
rollDeg = float(sys.argv[6])
gcps = json.loads(sys.argv[7])
width = float(sys.argv[8])
height = float(sys.argv[9])
locked = bool(int(sys.argv[10]))

# Get lat lng alt gcps
gcpLatLngAlt = gu.getLatLngAltGcps(gcps)

# Get images gcps
gcpXy = gu.getImageGcps(gcps)
gcpXy = gu.centerImageGcps(gcpXy, width, height)

# Compute normal focal
focal = gu.computeDiagonal(width, height)

if locked == False:
    (
        lngComp,
        latComp,
        altComp,
        azimuthComp,
        tiltComp,
        rollComp,
        focalComp,
        pComp,
        gcps,
        imageCoordinates,
        method,
    ) = gu.georeferencer(
        lng,
        lat,
        alt,
        azimuthDeg,
        tiltDeg,
        rollDeg,
        focal,
        width,
        height,
        gcps,
        plotBool=False,
    )
else:
    (
        lngComp,
        latComp,
        altComp,
        azimuthComp,
        tiltComp,
        rollComp,
        focalComp,
        pComp,
        gcps,
        imageCoordinates,
        method,
    ) = gu.georeferencerLocked(
        lng,
        lat,
        alt,
        azimuthDeg,
        tiltDeg,
        rollDeg,
        focal,
        width,
        height,
        gcps,
        plotBool=False,
    )

imageCoordinatesForGltf = gu.computeImageCoordinatesForGltf(pComp, width, height)
# Compute gcp errors
# Store results
result = {}
result["latitude"] = latComp
result["longitude"] = lngComp
result["altitude"] = altComp
result["focal"] = focalComp
result["tilt"] = tiltComp
result["roll"] = rollComp
result["azimuth"] = azimuthComp
result["GCPs"] = gcps
result["imageCoordinatesForGltf"] = imageCoordinatesForGltf
result["p"] = pComp

# Print results (node capture the printed parameters)
print(json.dumps(result))

sys.stdout.flush()
