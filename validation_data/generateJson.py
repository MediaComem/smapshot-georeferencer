#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 13:32:30 2020

@author: tproduit

This script is used to convert old gcps containing swiss coordinates in 
WGS 84 coordinates
"""
import json
from pyproj import CRS
from pyproj import Transformer

crs_4326 = CRS.from_epsg(4326)
crs_2056 = CRS.from_epsg(2056)
transformer = Transformer.from_crs(crs_2056, crs_4326)

with open("imagesObliqueIni.json") as json_file:
    images = json.load(json_file)


for image in images:

    gcps = image["gcp_json"]
    newGcps = []
    for gcp in gcps:
        gcp.pop("dxy", None)
        gcp.pop("errorPct", None)
        gcp.pop("xReproj", None)
        gcp.pop("yReproj", None)
        gcp.pop("errorClass", None)
        gcp.pop("eM", None)
        gcp.pop("ePix", None)
        gcp.pop("eClass", None)
        gcp["x"] = float(gcp["x"])
        gcp["y"] = float(gcp["y"])

        if "latitude" not in gcp:
            coord = transformer.transform(gcp["X"], gcp["Y"])
            gcp["longitude"] = coord[1]
            gcp["latitude"] = coord[0]

        if "altitude" not in gcp:
            gcp["altitude"] = float(gcp["Z"])

        gcp.pop("X", None)
        gcp.pop("Y", None)
        gcp.pop("Z", None)

        newGcps.append(gcp)
    image["gcp_json"] = newGcps

with open("imagesOblique.json", "w") as outfile:
    json.dump(images, outfile)
