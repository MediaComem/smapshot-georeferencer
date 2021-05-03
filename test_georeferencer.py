#!/usr/bin/env python
# coding: utf-8
import json
import numpy as np
import georeferencerUtils as gu
import matplotlib.pyplot as plt
import pymap3d as pm

plot = False
with open('./validation_data/imagesOblique.json') as json_file:
    images = json.load(json_file)

deltas = []
errors = []
georeferencerFailedCounter = 0
counter = 0
for image in images:
    print('--------')
    print ('image ID: ', image['id'])
    print('https://smapshot.heig-vd.ch/visit/{imageId}'.format(imageId=image['id']))

    counter+=1
    
    lng, lat, alt, azimuthDeg, tiltDeg, rollDeg, focal, gcps, width, height = gu.getGeolocalisationFromJson(image)
    
    # Get lat lng alt gcps
    gcpLatLngAlt = gu.getLatLngAltGcps(gcps)
    
    # Get images gcps
    gcpXy = gu.getImageGcps(gcps)
    gcpXy = gu.centerImageGcps(gcpXy, width, height)
    
    if plot == True:
        # Check values in DB
        ###
        # Convert cesium angles to LM angles
        azimuth, tilt, roll = gu.cesiumToLmAngles(azimuthDeg, tiltDeg, rollDeg)
       
        # Get GCP in ENU at the recorded location
        gcpEnuCheck = gu.convertEnu(gcpLatLngAlt, lat, lng, alt)
        p = [0., 0., 0., azimuth, tilt, roll, focal, 0, 0]
        xyForward = gu.project3Dto2D(gcpEnuCheck,p)
        plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
        plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
        plt.axis('square')
        plt.title('GCP projected with pose recorded in database')
        plt.show()

    
    # Run 10 time for each image
    for i in range(10):
        # Generate simulated apriori values
        lng0, lat0, alt0, azimuthDeg0, tiltDeg0, rollDeg0, focal0 = gu.generateSimulatedApriori(lng, lat, alt, azimuthDeg, tiltDeg, rollDeg, focal)
    
        # Call georeferencer with apriori values
        georeferencerSucceed = True
        try:
            lngComp, latComp, altComp, azimuthComp, tiltComp, rollComp, focalComp, pComp, gcps, imageCoordinates, method = gu.georeferencer(lng0, lat0, alt0, azimuthDeg0, tiltDeg0, rollDeg0, focal0, width, height, gcps, plotBool=False)
        except:
            georeferencerSucceed = False
            georeferencerFailedCounter += 1
            

        if georeferencerSucceed:
                
            # Store gcp errors
            for gcp in gcps:
                errors.append(gcp['dxy'])
                
            # Compute deltas in meter
            dE, dN, dU = pm.geodetic2enu(latComp, lngComp, altComp, lat, lng, alt, ell=None, deg=True)
            delta = np.sqrt(dE*dE + dN*dN + dU*dU)
            print('Test {i} delta:'.format(i=i+1), np.round(delta,1), method)
            if delta > 100:
                print('dE', np.round(dE,1))
                print('dN', np.round(dN,1))
                print('dU', np.round(dU,1))
                print ('Altitude DB: ', alt)
                print ('Altitude Calc: ', altComp)
                
                # Check result
                ###
                imagePointsCamera = gu.createImageMains(pComp[6], width/2, height/2).T
                imagePointsWorld = gu.camera2world(imagePointsCamera, pComp)
                imagePointsImage = gu.project3Dto2D(imagePointsWorld.T, pComp)
                
                # Check gcps
                gcpComp = gu.convertEnu(gcpLatLngAlt, latComp, lngComp, altComp)
                xyForward = gu.project3Dto2D(gcpComp, pComp)
                
                plt.plot([imagePointsImage[1,0], imagePointsImage[1,1], imagePointsImage[1,3], imagePointsImage[1,2], imagePointsImage[1,0]], [imagePointsImage[0,0], imagePointsImage[0,1], imagePointsImage[0,3], imagePointsImage[0,2], imagePointsImage[0,0]], '-k')
                plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
                plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
                plt.axis('square')
                plt.title('GCPs projected with computed pose')
                plt.show()
                

            dLat = lat-latComp
            dLng = lng-lngComp
            dAlt = alt-altComp
            if azimuthDeg < 0:
                azimuthDeg += 360
            if azimuthComp < 0:
                azimuthComp += 360
            if tiltDeg < 0:
                tiltDeg += 360
            if tiltComp < 0:
                tiltComp += 360
            if rollDeg < 0:
                rollDeg += 360
            if rollComp < 0:
                rollComp += 360
            dAz = (azimuthDeg-azimuthComp)
            dTilt = (tiltDeg-tiltComp)
            dRoll = (rollDeg-rollComp)
            dFocal = focal-focalComp
            deltas.append([dLat, dLng, dAlt, dAz, dTilt, dRoll, dFocal])
            
        else:
            print('Georeferencer failed')
            # Check values in DB
            ###
            # Convert cesium angles to LM angles
            azimuth, tilt, roll = gu.cesiumToLmAngles(azimuthDeg, tiltDeg, rollDeg)
            # Get lat lng alt gcps
            gcpLatLngAlt = gu.getLatLngAltGcps(gcps)
            # Get images gcps
            gcpXy = gu.getImageGcps(gcps)
            gcpXy = gu.centerImageGcps(gcpXy, width, height)
            # Get GCP in ENU at the recorded location
            gcpEnuCheck = gu.convertEnu(gcpLatLngAlt, lat, lng, alt)
            p = [0., 0., 0., azimuth, tilt, roll, focal, 0, 0]
            xyForward = gu.project3Dto2D(gcpEnuCheck,p)
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
            plt.axis('square')
            plt.title('GCP projected with pose recorded in database')
            plt.show()
            
    print (counter, '/', len(images))

print ('Number of failure:', georeferencerFailedCounter)