#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 10:52:38 2020

@author: tproduit
"""
import numpy as np
import cv2
from scipy.spatial.transform import Rotation as R
from random import random
import pymap3d as pm
from lmfit import minimize, Parameters
# Trick to avoid the installation of matplotlib in production server
try:
    import matplotlib.pyplot as plt
except:
    pass
try:
    import psycopg2
except:
    pass


def georeferencer(lng0, lat0, alt0, azimuth0, tilt0, roll0, focal0, width, height,
                  gcps, plotBool = False):

    """Compute camera pose given apriori values and gcps

    Parameters:
    lng0 (float): apriori longitude in degrees
    lat0 (float): apriori latitude in degrees
    alt0 (float): apriori altitude in meter
    azimuth0 (float): apriori azimuth in degrees (from cesium)
    tilt0 (float): apriori tilt in degrees (from cesium)
    roll0 (float): apriori roll in degrees (from cesium)
    focal0 (float): apriori focal in pixels
    width (float): image width
    height (float): image height
    gcps (dictionnary): gcps in smapshot format
    plotBool (boolean): if the plot must be shown

    Returns:
    float: longitude in degrees
    float: latitude in degrees
    float: altitude in meter
    float: azimuth in degrees
    float: tilt in degrees
    float: roll in degrees
    float: focal in pixels
    string: method for initialisation LM or PnP
    """
    method = 'LM'
    # Extract GCPs
    # ###
    # Get lat lng alt gcps
    gcpLatLngAlt = getLatLngAltGcps(gcps)
    # Get images gcps
    gcpXy = getImageGcps(gcps)
    # Center images gcps
    gcpXy = centerImageGcps(gcpXy, width, height)
    # Convert cesium angles to LM angles
    azimuth0, tilt0, roll0 = cesiumToLmAngles(azimuth0, tilt0, roll0)

    # Conversion to ENU
    # ###
    # Conversion ENU0: ENU at the apriori location
    gcpEnu0 = convertEnu(gcpLatLngAlt, lat0, lng0, alt0)

    if plotBool:
        # Plot inital state
        p = [0., 0., 0., azimuth0, tilt0, roll0, focal0, 0, 0]
        xyForward = project3Dto2D(gcpEnu0, p)
        plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
        plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
        plt.axis('square')
        plt.title('Initial state')
        plt.show()

    # Pose estimation with Levenberg-Marquardt
    # ###

    # Location of the image center
    cx = 0
    cy = 0

    # Pose parameters
    pApriori = [0., 0., 0., azimuth0, tilt0, roll0, focal0, cx, cy]
    pBool = [True,True,True,True,True,True,True,False,False]
    maxIterations = 100

    # LM
    isFeasible = True

    try:
        lng1, lat1, alt1, pCompEnu1 = poseLmEnu(lng0, lat0, alt0, pApriori, pBool, maxIterations, gcpXy, gcpEnu0, width, height, plotBool = False)
        gcpEnu1 = convertEnu(gcpLatLngAlt, lat1, lng1, alt1)

        # Compute error on each gcp
        xyComputed = project3Dto2D(gcpEnu1, pCompEnu1)
        gcps, medianError = computeGcpErrors(gcps, xyComputed.T, height, width)

        # Show result
        if plotBool:
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyComputed[0,:], xyComputed[1,:], 'ok')
            plt.title('LM Results in ENU comp')
            plt.axis('square')
            plt.show()

    except:
        isFeasible = False

    if isFeasible == False:
        method = 'PnP'
        # LM didn't converge: try OpenCV PnP
        ###
        # Compute pose with Pnp: angles are provided in ENU0
        EnuOpencv, eulersEnu0, focal = computePoseOpenCv(gcpEnu0, gcpXy, height, width)
        latOpencv, lngOpencv, altOpencv = pm.enu2geodetic(EnuOpencv[0], EnuOpencv[1], EnuOpencv[2], lat0, lng0, alt0, ell=None, deg=True)

        # Convert OpenCv angles to LM angles
        eulersEnu0 = eulersOpencvToLm(eulersEnu0)

        if plotBool:
            # Check OpenCV result
            p = [EnuOpencv[0], EnuOpencv[1], EnuOpencv[2], eulersEnu0[0], eulersEnu0[1], eulersEnu0[2], focal, 0, 0]
            xyForward = project3Dto2D(gcpEnu0,p)
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyForward[0,:], xyForward[1,:], 'og')
            plt.axis('square')
            plt.title ('OpenCV result')
            plt.show()

        # Convert euler angles to ENU at computed location
        Rot_z2 = R.from_euler('z', eulersEnu0[2]%(np.pi*2), degrees=False)
        Rot_x = R.from_euler('x', eulersEnu0[1]%(np.pi*2), degrees=False)
        Rot_z = R.from_euler('z', eulersEnu0[0]%(np.pi*2), degrees=False)
        Rot_enu0 = Rot_z2 * Rot_x * Rot_z
        R_ecef_enu0 = REcefToEnu(lng0*np.pi/180., lat0*np.pi/180.)
        R_ecef_enu1 = REcefToEnu(lngOpencv*np.pi/180., latOpencv*np.pi/180.)
        Rot_enu_opencv =  Rot_enu0 * R.from_matrix(R_ecef_enu0) * R.from_matrix(R_ecef_enu1.T)
        eulers_enu_opencv = Rot_enu_opencv.as_euler('zxz', degrees=False)

        # Convert GCPs to ENU located at the computed location
        gcpEnuOpencv = convertEnu(gcpLatLngAlt, latOpencv, lngOpencv, altOpencv)

        # Launch LM (using OpenCv as apriori value, get more accurate pose and angles in ENU)
        pApriori = [0., 0., 0., eulers_enu_opencv[0], eulers_enu_opencv[1], eulers_enu_opencv[2], focal, 0, 0]
        if plotBool:
            xyForward = project3Dto2D(gcpEnuOpencv,pApriori)
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyForward[0,:], xyForward[1,:], 'og')
            plt.axis('square')
            plt.title ('Initial state after OpenCV')
            plt.show()
        try:
            lng1, lat1, alt1, pCompEnu1 = poseLmEnu(lngOpencv, latOpencv, altOpencv, pApriori, pBool, maxIterations, gcpXy, gcpEnuOpencv, width, height, plotBool = False)
            gcpEnu1 = convertEnu(gcpLatLngAlt, lat1, lng1, alt1)
            # Show result
            if plotBool:
                xyForward = project3Dto2D(gcpEnu1, pCompEnu1)
                plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
                plt.plot(xyForward[0,:], xyForward[1,:], 'ok')
                plt.title('LM Results in ENU comp')
                plt.axis('square')
                plt.show()
        except:
            raise Exception('Georeferencer fails')

    # Extract Cesium angles from euler angles
    azimuthCompEnu1, tiltCompEnu1, rollCompEnu1 = lmToCesiumAngles(pCompEnu1[3], pCompEnu1[4], pCompEnu1[5])

    # Compute error on each gcp
    xyComputed = project3Dto2D(gcpEnu1, pCompEnu1)
    gcps, medianError = computeGcpErrors(gcps, xyComputed.T, height, width)

    # Compute image coordinates
    imageCoordinates = computeImageCoordinates(pCompEnu1, width, height)

    return lng1, lat1, alt1, azimuthCompEnu1%360, tiltCompEnu1%360, rollCompEnu1%360, pCompEnu1[6], pCompEnu1, gcps, imageCoordinates, method

def georeferencerLocked(lng0, lat0, alt0, azimuth0, tilt0, roll0, focal0, width, height,
                  gcps, plotBool = False):

    """Compute camera angles given apriori angle and fixed location

    Parameters:
    lng0 (float): exact longitude in degrees
    lat0 (float): exact latitude in degrees
    alt0 (float): exact altitude in meter
    azimuth0 (float): apriori azimuth in degrees (from cesium)
    tilt0 (float): apriori tilt in degrees (from cesium)
    roll0 (float): apriori roll in degrees (from cesium)
    focal0 (float): apriori focal in pixels
    width (float): image width
    height (float): image height
    gcps (dictionnary): gcps in smapshot format
    plotBool (boolean): if the plot must be shown

    Returns:
    float: longitude in degrees
    float: latitude in degrees
    float: altitude in meter
    float: azimuth in degrees
    float: tilt in degrees
    float: roll in degrees
    float: focal in pixels
    """
    method = 'LM'
    # Extract GCPs
    # ###
    # Get lat lng alt gcps
    gcpLatLngAlt = getLatLngAltGcps(gcps)
    # Get images gcps
    gcpXy = getImageGcps(gcps)
    # Center images gcps
    gcpXy = centerImageGcps(gcpXy, width, height)
    # Convert cesium angles to LM angles
    azimuth0, tilt0, roll0 = cesiumToLmAngles(azimuth0, tilt0, roll0)

    # Conversion to ENU
    # ###
    # Conversion ENU0: ENU at the apriori location
    gcpEnu0 = convertEnu(gcpLatLngAlt, lat0, lng0, alt0)

    if plotBool:
        # Plot inital state
        p = [0., 0., 0., azimuth0, tilt0, roll0, focal0, 0, 0]
        xyForward = project3Dto2D(gcpEnu0, p)
        plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
        plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
        plt.axis('square')
        plt.title('Initial state')
        plt.show()

    # Pose estimation with Levenberg-Marquardt
    # ###

    # Location of the image center
    cx = 0
    cy = 0

    # Pose parameters
    pApriori = [0., 0., 0., azimuth0, tilt0, roll0, focal0, cx, cy]
    pBool = [False,False,False,True,True,True,True,False,False]

    # LM
    try:
        pCompEnu0 = estimatePoseLmfit(pApriori, pBool, gcpXy, gcpEnu0, width, height)
        # Show result
        if plotBool:
            xyForward = project3Dto2D(gcpEnu0, pCompEnu0)
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyForward[0,:], xyForward[1,:], 'ok')
            plt.title('LM Results in ENU fixed')
            plt.axis('square')
            plt.show()

    except:
        raise Exception('Georeferencer fails')

    # Extract Cesium angles from euler angles
    azimuthCompEnu0, tiltCompEnu0, rollCompEnu0 = lmToCesiumAngles(pCompEnu0[3], pCompEnu0[4], pCompEnu0[5])

    # Compute error on each gcp
    xyComputed = project3Dto2D(gcpEnu0, pCompEnu0)
    gcps, medianError = computeGcpErrors(gcps, xyComputed.T, height, width)

    # Compute image coordinates
    imageCoordinates = computeImageCoordinates(pCompEnu0, width, height)

    return lng0, lat0, alt0, azimuthCompEnu0%360, tiltCompEnu0%360, rollCompEnu0%360, pCompEnu0[6], pCompEnu0.tolist(), gcps, imageCoordinates, method

def poseLmEnu(lng0, lat0, alt0, pAprioriEnu0, pBool, maxIt, gcpXy, gcpEnu0, imageWidth, imageHeight, plotBool = False):
    """Pose estimation with LM given an apriori ENU coordinate system, returns
    the pose in the computed ENU.

    Parameters:
    lng0 (float): apriori longitude in degrees
    lat0 (float): apriori longitude in degrees
    alt0 (float): apriori altitude in meter
    pAprioriEnu0 (array): apriori pose parameters
    pBool (array): pose parameters to be computed
    maxIt (integer): maximal number of LM iterations
    gcpXy (matrix): image coordinates
    gcpEnu0 (matrix): world coordinates in apriori ENU
    imageWidth (float): image width
    imageHeight (float): image height
    plotBool (boolean): if the plot must be shown


    Returns:
    float: longitude in degrees
    float: latitude in degrees
    float: altitude in meter
    array: pose parameters
    """

    failMessage = "LM don't converge"
    isFeasible = True
    try:
        pCompEnu0 = estimatePoseLmfit(pAprioriEnu0, pBool, gcpXy, gcpEnu0, imageWidth, imageHeight)
    except:
        raise Exception(failMessage)

    isFeasible = poseFeasible(pAprioriEnu0[0:3], pCompEnu0[0:3], pCompEnu0[6], imageWidth, imageHeight)

    if isFeasible:

        lat1, lng1, alt1 = pm.enu2geodetic(pCompEnu0[0], pCompEnu0[1], pCompEnu0[2], lat0, lng0, alt0, ell=None, deg=True)

        # Convert eulers angle to rotation matrix
        Rot_z2 = R.from_euler('z', pCompEnu0[5], degrees=False)
        Rot_x = R.from_euler('x', pCompEnu0[4], degrees=False)
        Rot_z = R.from_euler('z', pCompEnu0[3], degrees=False)
        Rot_enu0 = Rot_z2 * Rot_x * Rot_z

        # Show computed result in apriori coordinate system and check euler extraction
        if plotBool:
            eulers_enu0 = Rot_enu0.as_euler('zxz', degrees=False)
            pCheck = [pCompEnu0[0], pCompEnu0[1], pCompEnu0[2], eulers_enu0[0], eulers_enu0[1], eulers_enu0[2], pCompEnu0[6], 0, 0]
            xyForward = project3Dto2D(gcpEnu0, pCheck)
            plt.plot(gcpXy[:,0], gcpXy[:,1], 'or')
            plt.plot(xyForward[0,:], xyForward[1,:], 'ob')
            plt.title('LM Results in ENU apriori')
            plt.axis('square')
            plt.show()


        # Convert computed angles in ENU computed (enu1)
        R_ecef_enu0 = REcefToEnu(lng0*np.pi/180., lat0*np.pi/180.)
        R_ecef_enu1 = REcefToEnu(lng1*np.pi/180., lat1*np.pi/180.)
        REnu1 =  Rot_enu0 * R.from_matrix(R_ecef_enu0) * R.from_matrix(R_ecef_enu1.T)
        eulersEnu1 = REnu1.as_euler('zxz', degrees=False)
        pCompEnu1 = [0, 0, 0, eulersEnu1[0], eulersEnu1[1], eulersEnu1[2], pCompEnu0[6], 0, 0]

    else:
        raise Exception(failMessage)

    return  lng1, lat1, alt1, pCompEnu1

def createImageMains (focal, cu, cv):
    #Create matrix of mains image points
    imageMains3D = np.zeros((3,6))

    #Create four image corners
    imageMains3D[:,0] = (-cu,-cv,focal)
    imageMains3D[:,1] = (cu,-cv,focal)
    imageMains3D[:,2] = (-cu,cv,focal)
    imageMains3D[:,3] = (cu,cv,focal)

    #Create principal point
    imageMains3D[:,4] = (0, 0, focal)

    #Create center of projection
    imageMains3D[:,5] = (0, 0, 0)

    return imageMains3D

def computeImageCoordinatesForGltf (p, width, height):
    """Compute image corner and projection center coordinates for GLTF

    Parameters:
    p (array): pose parameters
    width (float): image width
    height (float): image height

    Returns:
    array: image points in world coordinates
    """
    # Related to this issue the image coordinates in the gltf must be rotated
    # https://github.com/CesiumGS/cesium/issues/6713
    p[3] = p[3] + np.pi/2
    # Not sure why
    p[4] = p[4] + np.pi
    p[5] = -p[5]

    # Create images points in camera coordinates
    imagePointsCamera = createImageMains(p[6], height/2, width/2).T
    # Convert to ENU
    imagePointsWorld = camera2world(imagePointsCamera, p)

    imageCoordinates = {}
    imageCoordinates['ul'] = imagePointsWorld.T[0,:].tolist()
    imageCoordinates['ll'] = imagePointsWorld.T[1,:].tolist()
    imageCoordinates['ur'] = imagePointsWorld.T[2,:].tolist()
    imageCoordinates['lr'] = imagePointsWorld.T[3,:].tolist()

    return imageCoordinates

def computeImageCoordinates (p, width, height):
    """Compute image corner and projection center coordinates

    Parameters:
    p (array): pose parameters
    width (float): image width
    height (float): image height

    Returns:
    array: image points in world coordinates
    """
    # Create images points in camera coordinates
    imagePointsCamera = createImageMains(p[6], width/2, height/2).T
    # Convert to ENU
    imagePointsWorld = camera2world(imagePointsCamera, p)

    imageCoordinates = {}
    imageCoordinates['ur'] = imagePointsWorld.T[0,:].tolist()
    imageCoordinates['ul'] = imagePointsWorld.T[1,:].tolist()
    imageCoordinates['lr'] = imagePointsWorld.T[2,:].tolist()
    imageCoordinates['ll'] = imagePointsWorld.T[3,:].tolist()

    return imageCoordinates

def getGeolocalisationFromDB(dbParams, imageId):
    """Get geolocalisation data of a image in the database

    Parameters:
    dbParams (dict): Connection parameters
    imageId (int): image identifier

    Returns:
    float: lng (degree)
    float: lat (degree)
    float: alt (degree)
    float: azimuth (degree)
    float: tilt (degree)
    float: roll (degree)
    float: focal (pixel)
    array: gcps
    int: image width
    int: image height

   """

    conn = psycopg2.connect(
        host=dbParams['host'],
        database=dbParams['database'],
        user=dbParams['user'],
        password=dbParams['password'],
        port=dbParams['port'])

    # create a cursor
    cur = conn.cursor()

    # execute a statement
    cur.execute("""
    select st_x(geolocalisations.location), st_y(geolocalisations.location), 
    st_z(geolocalisations.location),
    geolocalisations.azimuth, geolocalisations.tilt, geolocalisations.roll, 
    geolocalisations.focal,
    gcp_json, width, height, images.id from images
    left join geolocalisations on images.geolocalisation_id = geolocalisations.id
    where images.id = {}""".format(imageId)) #

    resall = cur.fetchall()
    res = resall[0]
    lng = float(res[0])
    lat = float(res[1])
    alt = float(res[2])
    azimuth = float(res[3])
    tilt= float(res[4])
    roll = float(res[5])
    focal = float(res[6])
    gcps = res[7]
    width = float(res[8])
    height = float(res[9])

    return lng, lat, alt, azimuth, tilt, roll, focal, gcps, width, height

def getGeolocalisationFromJson(image):
    """Get geolocalisation data of a image in the database

    Parameters:
    image (json): image from the validation data

    Returns:
    float: lng (degree)
    float: lat (degree)
    float: alt (degree)
    float: azimuth (degree)
    float: tilt (degree)
    float: roll (degree)
    float: focal (pixel)
    array: gcps
    int: image width
    int: image height

   """

    lng = float(image['lng'])
    lat = float(image['lat'])
    alt = float(image['alt'])
    azimuth = float(image['azimuth'])%360
    tilt= float(image['tilt'])%360
    roll = float(image['roll'])%360
    focal = float(image['focal'])
    gcps = image['gcp_json']
    width = float(image['width'])
    height = float(image['height'])

    return lng, lat, alt, azimuth, tilt, roll, focal, gcps, width, height

def randomInRange(min, max):

    return min + (random() * (max - min))

def REcefToEnu(lamba, phi):
    """Returns the rotation matrix from ecef to enu
    https://gis.stackexchange.com/questions/82998/trasformation-from-ecef-to-enu
    http://www.hydrometronics.com/downloads/Ellipsoidal%20Orthographic%20Projection.pdf

    Parameters:
    lamba (float): longitude in radian
    phi (float): latitude in radian

    Returns:
    matrix: rotation matrix

    """
    R = np.asarray([
            [-np.sin(lamba), np.cos(lamba), 0],
            [-np.sin(phi)*np.cos(lamba), -np.sin(phi)*np.sin(lamba), np.cos(phi)],
            [np.cos(phi)*np.cos(lamba), np.cos(phi)*np.sin(lamba), np.sin(phi)]
        ])

    return R

def generateSimulatedApriori(lng, lat, alt, azimuth, tilt, roll, focal):
    """Returns provided pose parameter with a small perturbation

    Parameters:
    lng (float): longitude in degrees
    lat (float): longitude in degrees
    alt (float): altitude in meter
    azimuth (float): azimuth in degrees
    tilt (float): tilt in degrees
    roll (float): roll in degrees
    focal (float): focal in pixels

    Returns:
    float: longitude in degrees
    float: longitude in degrees
    float: altitude in meter
    float: azimuth in degrees
    float: tilt in degrees
    float: roll in degrees
    float: focal in pixels

   """

    lng = lng + randomInRange(-0.02, 0.02)
    lat = lat + randomInRange(-0.02, 0.02)
    alt = alt + randomInRange(-500, 500)
    focal = focal + focal * randomInRange(-0.1, 0.1)
    azimuth = azimuth + randomInRange(-5, 5)
    tilt = tilt + randomInRange(-2, 2)
    roll = roll + randomInRange(-2, 2)

    return lng, lat, alt, azimuth, tilt, roll, focal


def getImageGcps(gcps):
    """
    Extract the gcps from the dict
    """
    gcp_xy = []
    for gcp in gcps:
        gcp_xy.append([float(gcp['x']), float(gcp['y'])])
    # List to numpy array
    gcp_xy = np.asarray(gcp_xy)
    return gcp_xy

def getLatLngAltGcps(gcps):
    """
    Extract the gcps from the database json
    """
    gcpLatLngAlt = []
    for gcp in gcps:
        gcpLatLngAlt.append([float(gcp['latitude']), float(gcp['longitude']), float(gcp['altitude'])])
    # List to numpy array
    gcpLatLngAlt = np.asarray(gcpLatLngAlt)
    return gcpLatLngAlt

def convertEnu(gcpLatLngAlt, lat, lng, alt):
    """Convert Lat Lng Alt to East North Up

    Parameters:
    gcpLatLngAlt (matrix): nx3 matrix with Lat, Lng, Height coordinates
    lat (float): latitude of the ENU origin in degrees
    lng (float): longitude of the ENU origin in degrees
    alt (float): altitude of the ENU origin in meter

    Returns:
    matrix: nx3 matrix with E, N, U coordinates

   """
    gcpEnu =  []
    for i in range(gcpLatLngAlt.shape[0]):
        E0, N0, U0 = pm.geodetic2enu(gcpLatLngAlt[i,0], gcpLatLngAlt[i,1], gcpLatLngAlt[i,2], lat, lng, alt, ell=None, deg=True)
        gcpEnu.append([E0, N0, U0])
    gcpEnu = np.asarray(gcpEnu)
    return gcpEnu

def centerImageGcps(gcps, width, height):
    """Put the origin of the images GCPs at the image center

    Parameters:
    gcps (matrix): nx2 matrix with x, y in pixel
    width (float): image width in pixel
    height (float): image height in pixel

    Returns:
    matrix: nx2 matrix with x, y in pixel

   """
    cx = width/2
    cy = height/2
    gcps = gcps-np.asarray([cx, cy])
    return gcps

def centerXYZ(XYZ):

    cg_XYZ = np.mean(XYZ,0)
    XYZ = np.subtract(XYZ, cg_XYZ)

    return (XYZ, cg_XYZ)

def createCameraMatrix(fx,fy,cx,cy):
    #Create matrix of internal orientation
    cameraMatrix = np.array([[fx,0,cx],[0,fy, cy],[0,0,1]])
    return cameraMatrix

def computePoseOpenCv (GcpXYZ, Gcpxy, width, height):

    """Compute pose with OpenCv:

    Parameters:
    GcpXYZ (array nx3): GCP in world coordinates
    Gcpxy(array nx2): GCP in image coordinates
    height (float): image height
    width (float): image width

    Returns:
    XYZ (array): XYZ coordinate of the camera
    eulers (array): three eulers angle (radian)
    focal (float): focal in pixel

    """
    # Center XYZ coordinates around the gravity center
    XYZ, cgXYZ = centerXYZ(GcpXYZ)
    XYZ = GcpXYZ

    # Initialise the error to find the minimum
    minError = np.Inf
    # Loop on the possible half angle of view
    for angle in range(10,70,2):

        # Compute corresponding focal and camera matrix
        diag = computeDiagonal(width, height)
        diag = diag/2
        focal = diag/np.tan(angle/180*np.pi)
        cameraMatrix = createCameraMatrix(focal, focal, 0, 0)

        # Pose estimation with OpenCV
        N, M = Gcpxy.shape
        imagePoints = np.ascontiguousarray(Gcpxy[:,:2]).reshape((N,1,2)) # Trick to avoid opencv error with some pnp methods
        retval, rvec, tvec = cv2.solvePnP(XYZ, imagePoints,cameraMatrix, None, None, None, False, cv2.SOLVEPNP_EPNP) #cv2.SOLVEPNP_UPNP

        # Compute error
        xyComputed, jacobian = cv2.projectPoints(XYZ, rvec, tvec, cameraMatrix, np.zeros(4))
        delta_xy = np.subtract(Gcpxy, xyComputed.squeeze())
        delta_x2 = np.multiply(delta_xy[:,0],delta_xy[:,0])
        delta_y2 = np.multiply(delta_xy[:,1],delta_xy[:,1])
        #delta = np.add(delta_x2, delta_y2)
        error = np.sum(np.add(delta_x2, delta_y2))

        # Check if the current focal provide better results
        if error < minError:
            #bestAngle = angle
            bestFocal = focal
            minError = error
            best_rvec = rvec
            best_tvec = tvec
            #best_xyComputed = xyComputed

    # Get best parameters
    rvec = best_rvec
    tvec = best_tvec
    focal = bestFocal

    ## Extract rotation matrix
    rmat = cv2.Rodrigues(rvec)[0]
    r = R.from_matrix(rmat.T)
    eulers = r.as_euler('zxz', degrees=False)

    ## Extract world pose
    Rt = cv2.Rodrigues(rvec)[0]
    Rmat = Rt.transpose()
    pos = -np.matrix(Rmat)*np.matrix(tvec)
    pos = np.squeeze(pos)
    r = R.from_matrix(Rmat)
    eulers = r.as_euler('zxz', degrees=False)

    return  [pos[0,0], pos[0,1], pos[0,2]], eulers, focal

def eulersOpencvToLm(eulers):
    """Convert eulers angles extracted with OpenCv to Eulers angles used by the
    LM algorithm

    Parameters:
    eulers (array): tree eulers angles

    Returns:
    array: tree eulers angles

   """

    return [np.pi-eulers[2], eulers[1]+np.pi, eulers[0]-np.pi/2]

def camera2image(xyz, p):
    """Project camera coordinates in the image plane (perspective transform).

    Parameters:
    xyz (matrix 3xn): n points in camera coordinates
    p (array or Parameters): 8 pose parameters

    Returns:
    xy (array 2xn): n points in image coordinates

    """
    # Get internal parameters
    if (type(p) is np.ndarray) or (type(p) is list):
        f = p[6]
        cx = p[7]
        cy = p[8]
    if type(p) is Parameters:
        f = p['f'].value
        cx = p['cx'].value
        cy = p['cy'].value

    # Perspective transform
    xs = xyz[0,:]
    ys = xyz[1,:]
    zs = xyz[2,:]

    x = -f*ys/zs
    y = -f*xs/zs

    xoff = np.add(x, cx)
    yoff = np.add(y, cy)

    return np.asarray([xoff, yoff])

def world2camera(XYZ, p):

    """Transform world coordinates in camera coordinates.

    Parameters:
    XYZ (matrix 3xn): n points in world coordinates
    p (array or Parameters): 8 pose parameters

    Returns:
    xyz (array 2xn): n points in camera coordinates

    """
    # Rigid body transformation
    if (type(p) is np.ndarray) or (type(p) is list):
        X = p[0]
        Y = p[1]
        Z = p[2]
        alpha = p[3]
        beta = p[4]
        gamma = p[5]
    if type(p) is Parameters:
        X = p['X'].value
        Y = p['Y'].value
        Z = p['Z'].value
        alpha = p['alpha'].value,
        beta = p['beta'].value
        gamma = p['gamma'].value

    # Compute ZYZ rotation matrix
    Rot_z2 = R.from_euler('z', gamma, degrees=False)
    Rot_x = R.from_euler('x', beta, degrees=False)
    Rot_z = R.from_euler('z', alpha, degrees=False)
    Rot = Rot_z2 * Rot_x * Rot_z
    Rot = Rot.as_matrix()

    # Vector camera DEM
    T = np.asarray([X, Y, Z])

    # Translate
    XYZc = np.subtract(XYZ, T)

    # Rotate according to rotation matrix
    xyz = np.dot(Rot, np.transpose(XYZc))
    return np.squeeze(xyz)

def camera2world(XYZ, x):

    """Transform world coordinates in camera coordinates.

    Parameters:
    XYZ (matrix 3xn): n points in world coordinates
    x (array): 8 pose parameters

    Returns:
    xyz (array 2xn): n points in camera coordinates

    """
    # Rigid body transformation

    # Compute ZYZ rotation matrix
    Rot_z2 = R.from_euler('z', x[5], degrees=False)
    Rot_x = R.from_euler('x', x[4], degrees=False)
    Rot_z = R.from_euler('z', x[3], degrees=False)
    Rot = Rot_z2 * Rot_x * Rot_z
    Rot = Rot.as_matrix().T
    xyz = np.dot(Rot, np.transpose(XYZ))
    return xyz

def project3Dto2D(XYZ, p):
    """Project world coordinates in image plane.

    Parameters:
    XYZ (matrix 3xn): n points in world coordinates
    p (array): 8 pose parameters

    Returns:
    xy (array 2xn): n points in camera coordinates

    """

    # Rigid body transfomation world -> camera
    xyz = world2camera(XYZ, p)

    #perspective transform
    xy = camera2image(xyz, p)

    return xy

def azCesiumToLm (azimuth):
    """Transform azimuth provided by cesium to azimuth needed by LM.

    Parameters:
    azimuth (float): angle in degree north is null

    Returns:
    azimuth (float): angle in degree east is null

    """
    if (azimuth < 0):
        azimuth = azimuth + 360

    azimuth = azimuth % 360

    if (azimuth >= 0.0) and (azimuth <= 90.0):
        az = -180 + azimuth
    elif (azimuth > 90.0) and (azimuth <= 180.0):
        az = azimuth - 180
    elif (azimuth > 180.0) and (azimuth <= 270.0):
        az = azimuth + 180
    elif (azimuth > 270.0) and (azimuth <= 360.0):
        az = azimuth -180
    return az

def tiltCesiumToLm (tilt):
    return 90 + tilt

def rollCesiumToLm (roll):
    return roll - 90

def azLmToCesium(azimuth):
    return azimuth + 180

def tiltLmToCesium (tilt):
    return tilt -90

def rollLmToCesium(roll):
    return 90 + roll

def cesiumToLmAngles (az, tilt, roll):
    """Transform angles provided by cesium to angles needed by LM.

    Parameters:
    azimuth (float): angle in degree north is null
    tilt (float): angle in degree horizontal is null
    roll (float): angle in degree

    Returns:
    azimuth (float): angle in rad
    tilt (float): angle in rad
    roll (float): angle in rad

    """
    az = deg2rad(azCesiumToLm(az))
    tilt = deg2rad(tiltCesiumToLm(tilt))
    roll = deg2rad(rollCesiumToLm(roll))
    return az, tilt, roll

def deg2rad(a):
    return np.pi*a/180

def rad2deg(a):
    return a*180/np.pi

def lmToCesiumAngles(az, tilt, roll):
    """Convert eulers angles (computed with LM) to angles required by Cesium.

    Parameters:
    azimuth (float): azimuth (radian)
    tilt (float): tilt (radian)
    roll (float): roll (radian)

    Returns:
    azimuth (float): azimuth (degree)
    tilt (float): tilt (degree)
    roll (float): roll (degree)

    """
    az = azLmToCesium(rad2deg(az))%360
    tilt = tiltLmToCesium(rad2deg(tilt))%360
    roll = rollLmToCesium(rad2deg(roll))%360
    return az, tilt, roll


def residual(params, xy, XYZ):
    """Function which compares the measurements with the model
    params (Parameters): 8 perspective parameters
    xy (array nx2): image coordinates GCP
    XYZ (array nx3): world coordinates GCP

    Returns:
    r (array nx1): residuals
    """
    xyBulle = project3Dto2D(XYZ, params)
    xBulle = xyBulle[0, :]
    yBulle = xyBulle[1, :]
    nObs = xBulle.size

    # Computation of the residuals
    dx = np.subtract(xy[:,0], np.transpose(xBulle))
    dy = np.subtract(xy[:,1], np.transpose(yBulle))

    # 1D vector of the residual
    r = []
    for i in range(nObs):
        r.append(dx[i])
        r.append(dy[i])
    r = np.asarray(r)

    return r

def estimatePoseLmfit(p, pBool, xy, XYZ, imageWidth, imageHeight):

    """Estimate pose with Levenberg-Marquardt algorithm using lmfit library.

    Parameters:
    p (array): 8 pose parameters
    pBool (array): 8 boolean True if the parameter must be computed
    xy (array nx2): n images gcp coordinates
    XYZ (array nx3): n world gcp coordinates

    Returns:
    p (array): (8 optimized pose parameters

    """

    params = Parameters()
    params.add('X', value=p[0], vary=pBool[0])
    params.add('Y', value=p[1], vary=pBool[1])
    params.add('Z', value=p[2], vary=pBool[2])
    params.add('alpha', value=p[3], vary=pBool[3])
    params.add('beta', value=p[4], vary=pBool[4])
    params.add('gamma', value=p[5], vary=pBool[5])
    # focal with min and max constraints
    halfDiag = computeDiagonal(imageWidth, imageHeight)/2
    minFocal = halfDiag/np.tan(70/180*np.pi) # fov is 140 degrees
    maxFocal = halfDiag/np.tan(5/180*np.pi) # fov is 10 degrees
    params.add('f', value=p[6], vary=pBool[6], min=minFocal, max=maxFocal)
    params.add('cx', value=p[7], vary=pBool[7])
    params.add('cy', value=p[8], vary=pBool[8])

    r = minimize(residual, params, args=(xy, XYZ), max_nfev= 1000)
    
    if r.success == False:
        raise Exception('Minimization fails')

    return np.asarray([r.params['X'].value, r.params['Y'].value,
                       r.params['Z'].value, r.params['alpha'].value,
                       r.params['beta'].value, r.params['gamma'].value,
                       r.params['f'].value, r.params['cx'].value,
                       r.params['cy'].value])

def poseFeasible(cameraAprioriXYZ, cameraXYZ, focal, height, width):

    """Check if the pose computed is feasible:
        * if the location is not too far (30km)
        * if the field of view is neither to narrow nor to wide

    Parameters:
    cameraAprioriXYZ(array): apriori image location
    cameraXYZ(array): computed image location
    focal (float): focal in pixel
    height (float): image height
    width (float): image width

    Returns:
    isFeasible (boolean): if the pose is feasible

    """

    isFeasible = False

    # Check fov
    fovValid = False
    imSide = np.max([width, height])
    # Compute fov in degrees
    fov = 2 * np.arctan(imSide / 2 / focal) * 180 / np.pi
    if ((fov > 10) and (fov < 140)):
        fovValid = True

    # Check distance
    distValid = False
    dX = cameraAprioriXYZ[0] - cameraXYZ[0]
    dY = cameraAprioriXYZ[1] - cameraXYZ[1]
    dZ = cameraAprioriXYZ[2] - cameraXYZ[2]
    dist = np.sqrt(dX * dX + dY * dY + dZ * dZ)
    if (dist < 30000):
        distValid = True

    if ((distValid) and (fovValid)):
        isFeasible = True

    return isFeasible

def generateCollada(p, colladaTemplate, colladaOutput, texture):
    """" Generate a collada file from pose parameters
    This function is a python translation of the javascript function
    It must be validated

    Parameters:
    p (array): 8 pose parameters
    colladaTemplate (string): path to the collada template
    colladaOutput (string): path to the collada output
    texture (string): ex: ../images/1024/{image_id}.jpg

    """

    # change the roll (sign + 90): there is a problem with the roll definition
    #p[5] = p[5] + np.pi/2

    # Compute image coordinate in world coordinates
    urCorner = [p[7], p[8], p[6]]
    ulCorner = [-p[7], p[8], p[6]]
    lrCorner = [p[7], -p[8], p[6]]
    llCorner = [-p[7], -p[8], p[6]]

    xyzCam = np.matrix([urCorner, ulCorner, lrCorner, llCorner])

    # Transform the image plan from camera coordinates to world coordinates
    XYZCam = camera2world(xyzCam, p);

    # Substract the camera location to get a 3D model close to the origin
    # -------------------------------------------------------------------
    # Vector camera DEM
    # If the pose is computed in ENU T is [0,0,0]
    T = p[0:3]
    XYZCamOff = np.subtract(XYZCam, T);



    # Scale the model
    # ---------------
    urCorner = XYZCamOff[:,0]
    ulCorner = XYZCamOff[:,1]
    delta = np.subtract(urCorner, ulCorner)
    dist = np.norm(delta)

    # max size of the biggest side of the image
    maxSize = 100
    ratio = maxSize / dist
    scaledXYZ = np.multiply(XYZCamOff, ratio)

    # Insert the computed coordinates in the collada template file
    # ------------------------------------------------------------
    # create the coordinate string for collada
    coordString = ""
    dim0 = XYZCamOff.shape[0]
    dim1 = XYZCamOff.shape[1]
    for i in range(dim1):
        for j in range(dim0):
            curCoord = "{:.1f}".format(scaledXYZ[j,i])
            coordString = coordString + curCoord + " "

    templateFile = open(colladaTemplate, 'r')
    xml = templateFile.read()
    templateFile.close()

    xml = templateFile.read()
    newXml = xml.replace("#IMAGECOORDINATES#", coordString)
    newXml = newXml.replace("#PATH2IMAGE#", texture)

    outputFile = open(colladaOutput, 'w')
    outputFile.write(newXml)
    outputFile.close()

    return

def computeDiagonal(size_x, size_y):
    #Compute the normal focal => diagonal = focal
    focal = np.sqrt(size_x*size_x + size_y*size_y)
    return focal

def computeGcpErrors(gcps, xyComputed, height, width):
    # Computation of the residuals
    # diagonal length
    diag = np.sqrt(height*height+width*width)
    errors = []
    i=0
    for gcp in gcps:
        dx = gcp['x'] - (xyComputed[i,0]+width/2);
        dy = gcp['y'] - (xyComputed[i,1]+height/2);
        dxy = np.sqrt(dx*dx+dy*dy);
        gcp['dxy'] = dxy
        gcp['errorPct'] = (dxy / diag) * 100;
        gcp['xReproj'] = xyComputed[i,0]+width/2
        gcp['yReproj'] = xyComputed[i,1]+height/2
        errors.append((dxy / diag) * 100)
        i+=1

    medianError = np.median(errors)
    return gcps, medianError

def RotZ(alpha):

    rotMat = [
      [np.cos(alpha), -np.sin(alpha), 0.0],
      [np.sin(alpha), np.cos(alpha), 0.0],
      [0.0, 0.0, 1.0]
    ]
    return rotMat;
