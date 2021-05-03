# Python georeferencer
The georeferencer is used in the project smapshot

## 1. Dependencies
* Pyton OpenCv
* Numpy
* scipy
* pymap3d
* lmfit
* matplotlib (only for testing / developpment)
* psycopg2 (only for testing / developpment)

## 2. General Workflow
* The front-end send to the back-end the current Cesium pose and the GCPs. The Cesium pose is the apriori pose.
* The API send these values to `georeferencer.py`
* The georeferencer returns the pose and the 3D image coordinates
* Collada (XML) file is generated (API)
* Collada is converted to GLTF (API)

## 3. Georeferencer algorithm
* The GCPs are converted in ENU (East-North-Up) system located at the apriori pose => ENU_0
* Levenberg-Marquardt minimise GCP errors and optimise pose in ENU_0 apriori 
* If LM don’t converge, the apriori pose is computed with OpenCV PnP 
* The computed pose define ENU_1
* The angles obtained in ENU_0 are converted to ENU_1
* The image coordinates (3D coordinates in ENU_1 of the image corners) are computed

## 4. Notes about the angles conversions.

### Cesium angles
https://cesium.com/docs/cesiumjs-ref-doc/HeadingPitchRoll.html

* heading: Heading is the rotation about the negative z axis. Heading is the rotation from the local north direction where a positive angle is increasing eastward.
* pitch: Pitch is the rotation about the negative y axis. Pitch is the rotation from the local east-north plane. Positive pitch angles are above the plane. Negative pitch angles are below the plane.
* roll: Roll is the rotation about the positive x axis. The roll is applied around the viewing direction. Positive value are turning around the viewing direction to the right.

where x, y, z is a local ENU coordinate system. z is up, y is East, x is North. The rotation around z is applied, the rotation around y'  is applied (right direction) and the rotation around x'' is applied (viewing direction)

### LM angles
LM angles are in radians, the ZXZ convention is applied: turn around Z to fix the heading, turn around X' to fix the tilt, turn around Z '' (viewing direction to fix the roll)

* azimuth 0 is East, positive values are going in the clockwise direction
* tilt 0 is Down, positive value are going, up is 2*pi
* roll 0 is Left, positive values are going in the anti clockwise direction

### Conversions
* The conversion from LM angles to cesium angles is done with `lmToCesiumAngles`
* The conversion from Cesium angles to LM angles is done with `cesiumToLmAngles`


## 5. Algorithm testing

### Simulate a call by the api

Locked location:
`python3 georeferencer.py -43.1566 -22.9495 396.7 248.64 0 0 '[{"X": 4285809.93609858,"Y": -4019884.1671358882,"Z": -2472208.9449742106,"altitude": 0.9068280845509122,"latitude": -22.955758337032876,"longitude": -43.16617119330748,"x": 2508.9841481669496,"y": 5424.2640760774075},{"X": 4279908.715931221,"Y": -4024656.385306744,"Z": -2475976.549851101,"altitude": 518.8414615478749,"latitude": -22.990725967893944,"longitude": -43.23948898145271,"x": 3564.422586858862,"y": 2109.774232471245},{"X": 4283313.735196287,"Y": -4023692.9425000357,"Z": -2472078.3600718803,"altitude": 674.3865588361209,"latitude": -22.951902340734442,"longitude": -43.209904490510844,"x": 8092.207424733941,"y": 1505.2004316663856}, {"X": 4284942.698891078,"Y": -4021241.690559589,"Z": -2471509.602966951,"altitude": 1.1019456558963059,"latitude": -22.94889966745083,"longitude": -43.18161010520996,"x": 8785.813359371354,"y": 3627.1173920537162}]' 9483 7096 1`

Free location:
`python3 georeferencer.py -43.1566 -22.9495 396.7 248.64 0 0 '[{"X": 4285809.93609858,"Y": -4019884.1671358882,"Z": -2472208.9449742106,"altitude": 0.9068280845509122,"latitude": -22.955758337032876,"longitude": -43.16617119330748,"x": 2508.9841481669496,"y": 5424.2640760774075},{"X": 4279908.715931221,"Y": -4024656.385306744,"Z": -2475976.549851101,"altitude": 518.8414615478749,"latitude": -22.990725967893944,"longitude": -43.23948898145271,"x": 3564.422586858862,"y": 2109.774232471245},{"X": 4283313.735196287,"Y": -4023692.9425000357,"Z": -2472078.3600718803,"altitude": 674.3865588361209,"latitude": -22.951902340734442,"longitude": -43.209904490510844,"x": 8092.207424733941,"y": 1505.2004316663856}, {"X": 4284942.698891078,"Y": -4021241.690559589,"Z": -2471509.602966951,"altitude": 1.1019456558963059,"latitude": -22.94889966745083,"longitude": -43.18161010520996,"x": 8785.813359371354,"y": 3627.1173920537162}]' 9483 7096 0`

### Free location
There are two validation datasets:

*  `/validation_data/imagesOblique.json`:  A set of 100 oblique images

* `/validation_data/imagesNadir.json`: A set of 100 nadir images

which can be tested with the script : `test_georeferencer.py`. 

* The pose recorded in the validation datasets will be slightly perturbed to generate the apriori pose. 
*  If the variable plot is set to True, the clicked GCPs and computed reprojections are shown. The GCPs should overlap precisely.
* If the computed location is far (100m) from the database location, the computed values are printed in the console. This is when the computed location is far from the database location, that the georeferencer performs either better or worse than the georeferencer used to generate the values stored in the database

### Fixed location
The validation datasets `/validation_data/imagesOblique.json` can be tested with the script : `test_georeferencer_locked.py`. 

* The location recorded in the validation datasets is fixed
* The angles are slightly pertubated
*  If the variable plot is set to True, the clicked GCPs and computed reprojections are shown. The GCPs should overlap precisely.
* If the computed angles have a difference of more than 1 degree with the database angles, the computed values are printed in the console.

## 6. Warning
The pose values in the validation datasets are the values extracted from smapshot database. Hence, they can't be used as ground truth. They were computed with various versions of the georeferencer and computed from GCPs having some errors.

For instance, when computing the pose of nadir (top-down) images the focal is highly correlated with the altitude. Hence some altitudes and focal lenghts in the database might be wrong.

The algorithm returns three euler angles to be applied in ZXZ order. As a remainder, different values of euler angles once combined can bring to the same result. Hence it is possible that the angles in the validation datasets don't match the computed angles but both results are correct.