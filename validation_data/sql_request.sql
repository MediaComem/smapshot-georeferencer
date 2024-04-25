SELECT
  json_agg(foo)
FROM
  (
    (
      SELECT
        images.id,
        collection_id,
        st_x(images.location) AS lng,
        st_y(images.location) AS lat,
        st_z(images.location) AS alt,
        images.azimuth % 360 AS azimuth,
        images.tilt % 360 AS tilt,
        images.roll % 360 AS roll,
        images.focal,
        height,
        width,
        gcp_json,
        view_type
      FROM
        images
        LEFT JOIN geolocalisations ON images.geolocalisation_id = geolocalisations.id
      WHERE
        images.state = 'validated'
        AND view_type LIKE 'nadir'
      ORDER BY
        random()
      LIMIT
        100
    )
  ) AS foo


SELECT
  json_agg(foo)
FROM
  (
    (
      SELECT
        images.id, 
        collection_id, 
        st_x(images.location) AS lng,
        st_y(images.location) AS lat,
        st_z(images.location) AS alt,
        images.azimuth % 360 AS azimuth,
        images.tilt % 360 AS tilt,
        images.roll % 360 AS roll,
        images.focal,
        height,
        width,
        gcp_json,
        view_type
      FROM
        images
        LEFT JOIN geolocalisations ON images.geolocalisation_id = geolocalisations.id
      WHERE
        images.state = 'validated'
        AND view_type NOT IN ('nadir', 'terrestrial')
      ODER BY
        random()
      LIMIT
        100
    )
  ) AS foo
