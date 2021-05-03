select json_agg(foo) from(
	(select images.id, collection_id, st_x(images.location) as lng, st_y(images.location) as lat , st_z(images.location) as alt, 
	images.azimuth%360 as azimuth, images.tilt%360 as tilt, images.roll%360 as roll, 
	images.focal, height, width, gcp_json, view_type
	from images 
	left join geolocalisations on images.geolocalisation_id = geolocalisations.id
	where images.state = 'validated' and view_type like 'nadir'
	order by random()
	limit 100)
) as foo

select json_agg(foo) from(
	(select images.id, collection_id, st_x(images.location) as lng, st_y(images.location) as lat , st_z(images.location) as alt, 
	images.azimuth%360 as azimuth, images.tilt%360 as tilt, images.roll%360 as roll, 
	images.focal, height, width, gcp_json, view_type
	from images 
	left join geolocalisations on images.geolocalisation_id = geolocalisations.id
	where images.state = 'validated' and view_type not in ('nadir', 'terrestrial')
	order by random()
	limit 100)
) as foo
