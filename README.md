Confirmed running on Mac!
Simply run "make" in the main directory.
To run, enter ".\as2 -i [inputfile] -o [outputfile]". Replace the brackets with the relevant files.

To set the recursion depth for reflection, add the flag "-d x", where x is the max recursion limit, i.e. "-d 2".
Default is 1.

With the "-ambient" tag, the ray tracer will calculate the ambient contribution from point and directional lights.
Without the tag, the ray tracer will only use the ambient light sources added to the scene.

Be careful about the order you place the vertices in the input file to ensure normals face the correct direction.
The normals are calculated by (c-b)x(a-b), if the vertices are specified as a, b, c, in order for both polygons and triangles.

This ray tracer supports anti aliasing with Monte Carlo and jitter sampling, axis aligned(hierarchical) bounding boxes.
To use Monte Carlo, add the flag "-mc x", where x is the number of samples per pixel, i.e. "-mc 10".
To use jitter sampling, add the flag "-j x y", where x is the number of rows and y is the number of columns, i.e. "-j 5 5". 
Default is sampling the center of each pixel.

This ray tracer can also use axis aligned bounding boxes and bounding volume hierarchies for more efficient intersection. 
Add the flag "-bvh x", where x is the maximum recursion depth. You can specify -1 for infinite depth, i.e. "-bvh". Running
on my development computer results in 2 sec execution time vs 30 sec for teapot.inp.

You can accomplish softer shadows by using "ltr" instead of "ltp". "ltr" creates a light with radius.
Format is "ltr x y z r g b radius num_samples falloff", replace with relevant values.

You can also specify a default "sky". In the input file, add the line "sky r g b", replacing r,g, and b with the relevant values.
For a blue sky, you'd put "sky 0 0 1".

To create a refracting material, add "ax ay az n" to the end of the mat line in the input file, replacing with relevant files.
An example would be "1 1 1 10".