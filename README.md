## Epipolar-Geometry
Relative Orientation of an Image Pair using Fundamental Matrix
The relative orientation of an image pair is defined by the epipolar geometry. Using algebraic projective geometry, the epipolar geometry can be represented by the fundamental matrix.

## Image acquisition:
Take pictures of a spatially structured object from two different views. Use thereby a general convergent image arrangement.
#Image pair orientation:
a.) Manually select at least 8 homologous points x ↔ x' in the image pair and implement an Octave / MATLAB function for the linear computation of the fundamental matrix F (Hint: use the normalized 8-point-algorithm).
b.) Mark the used points in both images (see get_points from Assignment 2) and draw the associated epipolar lines in the corresponding image. For drawing lines in homogeneous coordinates l = (a, b, c)T use the auxiliary function hline.m
Evaluation:
a.) Show the image pair and comment the line characteristics in brief.
b.) Calculate the geometric image error (symmetric epipolar distance) of F for all points.
