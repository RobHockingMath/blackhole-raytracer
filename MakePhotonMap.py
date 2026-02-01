#All this does is compiles the code for a basic one frame thing, not split into parts

import os,sys
import math
from math import *


success = 1
while success!=0:

	command="g++ -fopenmp -Ofast Main_PhotonMap_Geola.cpp Shape.cpp Shape.h "\
	+"Sphere.h Sphere.cpp Texture.h Texture.cpp Light.h "\
	+"Light.cpp Camera.h Camera.cpp World.h World.cpp "\
	+"Chimera.h Chimera.cpp "\
	+"bmp.h bmp.cpp " \
	+"poly_root.h invisibility.h invisibility.cpp Distance_fractal.h Distance_fractal.cpp " \
	+"Menger_Sponge.h Menger_Sponge.cpp " \
	+"Menger_Sponge_Type2.h Menger_Sponge_Type2.cpp " \
	+"Menger_Sponge_4DSlice_Type2.h Menger_Sponge_4DSlice_Type2.cpp " \
	+"Menger_Sponge_4DSlice.h Menger_Sponge_4DSlice.cpp " \
	+"Quaternion_Julia_Set.h Quaternion_Julia_Set.cpp Juliabulb.h Juliabulb.cpp Mandelbulb.h Mandelbulb.cpp " \
	+"Spacetime.h Spacetime.cpp Minkowski.h Minkowski.cpp Schwarzschild.h Schwarzschild.cpp " \
	+"nanoflann.hpp PhotonMapLight.h PhotonMapLight.cpp " \
	+"photon.hpp PhotonMap.h PhotonMap.cpp -o " \
	+"CreatePhotonMap"


	#sucess=0 only if the compilation is successful
	success=os.system(command)

	print(success)