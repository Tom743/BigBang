import vpython as vpy
from vpython import vector
import numpy as np
import math
from random import random, uniform
from parameters import *


def ranvec():
      # Generates a randomly orientated unit vector as numpy array
      theta = 2.0*math.pi*random()
      phi = math.acos(2.0*random()-1.0)-0.5*math.pi
      z = math.sin(phi)
      x = math.cos(phi)*math.sin(theta)
      y = math.cos(phi)*math.cos(theta)
      vec = np.array([x,y,z])
      return vec


scene = vpy.canvas(
	title="Stars",
  	range=L,
  	width=1280,
  	height=800,
  	forward=vector(-1,-1,-1))


# Create the stars. Uses numpy arrays for faster operations
stars = np.array([])
positions = np.array([])
masses = np.array([])
velocities = np.array([])
for star in range(Nstars):

	# Position
	posVector = L*ranvec()*(random()**(dist0))
	x = posVector[0]
	y = posVector[1]
	z = posVector[2]
	pos = vector(x,y,z)

	# Radius
	r = Rsun

	# Mass
	m = Msun

	# Color
	col = vector(uniform(0.7,1.0),uniform(0.7,1.0), uniform(0.7,1.0))

	# Create the star
	stars = np.append(stars, vpy.sphere(pos=pos, radius=r, color=col))
	positions = np.append(positions, pos)
	masses = np.append(masses, m)
	velocities = np.append(velocities, vector(0,0,0))


# Main loop
while True:
	vpy.rate(10)

	# Star to star scalar distances
	rscalar = positions-positions[:,np.newaxis]
	for row in range(np.shape(rscalar)[0]):
		for col in range(np.shape(rscalar)[1]):
			rscalar[row, col] = math.sqrt(rscalar[row, col].dot(rscalar[row, col]))
			if (rscalar[row, col] == 0):
				rscalar[row, col] = 1

	# Star to star vectorial distances
	r = positions-positions[:,np.newaxis]

	# Star to star forces
	F = G*masses*masses[:,np.newaxis]/(rscalar**2) # Module
	for row in range(np.shape(F)[0]):
		for col in range(np.shape(F)[1]):
			F[row][col] = r[row][col].norm()*F[row][col] # Complete vector

	# No self-forces
	for n in range(np.shape(F)[0]):
		F[n][n] = vector(0,0,0)

	# Update star positions
	velocities += (dt*np.sum(F, 1))/masses
	positions += velocities*dt

	# Merge if they are close enough
	br = False
	ran = Nstars
	for r in range(Nstars):
		for c in range(Nstars):

			# rscalar stores the distance between each star, so rscalar[1][0] == rscalar[0][1]. This conditional avoids looking in every value, as half is the same as the other half
			if c+r<Nstars:
				c+=r
			else:
				break 

			# Force the loop to the end if it reached the part of stars merged in this iteration of the main loop
			if (c >= ran or r >= ran):
				r , c = Nstars, Nstars # Force the loop to end
				br = True # Go straight to the end

			# Merger algorithm core
			if (r!=c and not br):
				if (rscalar[r, c] <= stars[r].radius*4/3 or rscalar[r, c] <= stars[c].radius*4/3):

					# Make one bigger and more massive
					new_radius = (stars[c].radius**3+stars[r].radius**3)**(1/3) # New radius maintaining the density
					new_mass = masses[c] + masses[r]
					new_velocity = (masses[c]*velocities[c]+masses[r]*velocities[r])/new_mass # The new velocity comes from the momenta of the other two
					new_position = (positions[c]+positions[r])/2 # New position is the average of the two merging stars
					col = vector(uniform(0.7,1.0), uniform(0.7,1.0), uniform(0.7,1.0))

					# Add it to the control arrays
					stars = np.append(stars, vpy.sphere(pos=new_position, radius=new_radius, color=col))	
					positions = np.append(positions, new_position)
					masses = np.append(masses, new_mass)
					velocities = np.append(velocities, new_velocity)

					# Make the originals disappear
					stars[r].visible = False
					stars[c].visible = False

					# Delete those stars from the control arrays
					if(c>r):
						stars = np.delete(stars, c)
						stars = np.delete(stars, r)
						positions = np.delete(positions, c)
						positions = np.delete(positions, r)
						velocities = np.delete(velocities, c)
						velocities = np.delete(velocities, r)
						masses = np.delete(masses, c)
						masses = np.delete(masses, r)
					if(c<r):
						stars = np.delete(stars, r)
						stars = np.delete(stars, c)
						positions = np.delete(positions, r)
						positions = np.delete(positions, c)
						velocities = np.delete(velocities, r)
						velocities = np.delete(velocities, c)
						masses = np.delete(masses, r)
						masses = np.delete(masses, c)
					rscalar = np.delete(np.delete(rscalar, r, 0), c, 1)
					Nstars -= 1

					# Fix the loop
					r-=1
					c-=1
					ran-=2

	# Expand the universe
	positions *= h0*dt+1.0

	# Update the spheres positions
	for i in range(stars.size):
		stars[i].pos = positions[i]
