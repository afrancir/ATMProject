# Written by Alessandro Bombelli, 10 February 2017
#______________________________________________________________________________#
#
# Code that loads a Departure/Arrival matrix file for a specific airport, 
# with hourly Departure/Arrival pairs. The convex hull of the distribution is 
# computed, and at most three linear inequality constraints are determined:
# 1) horizontal constraint that limits the max number of hourly arrivals
# 2) vertical constraint that limits the max number of hourly departures
# 3) sloped constraint that limits the hourly linear combination of departures
# and arrivals
# While constraints 1) and 2) are always defined, constraint 3) will or will
# not be present according to the specific properties of the convex hull
# computed
# The three constraints define an approximation of the Gilbo envelope
# characterizing the hourly departures/arrivals of the airport considered.
# References:
# 1 - Gilbo, E., "Airport Capacity: Representation, Estimation, Optimization",
# IEEE Transactions of Control Systems Technology, Vol. 1., No. 3, Sept. 1993
# 2 - Ramanujam, V., Balakrishnan, H., "Estimation of Arrival-Departure Capacity
# Tradeoffs in Multi-Airport Systems
#
#______________________________________________________________________________#
#
# Reference for DBSCAN in Python:
# http://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html ...
# sphx-glr-auto-examples-cluster-plot-dbscan-py
#______________________________________________________________________________#

##########################################################
### Loading all the packages necessary to run the code ###
##########################################################
from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import rand
from numpy import loadtxt
import os
from sklearn.cluster import DBSCAN
from sklearn import metrics

# Choosing the airport that we want to characterize in terms of
# departure/arrival envelope
airport_ID = 55

# String defining the file we want to load
GILBO_this_airport_string = os.getcwd() + "/GILBO_envelopes/" + str(airport_ID) + "_Gilbo.txt"
# Loading .txt file
GILBO_this_airport = loadtxt(GILBO_this_airport_string)

# Loading all departures
Departures = GILBO_this_airport[:,0];
# Loading all arrivals 
Arrivals   = GILBO_this_airport[:,1];
# Departures and arrivals matrix
Dep_Arr = GILBO_this_airport[:,0:2]

# Use DBSCAN to detect Departure/Arrival data points that are labeled as 
# outliers. Parameters used: epsilon=4, min_samples=2, metric=manhattan
db = DBSCAN(eps=4, min_samples=2,metric='manhattan').fit(Dep_Arr)
# Define a False vector with as many elements as the size of the Dep_Arr matrix
core_samples_mask = np.zeros_like(db.labels_,dtype=bool)
# Set data points that were not labeled as outliers as True
core_samples_mask[db.core_sample_indices_] = True

# Augment the data points that passed the outlier detection test with three
# additional data points that are necessary to compute the conver hull as
# desired
Aug_data = np.vstack((Dep_Arr[core_samples_mask,:],[0,0],[0,
max(Dep_Arr[core_samples_mask,1])],[max(Dep_Arr[core_samples_mask,0]),0]))

# Compute the convex hull
hull = ConvexHull(Aug_data)
# Extract the indices of the vertices of the convex hull from the hull 
# structure, and repeat the first  index to that we have a closed shape for
# plotting purposes. Note that the list of indices "Vertices" is only used for
# plotting purposes
Vertices = np.hstack((hull.vertices,hull.vertices[0]))
# Get the convex hull selecting only the vertices that belong to the
# hull.vertices list
CVX_hull = Aug_data[hull.vertices,:]

# Find rows of the CVX_hull characterized by no departures
zero_dep  = np.where(CVX_hull[:,0] == 0)[0]
# Find rows of the CVX_hull characterized by no arrivals
zero_arr  = np.where(CVX_hull[:,1] == 0)[0]
# Find rows of the CVX_hull characterized by nonzero departures and arrivals
slope_idx = np.setdiff1d(np.arange(len(CVX_hull)-1),np.hstack((zero_dep,zero_arr)))

# Determine if we have at least 2 data points of the convex hull that do not lie 
# on the axes. If yes, compute the least square line that best approximates the
# trend of those points
if len(slope_idx)>=2:
	# Label that will be used later in the code for plotting purposes
	sloped_constraint = 1
	# Points that we want to approximate with the least square method
	slope_values      = CVX_hull[slope_idx,:]
	# Matrix with departures (which would be our x) as first column, 
	# and 1s as second column (rememeber that we need to compute a and b
	# such that y = ax+b minimizes the error)
	A                 = np.vstack((slope_values[:,0],np.ones(len(slope_values[:,0])))).T
	# Arrivals
	y                 = slope_values[:,1]
	# Determining (i) a = slope and (ii) b = y-intercept
	a, b              = np.linalg.lstsq(A, y)[0]

###################################
### Plotting the Gilbo Envelope ###
###################################

# Defining fonts for axes and text in the figure
color_text        = (5./255,5./255,5./255)
color_data        = (40./255,40./255,40./255)
color_outliers    = (70./255,70./255,70./255)
color_contour     = (20./255,20./255,20./255)
color_constraints = (110./255,110./255,110./255)

axis_font  = {'fontname':'Arial', 'size':20}
const_font = {'family': 'Arial','color':color_text,'weight': 'normal','size': 16}

fig,ax = plt.subplots(1)
#######################################
### Scatter plot of all data points ###
#######################################
h1  =plt.scatter(Departures[core_samples_mask],Arrivals[core_samples_mask],color=color_data,alpha=0.7, edgecolors='none',marker='o',s=80)
###############################################
### Plotting the contour of the convex hull ###
###############################################
h2, =plt.plot(Aug_data[Vertices,0], Aug_data[Vertices,1],color=color_contour,ls='--',lw=2)
################################################
### Plotting the vertices of the convex hull ###
################################################
h3, =plt.plot(Aug_data[hull.vertices,0], Aug_data[hull.vertices,1],color=color_contour,marker='o',ls='none')
#############################################################################
### If the sloped constraint has been computed, plot the associated line, ###
### plus the horizontal/vertical constraints                              ###
#############################################################################
if sloped_constraint:
	h4 ,=plt.plot([max(Departures[core_samples_mask]),max(Departures[core_samples_mask]),(max(Arrivals[core_samples_mask])-b)/a,0],[0,a*max(Departures[core_samples_mask])+b,max(Arrivals[core_samples_mask]),max(Arrivals[core_samples_mask])],color=color_constraints,ls='-',lw=5)
###############################################################
### Otherwise, plot the the horizontal/vertical constraints ###
###############################################################
else:
	h4 ,=plt.plot([max(Departures[core_samples_mask]),max(Departures[core_samples_mask]),0],[0,max(Arrivals[core_samples_mask]),max(Arrivals[core_samples_mask])],color=color_constraints,ls='-',lw=5)
########################################
### If there are outliers, plot them ###
########################################
if len(Departures[~core_samples_mask]):
	h5  =plt.scatter(Departures[~core_samples_mask],Arrivals[~core_samples_mask], c=color_outliers,alpha=0.7,edgecolors='none',marker='v',s=80)

################################################################################
### Display (as inequalities) the different constraints that charaterize the ###
### current departure/arrival envelope                                       ###
################################################################################

# Definition of two constants that locate where the constraint inequalities
# will be plotted in the figure
horiz_alpha = 0.05
vert_alpha  = 0.3

##########################################################################
### If the sloped constraint has been computed, we highlight all three ###
### inequality constraints                                             ###
##########################################################################
if sloped_constraint:
	if (a>0 and b>0):
		plt.text(horiz_alpha*max(Departures),vert_alpha*max(Arrivals),r"$x \leq$ " + str("{:3.0f}".format(max(Departures[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.0f}".format(max(Arrivals[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.2f}".format(np.abs(a)))+ "$x$ +" + str("{:3.2f}".format(np.abs(b))),fontdict=const_font)
	elif (a>0 and b<0):
		plt.text(horiz_alpha*max(Departures),vert_alpha*max(Arrivals),r"$x \leq$ " + str("{:3.0f}".format(max(Departures[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.0f}".format(max(Arrivals[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.2f}".format(np.abs(a)))+ "$x$ - " + str("{:3.2f}".format(np.abs(b))),fontdict=const_font)
	elif (a<0 and b>0):
		plt.text(horiz_alpha*max(Departures),vert_alpha*max(Arrivals),r"$x \leq$ " + str("{:3.0f}".format(max(Departures[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.0f}".format(max(Arrivals[core_samples_mask]))) + "\n" "$y \leq$ - " + str("{:3.2f}".format(np.abs(a)))+ "$x$ + " + str("{:3.2f}".format(np.abs(b))),fontdict=const_font)
	else:
		plt.text(horiz_alpha*max(Departures),vert_alpha*max(Arrivals),r"$x \leq$ " + str("{:3.0f}".format(max(Departures[core_samples_mask]))) + "\n" "$y \leq$ " + str("{:3.0f}".format(max(Arrivals[core_samples_mask]))) + "\n" "$y \leq$ - " + str("{:3.2f}".format(np.abs(a)))+ "$x$ - " + str("{:3.2f}".format(np.abs(b))),fontdict=const_font)
#############################################################################
### If the sloped constraint has not been computed, we only highlight the ###
### vertical and horizontal constraints                                   ###
#############################################################################
else:
	plt.text(horiz_alpha*max(Departures),vert_alpha*max(Arrivals),r"$x \leq$ " + str("{:3.0f}".format(max(Departures))) + "\n" "$y \leq$ " + str("{:3.0f}".format(max(Arrivals))),fontdict=const_font)

# Plot grid
ax.grid(True)
# Define labels
ax.set_xlabel('Departures [1/h]',**axis_font)
ax.set_ylabel('Arrivals [1/h]',**axis_font)
# Define legend
if len(Departures[~core_samples_mask]):
	ax.legend([h1,h5,h2,h4], ["Processed Data Points","Outliers","Convex Hull","Gilbo Envelope"],scatterpoints=3,loc='upper left',fancybox=True,framealpha=0.5)
else:
	ax.legend([h1,h2,h4], ["Processed Data Points","Convex Hull","Gilbo Envelope"],scatterpoints=3,loc='upper left',fancybox=True,framealpha=0.5)
# Display figure
plt.show()
# Save figure
fig.savefig("Gilbo_"+str(airport_ID)+".pdf", dpi=fig.dpi)


