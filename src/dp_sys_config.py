# --Status ('MakeObs' to invent observations, 'Analyse' to do 4dVar)
mode = "MakeObs"
# Dimensions of pendulum (l1,l2,l3, metres)
l1 = 0.1
l2 = 0.2
l3 = 0.15
# Masses of pendulum nodes (m1,m2,m3, kg)
m1 = 0.1
m2 = 0.1
m3 = 0.1
# Acceleration due to gravity (m/s/s)
g = 10.0
# Time step (seconds)
time_step = 0.001
#Frequency of observation output (used only for MakeObs)
obs_freq = 50
# Length of integration (time steps)
length_of_integration = 5000
# Model state (initial guess(Analyse) or initial conditions(MakeObs)) angle 1,angle 2, roc angle 1, roc angle 2 (deg and deg/s)
angle1 = 160.0
angle2 = 100.0
roc_angle1 = -350.0
roc_angle2 = 100.0

#Noise for angles (deg, used only for MakeObs)
noise_pos = 5.0
#Noise for roc angles (deg/s used only for MakeObs)
noise_vel = 0.0

# Size of each observation window (time steps)
obs_window = 500
