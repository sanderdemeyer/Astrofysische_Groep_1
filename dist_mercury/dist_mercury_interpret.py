import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

def linear(x, a, b):
    return a*x+b

def sinewave(x, a, w, phi, offset):
    return a*np.sin(w*x+phi) + offset

PERIHELION = True

data = np.loadtxt("dist_mercury/Solar-System-New_PEFRL_1000.000000_0.001000.txt",unpack=True)

columns = np.shape(data)[0]

planets = (columns-1)//2
print(planets)

times = data[0,:]
dist = data[1,:]
angle = data[2,:]

body_names = ['Mercury', 'Venus', 'Earth', 'Moon', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto']
planet_names = ['Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune']

avg = 1

perihelion_times = []
perihelion_angles = []
perihelion_angles_averaged = []
moving_perihelion = []
kepler_data = []
major_axes = []
orbital_periods = []

for planet in range(planets):
    if planet != 3 and planet != 9:
        angle_help = [True if data[2*planet+1,i] > data[2*planet+1,i-1] else False for i in range(len(dist))]
        angles_of_max = [data[2*planet+2,i-1] if (angle_help[i-1] and not angle_help[i]) else 'no' for i in range(len(dist))]
        times_of_max = [times[i-1] if (angle_help[i-1] and not angle_help[i]) else 'no' for i in range(len(dist))]
        angles_of_max = [e for e in angles_of_max if e != 'no']
        times_of_max = [e for e in times_of_max if e != 'no']
        running_avg = [np.mean(angles_of_max[i-avg:i]) for i in range(len(angles_of_max))]

        popt, pcov = opt.curve_fit(linear, times_of_max[avg+1:], running_avg[avg+1:])
        moving_perihelion.append(-popt[0]*100*180*3600/np.pi)

        # plt.plot(times_of_max, running_avg)
        # plt.plot(times_of_max, [linear(x, popt[0], popt[1]) for x in times_of_max])
        # plt.show()

        perihelion_times.append(times_of_max)
        perihelion_angles.append(angles_of_max)
        perihelion_angles_averaged.append(running_avg)

        maximum_distance = np.max(data[2*planet+1,:])
        minimum_distance = np.min(data[2*planet+1,:])
        major_axis = (maximum_distance+minimum_distance)/2
        orbital_period = (times[-1]*365.25/len(times_of_max))

        average_dist = np.mean(data[2*planet+1,:])
        switch_help = [True if ((data[2*planet+1,i] > average_dist) and (data[2*planet+1,i-1] < average_dist)) else False for i in range(len(dist))]
        number_of_switches = np.sum(switch_help)
        orbital_period_new = (times[-1]*365.25/number_of_switches)

        print(f'for {body_names[planet]}, the orbital period is {orbital_period}, and the semi-major axis is {major_axis}, and the number of times of max is {len(times_of_max)}')
        kepler_data.append(major_axis**3/orbital_period_new**2)
        major_axes.append(major_axis)
        orbital_periods.append(orbital_period_new)


        print(number_of_switches)
        print(f'number_of_switches for {body_names[planet]} is {number_of_switches}')


print(kepler_data)


plt.scatter(planet_names, kepler_data, label = 'Own data')
plt.scatter(planet_names, np.array([7.496, 7.496, 7.496, 7.495, 7.504, 7.498, 7.506, 7.504])*10**(-6), label = 'Literature')
plt.legend()
plt.show()

orbital_periods_literature = [87.97, 224.70, 365.26, 686.98, 4331.865, 10760.265, 30684.6525, 60189.5475]
major_axes_literature = [0.39, 0.725, 1.0, 1.525, 5.2, 9.54, 19.185, 30.06]

print(orbital_periods)
print(major_axes)
fig = plt.figure()
ax = plt.gca()
ax.scatter(orbital_periods, major_axes, label = 'Own data')
ax.scatter(orbital_periods_literature, major_axes_literature, label = 'Literature')

(popt, pcov) = opt.curve_fit(linear, [np.log10(e) for e in orbital_periods], [np.log10(e) for e in major_axes])

print(popt)
X_data_scaled = np.linspace(1.5, 5, 1000)
Y_data_scaled = [linear(x, popt[0], popt[1]) for x in X_data_scaled]

#ax.plot()

ax.set_xscale("log")
ax.set_yscale("log")
ax.plot([10**e for e in X_data_scaled], [10**e for e in Y_data_scaled], c = 'black', label = 'Linear fit')
plt.xlabel('Orbital period [days]', fontsize = 15)
plt.ylabel('Semi-major axis [AU]', fontsize = 15)
plt.title("Proof of Kepler's 3rd law", fontsize = 15)
plt.legend(fontsize = 15)
plt.show()

plt.scatter(perihelion_times[0], perihelion_angles_averaged[0])
plt.xlabel('Time [years]', fontsize = 15)
plt.ylabel('Angle of perihelion [rad]', fontsize = 15)
plt.title('Precession of perihelion of Mercury', fontsize = 15)
plt.show()
print(len(angles_of_max))
