import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sigfig import round

from Solar import Celestial, main

class Simulation:  # class to create simulation object
    def __init__(self, parameters, beeman):
        self.beeman = beeman  # boolean to indicate whether to use Beeman method or Euler method

        # name of simulation
        if beeman:
            self.sim_name = "(Beeman)"  
        else:
            self.sim_name = "(Euler)"
        
        # checking whether to simulate satellite launch or not
        self.launch = True if beeman and input("Launch satellite from Earth to Mars? (y/n): ") == "y" else False
        print()

        self.niter = int(parameters[0])  # number of iterations
        self.time_interval = float(parameters[1])  # time interval between iterations
        self.G = float(parameters[2])  # gravitational constant from input file
        self.celestials = dict()  # dictionary to refer to celestial objects
        self.total_time = 0.0  # initialising total time elapsed in simulation

        # creating celestial objects and adding to dictionary
        for i in range(3, len(parameters)-4, 4):
            self.celestials[parameters[i].strip()] = Celestial(parameters[i].strip(), float(parameters[i+1]), float(parameters[i+2]), "#" + parameters[i+3])

        # creating satellite object if required
        if self.launch:
            self.celestials["Satellite"] = Celestial("Satellite",(590000/5.97219e+24), self.celestials["Earth"].orbital_radius + 4.2635e-5, "white")
            self.reached_mars = False  # boolean to indicate whether satellite has reached Mars
            self.reached_earth = False  # boolean to indicate whether satellite has returned to Earth

        self.list_of_energies = []  # list to store energies of system
        self.list_of_times = []  # list to store times when energies are calculated
 
    # function to be run at start of simulation
    def init(self):
        # initialising velocity and acceleration of each celestial object
        for celestial in self.celestials.values():
            if celestial.name != "Sun"  and celestial.name != "Satellite":
                celestial.initialise(self.celestials["Sun"], self.G)  # setting initial velocity and acceleration of planets
            elif celestial.name == "Satellite":
                # setting an initial velocity for the satellite such that it will do a flyby of mars and return to earth
                celestial.vel[1] = 7.323399766855
                celestial.acc[1] = (-self.G * self.celestials["Sun"].mass * celestial.pos[1]) / (np.linalg.norm(celestial.pos[1]) ** 3)
                celestial.acc[0] = np.copy(celestial.acc[1])
        return self.patches

    # function to update plot
    def animate(self, i):
        # total time elapsed in simulation
        self.total_time += self.time_interval

        patches_index = 0  # index to keep track of which patch to update

        # updating position, velocity and acceleration of each celestial object
        for celestial in self.celestials.values():
            if self.beeman:
                celestial.beeman_update(self.celestials, self.time_interval, self.G, self.planets_gravity)
            else:
                celestial.euler_update(self.celestials, self.time_interval, self.G, self.planets_gravity)

            # updating position of each celestial object in plot
            self.patches[patches_index].center = celestial.pos[1]
            patches_index += 1
        
        # adding energy of system to list to be plotted later
        self.list_of_energies.append(self.energy())

        self.orbital_period()  # checking if orbit is completed for each celestial object

        # checking if all the planets are aligned
        if self.doomsday():
            print("Doomsday! Total time elapsed = " + str(round(self.total_time, sigfigs=4)) + " earth years.")
            print()

        # checking if satellite has reached Mars
        if self.launch:
            self.satellite_journey()

        return self.patches
    
    # Experiment 1
    def orbital_period(self):
        for celestial in self.celestials.values():
            # condition to check if orbit is completed for each celestial object
            if celestial.pos[0, 1] < 0 and celestial.pos[1, 1] >= 0.0:
                celestial.years += 1
                # outputting number of years elapsed for each celestial object
                celestial.orbital_periods.append(self.total_time/celestial.years)
                print(celestial.name + " " + str(celestial.years) + " years = " + str(round(self.total_time, sigfigs=4)) + " earth years.")
                if celestial.name == "Earth":
                    # outputting total energy of system if earth completes orbit
                    print()
                    print("Time = " + str(round(self.total_time, sigfigs=4)) + " earth years. " + "Total energy = " + str(round(self.list_of_energies[-1], sigfigs=4)) + " J.")
                    print()

    # Experiment 2
    def energy(self):
        kinetic_energy = 0.0  # initialising kinetic energy
        potential_energy = 0.0  # initialising potential energy

        for primary_celestial in self.celestials.values():
            # adding kinetic energy of each celestial object
            kinetic_energy += 0.5 * primary_celestial.mass * (np.linalg.norm(primary_celestial.vel) ** 2)

            # adding potential energy of each celestial object
            for secondary_celestial in self.celestials.values():
                if secondary_celestial.name != primary_celestial.name:
                    potential_energy -= self.G * primary_celestial.mass * secondary_celestial.mass / np.linalg.norm(primary_celestial.pos[1] - secondary_celestial.pos[1])

        potential_energy = potential_energy / 2  # dividing potential energy by 2 due to double counting
        
        return (kinetic_energy + potential_energy) * (5.97219e+24 * 1.496e+11 * 1.496e+11) / (3.154e+7 * 3.154e+7)  # returning total energy in Joules

    # Experiment 3
    def satellite_journey(self):
        # checking if satellite has reached mars
        if np.linalg.norm(self.celestials["Satellite"].pos[1] - self.celestials["Mars"].pos[1]) < 0.05 and not self.reached_mars:
            self.reached_mars = True
            print()
            print("Satellite has reached mars. Total time elapsed = " + str(round(self.total_time, sigfigs=4)) + " earth years.")
            print("Distance between satellite and mars = " + str(round(np.linalg.norm(self.celestials["Satellite"].pos[1] - self.celestials["Mars"].pos[1]), sigfigs=4)) + " AU.")
            print()
        # checking if satellite has returned to earth
        if np.linalg.norm(self.celestials["Satellite"].pos[1] - self.celestials["Earth"].pos[1]) < 0.02 and self.reached_mars and not self.reached_earth:
            self.reached_earth = True
            print()
            print("Satellite has returned to Earth. Total time elapsed = " + str(round(self.total_time, sigfigs=4)) + " earth years.")
            print("Distance between satellite and earth = " + str(round(np.linalg.norm(self.celestials["Satellite"].pos[1] - self.celestials["Earth"].pos[1]), sigfigs=4)) + " AU.")
            print()

    # Experiment 5
    def doomsday(self):
        list_of_positions = []  # list to store unit vectors of positions of planets
        for celestial in self.celestials.values():
            if celestial.name != "Sun" and celestial.name != "Satellite":
                list_of_positions.append(celestial.pos[1] / np.linalg.norm(celestial.pos[1]))
        
        mean_vector = np.mean(list_of_positions, axis=0)  # calculating mean vector of unit vectors

        # checking if all the planets are aligned within 5 degrees of the mean vector
        for position in list_of_positions:
            if np.arccos(np.dot(mean_vector, position)) > 5 * np.pi / 180:
                return False
        return True

    # function to set initial plot and run animation
    def run(self):
         # checking whether to include other planets or not (Experiment 6)
        self.planets_gravity = True if self.beeman and input("Include gravitational influence of planets on each other? (y/n): ") == "y" else False
        print()

        print("----------Running simulation " + self.sim_name + "...")
        
        # creating figure and axes for simulation
        fig_1 = plt.figure(figsize=(8.5, 8.5)) 
        fig_1.suptitle("Solar System Simulation " + self.sim_name, fontsize=20)
        ax = plt.axes()

        # finding maximum orbital radius of all celestial objects
        self.max_orbit = 0.0
        for celestial in self.celestials.values():
            if celestial.orbital_radius > self.max_orbit:
                self.max_orbit = celestial.orbital_radius
        
        # adding celestial objects to patches list
        self.patches = []  # creating patches list
        for celestial in self.celestials.values():
            if celestial.name == "Sun":
                size = 0.04 * self.max_orbit
            elif celestial.name == "Satellite":
                size = 0.007 * self.max_orbit
            else:  # planets
                size = 0.015 * self.max_orbit
            self.patches.append(ax.add_patch(plt.Circle(celestial.pos[1], size, color=celestial.colour, animated=True, label=celestial.name)))
        
        # setting plot properties
        ax.axis('scaled')
        ax.set_xlim(-1.2 * self.max_orbit, 1.2 * self.max_orbit)
        ax.set_ylim(-1.2 * self.max_orbit, 1.2 * self.max_orbit)
        ax.set_xlabel('X-Coordinate (AU)')
        ax.set_ylabel('Y-Coordinate (AU)')
        ax.set_facecolor('xkcd:black')  # setting background colour to black
        ax.legend()  # adding legend to plot
        
        # creating animation object
        self.anim = FuncAnimation(fig_1, self.animate, init_func=self.init, frames=self.niter, repeat=False, interval=1, blit=True)
        
        plt.show()  # showing animation

        # checking if satellite has reached mars and returned to earth (Experiment 3)
        if self.launch:
            if not self.reached_earth and self.reached_mars:
                print()
                print("Satellite has not returned to earth. Total time elapsed = " + str(round(self.total_time, sigfigs=4)) + " earth years.")
                print()

        literature_orbit_period = dict()  # dictionary to store literature orbital periods of planets

        # extracting literature orbital periods from Orbital.txt file
        with open("Orbital.txt", "r") as file:
            file_data = file.readlines()
            for line in file_data:
                line = line.split(';')
                literature_orbit_period[line[0]] = float(line[1])
        file.close()

        # checking how mean orbital periods of planets compare to literature values
        for celestial in self.celestials.values():
            if celestial.name != "Sun" and celestial.name != "Satellite":
                if len(celestial.orbital_periods) != 0:
                    mean_orbital_period = np.mean(celestial.orbital_periods)
                    difference = mean_orbital_period - literature_orbit_period[celestial.name]
                    if np.abs(difference) <= 0.005 * literature_orbit_period[celestial.name]:
                        print(celestial.name + " has an orbital period within 0.5%" + " of it's literature value.")
                    else:
                        if difference < 0:
                            print(celestial.name + " has a shorter orbital period than expected.")
                        else:
                            print(celestial.name + " has a longer orbital period than expected.")
                    print(celestial.name + " has an mean orbital period of " + str(round(mean_orbital_period, sigfigs=4)) + " earth years.")
                else:
                    print(celestial.name + " has not completed an orbit.")
                print("Expected orbital period = " + str(literature_orbit_period[celestial.name]) + " earth years.")
                print()

        # printing whether total energy of system is conserved (Experiment 2)
        if len(self.list_of_energies) != 0 and np.max(self.list_of_energies) - np.min(self.list_of_energies) <= 0.001 * np.mean(self.list_of_energies):
            print("Total energy of system is conserved within 0.1%" + " of the the mean energy of the system.")
        else:
            print("Total energy of system is not conserved within 0.1%" + " of the the mean energy of the system.")
        
        # creating figure for energy plot (Experiment 2)
        plt.figure(figsize=(8.5, 8.5))
        plt.plot(np.linspace(0,self.total_time,len(self.list_of_energies)), self.list_of_energies)
        plt.title('Total Energy of System' + self.sim_name, fontsize=20)
        plt.xlabel('Time (Earth Years)')
        plt.ylabel('Total Energy (J)')
        plt.show()