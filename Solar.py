import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from sigfig import round  # module to round to significant figures or decimal places
from numpy.linalg import norm


class Celestial:  # class to define celestial objects
    def __init__(self, name, mass, orbital_radius, colour):
        self.name = name
        self.mass = mass 
        self.colour = colour
        self.orbital_radius = orbital_radius 
        self.pos = np.array([[0.0, 0.0], [orbital_radius, 0.0]])  # initialising position array containing old and current position vectors respectively
        self.vel = np.array([0.0, 0.0]) # initialising velocity array containing current velocity vector
        self.acc = np.array([[0.0, 0.0], [0.0, 0.0]])  # initialising acceleration array containing old and current acceleration vectors respectively
        self.years = 0  # initialising number of years elapsed for given celestial object
        self.orbital_period = 0.0  # property to store simulated orbital period of celestial for comparison at the end of simulation

    # function to initialise velocity and acceleration if not sun
    def initialise(self, sun, G):
        self.vel[1] = 7.327 if self.name == "Satellite" else np.sqrt(G * sun.mass / self.orbital_radius)  # setting inital y component of velocity (number is appropriate launch velocity for satellite to get to mars)
        self.acc[1] = (-G * sun.mass * self.pos[1]) / (norm(self.pos[1]) ** 3) # setting initial acceleration
        self.acc[0] = np.copy(self.acc[1]) # setting old acceleration to current acceleration for Beeman method

    # fuction to update position of celestial object
    def update_pos(self, interval, beeman):
        self.pos[0] = np.copy(self.pos[1])  # updating old position to check if orbit is completed later on in program
        self.pos[1] += self.vel * interval if not beeman else self.vel * interval + (4 * self.acc[1] - self.acc[0]) * (interval ** 2)/6  # updating position using Beeman method or Euler method

    # function to update velocity of celestial object
    def update_vel(self, celestial, interval, G, beeman, planets_gravity):
        if celestial.name != self.name and (planets_gravity or celestial.name == "Sun"):  # condition to prevent calculation with self
            direction = celestial.pos[1] - self.pos[1]  # direction vector for gravitational force
            new_acc = G * celestial.mass * direction / (norm(direction) ** 3) # calculating new acceleration
            self.vel += new_acc * interval if not beeman else (2 * new_acc + 5 * self.acc[1] - self.acc[0]) * interval/6  # updating velocity using Beeman method or Euler method
            if beeman:  # updating acceleration if Beeman method is used
                self.acc[0] = np.copy(self.acc[1])
                self.acc[1] = np.copy(new_acc)


class Simulation:  # class to create simulation object
    def __init__(self, parameters, beeman, planets_gravity, launch):
        self.beeman = beeman  # boolean to indicate whether to use Beeman method or Euler method
        self.planets_gravity = planets_gravity  # boolean to indicate whether to include planets in gravitational calculations
        self.launch = launch  # boolean to indicate whether to include satellite in simulation
        self.sim_name = "(Beeman)" if beeman else "(Euler)"  # simulation name to be displayed in plot title
        self.niter = int(parameters[0])  # number of iterations
        self.time_interval = float(parameters[1])  # time interval between iterations
        self.G = float(parameters[2])  # gravitational constant from input file
        self.celestials = dict()  # dictionary to refer to celestial objects
        self.total_time = 0.0  # initialising total time elapsed in simulation
        self.list_of_energies = []  # list to store energy of system at each iteration

        # creating celestial objects and adding to dictionary
        for i in range(3, len(parameters)-4, 4):
            self.celestials[parameters[i].strip()] = Celestial(parameters[i].strip(), float(parameters[i+1]), float(parameters[i+2]), "#" + parameters[i+3].strip())
            if parameters[i].strip() != "Sun":
                self.celestials[parameters[i].strip()].initialise(self.celestials["Sun"], self.G)  # setting initial velocity and acceleration of planets

        # creating satellite object if required
        if self.launch:
            self.celestials["Satellite"] = Celestial("Satellite",1e-21, self.celestials["Earth"].orbital_radius + 4.2635e-5, "white")
            self.celestials["Satellite"].initialise(self.celestials["Sun"], self.G)  # setting initial velocity and acceleration of satellite
            self.reached_mars = False
            self.reached_earth = False  # boolean values to indicate whether satellite has reached Mars and returned to Earth
            
    # function to update plot
    def animate(self, i):
        self.total_time += self.time_interval  # updating total time elapsed in simulation
        patches_index = 0  # index to keep track of which patch to update

        # updating position of each celestial object in plot
        for celestial in self.celestials.values():
            celestial.update_pos(self.time_interval, self.beeman)
            self.patches[patches_index].center = celestial.pos[1]
            patches_index += 1
        
        # updating velocity and acceleration of each celestial object
        for celestial_1 in self.celestials.values():
            for celestial_2 in self.celestials.values():
                celestial_1.update_vel(celestial_2, self.time_interval, self.G, self.beeman, self.planets_gravity)

        self.list_of_energies.append(self.energy())  # calculating and storing energy of system
        self.orbital_period()  # checking if orbit is completed for each celestial object
        self.doomsday()  # checking if all the planets are aligned
        self.satellite_journey() if self.launch else None  # checking if satellite has reached Mars and returned to Earth
        return self.patches
    
    # Experiment 1
    def orbital_period(self):
        for celestial in self.celestials.values():
            if celestial.name != "Sun" and celestial.name != "Satellite" and celestial.pos[0, 1] < 0 and celestial.pos[1, 1] >= 0.0:  # condition to check if orbit is completed for each celestial object
                celestial.years += 1
                if self.planets_gravity:
                    celestial.orbital_period = self.total_time/celestial.years if celestial.years == 1 else celestial.orbital_period
                else:
                    celestial.orbital_period = self.total_time/celestial.years if self.total_time <= 50 else celestial.orbital_period
                print(celestial.name + " " + str(celestial.years) + " years = " + str(round(self.total_time, decimals=3)) + " earth years.")  # outputting time elapsed for each celestial object
                # outputting total energy of system if earth completes orbit
                print("\nTime = " + str(round(self.total_time, decimals=3)) + " earth years. " + "Total energy = " + str(round(self.list_of_energies[-1], sigfigs=4)) + " J.\n") if celestial.name == "Earth" else None

    # Experiment 2
    def energy(self):
        kinetic_energy = 0.0  # initialising kinetic energy
        potential_energy = 0.0  # initialising potential energy
        for celestial_1 in self.celestials.values():  # adding kinetic energy of each celestial object
            kinetic_energy += celestial_1.mass * np.dot(celestial_1.vel, celestial_1.vel)/2
            for celestial_2 in self.celestials.values():  # adding potential energy of each celestial object
                # calculating potential energy for each pair of celestial objects and dividing by 2 to avoid double counting
                potential_energy -= self.G * celestial_1.mass * celestial_2.mass / (2 * norm(celestial_1.pos[1] - celestial_2.pos[1])) if celestial_2.name != celestial_1.name else 0.0
        return (kinetic_energy + potential_energy) * (5.97219e+24 * (1.496e+11 ** 2)) / (3.154e+7 ** 2)  # returning total energy in Joules (SI units)

    # Experiment 3
    def satellite_journey(self):
        if norm(self.celestials["Satellite"].pos[1] - self.celestials["Mars"].pos[1]) < 0.01 and not self.reached_mars:  # checking if satellite has reached mars
            self.reached_mars = True
            print("\nSatellite has reached Mars. Total time elapsed = " + str(round(self.total_time, decimals=3)) + " earth years.")
            print("Distance between Satellite and Mars = " + str(round(norm(self.celestials["Satellite"].pos[1] - self.celestials["Mars"].pos[1]), sigfigs=4)) + " AU.\n")
        if norm(self.celestials["Satellite"].pos[1] - self.celestials["Earth"].pos[1]) < 0.01 and self.reached_mars and not self.reached_earth: # checking if satellite has returned to earth
            self.reached_earth = True
            print("\nSatellite has returned to Earth. Total time elapsed = " + str(round(self.total_time, decimals=3)) + " earth years.")
            print("Distance between Satellite and Earth = " + str(round(norm(self.celestials["Satellite"].pos[1] - self.celestials["Earth"].pos[1]), sigfigs=4)) + " AU.\n")

    # Experiment 5
    def doomsday(self):
        list_of_positions = []  # list to store unit vectors of positions of planets
        # appending unit vectors of positions of planets to list
        for celestial in self.celestials.values():
            if celestial.name != "Sun" and celestial.name != "Satellite":
                list_of_positions.append(celestial.pos[1] / norm(celestial.pos[1]))
        mean_vector = np.mean(list_of_positions, axis=0)  # calculating mean vector of unit vectors
        aligned = True
        # checking if all the planets are aligned within 5 degrees of the mean vector
        for position in list_of_positions:
            if np.arccos(np.dot(mean_vector, position)) > 5 * np.pi / 180:
                aligned = False
                break
        print("Doomsday! Total time elapsed = " + str(round(self.total_time, decimals=3)) + " earth years.", end="\n\n") if aligned else None

    # function to set initial plot and run animation
    def run(self):
        print("----------Running simulation " + self.sim_name + "...\n")
        # creating figure and axes for simulation
        fig_1 = plt.figure(figsize=(8.5, 8.5)) 
        fig_1.suptitle("Solar System Simulation " + self.sim_name, fontsize=20)
        ax = plt.axes()

        # finding maximum orbital radius of all celestial objects
        self.max_orbit = 0.0
        for celestial in self.celestials.values():
            self.max_orbit = celestial.orbital_radius if celestial.orbital_radius > self.max_orbit else self.max_orbit
        
        # adding celestial objects to patches list
        self.patches = []  # creating patches list
        for celestial in self.celestials.values():
            if celestial.name == "Sun":
                size = 0.04 * self.max_orbit
            elif celestial.name == "Satellite":
                size = 0.007 * self.max_orbit
            else:  # planets
                size = 0.015 * self.max_orbit
            self.patches.append(ax.add_patch(plt.Circle(celestial.pos[1], size, color=celestial.colour, animated=True, label=celestial.name))) # type: ignore
        
        # setting plot properties
        ax.axis('scaled')
        ax.set_xlim(-1.2 * self.max_orbit, 1.2 * self.max_orbit)
        ax.set_ylim(-1.2 * self.max_orbit, 1.2 * self.max_orbit)
        ax.set_xlabel('X-Coordinate (AU)')
        ax.set_ylabel('Y-Coordinate (AU)')
        ax.set_facecolor('xkcd:black')  # setting background colour to black
        ax.legend()  # adding legend to plot
        
        # creating animation object
        self.anim = FuncAnimation(fig_1, self.animate, frames=self.niter, repeat=False, interval=1, blit=True)
        plt.show()  # showing animation

        # checking if satellite has reached mars and returned to earth (Experiment 3)
        if self.launch and self.reached_mars and not self.reached_earth:
            print("\nSatellite has not returned to earth. Total time elapsed = " + str(round(self.total_time, decimals=3)) + " earth years.\n")

        comparison_orbit_period = dict()  # creating dictionary to store literature orbital periods of planets

        # choosing which orbital periods to extract from file for comparison depending on which experiment is being run (1 or 6)
        file_name = "Orbital Periods (Literature).txt" if self.planets_gravity else "Orbital Periods (Planet Gravity).txt"
        with open(file_name, "r") as file:
            file_data = file.readlines()
            for line in file_data:
                line = line.split(';')
                comparison_orbit_period[line[0]] = float(line[1])
        file.close()

        # checking how simulated orbital period  of planets compares to literature values
        for celestial in self.celestials.values():
            if celestial.name != "Sun" and celestial.name != "Satellite":
                if not self.planets_gravity and celestial.name == "Jupiter":
                    continue
                if celestial.orbital_period != 0.0:
                    if np.abs(comparison_orbit_period[celestial.name] - celestial.orbital_period) <= 0.001 * comparison_orbit_period[celestial.name]:
                        print(celestial.name + " has a simulated orbital period within 0.1% " + ("of it's literature value."  if self.planets_gravity else "of it's period with Jupiter's influence (after 50 earth years)."))
                    elif celestial.orbital_period < comparison_orbit_period[celestial.name]:
                        print(celestial.name + " has a shorter orbital period than " + ("expected." if self.planets_gravity else "it's period with Jupiter's influence (after 50 earth years)."))
                    else:
                        print(celestial.name + " has a longer orbital period than " + ("expected." if self.planets_gravity else "it's period with Jupiter's influence (after 50 earth years)."))
                    print(celestial.name + " has a simulated orbital period of " + str(celestial.orbital_period) + (" earth years." if self.planets_gravity else " earth years (without Jupiter's influence)."))
                else:
                    print(celestial.name + " has not completed an orbit.")
                print(("Expected orbital period = " if self.planets_gravity else "Period with Jupiter's influence for 50 years = ") + str(comparison_orbit_period[celestial.name]) + " earth years.\n")

        # printing whether total energy of system is conserved (Experiment 2)
        if len(self.list_of_energies) != 0 and np.max(self.list_of_energies) - np.min(self.list_of_energies) <= 0.01 * np.mean(self.list_of_energies):
            print("Total energy of system is conserved within 1%" + " of the the mean energy of the system.")
        else:
            print("Total energy of system is not conserved within 1%" + " of the the mean energy of the system.")
        
        # creating figure for energy plot (Experiment 2)
        plt.figure(figsize=(8.5, 8.5))
        plt.plot(np.linspace(0,self.total_time,len(self.list_of_energies)), self.list_of_energies)
        plt.title('Total Energy of System' + self.sim_name, fontsize=20)
        plt.xlabel('Time (Earth Years)')
        plt.ylabel('Total Energy (J)')
        plt.show()


# main function to run program
def main():
    with open(input("Enter Filename: ") + ".txt") as file:  # opening input file to get parameters
        input_parameters = []
        for line in file:
            if not line.startswith("#"):  # condition to ignore comments
                input_parameters.append(line)
    file.close()

    planets_gravity = True if input("Include gravitational influence of planets on each other? (y/n): ") == "y" else False  # conditon to conduct Experiment 6
    launch = True if input("Launch satellite from Earth to Mars? (y/n): ") == "y" else False  # condition to conduct Experiment 3
    print()

    # creating simulation objects and running simulations
    beeman = Simulation(input_parameters, True, planets_gravity, launch)
    euler = Simulation(input_parameters, False, planets_gravity, launch)
    beeman.run()
    euler.run()

    # trimming list of energies so that they match in length and creating list of times for energy comparison plot (Experiment 4)
    if len(beeman.list_of_energies) != len(euler.list_of_energies):
        if len(beeman.list_of_energies) > len(euler.list_of_energies):
            beeman.list_of_energies = beeman.list_of_energies[:len(euler.list_of_energies)] 
            list_of_times = np.linspace(0, euler.total_time, len(euler.list_of_energies))  
        else:
            euler.list_of_energies = euler.list_of_energies[:len(beeman.list_of_energies)]
            list_of_times = np.linspace(0, beeman.total_time, len(beeman.list_of_energies))
    else:
        list_of_times = np.linspace(0, beeman.total_time, len(beeman.list_of_energies))
    
    plt.figure(figsize=(8.5, 8.5))  # creating figure and axes for energy comparison plot
    plt.plot(list_of_times, beeman.list_of_energies, label="Beeman", linestyle='-', color='red')
    plt.plot(list_of_times, euler.list_of_energies, label="Euler", linestyle='-', color='blue')
    plt.title('Evolution of Energy (Comparison Plot)', fontsize=20)
    plt.xlabel('Time (Earth Years)')
    plt.ylabel('Total Energy (J)')
    plt.legend()
    plt.show()


main()  # calling main function
