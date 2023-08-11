# Solar

## Introduction

In 1543 Copernicus proposed that the known planets (at that time only the planets as far out as Saturn) were in circular orbit around the sun. This theory was revised by Kepler and Brahe who showed that the orbits were in fact elliptical. Since then the orbits of the planets have been very accurately calculated and can be precisely predicted far into the future. NASA uses powerful computers to determine both the position of the planets and the path of its space missions, including the Mars landing probes.

This project aims to simulate the motion of the planets in two-dimensions and to predict launch conditions for satellites to successfully reach planets.

## Numerical integration

The solar system represents a gravitational “many body problem”, i.e. one that cannot be solved analytically, which means that a numerical integration scheme must be used. In this project we will use the three-step Beeman scheme. The Beeman algorithm is a stable method which predicts the position at the next time step by combining the current acceleration with the acceleration from the previous time step. This new position can then be used to calculate the new acceleration which, in turn, predicts the new velocity. The algorithm is given by:
###
(1.1)
$\vec{r}(t+Δt)=\vec{r}(t)+\vec{v}(t)Δt+\dfrac{1}{6}[4\vec{a}(t)−\vec{a}(t−Δt)]Δt^2$
###
and
###
(1.2)
$\vec{v}(t+Δt)=\vec{v}(t)+\dfrac{1}{6}[2\vec{a}(t+Δt)+5\vec{a}(t)−\vec{a}(t−Δt)]Δt$.
Where $Δt$ is the small time-step, $t$ is the current time and $\vec{r}$, $\vec{v}$ and $\vec{a}$ represent the position, velocity and acceleration vectors for the motion.

## Project task: Simulation

Write an object-oriented Python program to simulate the orbit of the inner planets (Mercury, Venus, Earth, Mars) around the sun. You should treat the simulation as a “many body problem”, so that adding new planets would not involve changing the design of the code. Your code should:

Read the planet details and simulation parameters from file.
Implement the Beeman integration scheme to update the position and velocity of the planets and the sun at each time step.
Show the orbit of the planets as they move around the sun in a graphical display.
Calculate and print the orbital periods of the planets in Earth years (you may also write them to file if you wish, but this is not required).
Regularly write out to a file the total energy of the system, i.e. the sum of kinetic and gravitational potential energy.
Hide  
### Total energy

- The total kinetic energy $K(t)$ of the system is the sum of the kinetic energies of all of the bodies in the system:

    $K(t)=∑_i\dfrac{1}{2}m_iv_i(t)^2$

    where $m_i$ and $v_i(t)$ are the mass and speed of the $i^{th}$ body. Note that the kinetic energy is positive (or zero for a stationary body).
- The (gravitational) potential energy $U(t)$ of the system is the sum of the potential energies between each pair of bodies:

    $U(t)=−\dfrac{1}{2}∑_{i≠j}\dfrac{Gm_im_j}{r_{ij}(t)}$

    where $mi$ and $m_j$ are the masses of bodies $i$ and $j$, $r_{ij}(t)$ is the distance between them and $G$ is the gravitational constant. Note that potential energy is defined between pairs of bodies - we can’t define the potential energy of a single body (unlike kinetic energy). Note also that we have to divide the total by two to avoid “double counting” pairs of bodies. Finally, note that potential energy is negative.
- The total energy of the system is simply the sum of the kinetic and potential energies:

    $E(t)=K(t)+U(t)$
###
For information on how to set the initial positions and velocities and how to determine the acceleration due to gravitational attraction, see CP5 Orbital Motion.

## Project task: Experiments

Once you have a working code, you should work through the experiments below.

### 1. Orbital periods

	Experiment 1 : Check how closely the orbital periods of the planets in your simulation match their actual orbital periods.
### 2. Energy conservation

The Beeman algorithm is a symplectic integrator, which means that the total energy of the system should be conserved over time.

	Experiment 2 : Check whether (or not) energy is conserved during your simulation. Illustrate your results graphically.
### 3. Satellite to Mars

Suppose you wish to launch a satellite from Earth to perform a fly-past of Mars.

	Experiment 3 : Search for an initial velocity (or range of velocities) that enables a satellite to get close to Mars. You should check:
How close to Mars your satellite gets.
What your satellite’s journey time is and how it compares to that of NASA’s Perseverance mission.
Whether your satellite ever gets back to Earth.
To perform this experiment you will need to extend your code.

You should launch the satellite at the start of your simulation. Take care to start the satellite from just off Earth (so that there is no division by zero) and also to ensure that its mass and starting velocity are realistic - a probe could get to Mars in a few days if it were fired fast enough, but carrying that much fuel would be prohibitive!

### Additional Experiments
In addition to the original three experiments, you must choose two of the three experiments below and work through them as well.

#### 4. Integrating with Direct Euler
Your simulation code implements the Beeman scheme to solve for the motions in the solar system. We use this method because it conserves total energy with time. A much simpler method, known as Direct Euler or Forward Euler, is used to solve many systems of ordinary differential equations. In the Direct Euler method, the position and velocity are updated as follows:

(1):    $\vec{r}(t + ∆t) = \vec{r}(t) + \vec{v}(t)∆t$

(2):    $\vec{v}(t + ∆t) = \vec{v}(t) + \vec{a}(t)∆t$

Modify your simulation code to use the Direct Euler method for updating position and velocity. This can be done through inheritance, but is not necessary. Run this version of your code for a few hundred years with the same timestep and compare graphically the evolution of total energy vs. time with the Beeman method.

#### 5. Doomsday Planetary Alignment
How often do planetary alignments occur? Modify your code to detect planetary alignments based on all planets being within some threshold (say $5\degree{}$) of the mean angle. Do this for the five innermost planets (i.e., out to and including Jupiter). You may optionally add some or all of the three outermost planets.

#### 6. Jupiter’s Influence on Planetary Orbits
Modify your code such that the planets only feel the gravity of the Sun and not any of the other planets. How does this affect the periods of the four inner planets?

## Project task: Report

Write up the work you have done as a project report.

### Guideline - Writing your project report

Important information on the required structure and length of your report together with advice on what should be included, style etc. is given in the section Guidance on Writing the Project Report.